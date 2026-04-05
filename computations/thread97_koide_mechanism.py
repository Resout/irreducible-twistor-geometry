"""
Thread 97: Can composition + Coxeter rotation produce θ = 2/9 directly?

The algebraic derivation gives θ_phys = θ_crystal - 1/h(E₈).
But the crystal COMPUTES eigenvalues dynamically. Does iterated
composition naturally shift the angle from 23/90 to 20/90?

Approach:
1. Build sign-sector eigenvalues under N compositions
2. Track the Koide angle as N grows
3. Look for the Coxeter correction emerging dynamically
4. Check if composition orbits have period 30 (= h(E₈))
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import numpy as np
from solver.algebra import H, MASS_DIM, K_STAR, BORN_FLOOR, DELTA
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

dim_sq = MASS_DIM ** 2

# S₃ structure
S3_PERMS = {
    "e":     [0, 1, 2],
    "(01)":  [1, 0, 2],
    "(02)":  [2, 1, 0],
    "(12)":  [0, 2, 1],
    "(012)": [1, 2, 0],
    "(021)": [2, 0, 1],
}

S3_CORR = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

def perm_mat_4x4(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return P

reps = {n: torch.kron(perm_mat_4x4(p), perm_mat_4x4(p)) for n, p in S3_PERMS.items()}

chi_s = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}
chi_t = {"e": 1, "(01)": 1, "(02)": 1, "(12)": 1, "(012)": 1, "(021)": 1}
chi_d = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}

def proj(chi, dim_chi):
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for n in S3_PERMS:
        P += chi[n] * reps[n]
    return P * dim_chi / 6.0

P_sign = proj(chi_s, 1)
P_triv = proj(chi_t, 1)
P_std = proj(chi_d, 2)

def build_T_B(B_joint):
    T = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for col in range(dim_sq):
        A = torch.zeros(dim_sq, dtype=torch.cfloat)
        A[col] = 1.0
        result = compose(A.reshape(MASS_DIM, MASS_DIM), B_joint)
        T[:, col] = result.reshape(dim_sq)
    return T

def sign_coupling(T):
    restricted = P_sign @ T @ P_sign
    evals = torch.linalg.eigvals(restricted)
    idx = evals.abs().argmax()
    return evals[idx]

def koide_angle_from_masses(m1, m2, m3):
    """Koide angle from three positive quantities."""
    m = sorted([m1, m2, m3], reverse=True)
    sqM = (math.sqrt(m[0]) + math.sqrt(m[1]) + math.sqrt(m[2])) / 3
    if sqM < 1e-30:
        return float('nan')
    ratio = math.sqrt(m[0]) / sqM
    cos_theta = (ratio - 1) / math.sqrt(2)
    if abs(cos_theta) > 1:
        return float('nan')
    return math.acos(cos_theta)

def koide_Q(m1, m2, m3):
    s = math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    return (m1 + m2 + m3) / s**2 if s > 1e-30 else float('nan')


# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("PART 1: SIGN-SECTOR EIGENVALUES UNDER ITERATED COMPOSITION")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Build base crystals for each conjugacy class
N_SEEDS = 100

print("\n  Building seed-averaged composition operators...")
class_T = {}
class_joints = {}
for class_name, members in [("identity", ["e"]),
                             ("3-cycles", ["(012)", "(021)"]),
                             ("transpositions", ["(01)", "(02)", "(12)"])]:
    T_avg = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    joint_avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    count = 0
    for name in members:
        for seed in range(N_SEEDS):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T_avg += build_T_B(ent.joint)
            joint_avg += ent.joint
            count += 1
    T_avg /= count
    joint_avg /= count
    class_T[class_name] = T_avg
    class_joints[class_name] = joint_avg

# Now: compose identity crystal with itself N times, track sign eigenvalue
print("\n  Iterated self-composition of each conjugacy class crystal:")
print(f"  {'N':>4s}  {'|λ_s(id)|':>12s}  {'|λ_s(3c)|':>12s}  {'|λ_s(tr)|':>12s}  {'θ_Koide':>10s}  {'Q':>8s}")
print(f"  {'-'*70}")

for N_comp in [1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30]:
    # T_B^N = T_B composed N times
    lambdas = {}
    for cn in ["identity", "3-cycles", "transpositions"]:
        T_N = class_T[cn].clone()
        for _ in range(N_comp - 1):
            T_N = T_N @ class_T[cn]
        lam = sign_coupling(T_N)
        lambdas[cn] = lam.abs().item()

    l1, l2, l3 = lambdas["identity"], lambdas["3-cycles"], lambdas["transpositions"]

    if l1 > 1e-30 and l2 > 1e-30 and l3 > 1e-30:
        Q = koide_Q(l1, l2, l3)
        theta = koide_angle_from_masses(l1, l2, l3)
        print(f"  {N_comp:>4d}  {l1:>12.6f}  {l2:>12.6f}  {l3:>12.6f}  {theta:>10.6f}  {Q:>8.5f}")
    else:
        print(f"  {N_comp:>4d}  {l1:>12.2e}  {l2:>12.2e}  {l3:>12.2e}  {'---':>10s}  {'---':>8s}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 2: COMPOSITION ORBIT STRUCTURE — PERIOD DETECTION")
print("Track λ_sign under T^N for each class, look for periodicity")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print("\n  Phase of sign eigenvalue under iteration (mod π):")
print(f"  {'N':>4s}  {'φ_id/π':>10s}  {'φ_3c/π':>10s}  {'φ_tr/π':>10s}  {'Δφ_id/π':>10s}")
print(f"  {'-'*55}")

prev_phase_id = None
for N_comp in range(1, 31):
    phases = {}
    for cn in ["identity", "3-cycles", "transpositions"]:
        T_N = class_T[cn].clone()
        for _ in range(N_comp - 1):
            T_N = T_N @ class_T[cn]
        lam = sign_coupling(T_N)
        phase = math.atan2(lam.imag.item(), lam.real.item()) / math.pi
        phases[cn] = phase

    delta_id = ""
    if prev_phase_id is not None:
        d = phases["identity"] - prev_phase_id
        # Wrap to [-1, 1]
        while d > 1: d -= 2
        while d < -1: d += 2
        delta_id = f"{d:>10.6f}"
    prev_phase_id = phases["identity"]

    if N_comp <= 10 or N_comp % 5 == 0:
        print(f"  {N_comp:>4d}  {phases['identity']:>10.6f}  {phases['3-cycles']:>10.6f}  {phases['transpositions']:>10.6f}  {delta_id}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 3: CROSS-CLASS COMPOSITION — THE MIXING EXPERIMENT")
print("Compose identity × 3-cycle and identity × transposition")
print("Does the mixed composition operator produce different θ?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Build cross-composition operators
# T_{B₁∘B₂}: what happens when you compose with B₁ then B₂?
print("\n  Cross-class composition: T_class1 @ T_class2")
print(f"  {'Composition':>30s}  {'|λ_sign|':>10s}  {'exponent':>10s}")
print(f"  {'-'*55}")

for cn1 in ["identity", "3-cycles", "transpositions"]:
    for cn2 in ["identity", "3-cycles", "transpositions"]:
        T_cross = class_T[cn1] @ class_T[cn2]
        lam = sign_coupling(T_cross)
        mag = lam.abs().item()
        # Express as 6^{-a}
        if mag > 1e-30:
            a = -math.log(mag) / math.log(6)
            print(f"  {cn1+' × '+cn2:>30s}  {mag:>10.6f}  6^{{-{a:.4f}}}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 4: THE COXETER ROTATION")
print("Does the sign eigenvalue rotate by 2π/30 per composition step?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# For the identity crystal, track the full complex sign eigenvalue
print("\n  Identity crystal: sign eigenvalue orbit")
print(f"  {'N':>4s}  {'|λ|':>10s}  {'φ/π':>10s}  {'φ/(2π/30)':>12s}  {'log|λ|/logH!':>14s}")
print(f"  {'-'*60}")

T_id = class_T["identity"]
T_N = torch.eye(dim_sq, dtype=torch.cfloat)

for N in range(31):
    if N > 0:
        T_N = T_N @ T_id
    lam = sign_coupling(T_N)
    mag = lam.abs().item()
    phase = math.atan2(lam.imag.item(), lam.real.item()) / math.pi

    coxeter_units = phase * math.pi / (2 * math.pi / 30) if N > 0 else 0
    log_ratio = math.log(mag) / math.log(math.factorial(H)) if mag > 1e-30 else float('inf')

    if N <= 10 or N % 5 == 0:
        print(f"  {N:>4d}  {mag:>10.6f}  {phase:>10.6f}  {coxeter_units:>12.4f}  {log_ratio:>14.4f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 5: KOIDE ANGLE FROM SIGN-SECTOR EIGENVALUES")
print("Use |λ_sign|^n at the Koide exponent n=35/H² and check θ")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Get precise sign couplings
l_id = sign_coupling(class_T["identity"]).abs().item()
l_3c = sign_coupling(class_T["3-cycles"]).abs().item()
l_tr = sign_coupling(class_T["transpositions"]).abs().item()

print(f"\n  Sign-sector eigenvalues (seed-averaged T):")
print(f"    identity:      |λ| = {l_id:.6f}  (6^{{-1/3}} = {6**(-1/3):.6f})")
print(f"    3-cycles:      |λ| = {l_3c:.6f}  (6^{{-a}}, a = {-math.log(l_3c)/math.log(6):.6f})")
print(f"    transpositions: |λ| = {l_tr:.6f}  (6^{{-59/30}} = {6**(-59/30):.6f})")

# The exponents
a0 = -math.log(l_id) / math.log(6)
a1 = -math.log(l_3c) / math.log(6) if l_3c > 0 else float('inf')
a2 = -math.log(l_tr) / math.log(6) if l_tr > 0 else float('inf')

print(f"\n  Exponents (|λ| = 6^{{-a}}):")
print(f"    a₀(identity)       = {a0:.6f}  (pred: 1/3 = {1/3:.6f})")
print(f"    a₁(3-cycles)       = {a1:.6f}")
print(f"    a₂(transpositions) = {a2:.6f}  (pred: 59/30 = {59/30:.6f})")

# For various candidate n values, compute Koide angle and Q
print(f"\n  Koide angle scan over exponent n:")
print(f"  {'n':>8s}  {'θ':>10s}  {'θ×H²':>8s}  {'Q':>8s}  {'note':>20s}")
print(f"  {'-'*60}")

n_koide = 35 / H**2  # The discovered Koide exponent

for n, note in [(1, "raw"),
                (2, "squared"),
                (n_koide, f"35/H²={n_koide:.4f}"),
                (n_koide - 1/30, f"35/H²-1/h(E₈)"),
                (n_koide + 1/30, f"35/H²+1/h(E₈)"),
                (H, "H"),
                (2*H, "2H"),
                (4, "MASS_DIM"),
                (math.pi, "π"),
                (math.e, "e"),
                (35/10, "35/dim(std)"),
                (35/16, "35/dim²"),
                (1/DELTA, "1/Δ"),
                (2/DELTA, "2/Δ")]:
    m1 = l_id ** n
    m2 = l_3c ** n
    m3 = l_tr ** n
    if m1 > 1e-300 and m2 > 1e-300 and m3 > 1e-300:
        Q = koide_Q(m1, m2, m3)
        theta = koide_angle_from_masses(m1, m2, m3)
        marker = " ←←" if abs(theta * H**2 - 2) < 0.05 else (" ←" if abs(Q - 2/3) < 0.01 else "")
        print(f"  {n:>8.4f}  {theta:>10.6f}  {theta*H**2:>8.4f}  {Q:>8.5f}  {note:>20s}{marker}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 6: THE 3-CYCLE EXPONENT — SEARCHING FOR a₂")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print(f"\n  Measured: a₁(3-cycles) = {a1:.8f}")
print(f"\n  Candidate expressions:")

candidates = {
    "5/7": 5/7,
    "1/√2": 1/math.sqrt(2),
    "ln(2)": math.log(2),
    "2/e": 2/math.e,
    "(H-1)/H!": (H-1)/math.factorial(H),
    "1-1/3": 2/3,
    "(1+1/H)/H": (1+1/H)/H,
    "K*+1/2": K_STAR + 0.5,
    "Δ/ln(H!)": DELTA / math.log(math.factorial(H)),
    "1-K*/(H-1)": 1 - K_STAR/(H-1),
    "(H²-1)/(H²+H+1)": (H**2-1)/(H**2+H+1),
    "a₀+a₂/Coxeter": a0 + a2/30,
    "√(a₀×a₂)": math.sqrt(a0 * a2),
    "(a₀+a₂)/H": (a0 + a2) / H,
    "(a₀×a₂)^(1/3)": (a0 * a2) ** (1/3),
    "a₀+Δ/ln6": a0 + DELTA / math.log(6),
    "1/(H-K*)": 1/(H - K_STAR),
    "(H-1)/(H+1/H)": (H-1)/(H+1/H),
    "2a₀+1/h(E₈)": 2*a0 + 1/30,
    "a₀×(a₂/a₀)^(1/2)": a0 * (a2/a0)**0.5,
    "H/(H²+1/H)": H/(H**2+1/H),
    "(2H-1)/(2H+1)": (2*H-1)/(2*H+1),
    "7/10": 7/10,
    "5/7": 5/7,
    "43/60": 43/60,
    "(a₀+a₂)/3": (a0 + a2) / 3,
    "a₂-a₀": a2 - a0,
    "1-a₀": 1 - a0,
    "(h(E₈)-a₂×h(E₈))/h(E₈)": 1 - a2,  # same as above conceptually
    "K*×H": K_STAR * H,
    "1-K*×H": 1 - K_STAR * H,
    "(2-a₂)/(2+a₀)": (2-a2)/(2+a0) if abs(2+a0) > 1e-10 else float('inf'),
    "ln3/ln6": math.log(3)/math.log(6),
    "ln(H)/ln(H!)": math.log(H)/math.log(math.factorial(H)),
    "1/ln6": 1/math.log(6),
    "(H²-2)/H²": (H**2-2)/H**2,
}

# Sort by closeness
ranked = sorted(candidates.items(), key=lambda x: abs(x[1] - a1))

for name, val in ranked[:15]:
    pct = abs(val - a1) / a1 * 100
    print(f"    {name:>25s} = {val:.8f}  (dev = {pct:.4f}%)")

# Also check rational approximation
print(f"\n  Best rational approximations p/q (q ≤ 100):")
best_rats = []
for q in range(1, 101):
    p = round(a1 * q)
    if p > 0:
        dev = abs(p/q - a1)
        best_rats.append((p, q, dev))
best_rats.sort(key=lambda x: x[2])
for p, q, dev in best_rats[:10]:
    print(f"    {p}/{q} = {p/q:.8f}  (dev = {dev:.6e}, {dev/a1*100:.4f}%)")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 7: RELATIONSHIP BETWEEN THE THREE EXPONENTS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print(f"\n  a₀ = {a0:.8f}")
print(f"  a₁ = {a1:.8f}")
print(f"  a₂ = {a2:.8f}")
print(f"\n  Sum: a₀+a₁+a₂ = {a0+a1+a2:.8f}")
print(f"  Product: a₀×a₁×a₂ = {a0*a1*a2:.8f}")
print(f"  Ratios: a₂/a₀ = {a2/a0:.6f}, a₁/a₀ = {a1/a0:.6f}, a₂/a₁ = {a2/a1:.6f}")

# Koide on the exponents themselves
Q_exp = (a0 + a1 + a2) / (math.sqrt(a0) + math.sqrt(a1) + math.sqrt(a2))**2
theta_exp = koide_angle_from_masses(a0, a1, a2)
print(f"\n  Koide Q(a₀, a₁, a₂) = {Q_exp:.8f}  (2/3 = {2/3:.8f})")
print(f"  Koide θ(a₀, a₁, a₂) = {theta_exp:.8f}")
print(f"  θ × H² = {theta_exp * H**2:.6f}")
print(f"  θ / (2/H²) = {theta_exp / (2/H**2):.6f}")

# Check if a₁ is determined by a₀, a₂ and the Koide relation
# If Q = 2/3, then a₁ is determined
print(f"\n  If Q(a) = 2/3 with a₀={a0:.6f}, a₂={a2:.6f}:")
# Q = (a₀+a₁+a₂) / (√a₀+√a₁+√a₂)² = 2/3
# Let x = √a₁. Then: (a₀+x²+a₂)/(√a₀+x+√a₂)² = 2/3
# 3(a₀+x²+a₂) = 2(√a₀+x+√a₂)²
# 3a₀+3x²+3a₂ = 2(a₀+x²+a₂+2x√a₀+2x√a₂+2√(a₀a₂))
# 3a₀+3x²+3a₂ = 2a₀+2x²+2a₂+4x√a₀+4x√a₂+4√(a₀a₂)
# x² + (a₀+a₂) - 4x(√a₀+√a₂) - 4√(a₀a₂) = 0
# x² - 4(√a₀+√a₂)x + (a₀+a₂-4√(a₀a₂)) = 0
sa0, sa2 = math.sqrt(a0), math.sqrt(a2)
A = 1
B = -4*(sa0 + sa2)
C = a0 + a2 - 4*sa0*sa2
disc = B**2 - 4*A*C
if disc >= 0:
    x1 = (-B + math.sqrt(disc)) / (2*A)
    x2 = (-B - math.sqrt(disc)) / (2*A)
    a1_pred_1 = x1**2
    a1_pred_2 = x2**2
    print(f"  Solution 1: a₁ = {a1_pred_1:.8f}  (measured: {a1:.8f}, dev: {abs(a1_pred_1-a1)/a1*100:.4f}%)")
    print(f"  Solution 2: a₁ = {a1_pred_2:.8f}  (measured: {a1:.8f}, dev: {abs(a1_pred_2-a1)/a1*100:.4f}%)")


print(f"\n\n{'='*80}")
print("SUMMARY")
print("=" * 80)
print(f"""
  Thread 97 findings will go here after computation.
""")
