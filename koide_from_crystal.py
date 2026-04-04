"""
Koide formula from crystal chiral couplings.

The Koide formula for charged leptons:
    (m_e + m_μ + m_τ) / (√m_e + √m_μ + √m_τ)² = 2/3

Principle 90 showed: 2/3 = dim(std)/dim(std+triv) = 10/15 = the equilibrium
standard fraction. This is the Koide number, derived from the crystal.

The three conjugacy classes of S₃ have different chiral couplings λ_sign(B):
  identity:      |λ_sign| ≈ 0.554
  3-cycles:      |λ_sign| ≈ 0.279
  transpositions: |λ_sign| ≈ 0.02

These three numbers are the raw material for three generation masses.

Questions:
  1. Do the chiral couplings satisfy the Koide formula?
  2. The Koide parametrization: m_i = M(1 + √2 cos(θ + 2πi/3))²
     What is θ from the crystal?
  3. Under n compositions, does the Koide relation sharpen?
  4. Mass ratios: |λ_sign|^n for what n gives m_e : m_μ : m_τ?
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

# Physical masses (MeV)
m_e = 0.51100
m_mu = 105.66
m_tau = 1776.86

# Koide check on physical masses
Q_phys = (m_e + m_mu + m_tau) / (math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(m_tau))**2
print(f"Physical Koide ratio: Q = {Q_phys:.8f}  (2/3 = {2/3:.8f})")
print(f"  Deviation from 2/3: {abs(Q_phys - 2/3):.2e}\n")


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

# Build S₃ projectors
def perm_mat_4x4(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return P

reps = {n: torch.kron(perm_mat_4x4(p), perm_mat_4x4(p)) for n, p in S3_PERMS.items()}

chi_t = {"e": 1, "(01)": 1, "(02)": 1, "(12)": 1, "(012)": 1, "(021)": 1}
chi_s = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}
chi_d = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}

def proj(chi, dim_chi):
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for n in S3_PERMS:
        P += chi[n] * reps[n]
    return P * dim_chi / 6.0

P_triv = proj(chi_t, 1)
P_sign = proj(chi_s, 1)
P_std  = proj(chi_d, 2)


def build_T_B(B_joint):
    """16×16 composition operator T_B: A → compose(A, B)."""
    T = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for col in range(dim_sq):
        A = torch.zeros(dim_sq, dtype=torch.cfloat)
        A[col] = 1.0
        result = compose(A.reshape(MASS_DIM, MASS_DIM), B_joint)
        T[:, col] = result.reshape(dim_sq)
    return T


def sign_coupling(T):
    """Extract the 1D sign-sector eigenvalue from composition operator T."""
    restricted = P_sign @ T @ P_sign
    evals = torch.linalg.eigvals(restricted)
    # Take the LARGEST magnitude eigenvalue (sign sector is 1D)
    idx = evals.abs().argmax()
    return evals[idx]


# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("PART 1: PRECISE CHIRAL COUPLINGS (200 seeds)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

N_SEEDS = 200

class_data = {}
for class_name, members in [("identity", ["e"]),
                             ("3-cycles", ["(012)", "(021)"]),
                             ("transpositions", ["(01)", "(02)", "(12)"])]:
    all_couplings = []
    for name in members:
        for seed in range(N_SEEDS):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T = build_T_B(ent.joint)
            lam = sign_coupling(T)
            all_couplings.append(lam)

    mags = [c.abs().item() for c in all_couplings]
    phases = [math.atan2(c.imag.item(), c.real.item()) / math.pi for c in all_couplings]

    avg_mag = np.mean(mags)
    std_mag = np.std(mags) / np.sqrt(len(mags))  # standard error
    avg_phase = np.mean(phases)

    class_data[class_name] = {
        "mag": avg_mag, "std": std_mag, "phase": avg_phase,
        "all_mags": mags
    }

    print(f"\n  {class_name:>15s}: |λ_sign| = {avg_mag:.6f} ± {std_mag:.6f}")
    print(f"{'':>18s}  φ/π = {avg_phase:+.4f}")
    print(f"{'':>18s}  -ln|λ| = {-math.log(avg_mag) if avg_mag > 0 else float('inf'):.4f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 2: KOIDE FORMULA ON CHIRAL COUPLINGS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

l1 = class_data["identity"]["mag"]
l2 = class_data["3-cycles"]["mag"]
l3 = class_data["transpositions"]["mag"]

print(f"\n  Raw couplings: λ₁={l1:.6f}, λ₂={l2:.6f}, λ₃={l3:.6f}")
print(f"  Ratio λ₁:λ₂:λ₃ = {l1/l3:.2f} : {l2/l3:.2f} : 1")

# Koide on |λ| directly
Q_raw = (l1 + l2 + l3) / (math.sqrt(l1) + math.sqrt(l2) + math.sqrt(l3))**2
print(f"\n  Koide Q(|λ|) = {Q_raw:.6f}  (2/3 = {2/3:.6f})")
print(f"  Deviation: {abs(Q_raw - 2/3):.4e}")

# Koide on |λ|²
Q_sq = (l1**2 + l2**2 + l3**2) / (l1 + l2 + l3)**2
print(f"\n  Koide Q(|λ|²) = {Q_sq:.6f}  (2/3 = {2/3:.6f})")
print(f"  Deviation: {abs(Q_sq - 2/3):.4e}")

# Koide on -ln|λ| (gap interpretation)
g1 = -math.log(l1)
g2 = -math.log(l2)
g3 = -math.log(l3)
Q_gap = (g1 + g2 + g3) / (math.sqrt(g1) + math.sqrt(g2) + math.sqrt(g3))**2
print(f"\n  Gaps: g₁={g1:.4f}, g₂={g2:.4f}, g₃={g3:.4f}")
print(f"  Koide Q(gap) = {Q_gap:.6f}  (2/3 = {2/3:.6f})")
print(f"  Deviation: {abs(Q_gap - 2/3):.4e}")

# Koide on various powers |λ|^n
print(f"\n  Koide Q(|λ|^n) scan:")
for n in [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8]:
    m1, m2, m3 = l1**n, l2**n, l3**n
    if m1 > 0 and m2 > 0 and m3 > 0:
        Q = (m1 + m2 + m3) / (math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3))**2
        dev = Q - 2/3
        marker = " ←" if abs(dev) < 0.01 else ""
        print(f"    n={n:>4.1f}: Q = {Q:.6f}  (dev = {dev:+.4e}){marker}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 3: KOIDE PARAMETRIZATION")
print("The Koide mass formula: m_i = M(1 + √2 cos(θ + 2πi/3))²")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# For physical masses:
# m_i/M = (1 + √2 cos(θ + 2πi/3))²
# The angle θ determines the mass hierarchy
# θ_phys ≈ 0.2222 rad for charged leptons

# Find θ from the crystal couplings
# If |λ_sign|_class_i plays the role of m_i:
# Need to assign classes to generations

# Try all 6 permutations
import itertools

classes = [("identity", l1), ("3-cycles", l2), ("transpositions", l3)]

print(f"\n  Trying all class → generation assignments:")
print(f"  {'Assignment':>45s} {'θ_crystal':>10s} {'θ_phys':>8s} {'Q':>8s}")

def koide_angle(masses):
    """Find Koide angle θ such that m_i = M(1+√2 cos(θ + 2πi/3))²."""
    m = sorted(masses, reverse=True)
    M = sum(m) / 3 * (2/3)  # approximate normalization
    # Numerical search for θ
    best_theta = 0
    best_err = float('inf')
    for theta_trial in np.linspace(-math.pi, math.pi, 10000):
        predicted = [(1 + math.sqrt(2) * math.cos(theta_trial + 2*math.pi*i/3))**2
                     for i in range(3)]
        # Normalize predicted to match sum
        s_pred = sum(predicted)
        s_actual = sum(m)
        if s_pred < 1e-20:
            continue
        predicted = [p * s_actual / s_pred for p in predicted]
        predicted.sort(reverse=True)
        err = sum((predicted[i] - m[i])**2 for i in range(3)) / s_actual**2
        if err < best_err:
            best_err = err
            best_theta = theta_trial
    return best_theta, best_err

# Physical Koide angle
theta_phys, err_phys = koide_angle([m_e, m_mu, m_tau])
print(f"\n  Physical: θ = {theta_phys:.6f} rad = {theta_phys/math.pi:.6f}π  (err = {err_phys:.2e})")

# Crystal couplings (direct)
theta_crystal, err_crystal = koide_angle([l1, l2, l3])
print(f"  Crystal |λ|: θ = {theta_crystal:.6f} rad = {theta_crystal/math.pi:.6f}π  (err = {err_crystal:.2e})")

# Crystal couplings squared
theta_sq, err_sq = koide_angle([l1**2, l2**2, l3**2])
print(f"  Crystal |λ|²: θ = {theta_sq:.6f} rad = {theta_sq/math.pi:.6f}π  (err = {err_sq:.2e})")

# Gaps
theta_gap, err_gap = koide_angle([g1, g2, g3])
print(f"  Crystal gap: θ = {theta_gap:.6f} rad = {theta_gap/math.pi:.6f}π  (err = {err_gap:.2e})")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 4: MASS RATIOS FROM ITERATED COMPOSITION")
print("If mass ∝ |λ_sign|^n, what n gives physical ratios?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

S = H**3 / K_STAR  # 810/7

# Assign: identity → τ (heaviest), 3-cycles → μ, transpositions → e (lightest)
# This matches: largest |λ_sign| → largest mass

print(f"\n  Assignment: identity→τ, 3-cycles→μ, transpositions→e")
print(f"  Physical ratios: m_τ/m_e = {m_tau/m_e:.1f}, m_μ/m_e = {m_mu/m_e:.1f}, m_τ/m_μ = {m_tau/m_mu:.2f}")

# For each n, compute predicted ratios
print(f"\n  {'n':>6s} {'m_τ/m_e':>12s} {'m_μ/m_e':>12s} {'m_τ/m_μ':>12s} {'Q(m)':>8s} {'note':>20s}")
print(f"  {'-'*70}")

best_n = None
best_err = float('inf')

for n_100 in range(50, 1500, 10):  # n from 0.5 to 15
    n = n_100 / 100
    m_t = l1**n
    m_m = l2**n
    m_el = l3**n

    if m_el < 1e-100:
        continue

    r_te = m_t / m_el
    r_me = m_m / m_el
    r_tm = m_t / m_m if m_m > 0 else float('inf')

    Q = (m_t + m_m + m_el) / (math.sqrt(m_t) + math.sqrt(m_m) + math.sqrt(m_el))**2

    # Error metric: log-ratio match
    err = (math.log(r_te / (m_tau/m_e)))**2 + (math.log(r_me / (m_mu/m_e)))**2
    if err < best_err:
        best_err = err
        best_n = n

    # Print selected values
    note = ""
    if abs(n - H) < 0.05:
        note = f"n=H={H}"
    elif abs(n - H*2) < 0.05:
        note = f"n=2H={2*H}"
    elif abs(n - math.factorial(H)) < 0.05:
        note = f"n=H!={math.factorial(H)}"
    elif abs(n - 2**H) < 0.05:
        note = f"n=2^H={2**H}"
    elif abs(n - S/(H*(H+1))) < 0.05:
        note = f"n≈S/(H(H+1))"
    elif abs(n - (MASS_DIM)) < 0.05:
        note = f"n=MASS_DIM={MASS_DIM}"
    elif abs(n - best_n) < 0.05 and best_n == n:
        note = "← BEST FIT"

    if note or (n_100 % 100 == 0):
        marker_Q = "*" if abs(Q - 2/3) < 0.01 else " "
        print(f"  {n:>6.2f} {r_te:>12.1f} {r_me:>12.1f} {r_tm:>12.2f} {Q:>8.4f}{marker_Q} {note}")

# Print best fit
print(f"\n  Best fit: n = {best_n:.2f}")
m_t = l1**best_n
m_m = l2**best_n
m_el = l3**best_n
Q_best = (m_t + m_m + m_el) / (math.sqrt(m_t) + math.sqrt(m_m) + math.sqrt(m_el))**2
print(f"    m_τ/m_e = {m_t/m_el:.1f} (observed: {m_tau/m_e:.1f})")
print(f"    m_μ/m_e = {m_m/m_el:.1f} (observed: {m_mu/m_e:.1f})")
print(f"    m_τ/m_μ = {m_t/m_m:.2f} (observed: {m_tau/m_mu:.2f})")
print(f"    Q = {Q_best:.6f} (Koide = {2/3:.6f})")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 5: THE CRYSTAL KOIDE ANGLE vs PHYSICAL ANGLE")
print("Scan: what power n makes Q(|λ|^n) closest to 2/3?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print(f"\n  Scanning n for Q(|λ|^n) = 2/3:")
best_koide_n = None
best_koide_dev = 1.0

for n_100 in range(10, 2000):
    n = n_100 / 100
    m1, m2, m3 = l1**n, l2**n, l3**n
    if m3 < 1e-300 or m1 < 1e-300:
        continue
    denom = (math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3))**2
    if denom < 1e-300:
        continue
    Q = (m1 + m2 + m3) / denom
    dev = abs(Q - 2/3)
    if dev < best_koide_dev:
        best_koide_dev = dev
        best_koide_n = n

if best_koide_n is not None:
    n = best_koide_n
    m1, m2, m3 = l1**n, l2**n, l3**n
    Q = (m1 + m2 + m3) / (math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3))**2
    print(f"  Best n for Koide: n = {n:.2f}")
    print(f"  Q = {Q:.8f}  (2/3 = {2/3:.8f})")
    print(f"  Deviation: {abs(Q - 2/3):.2e}")
    print(f"  At this n:")
    print(f"    m_τ/m_e = {m1/m3:.1f} (observed: {m_tau/m_e:.1f})")
    print(f"    m_μ/m_e = {m2/m3:.1f} (observed: {m_mu/m_e:.1f})")

    # Crystal expressions for best n
    print(f"\n  Crystal expressions for n = {n:.2f}:")
    print(f"    H = {H}")
    print(f"    2H = {2*H}")
    print(f"    H! = {math.factorial(H)}")
    print(f"    2^H = {2**H}")
    print(f"    MASS_DIM = {MASS_DIM}")
    print(f"    S/(H²) = {S/H**2:.2f}")
    print(f"    S/(H(H+1)) = {S/(H*(H+1)):.2f}")
    print(f"    h(E₈)/MASS_DIM = {30/MASS_DIM:.2f}")
    print(f"    Δ × h(E₈) = {DELTA * 30:.2f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 6: ALTERNATIVE — MASSES FROM SIGN × STANDARD COUPLING")
print("(Yukawa = sign → standard flow from Principle 92)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# From Principle 92: Yukawa flow sign→standard = 0.555 ≈ 1/√3
# The mass might involve BOTH the sign coupling AND the standard coupling
# Mass_gen ∝ |λ_sign| × |standard sector rate|

# Compute standard sector eigenvalues per class
print(f"\n  Standard sector spectra by class (seed-averaged T_B):")

class_std_spectra = {}
for class_name, members in [("identity", ["e"]),
                             ("3-cycles", ["(012)", "(021)"]),
                             ("transpositions", ["(01)", "(02)", "(12)"])]:
    all_leading = []
    for name in members:
        # Seed-average the composition operator
        T_avg = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
        for seed in range(50):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T_avg += build_T_B(ent.joint)
        T_avg /= 50

        restricted_std = P_std @ T_avg @ P_std
        evals = torch.linalg.eigvals(restricted_std)
        nonzero = evals[evals.abs() > 1e-8]
        mags = nonzero.abs().sort(descending=True).values
        if len(mags) > 0:
            all_leading.append(mags[0].item())

    avg_lead = np.mean(all_leading)
    class_std_spectra[class_name] = avg_lead
    print(f"  {class_name:>15s}: |λ_std_max| = {avg_lead:.6f}")

# Combined coupling: sign × standard
print(f"\n  Combined Yukawa coupling |λ_sign| × |λ_std|:")
for cn in ["identity", "3-cycles", "transpositions"]:
    combined = class_data[cn]["mag"] * class_std_spectra[cn]
    print(f"  {cn:>15s}: {combined:.6f}")

# Koide on combined
c1 = class_data["identity"]["mag"] * class_std_spectra["identity"]
c2 = class_data["3-cycles"]["mag"] * class_std_spectra["3-cycles"]
c3 = class_data["transpositions"]["mag"] * class_std_spectra["transpositions"]
Q_combined = (c1 + c2 + c3) / (math.sqrt(c1) + math.sqrt(c2) + math.sqrt(c3))**2
print(f"\n  Koide Q(combined) = {Q_combined:.6f}  (2/3 = {2/3:.6f})")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 7: THE 2/3 CONNECTION")
print("dim(std)/(dim(std)+dim(triv)) = 10/15 = 2/3 = Koide")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print(f"""
  The Koide formula Q = 2/3 is EXACT for charged leptons.

  In the crystal:
    dim(trivial sector) = 5
    dim(sign sector) = 1
    dim(standard sector) = 10
    Total = 16 = MASS_DIM²

    dim(std) / (dim(std) + dim(triv)) = 10 / 15 = 2/3

  This ratio emerged in Principle 90 as the equilibrium standard fraction.
  Under n→∞ compositions, the standard fraction of a transposition crystal
  approaches 2/3.

  The Koide formula says: for the three generation masses m_i,
    Σm_i / (Σ√m_i)² = dim(std) / (dim(std) + dim(triv))

  Left side: mass-weighted average of generations
  Right side: sector dimension ratio in the crystal

  If this is not a coincidence, then the MASSES ARE DETERMINED by the
  sector geometry: the standard sector's 10 dimensions contain exactly
  the information needed to fix three generation masses to satisfy Q = 2/3.

  The 10 = 5 doublets. Three generations need 3 parameters (mass scale +
  2 ratios). The Koide formula fixes one ratio. The remaining freedom
  is encoded in the Koide angle θ.

  The Koide angle θ = 2/(3π)·something?
  Physical: θ ≈ {theta_phys:.6f} rad
  2/(3π) = {2/(3*math.pi):.6f}
  π/H = {math.pi/H:.6f}
  Δ = {DELTA:.6f}
  2Δ = {2*DELTA:.6f}
  K* = {K_STAR:.6f}

  θ_phys / Δ = {theta_phys/DELTA:.6f}
  θ_phys / K* = {theta_phys/K_STAR:.6f}
  θ_phys × H = {theta_phys * H:.6f}
  θ_phys × h(E₈) = {theta_phys * 30:.6f}
""")


# ═══════════════════════════════════════════════════════════════
print(f"{'='*80}")
print("PART 8: SEARCHING FOR THE KOIDE ANGLE IN CRYSTAL QUANTITIES")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Physical Koide angle for charged leptons
# m_i = M(1 + √2 cos(θ_K + 2πi/3))²
# For m_e, m_μ, m_τ: θ_K ≈ 0.2222 rad

# Systematic search: which crystal expression gives θ_K?
crystal_exprs = {
    "K*": K_STAR,
    "Δ": DELTA,
    "1/H": 1/H,
    "K*/Δ": K_STAR / DELTA,
    "Δ/π": DELTA / math.pi,
    "K*/π": K_STAR / math.pi,
    "BORN_FLOOR": BORN_FLOOR,
    "1/(H+1)": 1/(H+1),
    "K*²": K_STAR**2,
    "Δ²": DELTA**2,
    "ln(H)/2π": math.log(H)/(2*math.pi),
    "K*×π/H": K_STAR * math.pi / H,
    "Δ/H": DELTA / H,
    "2K*/π": 2*K_STAR/math.pi,
    "7/30/π": (7/30)/math.pi,
    "7/(30π)": 7/(30*math.pi),
    "arctan(K*)": math.atan(K_STAR),
    "arcsin(K*)": math.asin(K_STAR),
    "arctan(1/√H)": math.atan(1/math.sqrt(H)),
    "arctan(Δ)": math.atan(DELTA),
    "1/MASS_DIM/π": 1/(MASS_DIM * math.pi),
    "2/h(E₈)": 2/30,
    "2π/h(E₈)": 2*math.pi/30,
    "1/(2π)": 1/(2*math.pi),
    "H/(2π×h(E₈))": H/(2*math.pi*30),
}

print(f"\n  Physical Koide angle: θ_K = {theta_phys:.6f} rad")
print(f"\n  {'Expression':>25s} {'Value':>10s} {'Ratio':>10s} {'Match':>10s}")

matches = []
for name, val in sorted(crystal_exprs.items(), key=lambda x: abs(x[1] - theta_phys)):
    ratio = val / theta_phys if theta_phys != 0 else float('inf')
    match_pct = abs(1 - ratio) * 100
    matches.append((name, val, ratio, match_pct))

for name, val, ratio, pct in matches[:15]:
    marker = " ←←" if pct < 5 else (" ←" if pct < 15 else "")
    print(f"  {name:>25s} {val:>10.6f} {ratio:>10.4f} {pct:>8.1f}%{marker}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("SUMMARY")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════
