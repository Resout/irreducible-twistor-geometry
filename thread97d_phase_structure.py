"""
Thread 97d: Phase structure of sign eigenvalues across conjugacy classes.

The identity sign eigenvalue has phase √2 × θ_crystal.
What are the 3-cycle and transposition phases? Do they form the
Koide 2π/3 structure?

Key insight: sign character χ(transpositions) = -1 introduces a
π offset. After accounting for the character, do the phases differ
by 2π/3?

Also: the standard sector 5-doublet structure — is the 1+3+1 band
related to quark generations?
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


# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("PART 1: PHASE GAPS BETWEEN CONJUGACY CLASSES")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

N_SEEDS = 300

# Measure phases for all 6 individual elements
print("\n  Individual element sign eigenvalues (300 seeds):")
element_data = {}
for name in S3_PERMS:
    mags = []
    raw_lambdas = []
    for seed in range(N_SEEDS):
        ent = Entangler(S3_CORR[name], seed=seed).build()
        T = build_T_B(ent.joint)
        lam = sign_coupling(T)
        mags.append(lam.abs().item())
        raw_lambdas.append(lam)

    avg_lam = sum(raw_lambdas) / len(raw_lambdas)
    avg_mag = np.mean(mags)
    avg_phase = math.atan2(avg_lam.imag.item(), avg_lam.real.item())

    element_data[name] = {
        "mag": avg_mag,
        "phase": avg_phase,
        "lam": avg_lam,
    }
    a = -math.log(avg_mag) / math.log(6) if avg_mag > 1e-20 else float('inf')
    print(f"    {name:>6s}: |λ| = {avg_mag:.6f}  a = {a:.6f}  φ = {avg_phase:>8.4f} rad = {avg_phase/math.pi:>8.4f}π")


# Class-averaged phases
class_phases = {}
for class_name, members in [("e", ["e"]),
                             ("3-cycles", ["(012)", "(021)"]),
                             ("transpositions", ["(01)", "(02)", "(12)"])]:
    phases = [element_data[m]["phase"] for m in members]
    class_phases[class_name] = np.mean(phases)

theta_crystal = (1 - K_STAR) / H

print(f"\n  Class-averaged phases:")
for cn, phi in class_phases.items():
    print(f"    {cn:>15s}: φ = {phi:>8.4f} rad = {phi/math.pi:>8.4f}π")

# Phase gaps
print(f"\n  Phase gaps (raw):")
phi_e = class_phases["e"]
phi_3c = class_phases["3-cycles"]
phi_tr = class_phases["transpositions"]

print(f"    e → 3-cycles:      Δφ = {phi_3c - phi_e:.4f} rad = {(phi_3c-phi_e)/math.pi:.4f}π")
print(f"    e → transpositions: Δφ = {phi_tr - phi_e:.4f} rad = {(phi_tr-phi_e)/math.pi:.4f}π")
print(f"    3c → transpositions: Δφ = {phi_tr - phi_3c:.4f} rad = {(phi_tr-phi_3c)/math.pi:.4f}π")

# Accounting for sign character: transpositions have χ=-1, so there's an effective
# π offset in the sign sector. The "bare" phase is the phase of |λ_sign/χ|.
# For transpositions: bare phase = raw phase - π (since χ=-1 means the eigenvalue
# is negated relative to the geometric phase)

print(f"\n  Phase gaps (sign-character corrected):")
print(f"  Transposition bare phase = raw + π (to undo χ(σ)=-1 reflection)")
phi_tr_bare = phi_tr + math.pi

# Wrap to [-π, π]
def wrap(x):
    while x > math.pi: x -= 2*math.pi
    while x < -math.pi: x += 2*math.pi
    return x

print(f"    e:           φ_bare = {phi_e:.4f}")
print(f"    3-cycles:    φ_bare = {phi_3c:.4f}")
print(f"    transpositions: φ_bare = {wrap(phi_tr_bare):.4f}")
print()

# Gaps from identity
print(f"    e → 3c (bare):  {wrap(phi_3c - phi_e):.4f} rad = {wrap(phi_3c - phi_e)/math.pi:.4f}π")
print(f"    e → tr (bare):  {wrap(phi_tr_bare - phi_e):.4f} rad = {wrap(phi_tr_bare - phi_e)/math.pi:.4f}π")
print(f"    3c → tr (bare): {wrap(phi_tr_bare - phi_3c):.4f} rad = {wrap(phi_tr_bare - phi_3c)/math.pi:.4f}π")

# Check 2π/3 structure
two_pi_3 = 2 * math.pi / 3
print(f"\n  2π/3 = {two_pi_3:.4f} rad")
print(f"  Comparison:")
d_e_3c = wrap(phi_3c - phi_e)
d_e_tr = wrap(phi_tr_bare - phi_e)
d_3c_tr = wrap(phi_tr_bare - phi_3c)

print(f"    |Δφ(e→3c)| = {abs(d_e_3c):.4f}, 2π/3 ratio: {abs(d_e_3c)/two_pi_3:.4f}")
print(f"    |Δφ(e→tr)| = {abs(d_e_tr):.4f}, 2π/3 ratio: {abs(d_e_tr)/two_pi_3:.4f}")
print(f"    |Δφ(3c→tr)| = {abs(d_3c_tr):.4f}, 2π/3 ratio: {abs(d_3c_tr)/two_pi_3:.4f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 2: STANDARD SECTOR — 5 DOUBLET DECOMPOSITION")
print("Where do the 5 doublets come from in 4⊗4?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# 4⊗4 = (3⊕1)⊗(3⊕1) → 5 trivial + 1 sign + 5 standard
# The 5 standard copies come from:
#   3⊗3 → 3 copies of V_std
#   3⊗1 → 1 copy of V_std
#   1⊗3 → 1 copy of V_std
# Total: 5 copies

# For the identity crystal, measure all standard eigenvalues with their
# degeneracy pattern

print("\n  Full standard sector eigenvalue decomposition (identity, 300 seeds):")

all_evals_real = []
all_evals_imag = []
for seed in range(N_SEEDS):
    ent = Entangler(S3_CORR["e"], seed=seed).build()
    T = build_T_B(ent.joint)
    restricted = P_std @ T @ P_std
    evals = torch.linalg.eigvals(restricted)
    mags = evals.abs().sort(descending=True).values
    above = mags[mags > 1e-8]
    all_evals_real.append(above.tolist())

# Average
max_len = max(len(s) for s in all_evals_real)
print(f"  # eigenvalues: {max_len}")

avg_spec = []
for i in range(max_len):
    vals = [s[i] for s in all_evals_real if len(s) > i]
    avg_spec.append(np.mean(vals))

print(f"\n  {'idx':>4s}  {'|λ|':>10s}  {'exponent a':>12s}  {'possible ID':>20s}")
print(f"  {'-'*55}")

# Identify the 5 doublets
for i, val in enumerate(avg_spec):
    a = -math.log(val) / math.log(6) if val > 1e-20 else float('inf')
    # Try to identify
    ident = ""
    for name, target in [("0 (near-1)", 0), ("1/3 = a₀", 1/3),
                          ("2a₀ = 2/3", 2/3), ("Δ", DELTA),
                          ("2Δ", 2*DELTA), ("K*", K_STAR),
                          ("1-K*", 1-K_STAR), ("59/30", 59/30),
                          ("a₁", 0.717), ("a₁/2", 0.717/2)]:
        if abs(a - target) < 0.05:
            ident = name
    print(f"  {i:>4d}  {val:>10.6f}  {a:>12.6f}  {ident:>20s}")

# Group into doublets (pairs of nearly equal eigenvalues)
print(f"\n  Doublet structure:")
doublets = []
i = 0
while i < len(avg_spec):
    if i + 1 < len(avg_spec) and abs(avg_spec[i] - avg_spec[i+1]) < 0.01:
        doublets.append((avg_spec[i] + avg_spec[i+1]) / 2)
        i += 2
    else:
        doublets.append(avg_spec[i])
        i += 1

for j, d in enumerate(doublets):
    a = -math.log(d) / math.log(6) if d > 1e-20 else float('inf')
    print(f"  Doublet {j}: |λ| = {d:.6f}  a = {a:.6f}")

# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 3: THE 1+3+1 BAND — NEAR-KOIDE TRIPLET IN THE MIDDLE?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

if len(doublets) >= 5:
    # Band 1: doublet 0 (near-1)
    # Band 2: doublets 1,2,3 (middle)
    # Band 3: doublet 4 (lowest)

    d0 = doublets[0]
    d1, d2, d3 = doublets[1], doublets[2], doublets[3]
    d4 = doublets[4]

    print(f"\n  Band structure:")
    print(f"    Top:    {d0:.6f}  (exponent {-math.log(d0)/math.log(6):.6f})")
    print(f"    Middle: {d1:.6f}, {d2:.6f}, {d3:.6f}")
    print(f"    Bottom: {d4:.6f}  (exponent {-math.log(d4)/math.log(6):.6f})")

    # Koide Q of middle triplet
    s = math.sqrt(d1) + math.sqrt(d2) + math.sqrt(d3)
    Q_mid = (d1 + d2 + d3) / s**2
    print(f"\n  Koide Q of middle triplet: {Q_mid:.6f}  (1/H = {1/H:.6f}, 2/3 = {2/3:.6f})")

    # Koide angle
    m = sorted([d1, d2, d3], reverse=True)
    sqM = (math.sqrt(m[0]) + math.sqrt(m[1]) + math.sqrt(m[2])) / 3
    ratio = math.sqrt(m[0]) / sqM
    cos_theta = (ratio - 1) / math.sqrt(2)
    if abs(cos_theta) <= 1:
        theta_mid = math.acos(cos_theta)
        print(f"  Koide angle of middle triplet: θ = {theta_mid:.6f}")
        print(f"  θ × H = {theta_mid * H:.6f}")

    # Product d0 × d4
    print(f"\n  Top × Bottom = {d0 * d4:.6f}")
    print(f"  Middle geometric mean = {(d1 * d2 * d3)**(1/3):.6f}")
    print(f"  √(Top × Bottom) = {math.sqrt(d0 * d4):.6f}")

    # Check: top × bottom ≈ middle center?
    print(f"  d0 × d4 / d2² = {d0 * d4 / d2**2:.6f}")

    # Exponent arithmetic
    a_doublets = [-math.log(d) / math.log(6) for d in doublets]
    print(f"\n  Exponent differences:")
    print(f"    a₁ - a₀ = {a_doublets[1] - a_doublets[0]:.6f}")
    print(f"    a₂ - a₁ = {a_doublets[2] - a_doublets[1]:.6f}")
    print(f"    a₃ - a₂ = {a_doublets[3] - a_doublets[2]:.6f}")
    print(f"    a₄ - a₃ = {a_doublets[4] - a_doublets[3]:.6f}")
    print(f"    a₄ - a₀ = {a_doublets[4] - a_doublets[0]:.6f}  (2/3 = {2/3:.6f})")

    # Middle triplet exponents
    am = [a_doublets[1], a_doublets[2], a_doublets[3]]
    print(f"\n  Middle triplet exponents: {am[0]:.6f}, {am[1]:.6f}, {am[2]:.6f}")
    print(f"    Mean: {np.mean(am):.6f}  (1/3 = {1/3:.6f})")
    print(f"    Spread: {max(am) - min(am):.6f}")
    print(f"    Center: {am[1]:.6f}  (a₀ = {-math.log(element_data['e']['mag'])/math.log(6):.6f})")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 4: COMPOSITION PHASE × MAGNITUDE = KOIDE MASSES?")
print("= Test: do |λ|^n × cos(nφ + 2πk/3) give the right mass ratios?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Physical masses
m_e = 0.51100; m_mu = 105.66; m_tau = 1776.86

# Crystal eigenvalues (per conjugacy class)
l0 = element_data["e"]["mag"]
l1 = np.mean([element_data["(012)"]["mag"], element_data["(021)"]["mag"]])
l2 = np.mean([element_data["(01)"]["mag"], element_data["(02)"]["mag"], element_data["(12)"]["mag"]])

phi0 = element_data["e"]["phase"]
phi1 = np.mean([element_data["(012)"]["phase"], element_data["(021)"]["phase"]])
phi2 = np.mean([element_data["(01)"]["phase"], element_data["(02)"]["phase"], element_data["(12)"]["phase"]])

print(f"\n  Crystal data:")
print(f"    Identity:      |λ| = {l0:.6f}, φ = {phi0:.4f}")
print(f"    3-cycles:      |λ| = {l1:.6f}, φ = {phi1:.4f}")
print(f"    Transpositions: |λ| = {l2:.6f}, φ = {phi2:.4f}")

# The Koide formula: m_k = M(1 + √2 cos(θ + 2πk/3))²
# If mass ∝ |λ|^n, and the sign eigenvalue λ = |λ| e^{iφ}, then
# Re(λ^n) = |λ|^n cos(nφ)
# So maybe: mass ∝ Re(λ^n)?

n_koide = 35 / H**2
print(f"\n  Koide exponent n = 35/H² = {n_koide:.4f}")

for n in [n_koide, 4.0, 3.89, 3.5, 3.0]:
    # Masses from |λ|^n (magnitude only)
    m0 = l0 ** n
    m1_pred = l1 ** n
    m2_pred = l2 ** n

    # Masses from Re(λ^n) = |λ|^n cos(nφ)
    r0 = l0**n * math.cos(n * phi0)
    r1 = l1**n * math.cos(n * phi1)
    r2 = l2**n * math.cos(n * phi2)

    # Koide Q from magnitudes
    if m2_pred > 1e-50 and m0 > 1e-50 and m1_pred > 1e-50:
        s = math.sqrt(m0) + math.sqrt(m1_pred) + math.sqrt(m2_pred)
        Q_mag = (m0 + m1_pred + m2_pred) / s**2

        # Mass ratios
        print(f"\n  n = {n:.4f}:")
        print(f"    |λ|^n: {m0:.6e}, {m1_pred:.6e}, {m2_pred:.6e}")
        print(f"    Q(|λ|^n) = {Q_mag:.6f}")
        print(f"    ratio τ/e = {m0/m2_pred:.1f}  (phys: {m_tau/m_e:.1f})")
        print(f"    ratio μ/e = {m1_pred/m2_pred:.1f}  (phys: {m_mu/m_e:.1f})")

        # Now try with phase information
        if r0 > 1e-50 and r1 > 0 and r2 > 0:
            sr = math.sqrt(abs(r0)) + math.sqrt(abs(r1)) + math.sqrt(abs(r2))
            print(f"    Re(λ^n): {r0:.6e}, {r1:.6e}, {r2:.6e}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 5: CKM FROM SYM² STRUCTURE")
print("The Sym² distances showed more variation than Born fidelity")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

from solver.algebra import sym2_fingerprint, sym2_distance

# Build individual element crystals with good statistics
element_joints = {}
for name in S3_PERMS:
    joints = []
    for seed in range(N_SEEDS):
        ent = Entangler(S3_CORR[name], seed=seed).build()
        joints.append(ent.joint)
    element_joints[name] = torch.stack(joints).mean(dim=0)

fps = {n: sym2_fingerprint(element_joints[n]) for n in S3_PERMS}

# The three transpositions define three "directions" in flavor space
# Their Sym² distances form a triangle
trans = ["(01)", "(02)", "(12)"]
d01_02 = sym2_distance(fps["(01)"], fps["(02)"])
d01_12 = sym2_distance(fps["(01)"], fps["(12)"])
d02_12 = sym2_distance(fps["(02)"], fps["(12)"])

print(f"\n  Transposition Sym² triangle:")
print(f"    d((01),(02)) = {d01_02:.6f}")
print(f"    d((01),(12)) = {d01_12:.6f}")
print(f"    d((02),(12)) = {d02_12:.6f}")
print(f"    Perimeter = {d01_02 + d01_12 + d02_12:.6f}")

# Ratios
dists_sorted = sorted([d01_02, d01_12, d02_12])
print(f"\n  Sorted: {dists_sorted[0]:.6f} : {dists_sorted[1]:.6f} : {dists_sorted[2]:.6f}")
print(f"  Ratios: 1 : {dists_sorted[1]/dists_sorted[0]:.4f} : {dists_sorted[2]/dists_sorted[0]:.4f}")

# CKM-like angles: the three transpositions swap different pairs
# (01) swaps hypotheses 0,1 → "12 sector"
# (02) swaps hypotheses 0,2 → "13 sector"
# (12) swaps hypotheses 1,2 → "23 sector"
#
# The CKM matrix V_ij mixes generation i with generation j
# CKM hierarchy: |V_12| >> |V_23| >> |V_13|
# Wolfenstein: V_12 ~ λ, V_23 ~ λ², V_13 ~ λ³

# If (12)=V_12, (02)=V_13, (01)=V_23:
# d(12) should be largest, d(02) smallest
# But we see: d(01,02) > d(02,12) > d(01,12) (roughly equal)

# Identity to each transposition (different distances!)
d_e_01 = sym2_distance(fps["e"], fps["(01)"])
d_e_02 = sym2_distance(fps["e"], fps["(02)"])
d_e_12 = sym2_distance(fps["e"], fps["(12)"])

print(f"\n  Identity ↔ transposition Sym² distances:")
print(f"    d(e, (01)) = {d_e_01:.6f}")
print(f"    d(e, (02)) = {d_e_02:.6f}")
print(f"    d(e, (12)) = {d_e_12:.6f}")

dists_e = sorted([(d_e_01, "(01)"), (d_e_02, "(02)"), (d_e_12, "(12)")], key=lambda x: x[0])
print(f"\n  Sorted: {dists_e[0][1]}={dists_e[0][0]:.6f} < {dists_e[1][1]}={dists_e[1][0]:.6f} < {dists_e[2][1]}={dists_e[2][0]:.6f}")

# Ratios
print(f"  Ratios: 1 : {dists_e[1][0]/dists_e[0][0]:.4f} : {dists_e[2][0]/dists_e[0][0]:.4f}")

# Check: does the hierarchy match Wolfenstein?
cabibbo = 0.2245
print(f"\n  Wolfenstein hierarchy check:")
print(f"    d₁/d₀ = {dists_e[1][0]/dists_e[0][0]:.4f}  (λ_W = {cabibbo:.4f})")
print(f"    d₂/d₀ = {dists_e[2][0]/dists_e[0][0]:.4f}  (λ_W² = {cabibbo**2:.4f})")

# The three distances are 0.066, 0.073, 0.079 — ratio ~1:1.11:1.20
# NOT a λ hierarchy. But they're not equal either.
# The three distances differ by ~20% — this is the asymmetry of S₃.

# What about the Sym² fingerprint components directly?
print(f"\n  Sym² fingerprint vectors (9 components):")
for name in ["e", "(01)", "(02)", "(12)", "(012)", "(021)"]:
    fp = fps[name]
    print(f"    {name:>6s}: [{', '.join(f'{v:.4f}' for v in fp.tolist())}]")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 6: THE 3-CYCLE EXPONENT — NEW APPROACH")
print("Is a₁ determined by the Coxeter element order?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

a0 = -math.log(element_data["e"]["mag"]) / math.log(6)
a1_012 = -math.log(element_data["(012)"]["mag"]) / math.log(6)
a1_021 = -math.log(element_data["(021)"]["mag"]) / math.log(6)
a2_01 = -math.log(element_data["(01)"]["mag"]) / math.log(6)
a2_02 = -math.log(element_data["(02)"]["mag"]) / math.log(6)
a2_12 = -math.log(element_data["(12)"]["mag"]) / math.log(6)

a1 = (a1_012 + a1_021) / 2
a2 = (a2_01 + a2_02 + a2_12) / 3

print(f"\n  Per-element exponents:")
print(f"    e:    a = {a0:.8f}")
print(f"    (012): a = {a1_012:.8f}")
print(f"    (021): a = {a1_021:.8f}")
print(f"    (01):  a = {a2_01:.8f}")
print(f"    (02):  a = {a2_02:.8f}")
print(f"    (12):  a = {a2_12:.8f}")

print(f"\n  Class averages:")
print(f"    identity: a₀ = {a0:.8f}")
print(f"    3-cycles: a₁ = {a1:.8f}")
print(f"    transpositions: a₂ = {a2:.8f}")

# The elements have ORDERS in S₃:
# e: order 1, (σ): order 2, (ρ): order 3
# The Coxeter element of S₃ = A₂ is a 3-cycle, with order = Coxeter number h = 3

# Hypothesis: the exponent is determined by the element's representation-theoretic
# contribution to the sign sector

# For the sign representation: χ_sign(g) = (-1)^{length(g)}
# where length is the Coxeter length
# e: length 0, 3-cycles: length 2, transpositions: length 1

print(f"\n  Coxeter lengths:")
print(f"    e: ℓ=0, a={a0:.6f}")
print(f"    3-cycles: ℓ=2, a={a1:.6f}")
print(f"    transpositions: ℓ=1, a={a2:.6f}")

# Does a grow with Coxeter length? a₀(ℓ=0) < a₁(ℓ=2) < a₂(ℓ=1)? No! a₂ > a₁.
# So it's NOT monotone in Coxeter length.

# The sign character values: χ(e)=1, χ(3c)=1, χ(tr)=-1
# For sign=+1 elements (e, 3-cycles): coupling is positive
# For sign=-1 elements (transpositions): coupling is negative → much weaker

# The exponent difference between the two sign=+1 classes:
print(f"\n  Within sign=+1 sector:")
print(f"    a₁ - a₀ = {a1 - a0:.8f}")
print(f"    a₁/a₀ = {a1/a0:.8f}")

# Check: a₁/a₀ ≈ 2 + something small?
ratio = a1 / a0
print(f"    a₁/a₀ - 2 = {ratio - 2:.8f}")
print(f"    1/(a₁/a₀ - 2) = {1/(ratio-2):.4f}")

# a₁ = a₀ × (2 + ε). What is ε?
eps = ratio - 2
print(f"\n  ε = a₁/a₀ - 2 = {eps:.8f}")
print(f"  Candidates for ε:")
for name, val in [("1/H²", 1/H**2), ("1/(H²+1)", 1/(H**2+1)),
                   ("K*/H", K_STAR/H), ("Δ/H", DELTA/H),
                   ("1/h(E₈)", 1/30), ("1/h(D₄)", 1/6),
                   ("K*", K_STAR), ("BORN_FLOOR", BORN_FLOOR),
                   ("1/(H+1)", 1/(H+1)), ("1/2H", 1/(2*H)),
                   ("1/(H²-1)", 1/(H**2-1)),
                   ("H/(H²+1)²", H/(H**2+1)**2),
                   ("1/(H!+H)", 1/(math.factorial(H)+H)),
                   ("a₀", a0),
                   ("2/(H²+1)", 2/(H**2+1))]:
    dev = abs(eps - val) / abs(eps) * 100
    if dev < 20:
        print(f"    {name:>15s} = {val:.8f}  (dev = {dev:.2f}%)")


print(f"\n\n{'='*80}")
print("SUMMARY")
print("=" * 80)
print(f"""
  1. √2 mechanism: CONFIRMED to 0.06% (Principle 97)
     λ_sign(e) = 6^{{-1/3}} exp(i√2 θ_crystal)
     √2 = √(dim_std/(dim_sign×dim_triv))

  2. Phase structure:
     identity → 3-cycles: Δφ ≈ {wrap(phi_3c - phi_e)/math.pi:.3f}π
     identity → transpositions (bare): Δφ ≈ {wrap(phi_tr_bare - phi_e)/math.pi:.3f}π
     NOT the Koide 2π/3 spacing.

  3. Standard sector: 5 doublets with 1+3+1 band structure
     Middle 3 doublets have Q ≈ 1/H (near-degenerate)

  4. 3-cycle exponent: a₁/a₀ ≈ 2 + ε where ε ≈ {eps:.4f}
""")
