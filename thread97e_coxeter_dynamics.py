"""
Thread 97e: The Coxeter correction as S₃ averaging.

θ_phys = θ_crystal − 1/h(E₈). Can the correction arise from averaging
the composition operator over all S₃ elements?

If T_avg = (1/6)Σ T_g, does the sign eigenvalue of T_avg have phase
√2 × θ_phys instead of √2 × θ_crystal?

Also: does the GROUP ALGEBRA element (sum over S₃ with sign characters)
produce the corrected angle?
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

CLASS_SIZES = {"e": 1, "(01)": 3, "(02)": 3, "(12)": 3, "(012)": 2, "(021)": 2}
# Careful: (01), (02), (12) are all in the transposition class (size 3)
# (012), (021) are in the 3-cycle class (size 2)
# e is alone (size 1)

def perm_mat_4x4(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return P

reps = {n: torch.kron(perm_mat_4x4(p), perm_mat_4x4(p)) for n, p in S3_PERMS.items()}

chi_s = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}

def proj_sign():
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for n in S3_PERMS:
        P += chi_s[n] * reps[n]
    return P / 6.0

P_sign = proj_sign()

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

theta_crystal = (1 - K_STAR) / H  # 23/90
theta_physical = 2 / H**2  # 2/9 = 20/90
sqrt2_theta_c = math.sqrt(2) * theta_crystal
sqrt2_theta_p = math.sqrt(2) * theta_physical

N_SEEDS = 100

# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("PART 1: S₃-AVERAGED COMPOSITION OPERATOR")
print("T_avg = (1/|S₃|) Σ_g T_g")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Build element composition operators averaged over seeds
element_T = {}
for name in S3_PERMS:
    T_avg = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for seed in range(N_SEEDS):
        ent = Entangler(S3_CORR[name], seed=seed).build()
        T_avg += build_T_B(ent.joint)
    T_avg /= N_SEEDS
    element_T[name] = T_avg

# S₃ uniform average: each element weighted by 1/6
T_S3_avg = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
for name in S3_PERMS:
    T_S3_avg += element_T[name]
T_S3_avg /= 6

lam_avg = sign_coupling(T_S3_avg)
print(f"\n  T_S3_avg sign eigenvalue:")
print(f"    |λ| = {lam_avg.abs().item():.8f}")
print(f"    phase = {math.atan2(lam_avg.imag.item(), lam_avg.real.item()):.8f} rad")
print(f"    √2 × θ_crystal = {sqrt2_theta_c:.8f}")
print(f"    √2 × θ_physical = {sqrt2_theta_p:.8f}")

phase_avg = math.atan2(lam_avg.imag.item(), lam_avg.real.item())
print(f"\n    Phase / √2 = {phase_avg / math.sqrt(2):.8f}")
print(f"    θ_crystal = {theta_crystal:.8f}")
print(f"    θ_physical = {theta_physical:.8f}")
print(f"    Phase / (√2 × θ_c) = {phase_avg / sqrt2_theta_c:.6f}")
print(f"    Phase / (√2 × θ_p) = {phase_avg / sqrt2_theta_p:.6f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 2: CLASS-WEIGHTED AVERAGES")
print("Weight by conjugacy class size: 1, 3, 2 → normalized")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Average composition operators per class
class_T = {
    "identity": element_T["e"],
    "transpositions": (element_T["(01)"] + element_T["(02)"] + element_T["(12)"]) / 3,
    "3-cycles": (element_T["(012)"] + element_T["(021)"]) / 2,
}

# Various weighted averages
print("\n  Various weighted averages of class operators:")

weightings = {
    "uniform (1,1,1)/3": (1/3, 1/3, 1/3),
    "class size (1,3,2)/6": (1/6, 3/6, 2/6),
    "sign char (1,-3,2)/6": (1/6, -3/6, 2/6),  # sign-weighted
    "sign=+1 only (1,0,2)/3": (1/3, 0, 2/3),
    "sign=+1 equal (1,0,1)/2": (1/2, 0, 1/2),
    "(1,0,0)": (1, 0, 0),  # identity only
    "dim-weighted (1,2,3)/6": (1/6, 2/6, 3/6),  # weight by rep dim
    "(1,-1,1)/3": (1/3, -1/3, 1/3),
    "Plancherel (1,9,4)/14": (1/14, 9/14, 4/14),  # dim² weighting
}

for name, (w_id, w_tr, w_3c) in weightings.items():
    T_w = w_id * class_T["identity"] + w_tr * class_T["transpositions"] + w_3c * class_T["3-cycles"]
    lam = sign_coupling(T_w)
    mag = lam.abs().item()
    phase = math.atan2(lam.imag.item(), lam.real.item())
    phase_over_sqrt2 = phase / math.sqrt(2) if abs(phase) > 1e-10 else 0

    match_c = abs(1 - phase / sqrt2_theta_c) * 100 if abs(sqrt2_theta_c) > 1e-10 else float('inf')
    match_p = abs(1 - phase / sqrt2_theta_p) * 100 if abs(sqrt2_theta_p) > 1e-10 else float('inf')

    marker = ""
    if match_p < 1: marker = " ← θ_phys!"
    elif match_c < 1: marker = " ← θ_crystal"

    print(f"    {name:>30s}: |λ|={mag:.6f}, φ={phase:>8.4f}, φ/√2={phase_over_sqrt2:.6f} "
          f"(θ_c={match_c:.2f}%, θ_p={match_p:.2f}%){marker}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 3: SIGN-PROJECTED COMPOSITION")
print("T_sign = Σ χ_sign(g) T_g  (the group algebra projection)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# The group algebra element for the sign representation:
# e_sign = (1/6) Σ_g χ_sign(g) g
# Applied to composition: T_sign = (1/6)(T_e - T_(01) - T_(02) - T_(12) + T_(012) + T_(021))

T_sign_proj = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
for name, chi in chi_s.items():
    T_sign_proj += chi * element_T[name]
T_sign_proj /= 6

lam_sign_proj = sign_coupling(T_sign_proj)
mag = lam_sign_proj.abs().item()
phase = math.atan2(lam_sign_proj.imag.item(), lam_sign_proj.real.item())

print(f"\n  Sign-projected composition operator:")
print(f"    |λ| = {mag:.8f}")
print(f"    phase = {phase:.8f} rad")
print(f"    phase / √2 = {phase/math.sqrt(2):.8f}")
print(f"    θ_crystal = {theta_crystal:.8f}")
print(f"    θ_physical = {theta_physical:.8f}")
print(f"    Match to θ_c: {abs(1 - phase/(sqrt2_theta_c))*100:.4f}%")
print(f"    Match to θ_p: {abs(1 - phase/(sqrt2_theta_p))*100:.4f}%")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 4: ITERATED S₃-AVERAGED COMPOSITION")
print("What happens when we compose with the S₃ average N times?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# T_avg^N
T_N = torch.eye(dim_sq, dtype=torch.cfloat)
T_uniform = T_S3_avg.clone()

print(f"\n  Iterated S₃-averaged composition:")
print(f"  {'N':>4s}  {'|λ_sign|':>12s}  {'phase':>10s}  {'phase/√2':>10s}  {'note':>20s}")
print(f"  {'-'*65}")

for N in range(1, 16):
    T_N = T_N @ T_uniform
    lam = sign_coupling(T_N)
    mag = lam.abs().item()
    phase = math.atan2(lam.imag.item(), lam.real.item())
    phase_norm = phase / math.sqrt(2)

    note = ""
    # Per-step phase increment
    if N == 1:
        phase_per_step = phase
    else:
        # phase should grow as N × δ
        phase_per_step = phase / N

    note = f"φ/N/√2={phase_per_step/math.sqrt(2):.6f}"

    if N <= 8 or N % 5 == 0:
        print(f"  {N:>4d}  {mag:>12.6e}  {phase:>10.6f}  {phase_norm:>10.6f}  {note}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 5: THE SIGN-ONLY AVERAGE (OPPOSITE SIGN REMOVED)")
print("T_sign_only = (T_e + T_(012) + T_(021)) / 3 (sign=+1 elements)")
print("= average over elements where χ_sign = +1")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

T_plus = (element_T["e"] + element_T["(012)"] + element_T["(021)"]) / 3

lam_plus = sign_coupling(T_plus)
mag = lam_plus.abs().item()
phase = math.atan2(lam_plus.imag.item(), lam_plus.real.item())

print(f"\n  Sign=+1 average composition operator:")
print(f"    |λ| = {mag:.8f}")
print(f"    phase = {phase:.8f} rad")
print(f"    phase / √2 = {phase/math.sqrt(2):.8f}")
print(f"    θ_crystal = {theta_crystal:.8f}")
print(f"    θ_physical = {theta_physical:.8f}")
print(f"    Match to θ_c: {abs(1 - phase/sqrt2_theta_c)*100:.4f}%")
print(f"    Match to θ_p: {abs(1 - phase/sqrt2_theta_p)*100:.4f}%")

# The sign=+1 average is interesting because it combines only the elements
# that preserve the sign sector. Does it produce the corrected angle?

# Also try: sign-weighted sum
T_sign_weighted = (1 * element_T["e"] + 1 * element_T["(012)"] + 1 * element_T["(021)"]
                   - 1 * element_T["(01)"] - 1 * element_T["(02)"] - 1 * element_T["(12)"]) / 6

lam_sw = sign_coupling(T_sign_weighted)
mag_sw = lam_sw.abs().item()
phase_sw = math.atan2(lam_sw.imag.item(), lam_sw.real.item())

print(f"\n  Sign-character-weighted composition (= P_sign applied to T):")
print(f"    |λ| = {mag_sw:.8f}")
print(f"    phase = {phase_sw:.8f} rad")
print(f"    phase / √2 = {phase_sw/math.sqrt(2):.8f}")
print(f"    Match to θ_c: {abs(1 - phase_sw/sqrt2_theta_c)*100:.4f}%")
print(f"    Match to θ_p: {abs(1 - phase_sw/sqrt2_theta_p)*100:.4f}%")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 6: THE KEY CALCULATION")
print("T_e − (1/h(E₈)) × T_Coxeter?")
print("The Coxeter element of A₂ is a 3-cycle (012)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# The Coxeter correction is 1/h(E₈) = 1/30.
# The Coxeter element of A₂ = S₃ is a product of simple reflections = (012).
# What if: T_corrected = T_e - (1/h(E₈)) × T_Coxeter?

h_E8 = 30
T_coxeter = element_T["(012)"]

for alpha_name, alpha in [("1/h(E₈)", 1/h_E8),
                           ("1/h(D₄)", 1/6),
                           ("K*", K_STAR),
                           ("1/(H²+1)", 1/(H**2+1)),
                           ("H/h(E₈)", H/h_E8),
                           ("(a₁-a₀)/h(E₈)", (0.717-0.334)/h_E8),
                           ("1/(H×h(D₄))", 1/(H*6)),
                           ("1/(2H)", 1/(2*H))]:
    T_corr = element_T["e"] - alpha * T_coxeter
    lam = sign_coupling(T_corr)
    mag = lam.abs().item()
    phase = math.atan2(lam.imag.item(), lam.real.item())

    match_c = abs(1 - phase/sqrt2_theta_c)*100
    match_p = abs(1 - phase/sqrt2_theta_p)*100
    marker = ""
    if match_p < 2: marker = " ← θ_phys"
    elif match_c < 2: marker = " ← θ_crystal"

    print(f"    T_e − {alpha_name:>15s}·T_Cox: |λ|={mag:.6f}, φ={phase:.6f}, "
          f"φ/√2={phase/math.sqrt(2):.6f} (θ_c:{match_c:.2f}%, θ_p:{match_p:.2f}%){marker}")

# Also try: T_e - α × T_avg(3-cycles)
print(f"\n  Using average 3-cycle operator:")
T_3c_avg = (element_T["(012)"] + element_T["(021)"]) / 2

for alpha_name, alpha in [("1/h(E₈)", 1/h_E8),
                           ("H/h(E₈)", H/h_E8),
                           ("1/h(D₄)", 1/6),
                           ("1/(H+1)", 1/(H+1))]:
    T_corr = element_T["e"] - alpha * T_3c_avg
    lam = sign_coupling(T_corr)
    mag = lam.abs().item()
    phase = math.atan2(lam.imag.item(), lam.real.item())

    match_c = abs(1 - phase/sqrt2_theta_c)*100
    match_p = abs(1 - phase/sqrt2_theta_p)*100
    marker = ""
    if match_p < 2: marker = " ← θ_phys"

    print(f"    T_e − {alpha_name:>15s}·T_3c: |λ|={mag:.6f}, φ={phase:.6f}, "
          f"φ/√2={phase/math.sqrt(2):.6f} (θ_c:{match_c:.2f}%, θ_p:{match_p:.2f}%){marker}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 7: COMPOSITION OF IDENTITY WITH 3-CYCLE")
print("T_e × T_3c (sequential composition, not sum)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

T_composed = element_T["e"] @ element_T["(012)"]
lam = sign_coupling(T_composed)
mag = lam.abs().item()
phase = math.atan2(lam.imag.item(), lam.real.item())

print(f"\n  T_e × T_(012):")
print(f"    |λ| = {mag:.8f}")
print(f"    phase = {phase:.8f}")
print(f"    phase / √2 = {phase/math.sqrt(2):.8f}")
print(f"    2×phase / √2 = {2*phase/math.sqrt(2):.8f}  (sum of identity and 3-cycle phases)")

# Check: does T_e × T_3c have phase = identity_phase + 3cycle_phase?
phi_e = math.atan2(sign_coupling(element_T["e"]).imag.item(),
                    sign_coupling(element_T["e"]).real.item())
phi_3c = math.atan2(sign_coupling(element_T["(012)"]).imag.item(),
                     sign_coupling(element_T["(012)"]).real.item())

print(f"\n  Phase additivity check:")
print(f"    φ(e) = {phi_e:.8f}")
print(f"    φ(3c) = {phi_3c:.8f}")
print(f"    φ(e) + φ(3c) = {phi_e + phi_3c:.8f}")
print(f"    φ(e×3c) = {phase:.8f}")
print(f"    Match: {abs(1 - phase/(phi_e + phi_3c))*100:.4f}%")

# Average of e and 3c phases
avg_phase = (phi_e + phi_3c) / 2
print(f"\n  Average phase: (φ(e)+φ(3c))/2 = {avg_phase:.8f}")
print(f"  avg / √2 = {avg_phase/math.sqrt(2):.8f}")
print(f"  θ_crystal = {theta_crystal:.8f}, θ_physical = {theta_physical:.8f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("SUMMARY")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print(f"""
  Key: θ_crystal = 23/90, θ_physical = 2/9 = 20/90.
  √2 × θ_crystal = {sqrt2_theta_c:.6f}
  √2 × θ_physical = {sqrt2_theta_p:.6f}

  The Coxeter correction θ_phys = θ_crystal − 1/h(E₈) was derived
  algebraically (Principle 96). This session tested whether the
  crystal dynamics can produce it. Results above.
""")
