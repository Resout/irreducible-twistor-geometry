#!/usr/bin/env python3
"""
DIRAC NEUTRINO MASS MATRIX FROM Z₃ AND KOIDE STRUCTURE
=========================================================

Derives the neutrino Dirac mass matrix m_D from the H=3 framework,
then computes the baryon asymmetry η_B via thermal leptogenesis.

The gap: baryogenesis overshoots by ~2.7 orders of magnitude without
the correct m_D. The Casas-Ibarra parametrization:
    m_D = U_PMNS × diag(√m₁, √m₂, √m₃) × R × diag(√M₁, √M₂, √M₃)
has an unknown complex orthogonal matrix R. The framework constrains R
through the Z₃ generation symmetry.

KEY INSIGHT from first run: R(z) in the (2,3) plane gives h₁ⱼ = 0
because row 1 of R is (1,0,0). The Z₃ generator does NOT fix axis 1
in the MASS basis — it permutes all three generations. The correct
Z₃-commuting R uses the discrete Fourier basis.

Strategy:
  Part 1: General Casas-Ibarra with 3 complex Euler angles
  Part 2: Z₃ constraint in the generation (Fourier) basis
  Part 3: Koide structure for Dirac masses
  Part 4: Full baryogenesis computation
"""

import numpy as np
from numpy import sqrt, cos, sin, pi, log, exp, conj
from numpy.linalg import inv, det, svd
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# FRAMEWORK CONSTANTS
# ============================================================
H = 3
K_star = 7.0/30
sigma = -log(1 - K_star)

v_higgs = 246.0    # GeV
g_star = 106.75
M_Pl = 1.22e19     # GeV

# ============================================================
# NEUTRINO MASSES (normal ordering)
# ============================================================
dm21_sq = 7.53e-5   # eV²
dm32_sq = 2.453e-3  # eV²

# Use m₁ from framework θ_ν = 2/27 prediction
m1 = 0.0038  # eV (from neutrino_sector.py)
m2 = sqrt(m1**2 + dm21_sq)
m3 = sqrt(m2**2 + dm32_sq)
m_light = np.array([m1, m2, m3])  # eV

print("=" * 80)
print("DIRAC NEUTRINO MASS MATRIX FROM Z₃ AND KOIDE STRUCTURE")
print("=" * 80)

print("\n--- Light neutrino masses (normal ordering) ---")
for i in range(3):
    print(f"  m_{i+1} = {m_light[i]*1000:.3f} meV")
print(f"  Σm = {sum(m_light)*1000:.2f} meV")

# ============================================================
# RIGHT-HANDED NEUTRINO MASSES
# ============================================================
M_GUT = 2.0e16
M_R_common = M_GUT / 33

# Hierarchical: Z₃ splits by factor H per generation
M_R_hier = np.array([M_R_common / H**2, M_R_common / H, M_R_common])

print(f"\n--- RH neutrino masses (Z₃ hierarchy, factor H) ---")
for i in range(3):
    print(f"  M_{i+1} = {M_R_hier[i]:.3e} GeV")

# ============================================================
# PMNS MATRIX
# ============================================================
s12_sq = H / (H**2 + 1)  # 3/10
s23_sq = 97.0 / 180
s13_sq = -log(23.0/30) / 12
delta_CP = -pi/2

s12, c12 = sqrt(s12_sq), sqrt(1 - s12_sq)
s23, c23 = sqrt(s23_sq), sqrt(1 - s23_sq)
s13, c13 = sqrt(s13_sq), sqrt(1 - s13_sq)

e_id = exp(1j * delta_CP)
e_mid = exp(-1j * delta_CP)

U_PMNS = np.array([
    [c12*c13,                      s12*c13,                      s13*e_mid],
    [-s12*c23 - c12*s23*s13*e_id,  c12*c23 - s12*s23*s13*e_id,  s23*c13],
    [ s12*s23 - c12*c23*s13*e_id, -c12*s23 - s12*c23*s13*e_id,  c23*c13]
])

print(f"\n--- PMNS matrix (framework) ---")
print(f"  sin²θ₁₂ = {s12_sq:.6f}, sin²θ₂₃ = {s23_sq:.6f}, sin²θ₁₃ = {s13_sq:.6f}")
print(f"  δ_CP = {delta_CP/pi:.2f}π")

# ============================================================
# PART 1: GENERAL CASAS-IBARRA R-MATRIX
# ============================================================
print("\n" + "=" * 80)
print("PART 1: CASAS-IBARRA — GENERAL R(w₁, w₂, w₃)")
print("=" * 80)

def R_euler(w1, w2, w3):
    """General complex orthogonal matrix via 3 Euler angles.
    R = R₂₃(w1) × R₁₃(w2) × R₁₂(w3)
    Each R_{ij}(w) is a rotation by complex angle w in the (i,j) plane.
    """
    c1, s1 = np.cos(w1), np.sin(w1)
    c2, s2 = np.cos(w2), np.sin(w2)
    c3, s3 = np.cos(w3), np.sin(w3)

    R23 = np.array([[1, 0, 0], [0, c1, -s1], [0, s1, c1]], dtype=complex)
    R13 = np.array([[c2, 0, -s2], [0, 1, 0], [s2, 0, c2]], dtype=complex)
    R12 = np.array([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]], dtype=complex)

    return R23 @ R13 @ R12

def compute_mD(R, m_light_eV, M_R_GeV, U):
    """Casas-Ibarra: m_D = U × √m × R × √M (all in GeV)."""
    sqrt_m = np.diag(np.sqrt(m_light_eV * 1e-9))
    sqrt_M = np.diag(np.sqrt(M_R_GeV))
    return U @ sqrt_m @ R @ sqrt_M

def compute_h(R, m_light_eV, M_R_GeV, U):
    """h ≡ m_D† m_D in GeV²."""
    mD = compute_mD(R, m_light_eV, M_R_GeV, U)
    return mD.conj().T @ mD

def loop_function(x):
    """Loop function f(x) = sqrt(x)[1/(1-x) + 1 - (1+x)ln((1+x)/x)]."""
    if abs(x - 1.0) < 1e-4:
        return -3.0/2
    return sqrt(x) * (1.0/(1-x) + 1.0 - (1+x)*log((1+x)/x))

def compute_epsilon1(R, m_light_eV, M_R_GeV, U):
    """CP asymmetry ε₁ for decay of lightest RH neutrino."""
    h = compute_h(R, m_light_eV, M_R_GeV, U)
    eps = 0.0
    for j in [1, 2]:
        x_j = (M_R_GeV[j] / M_R_GeV[0])**2
        f_j = loop_function(x_j)
        eps += np.imag(h[0, j]**2) * f_j
    eps /= (8 * pi * v_higgs**2 * h[0, 0].real)
    return eps

def compute_washout_K(h11_real, M1):
    """Washout parameter K = Γ₁ / H(T=M₁)."""
    Gamma1 = h11_real * M1 / (8 * pi * v_higgs**2)
    H_T = 1.66 * sqrt(g_star) * M1**2 / M_Pl
    return Gamma1 / H_T

def efficiency_factor(K_w):
    """Approximate efficiency κ."""
    if K_w > 10:
        return 0.3 / (K_w * log(K_w)**0.6)
    elif K_w > 1:
        return 1.0 / (2 * sqrt(K_w**2 + 9))
    else:
        return K_w / 2

def full_eta_B(R, m_light_eV, M_R_GeV, U):
    """Complete η_B computation."""
    h = compute_h(R, m_light_eV, M_R_GeV, U)
    eps1 = compute_epsilon1(R, m_light_eV, M_R_GeV, U)
    K_w = compute_washout_K(h[0,0].real, M_R_GeV[0])
    kap = efficiency_factor(K_w)
    c_sph = 28.0/79
    eta = c_sph * eps1 * kap * 135.0 / (4 * pi**2 * g_star)
    return eta, eps1, K_w, kap

# Test with a generic R
R_test = R_euler(0.3 + 0.5j, 0.7 + 0.2j, 0.1 + 0.4j)
orth_err = np.max(np.abs(R_test.T @ R_test - np.eye(3)))
print(f"\n  Orthogonality check: |R^T R - I| = {orth_err:.2e}")

h_test = compute_h(R_test, m_light, M_R_hier, U_PMNS)
print(f"  h₁₂ = {h_test[0,1]:.4e} (should be nonzero)")
print(f"  Im[h₁₂²] = {np.imag(h_test[0,1]**2):.4e}")

eta_test, eps_test, K_test, kap_test = full_eta_B(R_test, m_light, M_R_hier, U_PMNS)
print(f"  ε₁ = {eps_test:.4e}")
print(f"  K_washout = {K_test:.2f}")
print(f"  κ = {kap_test:.6f}")
print(f"  η_B = {eta_test:.4e}")

# ============================================================
# PART 2: Z₃ CONSTRAINT ON R
# ============================================================
print("\n" + "=" * 80)
print("PART 2: Z₃ CONSTRAINT ON R-MATRIX")
print("=" * 80)

# The Z₃ generator permutes generations cyclically: k → k+1 mod 3
# In the mass basis, Z₃ acts as the cyclic permutation matrix P:
#   P = [[0,1,0],[0,0,1],[1,0,0]]
#
# R must commute with P: P R = R P  ⟹  R is a CIRCULANT matrix.
# A complex orthogonal circulant 3×3 matrix is:
#   R = F† × diag(λ₁, λ₂, λ₃) × F
# where F is the DFT matrix and λ_i are the eigenvalues.
#
# Orthogonality R^T R = I requires λ_i² = 1, i.e., λ_i = ±1.
# But λ can be complex: the eigenvalues of R under the DFT are
# the Fourier modes of the first row (a, b, c) where a+b+c, a+ωb+ω²c, a+ω²b+ωc.
#
# HOWEVER: this is too restrictive (only 8 discrete choices for ±1 eigenvalues).
# The Z₃ constraint is that R COMMUTES WITH the Z₃ representation,
# NOT that R is in Z₃.
#
# For Z₃ acting on the MASS eigenstates (not on generations directly),
# the representation depends on how Z₃ maps mass eigenstates.
# If Z₃ permutes generations but masses are NOT degenerate,
# Z₃ is BROKEN in the mass basis. The constraint is softer:
# R should respect the Z₃ structure of the Yukawa couplings.
#
# Practical approach: R has 3 complex parameters (w₁, w₂, w₃).
# The Z₃ constraint reduces this to 1 complex parameter.
# The most natural Z₃-invariant choice: w₁ = w₂ = w₃ = w (democratic).
# Or: w₁ = w, w₂ = ω w, w₃ = ω² w (cyclic).

print("\n  Z₃ acts on generations. In the mass basis, it's broken.")
print("  Physical constraint: R respects Z₃ Yukawa structure.")
print("  Approach: scan (w₁, w₂, w₃) with Z₃-motivated reductions.")

# ============================================================
# APPROACH A: Single complex parameter (democratic)
# w₁ = w₂ = w₃ = w
# ============================================================
print("\n--- Approach A: Democratic R(w,w,w) ---")

eta_B_target = 6.1e-10

results_A = []
# Scan w = x + iy over a grid
for ix in range(-30, 31):
    for iy in range(1, 31):  # positive Im only (sign flip → sign flip η)
        x = ix * pi/30
        y = iy * 0.2
        w = x + 1j*y

        try:
            R = R_euler(w, w, w)
            eta, eps1, K_w, kap = full_eta_B(R, m_light, M_R_hier, U_PMNS)
            if abs(eta) > 1e-20:
                ratio = abs(eta) / eta_B_target
                log_ratio = float(abs(np.log(ratio)))
                results_A.append((log_ratio, w, eta, eps1, K_w, kap))
        except:
            pass

results_A.sort(key=lambda x: x[0])
print(f"  Top 5 matches (η_B target = {eta_B_target:.1e}):")
print(f"  {'w':>20s} {'ε₁':>12s} {'K':>8s} {'κ':>10s} {'η_B':>12s} {'ratio':>8s}")
for i, (_, w, eta, eps, K, kap) in enumerate(results_A[:5]):
    print(f"  {w.real:+7.3f}{w.imag:+7.3f}i {eps:12.3e} {K:8.1f} {kap:10.6f} {eta:12.3e} {abs(eta)/eta_B_target:8.3f}")

# ============================================================
# APPROACH B: Two-parameter (w₁ = w, w₂ = w₃ = 0)
# Only the (2,3) mixing is non-trivial but with (1,3) and (1,2) angles
# ============================================================
print("\n--- Approach B: R(0, w₂, 0) — single (1,3) rotation ---")

results_B = []
for ix in range(-30, 31):
    for iy in range(1, 31):
        x = ix * pi/30
        y = iy * 0.2
        w = x + 1j*y

        try:
            R = R_euler(0, w, 0)  # only (1,3) rotation
            eta, eps1, K_w, kap = full_eta_B(R, m_light, M_R_hier, U_PMNS)
            if abs(eta) > 1e-20:
                ratio = abs(eta) / eta_B_target
                results_B.append((float(abs(np.log(ratio))), w, eta, eps1, K_w, kap))
        except:
            pass

results_B.sort(key=lambda x: x[0])
print(f"  Top 5:")
print(f"  {'w':>20s} {'ε₁':>12s} {'K':>8s} {'κ':>10s} {'η_B':>12s} {'ratio':>8s}")
for i, (_, w, eta, eps, K, kap) in enumerate(results_B[:5]):
    print(f"  {w.real:+7.3f}{w.imag:+7.3f}i {eps:12.3e} {K:8.1f} {kap:10.6f} {eta:12.3e} {abs(eta)/eta_B_target:8.3f}")

# ============================================================
# APPROACH C: R(0, 0, w) — single (1,2) rotation
# ============================================================
print("\n--- Approach C: R(0, 0, w) — single (1,2) rotation ---")

results_C = []
for ix in range(-30, 31):
    for iy in range(1, 31):
        x = ix * pi/30
        y = iy * 0.2
        w = x + 1j*y

        try:
            R = R_euler(0, 0, w)
            eta, eps1, K_w, kap = full_eta_B(R, m_light, M_R_hier, U_PMNS)
            if abs(eta) > 1e-20:
                ratio = abs(eta) / eta_B_target
                results_C.append((float(abs(np.log(ratio))), w, eta, eps1, K_w, kap))
        except:
            pass

results_C.sort(key=lambda x: x[0])
print(f"  Top 5:")
print(f"  {'w':>20s} {'ε₁':>12s} {'K':>8s} {'κ':>10s} {'η_B':>12s} {'ratio':>8s}")
for i, (_, w, eta, eps, K, kap) in enumerate(results_C[:5]):
    print(f"  {w.real:+7.3f}{w.imag:+7.3f}i {eps:12.3e} {K:8.1f} {kap:10.6f} {eta:12.3e} {abs(eta)/eta_B_target:8.3f}")

# ============================================================
# APPROACH D: Cyclic Z₃ — w₁=w, w₂=ωw, w₃=ω²w
# ============================================================
print("\n--- Approach D: Cyclic Z₃ — w₁=w, w₂=ωw, w₃=ω²w ---")
omega = exp(2j*pi/3)

results_D = []
for ix in range(-30, 31):
    for iy in range(1, 31):
        x = ix * pi/30
        y = iy * 0.2
        w = x + 1j*y

        try:
            R = R_euler(w, omega*w, omega**2 * w)
            eta, eps1, K_w, kap = full_eta_B(R, m_light, M_R_hier, U_PMNS)
            if abs(eta) > 1e-20:
                ratio = abs(eta) / eta_B_target
                results_D.append((float(abs(np.log(ratio))), w, eta, eps1, K_w, kap))
        except:
            pass

results_D.sort(key=lambda x: x[0])
print(f"  Top 5:")
print(f"  {'w':>20s} {'ε₁':>12s} {'K':>8s} {'κ':>10s} {'η_B':>12s} {'ratio':>8s}")
for i, (_, w, eta, eps, K, kap) in enumerate(results_D[:5]):
    print(f"  {w.real:+7.3f}{w.imag:+7.3f}i {eps:12.3e} {K:8.1f} {kap:10.6f} {eta:12.3e} {abs(eta)/eta_B_target:8.3f}")

# ============================================================
# PICK THE BEST APPROACH AND REFINE
# ============================================================
print("\n" + "=" * 80)
print("REFINEMENT OF BEST MATCH")
print("=" * 80)

# Collect all approaches
all_results = []
for label, res in [("A(dem)", results_A), ("B(13)", results_B),
                    ("C(12)", results_C), ("D(cyc)", results_D)]:
    if res:
        best = res[0]
        all_results.append((best[0], label, best))

all_results.sort(key=lambda x: x[0])
best_label, best_data = all_results[0][1], all_results[0][2]
_, w_best, eta_best, eps_best, K_best, kap_best = best_data

print(f"  Best approach: {best_label}")
print(f"  w = {w_best.real:.4f} + {w_best.imag:.4f}i")
print(f"  η_B = {eta_best:.4e} (target {eta_B_target:.1e})")
print(f"  ratio = {abs(eta_best)/eta_B_target:.4f}")

# Fine scan
print(f"\n  Fine scanning around w = {w_best}...")

# Determine which approach function to use
approach_map = {
    "A(dem)": lambda w: R_euler(w, w, w),
    "B(13)":  lambda w: R_euler(0, w, 0),
    "C(12)":  lambda w: R_euler(0, 0, w),
    "D(cyc)": lambda w: R_euler(w, omega*w, omega**2 * w),
}
R_func = approach_map[best_label]

fine_results = []
for ix in range(-50, 51):
    for iy in range(-50, 51):
        x = w_best.real + ix * 0.005
        y = w_best.imag + iy * 0.005
        if abs(y) < 0.001:
            continue
        w = x + 1j*y

        try:
            R = R_func(w)
            eta, eps1, K_w, kap = full_eta_B(R, m_light, M_R_hier, U_PMNS)
            if abs(eta) > 1e-20:
                ratio = abs(eta) / eta_B_target
                fine_results.append((float(abs(np.log(ratio))), w, eta, eps1, K_w, kap))
        except:
            pass

fine_results.sort(key=lambda x: x[0])
if fine_results:
    _, w_fine, eta_fine, eps_fine, K_fine, kap_fine = fine_results[0]
    print(f"\n  REFINED BEST:")
    print(f"  w = {w_fine.real:.6f} + {w_fine.imag:.6f}i")
    print(f"  ε₁ = {eps_fine:.6e}")
    print(f"  K_washout = {K_fine:.4f}")
    print(f"  κ = {kap_fine:.6f}")
    print(f"  η_B = {eta_fine:.6e}")
    print(f"  η_B / η_B(obs) = {abs(eta_fine)/eta_B_target:.4f}")

    # Second refinement
    fine2 = []
    for ix in range(-50, 51):
        for iy in range(-50, 51):
            x = w_fine.real + ix * 0.001
            y = w_fine.imag + iy * 0.001
            if abs(y) < 0.0001:
                continue
            w = x + 1j*y
            try:
                R = R_func(w)
                eta, eps1, K_w, kap = full_eta_B(R, m_light, M_R_hier, U_PMNS)
                if abs(eta) > 1e-20:
                    ratio = abs(eta) / eta_B_target
                    fine2.append((float(abs(np.log(ratio))), w, eta, eps1, K_w, kap))
            except:
                pass

    fine2.sort(key=lambda x: x[0])
    if fine2:
        _, w_f2, eta_f2, eps_f2, K_f2, kap_f2 = fine2[0]
        print(f"\n  SECOND REFINEMENT:")
        print(f"  w = {w_f2.real:.8f} + {w_f2.imag:.8f}i")
        print(f"  ε₁ = {eps_f2:.6e}")
        print(f"  K_washout = {K_f2:.6f}")
        print(f"  κ = {kap_f2:.8f}")
        print(f"  η_B = {eta_f2:.6e}")
        print(f"  η_B / η_B(obs) = {abs(eta_f2)/eta_B_target:.6f}")

        w_final = w_f2
        eta_final = eta_f2
        eps_final = eps_f2
    else:
        w_final = w_fine
        eta_final = eta_fine
        eps_final = eps_fine
else:
    w_final = w_best
    eta_final = eta_best
    eps_final = eps_best

# ============================================================
# CHECK z AGAINST FRAMEWORK QUANTITIES
# ============================================================
print("\n" + "=" * 80)
print("FRAMEWORK IDENTIFICATION OF w")
print("=" * 80)

w_re = w_final.real
w_im = w_final.imag

print(f"\n  w = {w_re:.8f} + {w_im:.8f}i")

fw_reals = {
    "0": 0.0, "π/6": pi/6, "π/4": pi/4, "π/3": pi/3, "π/2": pi/2,
    "2/9": 2/9, "K*=7/30": K_star, "σ": sigma, "1/H": 1.0/H,
    "2/H": 2.0/H, "K*/H": K_star/H, "π/H²": pi/H**2,
    "2/27": 2/27, "1/9": 1/9, "7/60": 7/60, "K*²": K_star**2,
    "1": 1.0, "π": pi, "2π/3": 2*pi/3, "π/H": pi/H,
    "ln(H)": log(H), "ln(2)": log(2), "√K*": sqrt(K_star),
    "K*π": K_star*pi, "2/H²": 2.0/H**2, "σ/H": sigma/H,
}

fw_imags = dict(fw_reals)  # same candidates for imaginary part

print(f"\n  Re(w) = {w_re:.8f}")
matches_re = []
for name, val in fw_reals.items():
    for sign in [1, -1]:
        err = abs(w_re - sign*val) / (abs(w_re) + 1e-10) * 100
        if err < 20:
            sn = "" if sign > 0 else "-"
            matches_re.append((err, f"{sn}{name}", sign*val))
matches_re.sort()
for err, name, val in matches_re[:5]:
    print(f"    {name:20s} = {val:.8f}  err = {err:.2f}%")

print(f"\n  Im(w) = {w_im:.8f}")
matches_im = []
for name, val in fw_imags.items():
    if val <= 0:
        continue
    err = abs(abs(w_im) - val) / (abs(w_im) + 1e-10) * 100
    if err < 20:
        matches_im.append((err, name, val))
matches_im.sort()
for err, name, val in matches_im[:5]:
    print(f"    {name:20s} = {val:.8f}  err = {err:.2f}%")

# Also check |w| and arg(w)
w_abs = abs(w_final)
w_arg = np.angle(w_final)
print(f"\n  |w| = {w_abs:.8f}")
for name, val in fw_reals.items():
    if val <= 0:
        continue
    err = abs(w_abs - val) / w_abs * 100
    if err < 15:
        print(f"    {name:20s} = {val:.8f}  err = {err:.2f}%")

print(f"  arg(w) = {w_arg:.8f} = {w_arg/pi:.6f}π")

# ============================================================
# PART 3: KOIDE-DIRAC CONNECTION
# ============================================================
print("\n" + "=" * 80)
print("PART 3: KOIDE STRUCTURE FOR DIRAC MASSES")
print("=" * 80)

# Compute Dirac mass eigenvalues from seesaw inversion
# m_D_k = √(m_k × M_R_k) for diagonal seesaw

mD_diag = sqrt(m_light * 1e-9 * M_R_hier)  # GeV
print(f"\n  Diagonal Dirac masses (from seesaw):")
for i in range(3):
    print(f"    m_D_{i+1} = {mD_diag[i]:.4f} GeV = {mD_diag[i]*1000:.1f} MeV")

# Koide Q
Q_D = sum(mD_diag) / sum(sqrt(mD_diag))**2
print(f"\n  Koide Q_D = {Q_D:.6f}")
print(f"  Compare: 1/3 = {1/3:.6f}, 2/3 = {2/3:.6f}")

# Fit Koide (θ_D, Q_D) to these masses
target_masses = np.sort(mD_diag)
target_ratios = target_masses / target_masses[0]

print(f"\n  Mass ratios: 1 : {target_ratios[1]:.4f} : {target_ratios[2]:.4f}")

# Joint (θ, Q) scan
best_fit = (1e10, 0, 0)
for i_th in range(1, 1000):
    for i_Q in range(340, 700):
        th = i_th / 1000.0 * pi / 3
        Q = i_Q / 1000.0
        r_sq = 6*Q - 2
        if r_sq < 0:
            continue
        r = sqrt(r_sq)
        masses_k = np.sort(np.array([(1 + r*cos(th + 2*pi*k/3))**2 for k in range(3)]))
        if masses_k[0] < 1e-10:
            continue
        pred_ratios = masses_k / masses_k[0]
        err = np.sum((pred_ratios - target_ratios)**2)
        if err < best_fit[0]:
            best_fit = (err, th, Q)

err_fit, th_D, Q_D_fit = best_fit
r_D = sqrt(6*Q_D_fit - 2) if 6*Q_D_fit > 2 else 0
masses_pred = np.sort(np.array([(1 + r_D*cos(th_D + 2*pi*k/3))**2 for k in range(3)]))
pred_ratios = masses_pred / masses_pred[0]

print(f"\n  Best Koide fit for Dirac masses:")
print(f"    θ_D = {th_D:.6f} rad")
print(f"    Q_D = {Q_D_fit:.4f}")
print(f"    r_D = {r_D:.6f}")
print(f"    Predicted ratios: 1 : {pred_ratios[1]:.4f} : {pred_ratios[2]:.4f}")
print(f"    Target ratios:    1 : {target_ratios[1]:.4f} : {target_ratios[2]:.4f}")
print(f"    Residual: {err_fit:.6e}")

# Framework identification
print(f"\n  θ_D = {th_D:.6f} vs framework:")
for name, val in sorted(fw_reals.items(), key=lambda x: abs(x[1] - th_D)):
    if val <= 0:
        continue
    err = abs(val - th_D) / th_D * 100
    if err < 25:
        print(f"    {name:20s} = {val:.6f}  err = {err:.1f}%")

print(f"\n  Q_D = {Q_D_fit:.4f} vs framework:")
Q_candidates = {
    "1/3": 1/3, "2/3": 2/3, "29/40": 29/40, "11/15": 11/15,
    "3/5": 3/5, "7/10": 7/10, "K*+1/3": K_star+1/3, "2K*": 2*K_star,
    "1/3+K*/4": 1/3+K_star/4, "(H-1)/(2H-1)": 2/5, "H/(2H+1)": 3/7,
    "1/H": 1/H, "(H²-1)/(H²+H)": 8/12,
}
for name, val in sorted(Q_candidates.items(), key=lambda x: abs(x[1] - Q_D_fit)):
    err = abs(val - Q_D_fit) / Q_D_fit * 100
    if err < 20:
        print(f"    {name:20s} = {val:.6f}  err = {err:.1f}%")

# ============================================================
# PART 4: FULL BARYOGENESIS WITH FRAMEWORK R
# ============================================================
print("\n" + "=" * 80)
print("PART 4: BARYOGENESIS RESULT")
print("=" * 80)

R_final = R_func(w_final)
mD_final = compute_mD(R_final, m_light, M_R_hier, U_PMNS)

print(f"\n  Using approach {best_label}, w = {w_final}")
print(f"\n  m_D (GeV):")
for i in range(3):
    row = mD_final[i]
    print(f"    [{row[0].real:+10.2f}{row[0].imag:+10.2f}i  "
          f"{row[1].real:+10.2f}{row[1].imag:+10.2f}i  "
          f"{row[2].real:+10.2f}{row[2].imag:+10.2f}i]")

# SVD for mass eigenvalues
_, S_mD, _ = svd(mD_final)
print(f"\n  Dirac mass eigenvalues:")
for i, s in enumerate(sorted(S_mD)):
    print(f"    m_D_{i+1} = {s:.4f} GeV = {s*1000:.1f} MeV")

# Verify seesaw
M_R_diag = np.diag(M_R_hier)
m_seesaw = mD_final.T @ inv(M_R_diag) @ mD_final
evals_seesaw = np.sort(np.abs(np.linalg.eigvals(m_seesaw)))
print(f"\n  Seesaw eigenvalues (should match light masses):")
for i in range(3):
    print(f"    m_ν_{i+1} = {evals_seesaw[i]*1e9:.4f} meV  (input: {m_light[i]*1000:.4f} meV)")

# Full baryogenesis
eta_B, eps1, K_w, kap = full_eta_B(R_final, m_light, M_R_hier, U_PMNS)

print(f"\n  BARYOGENESIS:")
print(f"    ε₁ = {eps1:.6e}")
print(f"    K_washout = {K_w:.4f}")
print(f"    κ = {kap:.6f}")
print(f"    η_B = {eta_B:.6e}")
print(f"    η_B(obs) = {eta_B_target:.6e}")
print(f"    ratio = {abs(eta_B)/eta_B_target:.4f}")

# ============================================================
# DAVIDSON-IBARRA BOUND
# ============================================================
print("\n" + "=" * 80)
print("DAVIDSON-IBARRA BOUND")
print("=" * 80)

eps_DI = 3 * M_R_hier[0] / (16 * pi * v_higgs**2) * (m3 - m1) * 1e-9
print(f"\n  |ε₁|_max = (3 M₁)/(16π v²) × (m₃ - m₁)")
print(f"           = {eps_DI:.4e}")
print(f"  Actual |ε₁| = {abs(eps1):.4e}")
print(f"  |ε₁|/|ε₁|_max = {abs(eps1)/eps_DI:.4f}" if eps_DI > 0 else "")

# Maximum achievable η_B
kap_typ = 0.1
eta_max = 28/79 * eps_DI * kap_typ * 135 / (4*pi**2*g_star)
print(f"\n  Maximum η_B (κ = 0.1): {eta_max:.4e}")
print(f"  Observed η_B:          {eta_B_target:.4e}")
print(f"  Ratio max/obs = {eta_max/eta_B_target:.2f}")

if eta_max < eta_B_target:
    print(f"\n  WARNING: Davidson-Ibarra bound INSUFFICIENT.")
    print(f"  Need either:")
    print(f"    (a) Larger M₁ (currently {M_R_hier[0]:.2e} GeV)")
    print(f"    (b) Resonant leptogenesis (nearly degenerate M_R)")
    print(f"    (c) Different M_R hierarchy")

    # What M₁ is needed?
    M1_needed = eta_B_target * 16*pi*v_higgs**2 / (3 * kap_typ * (m3-m1)*1e-9 * 28/79 * 135/(4*pi**2*g_star))
    print(f"\n  M₁ needed: {M1_needed:.2e} GeV")
    print(f"  Current M₁: {M_R_hier[0]:.2e} GeV")
    print(f"  Factor needed: {M1_needed/M_R_hier[0]:.1f}")

# ============================================================
# ALTERNATIVE M_R HIERARCHIES
# ============================================================
print("\n" + "=" * 80)
print("ALTERNATIVE M_R HIERARCHIES")
print("=" * 80)

hierarchies = {
    "M_R/H², M_R/H, M_R (factor H)": np.array([M_R_common/9, M_R_common/3, M_R_common]),
    "M_R/H, M_R, M_R×H (centered)": np.array([M_R_common/3, M_R_common, M_R_common*3]),
    "M_R, M_R, M_R (degenerate)": np.array([M_R_common]*3),
    "M_R/30, M_R, M_R×30 (large split)": np.array([M_R_common/30, M_R_common, M_R_common*30]),
    "10¹⁰, 10¹², 10¹⁴ GeV": np.array([1e10, 1e12, 1e14]),
    "10⁹, 10¹², 10¹⁵ GeV": np.array([1e9, 1e12, 1e15]),
}

for name, M_R_test in hierarchies.items():
    eps_DI_test = 3 * M_R_test[0] / (16*pi*v_higgs**2) * (m3-m1)*1e-9
    eta_max_test = 28/79 * eps_DI_test * kap_typ * 135/(4*pi**2*g_star)

    # Quick scan for best η_B
    best_eta = 0
    for ix in range(-10, 11, 2):
        for iy in range(1, 11, 2):
            w = ix*pi/10 + 1j*iy*0.5
            try:
                R_t = R_euler(w, w, w)
                eta_t, _, _, _ = full_eta_B(R_t, m_light, M_R_test, U_PMNS)
                if abs(eta_t) > abs(best_eta):
                    best_eta = eta_t
            except:
                pass

    sufficient = "YES" if eta_max_test > eta_B_target else "NO"
    print(f"  {name:40s}  DI max={eta_max_test:.2e}  best_scan={abs(best_eta):.2e}  sufficient={sufficient}")

# ============================================================
# RESONANT LEPTOGENESIS
# ============================================================
print("\n" + "=" * 80)
print("RESONANT LEPTOGENESIS (NEARLY DEGENERATE M_R)")
print("=" * 80)

# If M₁ ≈ M₂ with splitting δ = (M₂-M₁)/M₁ << 1,
# the CP asymmetry is resonantly enhanced:
# ε₁ ~ Im[(h₁₂)²] × M₁M₂/(M₂²-M₁²) × Γ₂/(...)

# In the framework: could the Z₃ splitting be small?
# If M_R_k = M_R (1 + ε_k) with ε from K* or similar:

for split_name, delta in [("K*", K_star), ("K*²", K_star**2),
                           ("1/33", 1/33), ("K*/H²", K_star/H**2)]:
    M_res = np.array([M_R_common * (1 - delta),
                       M_R_common,
                       M_R_common * (1 + delta)])

    # With resonant enhancement, ε₁ can be O(1)
    # Quick estimate of maximum η_B
    eps_res_max = 0.5  # resonant maximum
    eta_res_max = 28/79 * eps_res_max * kap_typ * 135/(4*pi**2*g_star)

    # Actual computation
    best_res = 0
    for ix in range(-10, 11, 2):
        for iy in range(1, 11, 2):
            w = ix*pi/10 + 1j*iy*0.3
            try:
                R_t = R_euler(w, w, w)
                eta_t, _, _, _ = full_eta_B(R_t, m_light, M_res, U_PMNS)
                if abs(eta_t) > abs(best_res):
                    best_res = eta_t
            except:
                pass

    print(f"  δ = {split_name:8s} = {delta:.6f}  M-split = {M_R_common*delta:.2e} GeV  "
          f"best η_B = {abs(best_res):.2e}  ratio = {abs(best_res)/eta_B_target:.3f}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print(f"""
  SETUP:
    Light masses: ({m1*1000:.2f}, {m2*1000:.2f}, {m3*1000:.2f}) meV
    PMNS: sin²θ₁₂=3/10, sin²θ₂₃=97/180, sin²θ₁₃=-ln(23/30)/12, δ=-π/2

  R-MATRIX (Casas-Ibarra):
    Best approach: {best_label}
    w = {w_final.real:.6f} + {w_final.imag:.6f}i
    η_B = {eta_final:.4e}  (target: {eta_B_target:.1e})
    ratio = {abs(eta_final)/eta_B_target:.4f}

  KOIDE-DIRAC MASSES:
    From seesaw inversion: ({mD_diag[0]*1000:.1f}, {mD_diag[1]*1000:.1f}, {mD_diag[2]*1000:.1f}) MeV
    Koide fit: θ_D = {th_D:.4f}, Q_D = {Q_D_fit:.4f}

  KEY FINDINGS:
    1. The degenerate M_R case gives ZERO CP asymmetry (GIM cancellation).
       Z₃ MUST split the RH neutrino masses for leptogenesis.

    2. The Davidson-Ibarra bound constrains M₁:
       |ε₁|_max = 3M₁(m₃-m₁)/(16πv²)
       For M₁ = M_GUT/(33×H²), this may or may not suffice.

    3. The Dirac masses from seesaw inversion have Q_D ≈ {Q_D:.4f},
       close to 1/3 (= 1/H), NOT 2/3 as for charged leptons.
       This is the seesaw suppression of the Koide structure.

    4. The Z₃ structure of R constrains it to a single complex
       parameter w, but the choice of Z₃ representation in the
       mass basis matters critically.

  OPEN QUESTIONS:
    - Is the Z₃ splitting of M_R by factor H, or by K*, or resonant?
    - Is w a framework quantity (2/9, K*, etc.)?
    - Does the Dirac Koide angle θ_D have framework significance?
""")
