"""
Spectral gap as a Born probability of the transfer crystal.

The DS transfer operator at the fixed point IS a crystal —
a 4×4 relationship between "state before" and "state after."
Its eigenvalues are properties of that crystal.

The structural ingredients:
  - K* = 7/30 (from the conservation law)
  - Born floor: 26θ² = Σs² (the surface)
  - Sym²(C⁴): the combination is symmetric

The transfer crystal is built from these directly.
The eigenvalue λ₀ is a Born probability of the composed crystal.
"""

import numpy as np

# ============================================================
# The DS combination as a crystal
# ============================================================

# The DS combination of m with evidence e, for H=3:
#   s_i' = (s_i·e_i + s_i·φ + θ·e_i) / (1-K)
#   θ'   = θ·φ / (1-K)
#
# This is a LINEAR map in m (for fixed e).
# Linear map on C⁴ → C⁴ is a 4×4 matrix.
# That matrix IS the crystal.

def ds_crystal(e):
    """The DS combination as a 4×4 linear operator on m, for fixed evidence e.

    Returns the matrix T such that m_out = T @ m / (1-K(m,e)).
    But K depends on m, so this isn't quite linear...

    Actually: the pre-normalization output IS linear in m:
      s_i_raw = e_i * s_i + φ * s_i + e_i * θ = (e_i + φ) * s_i + e_i * θ
      θ_raw = φ * θ

    This is linear. The normalization 1/(1-K) is where it gets nonlinear
    (K depends on m). But at the FIXED POINT, K is constant (= 7/30).
    So the linearized operator IS the crystal.
    """
    a, b1, b2, phi = e

    # T[i,j] = coefficient of m[j] in the i-th output (pre-normalization)
    # s1_raw = (a + φ) * s1 + 0 * s2 + 0 * s3 + a * θ
    # s2_raw = 0 * s1 + (b1 + φ) * s2 + 0 * s3 + b1 * θ
    # s3_raw = 0 * s1 + 0 * s2 + (b2 + φ) * s3 + b2 * θ
    # θ_raw  = 0 * s1 + 0 * s2 + 0 * s3 + φ * θ

    T_pre = np.array([
        [a + phi,  0,        0,        a  ],
        [0,        b1 + phi, 0,        b1 ],
        [0,        0,        b2 + phi, b2 ],
        [0,        0,        0,        phi],
    ])

    return T_pre

# ============================================================
# The floor projection as a crystal
# ============================================================

def floor_crystal(m_pre):
    """Jacobian of the Born floor projection at a given pre-floor state.

    On the surface Born(θ)=1/27, the floor acts as a projection.
    The Jacobian of this projection maps perturbations in R⁴ to
    perturbations tangent to the surface.
    """
    s = m_pre[:3]
    sum_s = np.sum(s)
    L2_s = np.sqrt(np.sum(s**2))

    # The floor maps m_pre → (α*s₁, α*s₂, α*s₃, θ_new)
    # where α = (1-θ_new)/sum_s and 26θ_new² = α²·Σs²
    # θ_new = 1/(√26 · sum_s/L2_s + 1)

    R = sum_s / L2_s
    th_new = 1.0 / (np.sqrt(26) * R + 1)
    alpha = (1 - th_new) / sum_s

    # Jacobian of (s_out, θ_out) = (α(s)·s, θ_new(s)) w.r.t. s_pre
    # This is the projection crystal.
    eps = 1e-9
    J = np.zeros((4, 4))
    for j in range(4):
        mp = m_pre.copy()
        mp[j] += eps
        # Recompute floor
        s_p = mp[:3]
        sum_sp = np.sum(s_p)
        L2_sp = np.sqrt(np.sum(s_p**2))
        Rp = sum_sp / L2_sp
        th_p = 1.0 / (np.sqrt(26) * Rp + 1)
        alpha_p = (1 - th_p) / sum_sp
        out_p = np.array([alpha_p*s_p[0], alpha_p*s_p[1], alpha_p*s_p[2], th_p])

        out_0 = np.array([alpha*s[0], alpha*s[1], alpha*s[2], th_new])
        J[:, j] = (out_p - out_0) / eps

    return J

# ============================================================
# Build the transfer crystal from structural ingredients
# ============================================================

print("=" * 60)
print("THE TRANSFER CRYSTAL")
print("=" * 60)

# The known structural values
K_star = 7.0 / 30.0
one_minus_K = 1.0 - K_star  # = 23/30

# The fixed point (computed, but let's see if we even need it)
m_star = np.array([0.7868984462, 0.0292822600, 0.0292822600, 0.1545370339])
e_star = np.array([0.6312008879, 0.1202964949, 0.1202964949, 0.1282061222])

# The DS crystal (pre-normalization, linear in m)
T_ds = ds_crystal(e_star)
print("DS crystal (pre-normalization):")
print(T_ds.round(6))

# Normalize by 1/(1-K) at the fixed point
T_normalized = T_ds / one_minus_K
print(f"\nDS crystal (normalized by 1/(1-K) = 30/23):")
print(T_normalized.round(6))

# The pre-floor output
m_pre = T_normalized @ m_star
print(f"\nPre-floor output: {m_pre}")
print(f"L₁ = {np.sum(m_pre):.10f}")

# Floor projection crystal
P_floor = floor_crystal(m_pre)
print(f"\nFloor projection crystal:")
print(P_floor.round(6))

# The FULL transfer crystal: floor ∘ normalized_DS
T_full = P_floor @ T_normalized
print(f"\nFull transfer crystal (floor ∘ DS/(1-K)):")
print(T_full.round(6))

evals = np.sort(np.abs(np.linalg.eigvals(T_full)))[::-1]
print(f"\nEigenvalues: {evals}")
print(f"λ₀ = {evals[0]:.10f}")
print(f"Δ = -ln(λ₀) = {-np.log(evals[0]):.10f}")

# ============================================================
# Now: can we build this WITHOUT knowing m* and e*?
# ============================================================
print("\n" + "=" * 60)
print("STRUCTURAL CONSTRUCTION (no m*, no e*)")
print("=" * 60)

# What do we actually know?
# 1. K* = 7/30
# 2. Both m* and e* are on the Born floor surface
# 3. The fixed point has S₂ symmetry (s₂=s₃)
# 4. The DS combination is symmetric (Sym²)
#
# On the Born floor surface with S₂ symmetry, the state is
# determined by one parameter (s). Similarly the evidence.
#
# The DS crystal T_ds depends on e = (a, b, b, φ).
# The floor crystal P depends on the pre-floor output.
# The eigenvalue of T = P @ T_ds/(1-K) is what we want.
#
# KEY INSIGHT: the DS crystal T_ds is DIAGONAL in the
# singleton entries. The off-diagonal coupling comes ONLY
# from the θ column (commitment terms θ·e_i).
#
# In the S₂-symmetric sector, the 2×2 block that contains λ₀
# involves only the dominant singleton s and θ.
# The 2×2 DS crystal in this sector is:
#   [(a+φ)/(1-K),  a/(1-K)]
#   [φ_row...                ]
# Wait, let me think about this more carefully.

# The 4×4 T_normalized has structure:
# s₁ row: [(a+φ)/(1-K),  0,  0,  a/(1-K)]
# s₂ row: [0,  (b+φ)/(1-K),  0,  b/(1-K)]
# s₃ row: [0,  0,  (b+φ)/(1-K),  b/(1-K)]
# θ  row: [0,  0,  0,  φ/(1-K)]

print("T_normalized structure:")
print(f"  diagonal: {np.diag(T_normalized)}")
print(f"  θ column: {T_normalized[:, 3]}")

# The floor projection then mixes everything. But in the
# S₂-symmetric subspace, we only need the 2D block.

# Let me see what the eigenvalue depends on structurally.
# From Part 6 of the previous script, the eigenvalue is
# the TRACE of the 2×2 symmetric block (since det ≈ 0).
#
# trace = J_full_sym[0,0] + J_full_sym[1,1]
#       = (radial,radial) + (symmetric,symmetric) entries

# Let me parametrize. Define:
# r = s/w (the singleton ratio at the fixed point)
# At the fixed point, r = 26.87

# The DS crystal eigenvalues (before floor) are:
# a+φ, b+φ, b+φ, φ  (all divided by 1-K)
# The floor then projects, killing one direction.

# The SURVIVING eigenvalue is a weighted combination:
# λ₀ = α·(a+φ)/(1-K) + β·φ/(1-K)
# where α,β are the floor projection weights.

# Let me check: if the floor projects onto a specific direction,
# what's the effective eigenvalue?

# The floor has eigenvalues 0.868, 0.868, ~0, 0.
# Two eigenvalues equal = it preserves a 2D subspace.
# In that subspace, it scales by 0.868 (= α from the floor formula).

alpha_floor = float(np.sort(np.abs(np.linalg.eigvals(P_floor)))[-1])
print(f"\nFloor projection factor α = {alpha_floor:.10f}")

# The DS crystal (normalized) has diagonal values:
ds_diag = np.diag(T_normalized)
print(f"DS diagonal: {ds_diag}")

# Full eigenvalue ≈ α_floor × (DS eigenvalue in the surviving direction)
# The DS part in the S₂-symmetric sector has eigenvalue (a+φ)/(1-K)
ds_dominant = (e_star[0] + e_star[3]) / one_minus_K
ds_weak = (e_star[1] + e_star[3]) / one_minus_K
ds_theta = e_star[3] / one_minus_K
print(f"\nDS eigenvalues (normalized): dominant={(e_star[0]+e_star[3])/one_minus_K:.6f}, "
      f"weak={(e_star[1]+e_star[3])/one_minus_K:.6f}, θ={e_star[3]/one_minus_K:.6f}")

# Hmm, dominant DS eigenvalue = 0.99 * 0.868 ≈ 0.86, not 0.283.
# So it's not just α × DS_eigenvalue.

# Let me look at what the 2×2 block actually is in the
# (radial, θ) subspace BEFORE the floor kills one direction.

# ============================================================
# The structural formula
# ============================================================
print("\n" + "=" * 60)
print("SEARCHING FOR THE STRUCTURAL FORMULA")
print("=" * 60)

# The eigenvalue 0.2829 — what simple expressions give this?
target = 0.28291034

# Candidates involving K* = 7/30 and H = 3:
candidates = {
    "K*": K_star,
    "1-K*": 1 - K_star,
    "K*²": K_star**2,
    "(1-K*)²": (1-K_star)**2,
    "K*/(1+K*)": K_star/(1+K_star),
    "K*·(1-K*)": K_star*(1-K_star),
    "√(K*(1-K*))": np.sqrt(K_star*(1-K_star)),
    "1/(1+1/K*)": 1/(1+1/K_star),
    "K*^(1/K*)": K_star**(1/K_star),
    "(1-K*)^(1/(1-K*))": (1-K_star)**(1/(1-K_star)),
    "exp(-1/K*)·K*": np.exp(-1/K_star)*K_star,
    "K*·e^(-1)": K_star * np.exp(-1),
    "(7/30)·(23/30)": (7/30)*(23/30),
    "(23/30)^(30/23)": (23/30)**(30/23),
    "(23/30)^(7/23)": (23/30)**(7/23),
    "7/(3·(H²+1))": 7/(3*10),
    "7²/(H·(H²+1)²)": 49/(3*100),
    "η·K*": (4/27)*(7/30),
    "(H-1)²·K*/H²": 4*7/(9*30),
    "K*^(3/2)": K_star**1.5,
    "1/H - K*": 1/3 - 7/30,
    "(1/H - K*)·H": (1/3 - 7/30)*3,
    "7/23": 7/23,
    "7/(23+7)": 7/30,
    "7·23/30²": 7*23/900,
    "exp(-Δ_filter)·(1-something)": (23/30) * 0.37,  # wild guess
    "(23/30)^4": (23/30)**4,
    "(23/30)^(4/3)": (23/30)**(4/3),
    "(23/30)^(10/7)": (23/30)**(10/7),
    "K*/(1-K*²)": K_star/(1-K_star**2),
    "Born·(1-K*)": (1/27)*(23/30),
    "η/(1-K*)": (4/27)/(23/30),
    "η·(1+K*)": (4/27)*(1+7/30),
    "26/27·K*": (26/27)*(7/30),
    "(7/30)^(7/30)·(23/30)^(23/30)": (7/30)**(7/30) * (23/30)**(23/30),
}

print(f"Target: λ₀ = {target:.8f}\n")
ranked = sorted(candidates.items(), key=lambda x: abs(x[1] - target))
for name, val in ranked[:15]:
    print(f"  {name:40s} = {val:.8f}  (diff = {val-target:+.6f})")

# Check: is λ₀ related to the antisymmetric eigenvalue λ₁ = 0.28131?
lambda_1 = 0.28131300
print(f"\nλ₁ = {lambda_1:.8f}")
print(f"λ₀ - λ₁ = {target - lambda_1:.8f}")
print(f"λ₀/λ₁ = {target/lambda_1:.8f}")
print(f"λ₀² + λ₁² = {target**2 + lambda_1**2:.8f}")
print(f"λ₀ + λ₁ = {target + lambda_1:.8f}")
print(f"λ₀ · λ₁ = {target * lambda_1:.8f}")

# Check simple fractions near target
from fractions import Fraction
f = Fraction(target).limit_denominator(1000)
print(f"\nBest fraction: {f} = {float(f):.8f}")
f2 = Fraction(target).limit_denominator(100)
print(f"Simple fraction: {f2} = {float(f2):.8f}")

# What about λ₁?
f1 = Fraction(lambda_1).limit_denominator(1000)
print(f"λ₁ fraction: {f1} = {float(f1):.8f}")

# Sum
f_sum = Fraction(target + lambda_1).limit_denominator(1000)
print(f"λ₀+λ₁ fraction: {f_sum} = {float(f_sum):.8f}")
