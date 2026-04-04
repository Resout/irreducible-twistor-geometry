"""
Compute the Fubini-Study distance of one DS step.

The DS step takes m → DS(m, e*)/(1-K). At the fixed point,
this is the identity. For a perturbation δm, the step takes
m* + δm → m* + λ₀·δm (to leading order).

The Fubini-Study distance between m*+δm and m*+λ₀·δm is:
  d_FS = arccos(|<m*+δm, m*+λ₀·δm>| / (|m*+δm|·|m*+λ₀·δm|))

For small δm, this simplifies to something proportional to
(1-λ₀)|δm|.

The PHYSICAL step length should be related to the Fubini-Study
distance at some characteristic perturbation scale.
"""
import numpy as np

m_star = np.array([0.7868984462, 0.0292822600, 0.0292822600, 0.1545370339])
lambda_0 = 0.2829102944

# Fubini-Study distance on CP³:
# d_FS(ψ₁, ψ₂) = arccos(|<ψ₁,ψ₂>|/(|ψ₁||ψ₂|))

def fubini_study(m1, m2):
    """Fubini-Study distance between two mass functions on CP³."""
    inner = np.abs(np.dot(m1, np.conj(m2)))
    norm1 = np.linalg.norm(m1)
    norm2 = np.linalg.norm(m2)
    cos_d = inner / (norm1 * norm2)
    cos_d = min(cos_d, 1.0)  # numerical safety
    return np.arccos(cos_d)

# The characteristic scale on CP³:
# The maximum FS distance is π/2 (between orthogonal states).
# The Born floor constrains the state to a region of CP³.
# What's the FS diameter of the floor-constrained region?

# The floor says Born(θ) ≥ 1/27.
# The "most different" states consistent with the floor are:
# State 1: dominant s₁, floor θ
# State 2: dominant s₂, floor θ
# These differ in which hypothesis dominates.

# State 1: s₁ large, s₂=s₃ small, θ at floor
s_dom = 1 - 2*0.001 - m_star[3]  # nearly all mass on s₁
m1 = np.array([s_dom, 0.001, 0.001, m_star[3]])
m1 = m1 / np.sum(m1)

m2 = np.array([0.001, s_dom, 0.001, m_star[3]])
m2 = m2 / np.sum(m2)

m3 = np.array([0.001, 0.001, s_dom, m_star[3]])
m3 = m3 / np.sum(m3)

print("=" * 60)
print("FUBINI-STUDY DISTANCES ON CP³")
print("=" * 60)
print(f"  d(s₁-dom, s₂-dom) = {fubini_study(m1, m2):.6f}")
print(f"  d(s₁-dom, s₃-dom) = {fubini_study(m1, m3):.6f}")
print(f"  d(s₂-dom, s₃-dom) = {fubini_study(m2, m3):.6f}")
print(f"  Maximum FS distance: π/2 = {np.pi/2:.6f}")
print()

# The DS step at the fixed point contracts by factor λ₀:
# A perturbation at FS distance d from m* maps to FS distance λ₀·d.
# The "step size" is d·(1-λ₀).

# For a unit perturbation in the dominant direction:
delta = 0.01  # small perturbation
dm = np.array([delta, -delta/3, -delta/3, -delta/3])
m_pert = m_star + dm
m_pert_after = m_star + lambda_0 * dm  # leading order

d_before = fubini_study(m_star, m_pert)
d_after = fubini_study(m_star, m_pert_after)
d_step = fubini_study(m_pert, m_pert_after)

print(f"Unit perturbation test (δ={delta}):")
print(f"  d(m*, m*+δm) = {d_before:.6f}")
print(f"  d(m*, m*+λ₀δm) = {d_after:.6f}")
print(f"  d(m*+δm, m*+λ₀δm) = {d_step:.6f}")
print(f"  Ratio d_step/d_before = {d_step/d_before:.6f}")
print(f"  Expected: 1-λ₀ = {1-lambda_0:.6f}")
print()

# ============================================================
# The DS step length in physical units
# ============================================================
print("=" * 60)
print("DS STEP LENGTH")
print("=" * 60)
print()

# The FS diameter of CP³ is π/2.
# The DS step contracts perturbations by factor λ₀ per step.
# After n steps: perturbation ~ λ₀^n.
# The total FS path length traversed is:
# L = Σ_{k=0}^{n-1} d_step(k) = d₀ · Σ λ₀^k · (1-λ₀) = d₀ · (1-λ₀^n)
# For n→∞: L = d₀ (the initial FS distance).

# The "step length" per combination is the FS distance contracted:
# d_step = d · (1-λ₀) where d is the current perturbation distance.

# At the equilibrium, the characteristic perturbation scale is
# set by the evidence fluctuations. For σ_e ~ 0.01:
# d_characteristic ~ σ_e · (FS metric scale)

# The FS metric on CP³ at a point m:
# ds² = (|dm|²|m|² - |<dm,m>|²) / |m|⁴
# For L₁=1 masses with |m| ~ 1: ds² ~ |dm|² - |<dm,m>|²/|m|²

# The physical meaning: one DS step on the fibre corresponds to
# one "time step" in the transfer matrix formulation.
# The physical time step is a(β) (lattice spacing).
# The ratio c = a_DS/a(β) = 1.78 at β=2.3.

# Can we compute c from the geometry?
# The FS volume of CP³ is Vol(CP³) = π³/6.
# The FS volume of the floor-constrained region is smaller.
# The ratio might determine c.

print(f"Vol(CP³) = π³/6 = {np.pi**3/6:.6f}")
print()

# The floor-constrained region: Born(θ) ≥ 1/27.
# This constrains θ²/Σm² ≥ 1/27.
# On CP³ (projective), this is a specific subset.
# Its volume fraction determines how "much" of CP³ is accessible.

# Monte Carlo estimate of the floor-active fraction:
N_mc = 100000
np.random.seed(42)
n_active = 0
for _ in range(N_mc):
    m_rand = np.random.dirichlet([1,1,1,1])  # random L1=1 mass function
    born_theta = m_rand[3]**2 / np.sum(m_rand**2)
    if born_theta >= FLOOR:
        n_active += 1

frac = n_active / N_mc
print(f"Fraction of L1=1 simplex with Born(θ)≥1/27: {frac:.4f}")
print(f"  ({n_active}/{N_mc} samples)")
print()

# This is the fraction of the simplex, not of CP³.
# On CP³, the fraction would be different (projective measure).

# Let me try a different approach: the conversion factor c=1.78
# might be simply related to (1-K*) or 1/(1-K*).

print("Simple relationships:")
print(f"  1/(1-K*) = {1/(1-7/30):.4f}")
print(f"  c = 1.78")
print(f"  c × (1-K*) = {1.78 * (1-7/30):.4f}")
print(f"  c / (1/(1-K*)) = {1.78 / (1/(1-7/30)):.4f}")
print()
print(f"  H/(H-1) = {3/2:.4f}")
print(f"  c / (H/(H-1)) = {1.78/1.5:.4f}")
print()
print(f"  (H+1)/H = {4/3:.4f}")
print(f"  c × (H+1)/H = {1.78*4/3:.4f}")
print()

# c × (1-K*) = 1.78 × 23/30 = 1.365
# c / (30/23) = 1.78 × 23/30 = 1.365
# Not clean.

# c × H/(H-1) = 1.78 × 3/2 = 2.67 — this was the OLD ratio m_eff/Δ_filter!
print(f"  c × H/(H-1) = {1.78 * 3/2:.4f}")
print(f"  Old ratio m_eff/Δ_filter at β=2.3 = 2.67 ← matches!")
print()
print(f"  So c = (m_eff/Δ_filter) × (H-1)/H = 2.67 × 2/3 = 1.78")
print(f"  This is a consistency check, not a derivation.")
print()

# The conversion factor c doesn't have an obvious closed form
# in terms of the DS quantities. It depends on the lattice spacing a(β),
# which is a property of the lattice regularisation, not of the DS framework.
# Under Path (B), c is the ratio of two scales:
# a_DS (intrinsic to the DS framework) and a(β) (lattice artifact).
# a_DS should be derivable from the geometry of CP³.

print("CONCLUSION:")
print("The conversion factor c=1.78 at β=2.3 is the ratio a_DS/a(β).")
print("a_DS is an intrinsic geometric scale of the DS framework on CP³.")
print("Its derivation from first principles requires relating the")
print("Fubini-Study metric on CP³ to the physical metric on R⁴")
print("through the twistor fibration — an open computation.")
