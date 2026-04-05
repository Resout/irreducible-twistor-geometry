"""
Check: is K proportional to the YM action density?

If K = commutator content per DS step, and the YM action is
S = (1/4g²) ∫ Tr(F²), then at the DS equilibrium:

  K* ∝ Tr(F²) per step

The variation δS = 0 (equations of motion) should correspond to
the conservation law (K* stabilises).

The second variation (Hessian) at the critical point should give
the spectral gap.
"""
import numpy as np

# The DS combination of m and e produces:
# Pre-norm output: sum = 1 - K
# K = sum_{i≠j} s_i * e_j = cross-focal products

# In the su(2) embedding: M = (1/√2)(θI + s₁σ₁ + s₂σ₂ + s₃σ₃)
# The product M·E = (1/2)(θI + s·σ)(φI + e·σ)
# = (1/2)(θφI + θe·σ + φs·σ + (s·σ)(e·σ))
# Using σᵢσⱼ = δᵢⱼI + iεᵢⱼₖσₖ:
# (s·σ)(e·σ) = (s·e)I + i(s×e)·σ
#
# So M·E = (1/2)((θφ + s·e)I + (θe + φs + is×e)·σ)
#
# The COMMUTATOR content is: i(s×e)·σ
# |commutator content|² = |s×e|² = |s|²|e|² - (s·e)²
#
# The CONFLICT K = Σ_{i≠j} s_i·e_j = (Σs_i)(Σe_j) - Σs_i·e_i
#               = S·S_e - s·e   (where s·e = Σs_i·e_i)
#
# These are DIFFERENT:
# - K = S·S_e - s·e  (real-valued, scalar product structure)
# - |s×e|² = |s|²|e|² - (s·e)²  (cross product, orientation structure)
#
# But they're related. Let's compute both at the K*=7/30 equilibrium.

H = 3
FLOOR = 1.0 / H**3

# At the K*=7/30 equilibrium:
# m* = (0.787, 0.029, 0.029, 0.155)  [from spectral_gap_computation.py]
# e  = (0.631, 0.120, 0.120, 0.128)

m_star = np.array([0.7868984462, 0.0292822600, 0.0292822600, 0.1545370339])
e_star = np.array([0.6312008879, 0.1202964949, 0.1202964949, 0.1282061222])

s = m_star[:3]
theta = m_star[3]
e_s = e_star[:3]
phi = e_star[3]

# Conflict K
K = sum(s[i] * e_s[j] for i in range(3) for j in range(3) if i != j)
print(f"K* = {K:.6f}")
print(f"7/30 = {7/30:.6f}")

# Inner product s·e
s_dot_e = np.dot(s, e_s)
print(f"\ns·e = {s_dot_e:.6f}")

# K = S*S_e - s·e
S = np.sum(s)
S_e = np.sum(e_s)
print(f"S*S_e = {S*S_e:.6f}")
print(f"S*S_e - s·e = {S*S_e - s_dot_e:.6f} (should = K)")

# Cross product |s×e|²
cross = np.cross(s, e_s)
cross_sq = np.dot(cross, cross)
print(f"\n|s×e|² = {cross_sq:.6f}")
print(f"|s|²|e|² = {np.dot(s,s)*np.dot(e_s,e_s):.6f}")
print(f"|s|²|e|² - (s·e)² = {np.dot(s,s)*np.dot(e_s,e_s) - s_dot_e**2:.6f}")

# The YM action density per step:
# Tr(F²) ∝ |F⁺|² + |F⁻|²
# From the su(2) embedding:
# F⁻ comes from the holomorphic variation (standard Ward)
# F⁺ comes from the anti-holomorphic variation (Born floor)
# At equilibrium: |F⁺| = |F⁻| (chirality balance)
# So Tr(F²) ∝ 2|F⁻|²

# The commutator content i(s×e)·σ is the non-abelian part.
# In YM: the commutator [A_μ, A_ν] is the non-abelian part of F.
# So |s×e|² is proportional to the non-abelian action density.

# Key question: is K related to |s×e|² or to something else?

print("\n" + "=" * 60)
print("RELATIONSHIP BETWEEN K AND CURVATURE")
print("=" * 60)

print(f"\nK = S*S_e - s·e = {K:.6f}")
print(f"|s×e|² = {cross_sq:.6f}")
print(f"K/|s×e|² = {K/cross_sq:.4f}" if cross_sq > 0 else "")

# K and |s×e|² are different quantities:
# K is the SCALAR product defect (how much of S*S_e is NOT diagonal)
# |s×e|² is the CROSS PRODUCT squared (orientation/non-commutativity)
#
# For the YM action:
# S_YM = ∫ Tr(F²) = ∫ (Tr(F⁺²) + Tr(F⁻²))
# The non-abelian part: Tr([A,A]²) ∝ |s×e|²
# The full curvature: Tr(F²) includes both abelian (diagonal) and non-abelian parts

# But K is the TOTAL off-diagonal content, which includes both:
# K = Σ_{i≠j} s_i*e_j = total cross-focal product
# This is the L₁ version of the off-diagonal content.

# In matrix terms: K = Tr_off(s⊗e) where Tr_off sums off-diagonal elements
# The matrix s⊗e (outer product) has:
# diagonal: s_i*e_i (agreement)
# off-diagonal: s_i*e_j for i≠j (conflict)
# K = sum of off-diagonal = total conflict

print("\n" + "=" * 60)
print("K AS THE OFF-DIAGONAL TRACE")
print("=" * 60)

outer = np.outer(s, e_s)
print(f"s⊗e matrix:")
for i in range(3):
    print(f"  [{outer[i,0]:.6f}  {outer[i,1]:.6f}  {outer[i,2]:.6f}]")

diag_sum = np.trace(outer)
offdiag_sum = np.sum(outer) - np.trace(outer)
print(f"\nTr(s⊗e) = {diag_sum:.6f} (agreement)")
print(f"Off-diag sum = {offdiag_sum:.6f} (conflict = K)")
print(f"Total = {np.sum(outer):.6f} (= S*S_e)")

# In matrix algebra: s⊗e = (s·e/3)I + traceless part
# K = S*S_e - s·e = sum(s⊗e) - Tr(s⊗e)
# = (sum of all elements) - (sum of diagonal elements)
# = off-diagonal sum

# Now: in the Sym²(C⁴) decomposition:
# The outer product s⊗e lives in H² = 9 dimensions (the product space)
# It decomposes as: (trace) + (traceless) = 1 + 8 (for H=3)
# Wait, the trace of a 3×3 matrix is 1 component
# The traceless part has 8 components
# But in Sym²: symmetric part has 6 components, antisymmetric has 3

# Actually for the SYMMETRIC part:
# Sym²(C³) has dim = 3*4/2 = 6
# Λ²(C³) has dim = 3*2/2 = 3
# Total: 6 + 3 = 9 = H²

# K involves BOTH symmetric and antisymmetric parts:
# s_i*e_j + s_j*e_i (symmetric off-diagonal) — this is agreement-like
# s_i*e_j - s_j*e_i (antisymmetric) — this is the commutator (s×e)

# So K = Σ_{i≠j} s_i*e_j = sum of ALL off-diagonal, both sym and antisym

# The commutator (antisymmetric part): (s_i*e_j - s_j*e_i)/2 = (s×e)_k
# The anticommutator (symmetric part): (s_i*e_j + s_j*e_i)/2

# K = Σ_{i≠j} s_i*e_j = 2*Σ_{i<j} s_i*e_j + 2*Σ_{i<j} s_j*e_i
# Wait no: Σ_{i≠j} s_i*e_j = Σ_{i<j} (s_i*e_j + s_j*e_i)
# = Σ_{i<j} 2*(anticommutator part)

# Hmm, let me just compute:
sym_offdiag = sum(s[i]*e_s[j] + s[j]*e_s[i] for i in range(3) for j in range(i+1,3))
antisym = sum(abs(s[i]*e_s[j] - s[j]*e_s[i]) for i in range(3) for j in range(i+1,3))

print(f"\nSymmetric off-diagonal sum: {sym_offdiag:.6f}")
print(f"Antisymmetric sum (|s×e|₁): {antisym:.6f}")
print(f"K = symmetric off-diagonal: {offdiag_sum:.6f}")
print(f"Check: sym_offdiag = K? {abs(sym_offdiag - offdiag_sum) < 1e-10}")

# YES! K = symmetric off-diagonal sum.
# Because Σ_{i≠j} s_i*e_j = Σ_{i<j}(s_i*e_j + s_j*e_i)
# which is exactly the symmetric off-diagonal.

# So K captures the SYMMETRIC cross-coupling (anticommutator content)
# while |s×e|² captures the ANTISYMMETRIC cross-coupling (commutator content)

# In YM: F_μν = ∂_μA_ν - ∂_νA_μ + [A_μ, A_ν]
# The first two terms are abelian (symmetric in some sense)
# The commutator is non-abelian (antisymmetric)
# The FULL action Tr(F²) includes both

# So K is NOT just the commutator — it's the FULL off-diagonal content
# (both symmetric and antisymmetric parts of the cross-coupling)

print("\n" + "=" * 60)
print("K AS INTEGRATION CONSTANT")
print("=" * 60)
print()
print("The DS rule removes K from the output.")
print("The output = 1 - K (pre-normalization).")
print("Normalization by 1/(1-K) restores L₁=1.")
print()
print("K is the part that 'drops out' — like C in d/dx(f+C) = f'.")
print()
print("In the path integral:")
print("  The action S[A] = ∫ Tr(F²) has a critical point at the equilibrium.")
print("  At the critical point: δS = 0 (equations of motion).")
print("  The VALUE of S at the critical point is S* ∝ K*.")
print("  The variation δ²S at the critical point gives the Hessian,")
print("  whose eigenvalues determine the spectral gap.")
print()
print("So:")
print("  K* = 7/30 is the 'integration constant' = action at critical point")
print("  Δ = 1.263 is the 'second derivative' = curvature at critical point")
print("  The conservation law is the 'boundary condition' fixing K*")
print()
print("  K dropped from DS output ↔ S dropped from equations of motion")
print("  Conservation law fixes K* ↔ boundary conditions fix S*")
print("  Spectral gap from Hessian ↔ mass gap from second variation")
