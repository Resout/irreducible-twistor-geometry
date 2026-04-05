"""
Analytical spectral gap derivation
===================================

Goal: derive λ₀ = 0.2829 in closed form from the DS algebra at H=3.

The transfer operator at the K*=7/30 fixed point is:
  T(m) = floor(DS(m, e*) / (1-K))

where m* and e* are the fixed point mass function and evidence.

At the fixed point, Born(θ) = 1/27 exactly (floor saturated).
The Jacobian dT/dm splits into:
  1. The DS polynomial part (holomorphic)
  2. The floor enforcement (the only nonlinear part)

Strategy: compute the Jacobian symbolically.
"""

import numpy as np
from fractions import Fraction
import sympy as sp

# ============================================================
# Part 1: The fixed point in exact arithmetic
# ============================================================
print("=" * 60)
print("PART 1: Fixed point structure")
print("=" * 60)

# At the fixed point:
# - m* = (s, w, w, θ) with s + 2w + θ = 1 (L₁=1)
# - Born(θ) = θ²/(s² + 2w² + θ²) = 1/27
# - K(m*, e*) = 7/30
#
# The fixed-point equation m* = DS(m*, e*) / (1-K) gives:
# s = (s·a + s·φ + θ·a) / (1-K)
# w = (w·b + w·φ + θ·b) / (1-K)
# θ = θ·φ / (1-K)
#
# where e* = (a, b, b, φ), and K = 2sw·(1/a)(sums of cross terms)
#
# From the θ equation: θ = θφ/(1-K), so φ = 1-K = 23/30
# (This is exact!)

K_star = sp.Rational(7, 30)
phi = 1 - K_star  # = 23/30
print(f"K* = {K_star} = {float(K_star):.10f}")
print(f"φ (evidence ignorance) = 1 - K* = {phi} = {float(phi):.10f}")

# From θ equation: φ = 1-K = 23/30. This means:
# a + 2b + φ = 1 (evidence L₁)
# a + 2b = 1 - 23/30 = 7/30

print(f"a + 2b = {1 - phi} = {float(1-phi):.10f}")

# From s equation: s(1-K) = s·a + s·φ + θ·a
# s(1-K) = s(a+φ) + θ·a
# s(1-K-a-φ) = θ·a
# s(-K-a+1-K) = θ·a  ... let me redo this
# s·(1-K) = s·a + s·φ + θ·a
# s·(23/30) = s·a + s·(23/30) + θ·a
# 0 = s·a + θ·a
# 0 = a(s + θ)
#
# This gives a = 0, which can't be right. Let me re-examine.
#
# Wait — the φ in the DS rule is the evidence ignorance component,
# and a is the evidence dominant singleton. But at the fixed point
# under EXTERNAL evidence (not self-evidence), m* ≠ e*.
#
# Fixed point equation: m* = DS(m*, e*)/(1-K*)
#
# θ* = θ*·φ / (1-K*)  =>  1-K* = φ  =>  φ = 23/30
# s* = (s*·a + s*·φ + θ*·a) / (1-K*)
# w* = (w*·b + w*·φ + θ*·b) / (1-K*)
#
# From s: s*(1-K*) = s*a + s*φ + θ*a = s*(a+φ) + θ*a
#          s*(1-K*-a-φ) = θ*a
#          s*(-a) = θ*a   (since 1-K*=φ, so 1-K*-a-φ=-a)
#          Wait: 1-K* - a - φ = φ - a - φ = -a
#          So: -s*a = θ*a
#          If a ≠ 0: s* = -θ*  ... impossible for positive masses
#
# Something is wrong. Let me be more careful.
#
# Actually the conflict K depends on BOTH m and e:
# K = Σ_{i≠j} s_i * e_j  (m singletons × e singletons, cross terms)
#
# For m* = (s, w, w, θ) and e* = (a, b, b, φ):
# K = s·b + s·b + w·a + w·b + w·a + w·b
#   = 2sb + 2wa + 2wb
#   = 2sb + 2w(a+b)
#
# The DS output (pre-normalization):
# s_out = s·a + s·φ + θ·a  (agreement + commitment)
# w_out = w·b + w·φ + θ·b
# θ_out = θ·φ
#
# After normalization by 1/(1-K):
# s* = (sa + sφ + θa)/(1-K)
# w* = (wb + wφ + θb)/(1-K)
# θ* = θφ/(1-K)
#
# θ equation: θ*(1-K) = θ*φ  =>  1-K = φ  =>  φ = 23/30 ✓
#
# s equation: s*(1-K) = s*a + s*φ + θ*a
#             s*·23/30 = s*(a + 23/30) + θ*a
#             0 = s*a + θ*a
#             a(s* + θ*) = 0
#
# Since s* > 0 and θ* > 0, we need a = 0.
# But the evidence has a ≠ 0 numerically (a ≈ 0.631).
#
# I think the issue is that the fixed point is NOT under
# the condition m = DS(m,e)/(1-K). The floor enforcement
# CHANGES the output. At the fixed point, floor is active.
# So the fixed point equation is:
#   m* = floor(DS(m*, e*) / (1-K*))
# not just DS(m*,e*)/(1-K*).

print("\nThe fixed point includes floor enforcement!")
print("m* = floor(DS(m*, e*)/(1-K*))")
print("The floor modifies the pre-floor output to satisfy Born(θ)=1/27")

# ============================================================
# Part 2: Jacobian structure at the fixed point
# ============================================================
print("\n" + "=" * 60)
print("PART 2: Jacobian decomposition")
print("=" * 60)

# The full map is T = floor ∘ normalize ∘ DS_raw
# dT = d(floor) · d(normalize) · d(DS_raw)
#
# d(DS_raw) is polynomial — easy
# d(normalize) is division by (1-K(m,e)) — rational
# d(floor) is the projection onto {Born(θ)≥1/27}
#
# At the fixed point, floor is active with equality.
# The floor acts as a specific projection.

# Let me compute the Jacobian components numerically but understand
# the structure analytically.

# Known values (from spectral_gap_computation.py):
m_star = np.array([0.7868984462, 0.0292822600, 0.0292822600, 0.1545370339])
e_star = np.array([0.6312008879, 0.1202964949, 0.1202964949, 0.1282061222])

s, w, w2, th = m_star
a, b, b2, phi_val = e_star

print(f"m* = ({s:.10f}, {w:.10f}, {w:.10f}, {th:.10f})")
print(f"e* = ({a:.10f}, {b:.10f}, {b:.10f}, {phi_val:.10f})")
print(f"φ = {phi_val:.10f}, 23/30 = {23/30:.10f}")
print(f"φ ≈ 23/30? {abs(phi_val - 23/30) < 0.01}")

# Check: φ should NOT be exactly 23/30 because the floor changes things
K_actual = 2*s*b + 2*w*(a+b)
print(f"\nK(m*,e*) = {K_actual:.10f}")
print(f"7/30 = {7/30:.10f}")
print(f"1-K = {1-K_actual:.10f}")

# Pre-floor DS output
s_pre = (s*a + s*phi_val + th*a) / (1 - K_actual)
w_pre = (w*b + w*phi_val + th*b) / (1 - K_actual)
th_pre = th*phi_val / (1 - K_actual)

print(f"\nPre-floor output: ({s_pre:.10f}, {w_pre:.10f}, {w_pre:.10f}, {th_pre:.10f})")
print(f"L₁ pre-floor: {s_pre + 2*w_pre + th_pre:.10f}")

born_pre = th_pre**2 / (s_pre**2 + 2*w_pre**2 + th_pre**2)
print(f"Born(θ) pre-floor: {born_pre:.10f}")
print(f"1/27 = {1/27:.10f}")
print(f"Floor activates: {born_pre < 1/27 + 1e-10}")

# After floor enforcement → m*
print(f"\nAfter floor: m* = ({m_star[0]:.10f}, {m_star[1]:.10f}, {m_star[2]:.10f}, {m_star[3]:.10f})")
print(f"Born(θ) after floor: {m_star[3]**2/np.sum(m_star**2):.10f}")

# ============================================================
# Part 3: Decompose the Jacobian
# ============================================================
print("\n" + "=" * 60)
print("PART 3: Jacobian decomposition — DS part vs floor part")
print("=" * 60)

eps = 1e-8

def ds_raw(m, e):
    """DS combination without floor, with normalization."""
    s_m = m[:3]
    th_m = m[3]
    s_e = e[:3]
    th_e = e[3]

    s_new = s_m * s_e + s_m * th_e + th_m * s_e
    th_new = th_m * th_e
    K = sum(s_m[i]*s_e[j] for i in range(3) for j in range(3) if i!=j)

    out = np.zeros(4)
    out[:3] = s_new / (1-K)
    out[3] = th_new / (1-K)
    return out

def floor_enforce(m):
    """Enforce Born(θ)≥1/27, maintaining L₁=1."""
    s = m[:3].copy()
    th = m[3]
    s_sq = np.sum(s**2)
    born = th**2 / (s_sq + th**2)

    if born >= 1/27 - 1e-15:
        return m.copy()

    # Need to boost θ and rescale s to maintain L₁
    # Born(θ_new) = θ_new² / (α²Σs² + θ_new²) = 1/27
    # L₁: αΣs + θ_new = 1
    # From Born: 27θ_new² = α²Σs² + θ_new²  =>  26θ_new² = α²Σs²
    # From L₁: α = (1-θ_new)/Σs where Σs = sum(s)

    sum_s = np.sum(s)

    lo, hi = th, 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        alpha = (1 - mid) / sum_s if sum_s > 0 else 0
        s_trial = s * alpha
        born_trial = mid**2 / (np.sum(s_trial**2) + mid**2)
        if born_trial < 1/27:
            lo = mid
        else:
            hi = mid

    th_new = (lo + hi) / 2
    alpha = (1 - th_new) / sum_s
    out = np.zeros(4)
    out[:3] = s * alpha
    out[3] = th_new
    return out

def full_transfer(m, e):
    return floor_enforce(ds_raw(m, e))

# Jacobian of DS_raw (without floor)
J_ds = np.zeros((4, 4))
f0_ds = ds_raw(m_star, e_star)
for j in range(4):
    mp = m_star.copy()
    mp[j] += eps
    J_ds[:, j] = (ds_raw(mp, e_star) - f0_ds) / eps

# Jacobian of floor at the pre-floor output
J_floor = np.zeros((4, 4))
pre_floor_out = ds_raw(m_star, e_star)
for j in range(4):
    mp = pre_floor_out.copy()
    mp[j] += eps
    J_floor[:, j] = (floor_enforce(mp) - floor_enforce(pre_floor_out)) / eps

# Full Jacobian = J_floor @ J_ds
J_full = J_floor @ J_ds

print("Jacobian of DS_raw (no floor):")
for i in range(4):
    print(f"  [{J_ds[i,0]:+.8f}  {J_ds[i,1]:+.8f}  {J_ds[i,2]:+.8f}  {J_ds[i,3]:+.8f}]")

evals_ds = np.sort(np.abs(np.linalg.eigvals(J_ds)))[::-1]
print(f"  eigenvalues: {evals_ds}")

print("\nJacobian of floor projection:")
for i in range(4):
    print(f"  [{J_floor[i,0]:+.8f}  {J_floor[i,1]:+.8f}  {J_floor[i,2]:+.8f}  {J_floor[i,3]:+.8f}]")

evals_floor = np.sort(np.abs(np.linalg.eigvals(J_floor)))[::-1]
print(f"  eigenvalues: {evals_floor}")

print("\nFull Jacobian (J_floor @ J_ds):")
for i in range(4):
    print(f"  [{J_full[i,0]:+.8f}  {J_full[i,1]:+.8f}  {J_full[i,2]:+.8f}  {J_full[i,3]:+.8f}]")

evals_full = np.sort(np.abs(np.linalg.eigvals(J_full)))[::-1]
print(f"  eigenvalues: {evals_full}")

# Project to L₁=1 tangent space
V = np.array([[1,0,0,-1], [0,1,0,-1], [0,0,1,-1]], dtype=float).T
VTV_inv = np.linalg.inv(V.T @ V)

J_ds_proj = VTV_inv @ V.T @ J_ds @ V
J_full_proj = VTV_inv @ V.T @ J_full @ V

print("\nDS-only projected eigenvalues:", np.sort(np.abs(np.linalg.eigvals(J_ds_proj)))[::-1])
print("Full projected eigenvalues:", np.sort(np.abs(np.linalg.eigvals(J_full_proj)))[::-1])

# ============================================================
# Part 4: What does the DS-only Jacobian look like analytically?
# ============================================================
print("\n" + "=" * 60)
print("PART 4: Analytical DS Jacobian")
print("=" * 60)

# The DS step (pre-normalization) for m=(s₁,s₂,s₃,θ) with evidence e=(a,b,b,φ):
# s₁_out = (s₁a + s₁φ + θa) / (1-K)
# s₂_out = (s₂b + s₂φ + θb) / (1-K)  [same for s₃]
# θ_out  = θφ / (1-K)
# K = s₁(b+b) + s₂(a+b) + s₃(a+b)
#   = 2s₁b + (s₂+s₃)(a+b)
#
# ds₁_out/ds₁ = [(a+φ)(1-K) + (s₁a + s₁φ + θa)·2b] / (1-K)²
#             = (a+φ)/(1-K) + s₁_out · 2b/(1-K)
#
# At the fixed point s₁_out = s₁ = s, etc., so:
# ds₁/ds₁ = (a+φ)/(1-K) + 2bs/(1-K)

# Let me compute these derivatives analytically
K_val = K_actual
one_minus_K = 1 - K_val

# The four state variables
# m = (s, w, w, θ) evidence = (a, b, b, φ)

# ∂K/∂s = 2b
# ∂K/∂w = (a+b) for each w component
# ∂K/∂θ = 0

# ∂(s₁_raw)/∂s₁ = (a+φ) before normalization
# ∂(s₁_raw)/∂s₂ = 0
# ∂(s₁_raw)/∂θ = a
# etc.

# After normalization by 1/(1-K):
# ∂s₁_out/∂x = [∂(s₁_raw)/∂x · (1-K) + s₁_raw · ∂K/∂x] / (1-K)²
#             = ∂(s₁_raw)/∂x / (1-K) + s₁_out · ∂K/∂x / (1-K)

dK_ds = 2*b
dK_dw = a + b  # per w component
dK_dth = 0

# s₁_raw derivatives
ds1r_ds = a + phi_val
ds1r_dw1 = 0
ds1r_dw2 = 0
ds1r_dth = a

# s₂_raw derivatives (w component)
dw1r_ds = 0
dw1r_dw1 = b + phi_val
dw1r_dw2 = 0
dw1r_dth = b

# θ_raw derivatives
dthr_ds = 0
dthr_dw = 0
dthr_dth = phi_val

# After 1/(1-K) normalization, using quotient rule:
# d(f/(1-K))/dx = (df/dx)/(1-K) + f/(1-K) · (dK/dx)/(1-K)
#               = (df/dx)/(1-K) + output · dK/dx / (1-K)

J_analytic = np.zeros((4,4))

# Row for s (index 0)
s_out = s  # at fixed point
J_analytic[0,0] = ds1r_ds/one_minus_K + s_out * dK_ds/one_minus_K
J_analytic[0,1] = ds1r_dw1/one_minus_K + s_out * dK_dw/one_minus_K
J_analytic[0,2] = ds1r_dw2/one_minus_K + s_out * dK_dw/one_minus_K
J_analytic[0,3] = ds1r_dth/one_minus_K + s_out * dK_dth/one_minus_K

# Row for w₁ (index 1)
w_out = w  # at fixed point
J_analytic[1,0] = dw1r_ds/one_minus_K + w_out * dK_ds/one_minus_K
J_analytic[1,1] = dw1r_dw1/one_minus_K + w_out * dK_dw/one_minus_K
J_analytic[1,2] = dw1r_dw2/one_minus_K + w_out * dK_dw/one_minus_K
J_analytic[1,3] = dw1r_dth/one_minus_K + w_out * dK_dth/one_minus_K

# Row for w₂ (index 2) — same as w₁ by symmetry
J_analytic[2,:] = J_analytic[1,:]

# Row for θ (index 3)
th_out = th
J_analytic[3,0] = dthr_ds/one_minus_K + th_out * dK_ds/one_minus_K
J_analytic[3,1] = dthr_dw/one_minus_K + th_out * dK_dw/one_minus_K
J_analytic[3,2] = dthr_dw/one_minus_K + th_out * dK_dw/one_minus_K
J_analytic[3,3] = dthr_dth/one_minus_K + th_out * dK_dth/one_minus_K

print("Analytical DS Jacobian (no floor):")
for i in range(4):
    print(f"  [{J_analytic[i,0]:+.8f}  {J_analytic[i,1]:+.8f}  {J_analytic[i,2]:+.8f}  {J_analytic[i,3]:+.8f}]")

print("\nNumerical DS Jacobian (no floor):")
for i in range(4):
    print(f"  [{J_ds[i,0]:+.8f}  {J_ds[i,1]:+.8f}  {J_ds[i,2]:+.8f}  {J_ds[i,3]:+.8f}]")

print(f"\nMax difference: {np.max(np.abs(J_analytic - J_ds)):.2e}")

# ============================================================
# Part 5: Structure of the floor Jacobian
# ============================================================
print("\n" + "=" * 60)
print("PART 5: Floor Jacobian structure")
print("=" * 60)

# The floor enforcement at Born(θ)=1/27 exactly:
# Given pre-floor m_pre, compute m_out such that:
# 1. Born(θ_out) = 1/27
# 2. s_out_i ∝ s_pre_i (proportional rescaling of singletons)
# 3. L₁(m_out) = 1
#
# This means: s_out = α · s_pre, θ_out = θ_new
# with α = (1-θ_new)/Σs_pre and Born(θ_new) = 1/27
#
# Born: θ_new² / (α²Σs_pre² + θ_new²) = 1/27
# 26θ_new² = α²Σs_pre²
# α = √(26) · θ_new / √(Σs_pre²)
#
# L₁: α·Σs_pre + θ_new = 1
# √(26)·θ_new·Σs_pre/√(Σs_pre²) + θ_new = 1
# θ_new · (√(26)·Σs_pre/√(Σs_pre²) + 1) = 1
#
# Let R = Σs_pre / √(Σs_pre²) (ratio of L₁ to L₂ for singletons)
# θ_new = 1 / (√26·R + 1)
# α = √26·θ_new / √(Σs_pre²) = √26 / ((√26·R+1)·√(Σs_pre²))

# At the fixed point pre-floor values:
s_pre_vec = np.array([s_pre, w_pre, w_pre])
sum_s_pre = np.sum(s_pre_vec)
L2_s_pre = np.sqrt(np.sum(s_pre_vec**2))
R = sum_s_pre / L2_s_pre

th_new_formula = 1 / (np.sqrt(26)*R + 1)
alpha_formula = np.sqrt(26) * th_new_formula / L2_s_pre

print(f"Pre-floor output: s_pre=({s_pre:.8f}, {w_pre:.8f}, {w_pre:.8f}), θ_pre={th_pre:.8f}")
print(f"Σs_pre = {sum_s_pre:.8f}")
print(f"√(Σs²_pre) = {L2_s_pre:.8f}")
print(f"R = Σs/√(Σs²) = {R:.8f}")
print(f"θ_new = 1/(√26·R+1) = {th_new_formula:.8f}")
print(f"m*[3] = {m_star[3]:.8f}")
print(f"α = {alpha_formula:.8f}")
print(f"α·s_pre = {alpha_formula*s_pre:.8f}, m*[0] = {m_star[0]:.8f}")

# The floor projection Jacobian has a specific structure.
# It depends on s_pre through R and L2_s_pre.
# Since the dominant singleton s_pre >> w_pre, R is close to 1.

# Key insight: the floor projection is rank-deficient.
# It maps 4D input to the 3D surface {Born=1/27, L₁=1}.
# One eigenvalue should be 0 (the normal to the constraint surface).

print(f"\nFloor Jacobian eigenvalues: {evals_floor}")
print(f"Rank of floor Jacobian: {np.linalg.matrix_rank(J_floor, tol=1e-6)}")

# ============================================================
# Part 6: Can we get the eigenvalue analytically?
# ============================================================
print("\n" + "=" * 60)
print("PART 6: Eigenvalue structure")
print("=" * 60)

# The projected Jacobian on L₁=1 tangent space has:
# λ₀ ≈ 0.2829, λ₁ ≈ 0.2813, λ₂ ≈ 0
#
# The near-degeneracy λ₀ ≈ λ₁ comes from the S₂ symmetry
# (hypotheses 2 and 3 are exchangeable).
#
# λ₂ ≈ 0 corresponds to the floor-killed direction.
#
# Can we find λ₀ analytically?

# The 3×3 projected Jacobian has a block structure from S₂ symmetry:
# symmetric direction: (0, 1, 1, 0)/√2 (same perturbation to both w's)
# antisymmetric direction: (0, 1, -1, 0)/√2 (opposite perturbation)
# radial direction: (1, 0, 0, -1)/√2 (shift between s and θ)

# In the tangent space basis {e₁-e₄, e₂-e₄, e₃-e₄}:
# symmetric = (0, 1, 1)/√2, antisymmetric = (0, 1, -1)/√2, radial = (1, 0, 0)

# The S₂ symmetry means J_proj commutes with the exchange (e₂↔e₃).
# So it block-diagonalizes into:
# - 2×2 block in the S₂-symmetric subspace {radial, symmetric}
# - 1×1 block in the antisymmetric direction

# The antisymmetric eigenvalue IS λ₁.
# The 2×2 block contains λ₀ and λ₂.

J_proj = VTV_inv @ V.T @ J_full @ V

# Transform to the symmetry-adapted basis
P = np.array([
    [1, 0, 0],           # radial: e₁-e₄
    [0, 1/np.sqrt(2), 1/np.sqrt(2)],   # symmetric
    [0, 1/np.sqrt(2), -1/np.sqrt(2)]   # antisymmetric
]).T

J_sym = np.linalg.inv(P) @ J_proj @ P
print("Jacobian in symmetry-adapted basis:")
for i in range(3):
    print(f"  [{J_sym[i,0]:+.8f}  {J_sym[i,1]:+.8f}  {J_sym[i,2]:+.8f}]")

print(f"\n2×2 symmetric block:")
J_2x2 = J_sym[:2,:2]
for i in range(2):
    print(f"  [{J_2x2[i,0]:+.8f}  {J_2x2[i,1]:+.8f}]")
evals_2x2 = np.linalg.eigvals(J_2x2)
print(f"  eigenvalues: {sorted(np.abs(evals_2x2), reverse=True)}")

print(f"\nAntisymmetric eigenvalue: {J_sym[2,2]:.8f}")
print(f"Expected λ₁ = 0.2813")

print(f"\nOff-diagonal blocks (should be ~0 by symmetry):")
print(f"  (0,2): {J_sym[0,2]:.2e}")
print(f"  (1,2): {J_sym[1,2]:.2e}")
print(f"  (2,0): {J_sym[2,0]:.2e}")
print(f"  (2,1): {J_sym[2,1]:.2e}")
