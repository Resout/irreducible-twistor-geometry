"""
Spectral gap on the Born floor surface
=======================================

The Born floor doesn't act — it's the shape of the space.
The dynamics lives on the surface from the start.

On the Born floor surface with S₂ symmetry (s₂=s₃=w):
  θ² = (s² + 2w²)/26        ... Born(θ) = 1/27
  s + 2w + θ = 1             ... L₁ = 1

Two constraints on three variables (s, w, θ) → one free parameter.
Use s as the parameter. Everything is algebraic.
"""

import numpy as np
from sympy import *

# ============================================================
# Part 1: The surface parametrized by s
# ============================================================
print("=" * 60)
print("PART 1: Born floor surface, parametrized by s")
print("=" * 60)

s_sym, w_sym, th_sym = symbols('s w theta', positive=True)

# Constraints
eq1 = s_sym + 2*w_sym + th_sym - 1      # L₁ = 1
eq2 = 26*th_sym**2 - s_sym**2 - 2*w_sym**2  # Born(θ) = 1/27

# Eliminate θ using L₁: θ = 1 - s - 2w
th_expr = 1 - s_sym - 2*w_sym

# Substitute into Born floor equation
born_eq = 26*(1 - s_sym - 2*w_sym)**2 - s_sym**2 - 2*w_sym**2
born_expanded = expand(born_eq)
print(f"Born floor equation (in s, w):")
print(f"  {born_expanded} = 0")

# Solve for w in terms of s
w_solutions = solve(born_eq, w_sym)
print(f"\nSolutions for w(s):")
for i, sol in enumerate(w_solutions):
    print(f"  w_{i} = {sol}")
    print(f"       = {simplify(sol)}")

# The physically meaningful solution (w > 0, small)
# At the fixed point: s ≈ 0.787, w ≈ 0.029
# Let's check which solution gives this
for i, sol in enumerate(w_solutions):
    w_val = float(sol.subs(s_sym, Rational(787, 1000)))
    print(f"  w_{i}(0.787) = {w_val:.6f}")

# Take the solution that gives small positive w
# (the one with the minus sign in the quadratic formula)
w_of_s = w_solutions[0]  # Check which one
w_test = float(w_of_s.subs(s_sym, Rational(787, 1000)))
if w_test > 0.1:  # wrong branch
    w_of_s = w_solutions[1]

print(f"\nPhysical branch: w(s) = {w_of_s}")
print(f"  w(0.787) = {float(w_of_s.subs(s_sym, Rational(787, 1000))):.6f}")

th_of_s = 1 - s_sym - 2*w_of_s
th_of_s_simplified = simplify(th_of_s)
print(f"  θ(s) = {th_of_s_simplified}")
print(f"  θ(0.787) = {float(th_of_s_simplified.subs(s_sym, Rational(787, 1000))):.6f}")

# ============================================================
# Part 2: The DS step restricted to the surface
# ============================================================
print("\n" + "=" * 60)
print("PART 2: DS step on the surface")
print("=" * 60)

# Evidence e* = (a, b, b, φ) with a + 2b + φ = 1
a_sym, b_sym, phi_sym = symbols('a b phi', positive=True)

# DS combination of m=(s,w,w,θ) with e=(a,b,b,φ):
# Pre-normalization:
# s_raw = s*a + s*φ + θ*a = s(a+φ) + θa
# w_raw = w*b + w*φ + θ*b = w(b+φ) + θb
# θ_raw = θ*φ
# K = 2s*b + 2w*(a+b)
# Output = raw / (1-K)

K_expr = 2*s_sym*b_sym + 2*w_sym*(a_sym + b_sym)
s_raw = s_sym*(a_sym + phi_sym) + th_sym*a_sym
w_raw = w_sym*(b_sym + phi_sym) + th_sym*b_sym
th_raw = th_sym*phi_sym

s_out = s_raw / (1 - K_expr)
w_out = w_raw / (1 - K_expr)
th_out = th_raw / (1 - K_expr)

print(f"K = {K_expr}")
print(f"s_out = {s_out}")
print(f"w_out = {w_out}")
print(f"θ_out = {th_out}")

# Now substitute w = w(s) and θ = θ(s) to get everything in terms of s
# The output lives on the same surface (the dynamics preserves the surface,
# after floor enforcement — but the surface IS the constraint, so the
# output must be projected back).
#
# Actually: the DS step does NOT preserve the Born floor surface.
# It pushes θ_raw = θφ/(1-K) which has Born(θ) < 1/27.
# The floor then projects back.
#
# But the floor is the surface. So the effective map on the surface is:
# 1. DS step gives (s_out, w_out, θ_out) OFF the surface
# 2. Project back to surface: find s_new such that the projection of
#    (s_out, w_out, θ_out) onto the surface has parameter s_new

# The projection preserves singleton ratios: s/w is unchanged.
# (The floor rescales all singletons proportionally.)
# So s_new / w_new = s_out / w_out
# And (s_new, w_new, θ_new) is on the surface.

# The singleton ratio after DS step:
ratio_out = simplify(s_out / w_out)
print(f"\ns_out/w_out = {ratio_out}")

# On the surface: s_new/w_new = s_out/w_out = R(s)
# And s_new, w_new satisfy the surface equations.
# With s_new/w_new = R, we have s_new = R * w_new.
# Substituting into the surface equation:
# 26(1 - R*w - 2w)² = R²w² + 2w²
# 26(1 - w(R+2))² = w²(R²+2)
# Let's solve for w_new:

R_sym = symbols('R', positive=True)
w_new_eq = 26*(1 - w_sym*(R_sym + 2))**2 - w_sym**2*(R_sym**2 + 2)
print(f"\nSurface equation with ratio R:")
print(f"  {expand(w_new_eq)} = 0")

w_new_sols = solve(w_new_eq, w_sym)
print(f"\nw_new(R) solutions:")
for sol in w_new_sols:
    print(f"  {simplify(sol)}")

# ============================================================
# Part 3: The 1D map s → s_new
# ============================================================
print("\n" + "=" * 60)
print("PART 3: The effective 1D map on the surface")
print("=" * 60)

# The map is: given s on the surface,
# 1. Compute w(s), θ(s) from the surface equations
# 2. DS step gives s_out, w_out (off surface)
# 3. Ratio R = s_out/w_out is preserved by floor projection
# 4. Find s_new on the surface with ratio R: s_new = R * w_new(R)

# Let's do this numerically first to verify we get the right fixed point
# and eigenvalue.

def surface_state(s_val):
    """Given s, return (s, w, θ) on the Born floor surface."""
    # Solve 26(1-s-2w)² = s² + 2w² for w
    # 26 - 52s - 104w + 26s² + 104sw + 104w² = s² + 2w²
    # 102w² + (104s-104)w + (25s² - 52s + 26) = 0
    A = 102
    B = 104*s_val - 104
    C = 25*s_val**2 - 52*s_val + 26
    disc = B**2 - 4*A*C
    if disc < 0:
        return None
    w1 = (-B + np.sqrt(disc)) / (2*A)
    w2 = (-B - np.sqrt(disc)) / (2*A)
    # Pick the smaller positive one (the physical branch)
    if w2 > 0 and w2 < w1:
        w_val = w2
    else:
        w_val = w1
    th_val = 1 - s_val - 2*w_val
    if th_val < 0:
        return None
    return s_val, w_val, th_val

def ds_step_and_project(s_val, a_val, b_val, phi_val):
    """One DS step on the surface, return new s value."""
    state = surface_state(s_val)
    if state is None:
        return None
    s_v, w_v, th_v = state

    # DS combination
    K = 2*s_v*b_val + 2*w_v*(a_val + b_val)
    s_out = (s_v*(a_val + phi_val) + th_v*a_val) / (1 - K)
    w_out = (w_v*(b_val + phi_val) + th_v*b_val) / (1 - K)

    # Ratio preserved by floor projection
    R = s_out / w_out

    # Find s_new on surface with this ratio
    # s_new = R * w_new, where w_new satisfies surface eq with s = R*w
    # 102w² + (104*R*w - 104)w + (25*R²*w² - 52*R*w + 26) = 0
    # Actually let me substitute s = R*w directly:
    # 26(1 - Rw - 2w)² = R²w² + 2w²
    # 26(1 - w(R+2))² = w²(R²+2)
    # Let u = w(R+2), then w = u/(R+2):
    # 26(1-u)² = u²(R²+2)/(R+2)²
    # 26(R+2)²(1-u)² = u²(R²+2)
    # Let c = (R²+2)/(26(R+2)²)
    # (1-u)² = c*u²
    # 1 - 2u + u² = cu²
    # (1-c)u² - 2u + 1 = 0
    c = (R**2 + 2) / (26 * (R + 2)**2)
    # (1-c)u² - 2u + 1 = 0
    A_q = 1 - c
    B_q = -2
    C_q = 1
    disc = B_q**2 - 4*A_q*C_q
    if disc < 0:
        return None
    u1 = (-B_q + np.sqrt(disc)) / (2*A_q)
    u2 = (-B_q - np.sqrt(disc)) / (2*A_q)
    # u = w(R+2), want 0 < u < 1 (so that θ > 0)
    candidates = [u for u in [u1, u2] if 0 < u < 1]
    if not candidates:
        return None
    u = min(candidates)  # smaller w → larger θ → physical
    w_new = u / (R + 2)
    s_new = R * w_new
    return s_new

# Known evidence at K*=7/30
a_val = 0.6312008879
b_val = 0.1202964949
phi_val = 0.1282061222

# Test: find fixed point by iteration
s_val = 0.5
for i in range(100):
    s_new = ds_step_and_project(s_val, a_val, b_val, phi_val)
    if s_new is None:
        print(f"  Step {i}: failed")
        break
    if abs(s_new - s_val) < 1e-15:
        print(f"  Converged at step {i}")
        break
    s_val = s_new

s_star = s_val
state_star = surface_state(s_star)
print(f"\n  s* = {s_star:.10f}")
print(f"  Full state: s={state_star[0]:.10f}, w={state_star[1]:.10f}, θ={state_star[2]:.10f}")
print(f"  Compare m* from before: s=0.7868984462, w=0.0292822600, θ=0.1545370339")

# Compute derivative of the 1D map at the fixed point
ds = 1e-10
s_plus = ds_step_and_project(s_star + ds, a_val, b_val, phi_val)
s_minus = ds_step_and_project(s_star - ds, a_val, b_val, phi_val)
lambda_1d = (s_plus - s_minus) / (2*ds)

print(f"\n  *** 1D MAP EIGENVALUE ***")
print(f"  dF/ds at s* = {lambda_1d:.10f}")
print(f"  Expected λ₀ = 0.2829103480")
print(f"  Match: {abs(lambda_1d - 0.2829103480) < 1e-6}")

# ============================================================
# Part 4: Analytical derivative of the 1D map
# ============================================================
print("\n" + "=" * 60)
print("PART 4: Analytical derivative")
print("=" * 60)

# The 1D map F(s) = s_new works as:
# 1. (s, w(s), θ(s)) → DS step → (s_out, w_out)
# 2. R = s_out/w_out
# 3. s_new = R · w_new(R)  where w_new comes from the surface equation
#
# dF/ds = d(R·w_new(R))/ds = (dR/ds)·(w_new + R·dw_new/dR)
#       = (dR/ds) · d(R·w_new)/dR
#
# This is a chain: ds → dR → ds_new

# Compute dR/ds numerically
R_star = s_star / state_star[1]  # s/w at fixed point
print(f"  R* = s*/w* = {R_star:.10f}")

# dR/ds
state_plus = surface_state(s_star + ds)
state_minus = surface_state(s_star - ds)
K_plus = 2*(s_star+ds)*b_val + 2*state_plus[1]*(a_val + b_val)
s_out_plus = ((s_star+ds)*(a_val+phi_val) + state_plus[2]*a_val) / (1-K_plus)
w_out_plus = (state_plus[1]*(b_val+phi_val) + state_plus[2]*b_val) / (1-K_plus)
R_plus = s_out_plus / w_out_plus

K_minus = 2*(s_star-ds)*b_val + 2*state_minus[1]*(a_val + b_val)
s_out_minus = ((s_star-ds)*(a_val+phi_val) + state_minus[2]*a_val) / (1-K_minus)
w_out_minus = (state_minus[1]*(b_val+phi_val) + state_minus[2]*b_val) / (1-K_minus)
R_minus = s_out_minus / w_out_minus

dR_ds = (R_plus - R_minus) / (2*ds)
print(f"  dR/ds = {dR_ds:.10f}")

# d(R·w_new)/dR — the surface projection factor
dR = 1e-8
def s_from_R(R):
    c = (R**2 + 2) / (26 * (R + 2)**2)
    A_q = 1 - c
    B_q = -2
    C_q = 1
    disc = B_q**2 - 4*A_q*C_q
    u1 = (-B_q + np.sqrt(disc)) / (2*A_q)
    u2 = (-B_q - np.sqrt(disc)) / (2*A_q)
    candidates = [u for u in [u1, u2] if 0 < u < 1]
    u = min(candidates)
    w = u / (R + 2)
    return R * w

dsn_dR = (s_from_R(R_star + dR) - s_from_R(R_star - dR)) / (2*dR)
print(f"  d(s_new)/dR = {dsn_dR:.10f}")
print(f"  dF/ds = dR/ds · ds_new/dR = {dR_ds * dsn_dR:.10f}")
print(f"  Direct computation: {lambda_1d:.10f}")
print(f"  Match: {abs(dR_ds * dsn_dR - lambda_1d) < 1e-6}")

# ============================================================
# Part 5: Try to get R* and dR/ds in closed form
# ============================================================
print("\n" + "=" * 60)
print("PART 5: Structure of R at fixed point")
print("=" * 60)

# At the fixed point, R = s_out/w_out = s*/w* (because it maps to itself)
# s_out = (s(a+φ) + θa)/(1-K)
# w_out = (w(b+φ) + θb)/(1-K)
# R = (s(a+φ) + θa) / (w(b+φ) + θb)
# At fixed point: R = s/w
# So: s/w = (s(a+φ) + θa) / (w(b+φ) + θb)
# Cross multiply: s(w(b+φ) + θb) = w(s(a+φ) + θa)
# sw(b+φ) + sθb = ws(a+φ) + wθa
# sw(b+φ-a-φ) = wθa - sθb
# sw(b-a) = θ(wa - sb)
# sw(b-a) = -θ(sb - wa)
# sw(b-a) + θ(sb-wa) = 0
# (b-a)(sw) + θ(sb-wa) = 0
# (b-a)sw + θs·b - θw·a = 0
# s(b-a)w + sb·θ - wa·θ = 0
# s·b(w+θ) - a·w(s+θ) = 0
# sb(w+θ) = aw(s+θ)
# b/a = w(s+θ) / (s(w+θ))

s_v, w_v, th_v = state_star
lhs = b_val / a_val
rhs = w_v * (s_v + th_v) / (s_v * (w_v + th_v))
print(f"  Fixed point condition: b/a = w(s+θ)/(s(w+θ))")
print(f"  b/a = {lhs:.10f}")
print(f"  w(s+θ)/(s(w+θ)) = {rhs:.10f}")
print(f"  Match: {abs(lhs - rhs) < 1e-6}")

# This is a clean algebraic relation between evidence and fixed point!
# Combined with:
# - Born floor: 26θ² = s² + 2w²
# - L₁: s + 2w + θ = 1
# - K = 7/30
# - K = 2sb + 2w(a+b)
# - a + 2b + φ = 1
# - θ equation: φ/(1-K) gives θ_out/θ_in

# φ equation from θ fixed point:
# θ = θφ/(1-K)  =>  φ = 1-K = 23/30
# But we saw φ ≈ 0.128, not 23/30 ≈ 0.767!
# Because the floor CHANGES θ_out. The pre-floor θ_out = θφ/(1-K)
# is NOT equal to θ.
# On the surface: the map is s→s, not component-by-component.
# The θ equation doesn't give φ = 1-K anymore.

print(f"\n  φ = {phi_val:.10f}")
print(f"  1-K = {1-7/30:.10f}")
print(f"  φ ≠ 1-K because floor changes θ")

# But φ = 1 - a - 2b
print(f"  a + 2b = {a_val + 2*b_val:.10f}")
print(f"  φ = 1 - (a+2b) = {1 - a_val - 2*b_val:.10f}")

# K = 2sb + 2w(a+b) = 7/30
K_check = 2*s_v*b_val + 2*w_v*(a_val + b_val)
print(f"\n  K = 2sb + 2w(a+b) = {K_check:.10f}")
print(f"  7/30 = {7/30:.10f}")

# How many unknowns, how many equations?
# Unknowns: s, w, θ, a, b, φ (6)
# Equations:
# 1. s + 2w + θ = 1
# 2. a + 2b + φ = 1
# 3. 26θ² = s² + 2w²
# 4. K = 2sb + 2w(a+b) = 7/30
# 5. Fixed point: b/a = w(s+θ)/(s(w+θ))
# 6. Evidence also on floor?: 26φ² = a² + 2b²  (maybe)

born_e = phi_val**2 / (a_val**2 + 2*b_val**2 + phi_val**2)
print(f"\n  Born(φ) of evidence = {born_e:.10f}")
print(f"  1/27 = {1/27:.10f}")
print(f"  Evidence also on floor: {abs(born_e - 1/27) < 0.001}")
# Evidence is NOT necessarily on the floor.

# So we have 6 unknowns and 5 equations. One degree of freedom.
# That degree of freedom is the evidence strength — how concentrated
# the evidence is. K*=7/30 constrains it but doesn't fix it uniquely
# (different concentrations can give the same K*).
#
# Or does the self-consistency fully determine it?
# The evidence e* is determined by the gauge field at equilibrium.
# In the DS framework, e* is whatever produces K*=7/30.
# Given the symmetric structure (two weak hypotheses),
# the equations might close.

# Actually: there's another equation. The evidence must also
# satisfy its own Born floor. Let me check:
print(f"\n  Evidence Born(φ) = {born_e:.6f}")
print(f"  This is {'≥' if born_e >= 1/27 else '<'} 1/27")
# Evidence CAN have Born(φ) > 1/27. It doesn't have to be on the floor.
