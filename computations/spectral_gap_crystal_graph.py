"""
Build a crystal graph from the six fixed-point equations.

Variables: s, w, θ, a, b, φ
Equations:
  1. s + 2w + θ = 1
  2. a + 2b + φ = 1
  3. 26θ² = s² + 2w²
  4. 26φ² = a² + 2b²
  5. K = 2sb + 2w(a+b) = 7/30
  6. b/a = w(s+θ) / (s(w+θ))

Every equation is a relationship between variables.
Every relationship is a crystal.
Build the graph, compose, and see what's there.
"""

import numpy as np
from itertools import combinations

# ============================================================
# The six equations as relationships
# ============================================================

# Equation 1: s + 2w + θ = 1
# Relates s, w, θ. Budget constraint for state.
# If s goes up, w+θ must go down. Anticorrelated.

# Equation 2: a + 2b + φ = 1
# Same structure for evidence.

# Equation 3: 26θ² = s² + 2w²
# Born floor. θ is geometrically determined by s and w.
# θ = √((s² + 2w²)/26)

# Equation 4: 26φ² = a² + 2b²
# Same for evidence.

# Equation 5: 2sb + 2w(a+b) = 7/30
# The coupling equation. Mixes state and evidence variables.

# Equation 6: b/a = w(s+θ) / (s(w+θ))
# The fixed-point condition. The deepest coupling.

# ============================================================
# Solve the system numerically to high precision
# ============================================================

from scipy.optimize import fsolve

def system(x):
    s, w, th, a, b, phi = x
    return [
        s + 2*w + th - 1,                          # eq 1
        a + 2*b + phi - 1,                          # eq 2
        26*th**2 - s**2 - 2*w**2,                   # eq 3
        26*phi**2 - a**2 - 2*b**2,                  # eq 4
        2*s*b + 2*w*(a+b) - 7/30,                   # eq 5
        b*s*(w+th) - a*w*(s+th),                    # eq 6
    ]

# Initial guess from known values
x0 = [0.787, 0.029, 0.155, 0.631, 0.120, 0.128]
sol = fsolve(system, x0, full_output=True)
x_sol = sol[0]
s, w, th, a, b, phi = x_sol

print("=" * 60)
print("SOLUTION")
print("=" * 60)
print(f"s   = {s:.12f}")
print(f"w   = {w:.12f}")
print(f"θ   = {th:.12f}")
print(f"a   = {a:.12f}")
print(f"b   = {b:.12f}")
print(f"φ   = {phi:.12f}")
print(f"\nResiduals: {[f'{r:.2e}' for r in system(x_sol)]}")

# Verify
print(f"\nChecks:")
print(f"  s + 2w + θ = {s + 2*w + th:.12f}")
print(f"  a + 2b + φ = {a + 2*b + phi:.12f}")
print(f"  26θ² - s² - 2w² = {26*th**2 - s**2 - 2*w**2:.2e}")
print(f"  26φ² - a² - 2b² = {26*phi**2 - a**2 - 2*b**2:.2e}")
print(f"  K = 2sb + 2w(a+b) = {2*s*b + 2*w*(a+b):.12f}")
print(f"  b/a = {b/a:.12f}")
print(f"  w(s+θ)/(s(w+θ)) = {w*(s+th)/(s*(w+th)):.12f}")

# ============================================================
# All pairwise relationships
# ============================================================
print("\n" + "=" * 60)
print("PAIRWISE RELATIONSHIPS")
print("=" * 60)

vars_dict = {'s': s, 'w': w, 'θ': th, 'a': a, 'b': b, 'φ': phi}
var_names = list(vars_dict.keys())

# For each pair, which equations connect them?
connections = {
    ('s', 'w'):  ["eq1: s+2w+θ=1", "eq3: 26θ²=s²+2w²", "eq5: K coupling", "eq6: fixed point"],
    ('s', 'θ'):  ["eq1: s+2w+θ=1", "eq3: 26θ²=s²+2w²", "eq6: fixed point"],
    ('w', 'θ'):  ["eq1: s+2w+θ=1", "eq3: 26θ²=s²+2w²", "eq6: fixed point"],
    ('a', 'b'):  ["eq2: a+2b+φ=1", "eq4: 26φ²=a²+2b²", "eq5: K coupling", "eq6: fixed point"],
    ('a', 'φ'):  ["eq2: a+2b+φ=1", "eq4: 26φ²=a²+2b²"],
    ('b', 'φ'):  ["eq2: a+2b+φ=1", "eq4: 26φ²=a²+2b²"],
    ('s', 'a'):  ["eq5: K coupling", "eq6: fixed point"],
    ('s', 'b'):  ["eq5: K coupling", "eq6: fixed point"],
    ('w', 'a'):  ["eq5: K coupling", "eq6: fixed point"],
    ('w', 'b'):  ["eq5: K coupling", "eq6: fixed point"],
    ('s', 'φ'):  ["indirect via a,b"],
    ('w', 'φ'):  ["indirect via a,b"],
    ('θ', 'a'):  ["eq6: fixed point"],
    ('θ', 'b'):  ["indirect via s,w,a"],
    ('θ', 'φ'):  ["indirect"],
}

for pair in combinations(var_names, 2):
    key = pair if pair in connections else (pair[1], pair[0])
    conns = connections.get(key, connections.get((pair[1],pair[0]), ["none direct"]))
    print(f"  {pair[0]:2s} ↔ {pair[1]:2s}: {', '.join(conns)}")

# ============================================================
# Build 3×3 correlation matrices for each equation
# ============================================================
print("\n" + "=" * 60)
print("CORRELATION MATRICES (each equation)")
print("=" * 60)

# For each equation, perturb variables and see how outputs co-vary.
# This gives the correlation structure = the crystal.

def numerical_correlation(f, var_indices, x_base, eps=1e-6):
    """Compute correlation matrix between two variables under constraint f=0.

    var_indices: (i, j) — the two variable indices to correlate
    Returns 3×3 matrix: how perturbation of var_i affects var_j
    """
    i, j = var_indices
    n = len(x_base)

    # Sample perturbations of variable i, see response of variable j
    # Use implicit function theorem: df/dx_j * dx_j = -df/dx_i * dx_i

    # Jacobian of f at x_base
    J = np.zeros(n)
    f0 = f(x_base)
    for k in range(n):
        xp = x_base.copy()
        xp[k] += eps
        J[k] = (f(xp) - f0) / eps

    # Response: dx_j/dx_i = -J[i]/J[j] (implicit function theorem)
    if abs(J[j]) > 1e-15:
        response = -J[i] / J[j]
    else:
        response = 0

    return response

# The six constraint functions
def f1(x): return x[0] + 2*x[1] + x[2] - 1  # s + 2w + θ = 1
def f2(x): return x[3] + 2*x[4] + x[5] - 1  # a + 2b + φ = 1
def f3(x): return 26*x[2]**2 - x[0]**2 - 2*x[1]**2  # Born floor state
def f4(x): return 26*x[5]**2 - x[3]**2 - 2*x[4]**2  # Born floor evidence
def f5(x): return 2*x[0]*x[4] + 2*x[1]*(x[3]+x[4]) - 7/30  # K = 7/30
def f6(x): return x[4]*x[0]*(x[1]+x[2]) - x[3]*x[1]*(x[0]+x[2])  # fixed pt

# For each equation and each variable pair IN that equation,
# compute the response (how much j moves when i moves, on-constraint)

print("Response matrix (dx_j/dx_i on each constraint surface):\n")

equations = [
    ("eq1: s+2w+θ=1", f1, [0,1,2]),
    ("eq2: a+2b+φ=1", f2, [3,4,5]),
    ("eq3: Born state", f3, [0,1,2]),
    ("eq4: Born evidence", f4, [3,4,5]),
    ("eq5: K=7/30", f5, [0,1,3,4]),
    ("eq6: fixed point", f6, [0,1,2,3,4]),
]

idx_to_name = {0:'s', 1:'w', 2:'θ', 3:'a', 4:'b', 5:'φ'}

for name, f, var_idxs in equations:
    print(f"  {name}:")
    for i in var_idxs:
        for j in var_idxs:
            if i != j:
                r = numerical_correlation(f, (i,j), x_sol)
                print(f"    d{idx_to_name[j]}/d{idx_to_name[i]} = {r:+.6f}")
    print()

# ============================================================
# The Jacobian of the full system
# ============================================================
print("=" * 60)
print("FULL SYSTEM JACOBIAN")
print("=" * 60)

eps = 1e-9
n = 6
J_sys = np.zeros((n, n))
f_base = np.array(system(x_sol))
for j in range(n):
    xp = x_sol.copy()
    xp[j] += eps
    J_sys[:, j] = (np.array(system(xp)) - f_base) / eps

print("Jacobian of the constraint system:")
names = ['s', 'w', 'θ', 'a', 'b', 'φ']
print(f"{'':>4s}", end="")
for n in names:
    print(f"  {n:>10s}", end="")
print()
for i in range(6):
    print(f"{names[i]:>4s}", end="")
    for j in range(6):
        print(f"  {J_sys[i,j]:+10.5f}", end="")
    print()

evals_sys = np.linalg.eigvals(J_sys)
print(f"\nEigenvalues of constraint Jacobian: {np.sort(np.abs(evals_sys))[::-1].round(4)}")

# The null space of the Jacobian = the direction along the solution manifold
# If the system is fully determined (6 eq, 6 unkn), the Jacobian should be
# full rank at the solution.
rank = np.linalg.matrix_rank(J_sys, tol=1e-6)
print(f"Rank: {rank} (should be 6 for isolated solution)")

# ============================================================
# Compose: state-side crystal ∘ coupling crystal ∘ evidence-side crystal
# ============================================================
print("\n" + "=" * 60)
print("COMPOSITION: state ↔ coupling ↔ evidence")
print("=" * 60)

# The state-side (eqs 1,3): s,w,θ constrained to Born floor surface
# The evidence-side (eqs 2,4): a,b,φ constrained to Born floor surface
# The coupling (eqs 5,6): connects state to evidence
#
# On each surface, there's one degree of freedom.
# The coupling maps one to the other.
# The COMPOSED map: s → (via coupling) → a → (via evidence surface) → back to s?
# No — the map is s → s_next via the DS dynamics.
# The DS dynamics IS the composition of these crystals.

# Let me think about this differently.
# The transfer operator T takes m to m_next.
# m lives on the state surface (1 DOF, parametrized by s).
# The evidence e lives on the evidence surface (1 DOF, parametrized by a).
# The coupling fixes a in terms of s (via eqs 5,6).
# So T: s → s is a 1D map.

# What's dT/ds? It's the chain:
# ds → (surface) → dw, dθ → (DS combination) → ds_out, dw_out → (floor) → ds_new

# But the surface, DS, and floor are all encoded in the six equations.
# The implicit function theorem on the full system gives us all the
# partial derivatives we need.

# Specifically: at the fixed point, the system F(s,w,θ,a,b,φ)=0 determines
# all six variables. If we treat one as free (say s), the others are
# functions of s. But at the fixed point, s is also determined — there's
# no free parameter. The 1D map eigenvalue comes from the DYNAMICS, not
# the constraint.

# The dynamics says: given current (s,w,θ) on the state surface,
# and evidence (a,b,φ) on the evidence surface,
# the NEXT state is DS(m,e) projected to the state surface.
# At the fixed point, m_next = m. The eigenvalue is d(m_next)/d(m).

# This is already what we computed: λ₀ = 0.2829.
# The question is whether the six equations encode it algebraically.

# The constraint Jacobian J_sys has the structure to answer this.
# If we perturb s and ask "how much does s_next change?", the answer
# flows through all six equations simultaneously.

# ============================================================
# Direct: the 1D eigenvalue from the constraint Jacobian
# ============================================================
print("\n" + "=" * 60)
print("THE EIGENVALUE FROM THE CONSTRAINT STRUCTURE")
print("=" * 60)

# The six equations at the NEXT time step would be:
# s' + 2w' + θ' = 1
# (a,b,φ same — evidence is fixed)
# 26θ'² = s'² + 2w'²
# DS: s' = (s(a+φ) + θa)/(1-K), etc, then floor project
#
# But we already know the 1D derivative is 0.2829.
# Can we express it in terms of the ratios in the system?

# Key ratios at the fixed point:
print("Key ratios:")
r = s/w
print(f"  s/w = {r:.8f}")
print(f"  s/θ = {s/th:.8f}")
print(f"  a/b = {a/b:.8f}")
print(f"  a/φ = {a/phi:.8f}")
print(f"  s/a = {s/a:.8f}")
print(f"  w/b = {w/b:.8f}")
print(f"  θ/φ = {th/phi:.8f}")

# The transfer map derivative: ds_new/ds = 0.2829
# From the chain rule through the constraint system:
# The DS step gives ds_out/ds, and the floor projection gives ds_new/ds_out.
# Their product is 0.2829.

# DS step (before floor): s_out = (s(a+φ) + θa)/(1-K)
# At fixed point, s_out is a function of s (through the state surface).
# ds_out/ds = [(a+φ) + (dθ/ds)·a] / (1-K)  +  s_out · (dK/ds)/(1-K)

# On the state surface (eqs 1,3): θ and w are functions of s.
# dθ/ds from implicit differentiation of eqs 1,3:

# eq1: 1 + 2·dw/ds + dθ/ds = 0
# eq3: 52θ·dθ/ds = 2s + 4w·dw/ds

# From eq1: dθ/ds = -1 - 2·dw/ds
# Substitute into eq3: 52θ(-1-2dw/ds) = 2s + 4w·dw/ds
# -52θ - 104θ·dw/ds = 2s + 4w·dw/ds
# -52θ - 2s = (4w + 104θ)·dw/ds
# dw/ds = (-52θ - 2s) / (4w + 104θ)

dw_ds = (-52*th - 2*s) / (4*w + 104*th)
dth_ds = -1 - 2*dw_ds

print(f"\nOn-surface derivatives (state):")
print(f"  dw/ds  = {dw_ds:.10f}")
print(f"  dθ/ds  = {dth_ds:.10f}")

# Similarly for evidence surface:
# eq2: da + 2db + dφ = 0 → dφ/da = -1 - 2·db/da
# eq4: 52φ·dφ/da = 2a + 4b·db/da
# db/da = (-52φ - 2a) / (4b + 104φ)

db_da = (-52*phi - 2*a) / (4*b + 104*phi)
dphi_da = -1 - 2*db_da

print(f"\nOn-surface derivatives (evidence):")
print(f"  db/da  = {db_da:.10f}")
print(f"  dφ/da  = {dphi_da:.10f}")

# Now: K = 2sb + 2w(a+b) = 7/30
# At the fixed point, evidence is fixed (a,b,φ don't change with the dynamics).
# So dK/ds = 2b + 2(dw/ds)(a+b)

dK_ds = 2*b + 2*dw_ds*(a+b)
print(f"\n  dK/ds  = {dK_ds:.10f}")

# DS pre-floor output derivatives:
# s_out = (s(a+φ) + θa) / (1-K)
# At fixed K=7/30:
# ds_out/ds = [(a+φ) + dθ/ds · a] / (1-K)

# Wait — K depends on s through w, so 1-K also varies.
# ds_out/ds = d/ds [(s(a+φ) + θa) / (1-K)]
# = [(a+φ + dθ/ds·a)(1-K) + (s(a+φ)+θa)·dK/ds] / (1-K)²

numerator_ds = (a + phi + dth_ds*a) * (1 - 7/30) + (s*(a+phi) + th*a) * dK_ds
ds_out_ds = numerator_ds / (1 - 7/30)**2

# Similarly for w_out:
numerator_dw = (dw_ds*(b+phi) + dth_ds*b) * (1-7/30) + (w*(b+phi)+th*b) * dK_ds
dw_out_ds = numerator_dw / (1-7/30)**2

print(f"\nDS output derivatives (pre-floor):")
print(f"  ds_out/ds = {ds_out_ds:.10f}")
print(f"  dw_out/ds = {dw_out_ds:.10f}")

# The ratio R = s_out/w_out. Its derivative:
R_val = s/w  # at fixed point, s_out/w_out = s/w
dR_ds = (ds_out_ds * w - s * dw_out_ds) / w**2
# Wait, R = s_out/w_out, not s/w. But at fixed point after floor, s_new/w_new = s_out/w_out.
# And s_new = s, w_new = w. So s_out/w_out = s/w? No — floor changes the values.
# s_out ≠ s. s_out is the pre-floor value. After floor: s_new = s.
# The floor preserves the ratio: s_new/w_new = s_out/w_out.
# So at fixed point: s_out/w_out = s/w.

# Actually let me compute s_out and w_out:
s_out_val = (s*(a+phi) + th*a) / (1 - 7/30)
w_out_val = (w*(b+phi) + th*b) / (1 - 7/30)
R_out = s_out_val / w_out_val
print(f"\n  s_out = {s_out_val:.10f}")
print(f"  w_out = {w_out_val:.10f}")
print(f"  R_out = s_out/w_out = {R_out:.10f}")
print(f"  s/w = {s/w:.10f}")
print(f"  Match: {abs(R_out - s/w) < 1e-6}")

# Good. Now dR/ds:
dR_ds_val = (ds_out_ds * w_out_val - s_out_val * dw_out_ds) / w_out_val**2
print(f"  dR/ds = {dR_ds_val:.10f}")

# Floor projection: given R, find s_new on the state surface.
# s_new = R · w_new(R)
# where w_new satisfies 26(1-w(R+2))² = w²(R²+2)
# Let c = (R²+2)/(26(R+2)²)
# (1-c)u² - 2u + 1 = 0 where u = w(R+2)
# u = [2 - √(4-4(1-c))] / (2(1-c)) = [1-√c]/(1-c) = 1/(1+√c)
# (taking the smaller root for physical branch)

# So u = 1/(1+√c) where c = (R²+2)/(26(R+2)²)
# w = u/(R+2) = 1/((1+√c)(R+2))
# s_new = R·w = R/((1+√c)(R+2))

c_val = (R_val**2 + 2) / (26 * (R_val + 2)**2)
sqrt_c = np.sqrt(c_val)
s_new_formula = R_val / ((1 + sqrt_c) * (R_val + 2))
print(f"\n  c = (R²+2)/(26(R+2)²) = {c_val:.10f}")
print(f"  √c = {sqrt_c:.10f}")
print(f"  s_new = R/((1+√c)(R+2)) = {s_new_formula:.10f}")
print(f"  s = {s:.10f}")
print(f"  Match: {abs(s_new_formula - s) < 1e-6}")

# ds_new/dR:
# s_new = R / ((1+√c)(R+2))
# This is a function of R alone. Let's differentiate.
# dc/dR = d/dR [(R²+2)/(26(R+2)²)]
#       = [2R·26(R+2)² - (R²+2)·26·2(R+2)] / [26(R+2)²]²
#       = [2R(R+2) - 2(R²+2)] / [26(R+2)³]
#       = [2R²+4R - 2R²-4] / [26(R+2)³]
#       = [4R-4] / [26(R+2)³]
#       = 4(R-1) / [26(R+2)³]

dc_dR = 4*(R_val - 1) / (26 * (R_val + 2)**3)

# d(√c)/dR = dc/dR / (2√c)
dsqrtc_dR = dc_dR / (2 * sqrt_c)

# s_new = R / ((1+√c)(R+2)) = R / g(R) where g = (1+√c)(R+2)
# ds_new/dR = [g - R·dg/dR] / g²
# dg/dR = d√c/dR · (R+2) + (1+√c)

dg_dR = dsqrtc_dR * (R_val + 2) + (1 + sqrt_c)
g_val = (1 + sqrt_c) * (R_val + 2)
ds_new_dR = (g_val - R_val * dg_dR) / g_val**2

print(f"\n  ds_new/dR = {ds_new_dR:.10f}")

# THE EIGENVALUE:
lambda_0 = dR_ds_val * ds_new_dR
print(f"\n  *** λ₀ = dR/ds · ds_new/dR = {lambda_0:.10f} ***")
print(f"  Expected: 0.2829103480")
print(f"  Match: {abs(lambda_0 - 0.2829103) < 1e-4}")

# ============================================================
# Now: is this expressible in closed form?
# ============================================================
print("\n" + "=" * 60)
print("CLOSED FORM ANALYSIS")
print("=" * 60)

# λ₀ = (dR/ds) · (ds_new/dR)
#
# dR/ds comes from the DS combination (polynomial) differentiated
# along the Born floor surface (algebraic).
#
# ds_new/dR comes from the floor projection formula
# s_new = R/((1+√c)(R+2)) where c = (R²+2)/(26(R+2)²)
#
# Both pieces are algebraic functions of the fixed-point values.
# The fixed-point values are roots of a polynomial system.
#
# So λ₀ is an algebraic number. It's a root of some polynomial.
# The question is: which one?

# Let me compute the minimal polynomial numerically.
# If λ₀ is a root of a low-degree polynomial with small integer
# coefficients, we can find it.

from numpy.polynomial import polynomial as P

# Try to find integer relation: a₀ + a₁λ + a₂λ² + ... = 0
lam = lambda_0
powers = [lam**k for k in range(20)]

# Use numpy lstsq to find approximate integer coefficients
for degree in range(2, 15):
    A = np.array(powers[:degree+1]).reshape(1, -1)
    # This is underdetermined for one point. Need more structure.
    # Instead, compute many powers and look for linear dependence.
    pass

# Simpler: just print the value to many digits and see if we recognize it
print(f"λ₀ = {lambda_0:.15f}")
print(f"λ₀² = {lambda_0**2:.15f}")
print(f"1/λ₀ = {1/lambda_0:.15f}")
print(f"Δ = -ln(λ₀) = {-np.log(lambda_0):.15f}")

# Check if Δ has a nice form
Delta = -np.log(lambda_0)
print(f"\nΔ/ln(2) = {Delta/np.log(2):.10f}")
print(f"Δ/ln(3) = {Delta/np.log(3):.10f}")
print(f"Δ·30/7 = {Delta*30/7:.10f}")
print(f"Δ·23/7 = {Delta*23/7:.10f}")
print(f"Δ·30/23 = {Delta*30/23:.10f}")
print(f"e^Δ = 1/λ₀ = {np.exp(Delta):.10f}")
print(f"e^Δ - 1 = {np.exp(Delta)-1:.10f}")
print(f"30/23 = {30/23:.10f}")
print(f"(30/23)^(30/23) = {(30/23)**(30/23):.10f}")
