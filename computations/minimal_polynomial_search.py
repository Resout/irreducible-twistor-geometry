"""
Minimal polynomial search for lambda_0.

Strategy:
- Use the symmetric fixed-point structure (s2=s3) to reduce to fewer variables
- Solve the fixed-point equations directly using mpmath.findroot (Newton's method)
  starting from good initial values (known to 15 digits)
- Then compute Jacobian analytically, get characteristic polynomial, run PSLQ
"""

from mpmath import mp, mpf, sqrt, log, nstr, matrix, eye, det, fsum, pi
import mpmath
import sys

mp.dps = 300  # 300 digits for headroom

H = mpf(3)
FLOOR = mpf(1) / mpf(27)
KSTAR = mpf(7) / mpf(30)

print("="*70)
print("STEP 1: Solve fixed-point equations via Newton's method")
print("="*70)
sys.stdout.flush()

# At symmetric equilibrium: m* = (a, b, b, theta), e* = (w, v, v, phi)
# with s2=s3, e2=e3
#
# Constraints:
# (1) a + 2b + theta = 1
# (2) w + 2v + phi = 1
# (3) 26*theta^2 = a^2 + 2*b^2   [Born floor state]
# (4) 26*phi^2 = w^2 + 2*v^2     [Born floor evidence]
# (5) K = 2*a*v + 2*b*w + 2*b*v = 7/30
# (6) Fixed-point: Phi(m*, e*) = m*
#
# From (1): theta = 1 - a - 2b
# From (2): phi = 1 - w - 2v
# From (3): 26*(1-a-2b)^2 = a^2 + 2b^2
# From (4): 26*(1-w-2v)^2 = w^2 + 2v^2
#
# DS step outputs (before floor):
#   a_ds = (a*w + a*phi + theta*w) / (1-K)
#   b_ds = (b*v + b*phi + theta*v) / (1-K)
#   theta_ds = theta*phi / (1-K)
#
# After floor: need to solve floor quadratic.
# Fixed point: floor(DS(m*,e*)) = m*
# This means: a = alpha * a_ds, b = alpha * b_ds, theta = t(R_ds)
# where R_ds = (a_ds^2 + 2*b_ds^2)/(a_ds + 2*b_ds)^2
# and alpha = (1-t)/(a_ds + 2*b_ds)
#
# The key relation from fixed-point: the RATIO a/b = alpha*a_ds / (alpha*b_ds) = a_ds/b_ds
# So: a/b = (a*w + a*phi + theta*w) / (b*v + b*phi + theta*v)
# This gives one equation.
#
# Also: theta = t(R_ds) where R_ds depends on a_ds, b_ds

# Let's use 4 free variables: a, b, v, w
# (theta, phi determined by (1),(2))
# Equations:
# (3): 26*(1-a-2b)^2 = a^2 + 2b^2
# (4): 26*(1-w-2v)^2 = w^2 + 2v^2
# (5): K = 2av + 2bw + 2bv = 7/30
# (6): a/b = a_ds/b_ds  (ratio preservation at fixed point)

# Actually, let me just set up 4 equations in 4 unknowns and use findroot.
# Variables: a, b, w, v

def equations(a, b, w, v):
    theta = 1 - a - 2*b
    phi = 1 - w - 2*v

    # Eq 1: Born floor state
    eq1 = 26*theta**2 - a**2 - 2*b**2

    # Eq 2: Born floor evidence
    eq2 = 26*phi**2 - w**2 - 2*v**2

    # Eq 3: K = 7/30
    K = 2*a*v + 2*b*w + 2*b*v
    eq3 = K - KSTAR

    # Eq 4: ratio preservation (fixed-point condition)
    # a_ds/b_ds = a/b
    a_ds_num = a*w + a*phi + theta*w  # numerator (before 1/(1-K))
    b_ds_num = b*v + b*phi + theta*v
    # a_ds/b_ds = a/b means a_ds_num * b = b_ds_num * a
    eq4 = a_ds_num * b - b_ds_num * a

    return (eq1, eq2, eq3, eq4)

# Good starting point from double precision
a0 = mpf('0.78689845')
b0 = mpf('0.02928226')
w0 = mpf('0.63120089')
v0 = mpf('0.12029649')

print("Solving with Newton's method (mpmath.findroot)...")
sys.stdout.flush()

result = mpmath.findroot(equations, (a0, b0, w0, v0), tol=mpf(10)**(-mp.dps+20))
a, b, w, v = result
theta = 1 - a - 2*b
phi = 1 - w - 2*v

m_star = [a, b, b, theta]
e_star = [w, v, v, phi]

# Verify
K_check = 2*a*v + 2*b*w + 2*b*v
born_m = theta**2 / (a**2 + 2*b**2 + theta**2)
born_e = phi**2 / (w**2 + 2*v**2 + phi**2)

print(f"K* residual: {nstr(K_check - KSTAR, 15)}")
print(f"Born(m*) = {nstr(born_m, 30)}  (should be 1/27 = {nstr(FLOOR, 30)})")
print(f"Born(e*) = {nstr(born_e, 30)}")
print(f"Born floor eq residual (state): {nstr(26*theta**2 - a**2 - 2*b**2, 15)}")
print(f"Born floor eq residual (evid):  {nstr(26*phi**2 - w**2 - 2*v**2, 15)}")

# Check ratio preservation
a_ds_num = a*w + a*phi + theta*w
b_ds_num = b*v + b*phi + theta*v
print(f"Ratio preservation: {nstr(a_ds_num*b - b_ds_num*a, 15)}")

print(f"\na  = {nstr(a, 60)}")
print(f"b  = {nstr(b, 60)}")
print(f"th = {nstr(theta, 60)}")
print(f"w  = {nstr(w, 60)}")
print(f"v  = {nstr(v, 60)}")
print(f"ph = {nstr(phi, 60)}")
sys.stdout.flush()

print()
print("="*70)
print("STEP 2: Analytical Jacobian of Phi = floor o DS")
print("="*70)
sys.stdout.flush()

# DS step (using full 4-component vectors)
def ds_step(m, e):
    s1, s2, s3, th = m
    e1, e2, e3, ph = e
    sn1 = s1*e1 + s1*ph + th*e1
    sn2 = s2*e2 + s2*ph + th*e2
    sn3 = s3*e3 + s3*ph + th*e3
    tn = th*ph
    K = s1*(e2+e3) + s2*(e1+e3) + s3*(e1+e2)
    d = 1 - K
    return [sn1/d, sn2/d, sn3/d, tn/d], K

m_ds, K_val = ds_step(m_star, e_star)
print(f"DS output: [{nstr(m_ds[0],30)}, {nstr(m_ds[1],30)}, {nstr(m_ds[2],30)}, {nstr(m_ds[3],30)}]")
print(f"K = {nstr(K_val, 30)}")
born_ds = m_ds[3]**2 / fsum([x**2 for x in m_ds])
print(f"Born at DS output: {nstr(born_ds, 15)} (< 1/27 = {nstr(FLOOR, 15)}: floor fires)")
sys.stdout.flush()

# 4x4 DS Jacobian
ee = [w, v, v]
ph = phi
ss = [a, b, b]
th = theta
K = K_val
d = 1 - K

dKds = [(1 - ph) - ee[i] for i in range(3)]
N = [ss[i]*(ee[i]+ph) + th*ee[i] for i in range(3)]
N_th = th*ph

J_ds = matrix(4, 4)
for i in range(3):
    for j in range(3):
        dNi_dsj = (ee[i] + ph) if i == j else mpf(0)
        J_ds[i,j] = (dNi_dsj * d + N[i] * dKds[j]) / d**2
    J_ds[i,3] = (ee[i] * d + N[i] * 0) / d**2  # dK/dtheta = 0

for j in range(3):
    J_ds[3,j] = (0 * d + N_th * dKds[j]) / d**2
J_ds[3,3] = (ph * d) / d**2

# Floor Jacobian
S_ds = m_ds[0] + m_ds[1] + m_ds[2]
Q_ds = m_ds[0]**2 + m_ds[1]**2 + m_ds[2]**2
R = Q_ds / S_ds**2
sqrtR26 = sqrt(26 * R)
t = (sqrtR26 - R) / (26 - R)
alpha_fl = (1 - t) / S_ds

# dt/dR
f = sqrtR26 - R
fp = 13/sqrtR26 - 1
g = 26 - R
dt_dR = (fp * g + f) / g**2

# dR/ds_i^ds
s_ds = [m_ds[0], m_ds[1], m_ds[2]]
dRds_ds = [2*(s_ds[i]/S_ds - R)/S_ds for i in range(3)]

dt_ds_ds = [dt_dR * dRds_ds[i] for i in range(3)]
dalpha_ds_ds = [(-dt_ds_ds[i])/S_ds - alpha_fl/S_ds for i in range(3)]

J_fl = matrix(4, 4)
for i in range(3):
    for j in range(3):
        J_fl[i,j] = (alpha_fl if i==j else mpf(0)) + s_ds[i]*dalpha_ds_ds[j]
    J_fl[i,3] = mpf(0)

for j in range(3):
    J_fl[3,j] = dt_ds_ds[j]
J_fl[3,3] = mpf(0)

# Total Jacobian
J_total = J_fl * J_ds

# Project to L1=1 tangent space
J_proj = matrix(3, 3)
for i in range(3):
    for j in range(3):
        J_proj[i,j] = J_total[i,j] - J_total[i,3]

print(f"\nProjected 3x3 Jacobian:")
for i in range(3):
    print(f"  [{nstr(J_proj[i,0],30)}, {nstr(J_proj[i,1],30)}, {nstr(J_proj[i,2],30)}]")

# Characteristic polynomial
tr_J = J_proj[0,0] + J_proj[1,1] + J_proj[2,2]
cof = (J_proj[0,0]*J_proj[1,1] - J_proj[0,1]*J_proj[1,0]) + \
      (J_proj[0,0]*J_proj[2,2] - J_proj[0,2]*J_proj[2,0]) + \
      (J_proj[1,1]*J_proj[2,2] - J_proj[1,2]*J_proj[2,1])
det_J = det(J_proj)

print(f"\ntr(J_proj) = {nstr(tr_J, 70)}")
print(f"cofactor   = {nstr(cof, 70)}")
print(f"det(J_proj)= {nstr(det_J, 15)} (should be ~0)")

# Eigenvalues
disc = tr_J**2 - 4*cof
lam0 = (tr_J + sqrt(disc)) / 2
lam1 = (tr_J - sqrt(disc)) / 2

print(f"\nlambda_0 = {nstr(lam0, 70)}")
print(f"lambda_1 = {nstr(lam1, 70)}")
print(f"Delta    = {nstr(-log(lam0), 70)}")

det_IJ = (1-lam0)*(1-lam1)
print(f"det(I-J) = {nstr(det_IJ, 70)}")
sys.stdout.flush()

# Cross-check
lam0_ref = mpf('0.28291034660315143666181958723824224253434')
print(f"\nlam0 cross-check error: {nstr(abs(lam0 - lam0_ref), 15)}")

print()
print("="*70)
print("STEP 3: PSLQ searches")
print("="*70)
sys.stdout.flush()

def pslq_search(x, max_degree, label, maxcoeff=10**12):
    for deg in range(2, max_degree + 1):
        powers = [x**k for k in range(deg + 1)]
        try:
            rel = mpmath.pslq(powers, maxcoeff=maxcoeff, maxsteps=10000)
            if rel is not None:
                val = fsum([mpf(rel[k]) * x**k for k in range(deg + 1)])
                if abs(val) < mpf(10)**(-100):
                    print(f"  {label}: FOUND degree {deg}")
                    print(f"  Coefficients: {rel}")
                    print(f"  Residual: {nstr(abs(val), 15)}")
                    sys.stdout.flush()
                    return rel, deg
        except:
            pass
        if deg % 10 == 0:
            print(f"  {label}: searched up to degree {deg}...")
            sys.stdout.flush()
    print(f"  {label}: NOT FOUND up to degree {max_degree}")
    sys.stdout.flush()
    return None, None

print("\n--- tr(J_proj) ---")
pslq_search(tr_J, 80, "tr(J)")

print("\n--- cofactor sum ---")
pslq_search(cof, 80, "cofactor")

print("\n--- discriminant tr^2 - 4*cof ---")
pslq_search(disc, 80, "disc")

print("\n--- sqrt(disc) ---")
pslq_search(sqrt(disc), 80, "sqrt(disc)")

print("\n--- lambda_0 ---")
pslq_search(lam0, 100, "lambda_0")

print("\n--- lambda_0^2 ---")
pslq_search(lam0**2, 80, "lambda_0^2")

print("\n--- det(I-J) ---")
pslq_search(det_IJ, 80, "det(I-J)")

# Also try: are tr(J) or cof in Q(sqrt(3))?
# If so, try PSLQ on {1, sqrt(3), tr(J), sqrt(3)*tr(J)}
s3 = sqrt(mpf(3))
print("\n--- tr(J) in Q(sqrt(3))? ---")
pslq_search_qsqrt3 = mpmath.pslq([mpf(1), s3, tr_J, s3*tr_J], maxcoeff=10**12)
if pslq_search_qsqrt3:
    print(f"  tr(J) relation over Q(sqrt(3)): {pslq_search_qsqrt3}")
else:
    print(f"  No linear relation found")

print("\nDone.")
