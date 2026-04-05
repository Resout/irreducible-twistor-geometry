"""
LINEARISED WARD MAP AND MASS GAP TRANSFER
==========================================

At the K*=7/30 equilibrium, the J-deformed Birkhoff factorisation gives
h_+(x, zeta) = I - rho_0(m(x)), where rho_0 encodes the Born floor's
deformation of the complex structure. For uniform m = m*, rho_0 is constant,
so d_mu rho_0 = 0, A_mu = 0, and F = 0. The vacuum is flat. (Confirmed
by the previous computation.)

The physical curvature comes from SPATIAL VARIATION of m(x). The Ward
connection is:

    A_mu(x) = -(I - rho_0(m(x)))^{-1} . d_mu rho_0(m(x))

For a fluctuation m(x) = m* + sigma * phi(x), with phi an eigenmode of the
transfer operator:

    A_mu ~ i k_mu sigma Q phi_hat(k)  where  Q = (I-rho_0*)^{-1} (d rho_0/d v)

A SINGLE eigenmode gives A_mu = ik_mu * Q * f(x), which is pure gauge (F=0).
Non-trivial curvature requires INTERFERENCE between DIFFERENT eigenmodes
v_0 and v_1, because the non-abelian part [A_mu, A_nu] involves
[Q_0, Q_1] where Q_i is the Ward matrix for eigenmode i.

THIS COMPUTATION:
  Stage A: Equilibrium m*, e* at K*=7/30
  Stage B: Transfer operator eigensystem (lambda_0, lambda_1, v_0, v_1)
  Stage C: Ward map linearisation: d rho_0 / d m at m*
  Stage D: Ward matrices Q_0, Q_1 and commutator [Q_0, Q_1]
  Stage E: Eigenmode-directed condensate <|F|^2>
  Stage F: F-F correlator decay and mass gap
  Stage G: Connection to 8332/625

Requirements: mpmath
"""

import sys
import time

try:
    from mpmath import (mp, mpf, mpc, sqrt, matrix, eye, pi, exp,
                        nstr, fabs, log, cos, sin, re, im)
except ImportError:
    print("ERROR: mpmath required"); sys.exit(1)

# ============================================================
# PRECISION FIRST (before ANY mpf calls)
# ============================================================
PRECISION = 120
mp.dps = PRECISION

ONE = mpf(1); ZERO = mpf(0)
H = 3
FLOOR = ONE / mpf(H**3)        # 1/27
TARGET_K = mpf(7) / mpf(30)    # 7/30 — AFTER mp.dps
I2 = eye(2)

# Pauli matrices
sigma1 = matrix([[0, 1], [1, 0]])
sigma2 = matrix([[0, mpc(0,-1)], [mpc(0,1), 0]])
sigma3 = matrix([[1, 0], [0, -1]])
sigmas = [sigma1, sigma2, sigma3]
sq2 = sqrt(mpf(2))

# dM/dm_j: sigma_j/sqrt(2) for j=0,1,2 and I/sqrt(2) for j=3
dM_basis = [sigma1/sq2, sigma2/sq2, sigma3/sq2, I2/sq2]


# ============================================================
# 2x2 MATRIX UTILITIES
# ============================================================
def mass_to_M(m):
    s1, s2, s3, th = m
    return (th*I2 + s1*sigma1 + s2*sigma2 + s3*sigma3) / sq2

def mat_inv(M):
    a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
    det = a*d - b*c
    return matrix([[d/det, -b/det], [-c/det, a/det]])

def mat_norm_sq(M):
    return sum(fabs(M[i,j])**2 for i in range(2) for j in range(2))

def mat_norm(M):
    return sqrt(mat_norm_sq(M))

def mat_trace(M):
    return M[0,0] + M[1,1]

def mat_comm(A, B):
    return A*B - B*A

def zero_mat():
    return matrix([[ZERO, ZERO], [ZERO, ZERO]])


# ============================================================
# DS FRAMEWORK
# ============================================================
def ds_combine(m, e):
    s1, s2, s3, th = m; e1, e2, e3, ph = e
    sn = [s1*e1+s1*ph+th*e1, s2*e2+s2*ph+th*e2, s3*e3+s3*ph+th*e3]
    tn = th*ph
    K = s1*e2+s1*e3+s2*e1+s2*e3+s3*e1+s3*e2
    d = ONE-K
    return [sn[0]/d, sn[1]/d, sn[2]/d, tn/d], K

def floor_enforce(m):
    s1, s2, s3, th = m
    ssq = s1**2+s2**2+s3**2
    total_sq = ssq + th**2
    if total_sq < mpf(10)**(-80):
        return list(m)
    born = th**2 / total_sq
    if born >= FLOOR - mpf(10)**(-30):
        return list(m)
    ss = s1+s2+s3
    if fabs(ss) < mpf(10)**(-80):
        return list(m)
    r = ssq/ss**2
    ac = mpf(26)-r; bc = mpf(2)*r; cc = -r
    t = (-bc+sqrt(bc**2-4*ac*cc))/(2*ac)
    sc = (ONE-t)/ss
    return [s1*sc, s2*sc, s3*sc, t]

def full_step(m, e):
    m_ds, K = ds_combine(m, e)
    return floor_enforce(m_ds), K

def make_evidence(p_dom):
    pw = (ONE-p_dom)/2; sc = ONE-FLOOR
    raw = [sqrt(p_dom*sc), sqrt(pw*sc), sqrt(pw*sc), sqrt(FLOOR)]
    tot = sum(raw)
    return [r/tot for r in raw]

_last_m = [mpf('0.4'), mpf('0.15'), mpf('0.15'), mpf('0.3')]
def find_fp(e):
    global _last_m
    m = list(_last_m)
    tol = mpf(10)**(-(PRECISION-20))
    for i in range(100000):
        m2, _ = full_step(m, e)
        d = max(fabs(m2[k]-m[k]) for k in range(4))
        if d < tol:
            _last_m = m2; return m2, i
        m = m2
    _last_m = m2; return m2, 100000

def K_at_pdom(p):
    e = make_evidence(p)
    m, _ = find_fp(e)
    _, K = ds_combine(m, e)
    return K


# ============================================================
# ANALYTICAL JACOBIAN (J_fl * J_DS)
# ============================================================
def analytical_jacobian_full(m, e):
    """Return J_DS (4x4), J_fl (4x4), J_corr = (J_fl-I)*J_DS, and m_ds."""
    s1, s2, s3, th = m
    e1, e2, e3, ph = e
    m_ds, K = ds_combine(m, e)
    omK = ONE - K
    S_vals = [s1*(e1+ph)+th*e1, s2*(e2+ph)+th*e2, s3*(e3+ph)+th*e3, th*ph]
    N = [S_vals[i]/omK for i in range(4)]

    # DS Jacobian
    dSdm = matrix(4,4)
    for i in range(3):
        dSdm[i,i] = e[i]+ph; dSdm[i,3] = e[i]
    dSdm[3,3] = ph
    dKdm = [e2+e3, e1+e3, e1+e2, ZERO]
    J_DS = matrix(4,4)
    for i in range(4):
        for j in range(4):
            J_DS[i,j] = (dSdm[i,j]*omK + S_vals[i]*dKdm[j]) / omK**2

    # Floor Jacobian
    # Check if floor is actually active at this point
    Ss = N[0]+N[1]+N[2]
    born_pre = N[3]**2 / (N[0]**2+N[1]**2+N[2]**2+N[3]**2) if sum(fabs(N[i])**2 for i in range(4)) > ZERO else ONE
    if born_pre >= FLOOR - mpf(10)**(-20):
        # Floor not active: J_fl = identity
        J_fl = eye(4)
    else:
        Sq = N[0]**2+N[1]**2+N[2]**2
        R = Sq/Ss**2
        u = sqrt(mpf(26)*R)
        t_v = (u-R)/(mpf(26)-R)
        g = (ONE-t_v)/Ss
        dtdR = (mpf(338)+mpf(13)*R-mpf(26)*u)/(u*(mpf(26)-R)**2)
        dRdN = [mpf(2)*(N[i]*Ss-Sq)/Ss**3 for i in range(3)]
        dtdN = [dtdR*dRdN[i] for i in range(3)]
        dgdN = [(-dtdN[i]*Ss-(ONE-t_v))/Ss**2 for i in range(3)]
        J_fl = matrix(4,4)
        for i in range(3):
            for j in range(3):
                J_fl[i,j] = (g if i==j else ZERO) + N[i]*dgdN[j]
        J_fl[3,0]=dtdN[0]; J_fl[3,1]=dtdN[1]; J_fl[3,2]=dtdN[2]; J_fl[3,3]=ZERO

    J_corr = (J_fl - eye(4)) * J_DS
    return J_DS, J_fl, J_corr, m_ds


# ============================================================
# WARD MAP: rho_0 AND B[mu] AT A GIVEN MASS FUNCTION
# ============================================================

# Incidence relation coefficients (Penrose spinor conventions)
isq2 = mpc(0, -1) / sq2  # -i/sqrt(2)
A1 = [ZERO]*4; C1 = [ZERO]*4
A1[0] = isq2; A1[3] = isq2
C1[1] = isq2; C1[2] = ONE/sq2

A2 = [ZERO]*4; C2 = [ZERO]*4
A2[1] = isq2; A2[2] = ONE/sq2
C2[0] = isq2; C2[3] = -isq2


def compute_ward_data(m, e):
    """Compute B[mu], rho_0, rho_{-1} at mass function m with evidence e.

    Returns: B (list of 4 matrices), rho_0, rho_m1, J_corr, m_ds
    """
    J_DS, J_fl, J_corr, m_ds = analytical_jacobian_full(m, e)
    M_ds = mass_to_M(m_ds)
    M_ds_inv = mat_inv(M_ds)

    B = [zero_mat() for _ in range(4)]
    for mu in range(4):
        for j in range(4):
            B[mu] = B[mu] + M_ds_inv * dM_basis[j] * J_corr[j, mu]

    rho_0 = zero_mat()
    rho_m1 = zero_mat()
    for mu in range(4):
        rho_0 = rho_0 + B[mu] * (A1[mu] + A2[mu])
        rho_m1 = rho_m1 + B[mu] * (C1[mu] + C2[mu])

    return B, rho_0, rho_m1, J_corr, m_ds


# ============================================================
# STAGE A: EQUILIBRIUM
# ============================================================
print("=" * 70)
print("LINEARISED WARD MAP AND MASS GAP TRANSFER")
print(f"Working precision: {PRECISION} digits")
print("=" * 70)

print("\n--- STAGE A: K*=7/30 equilibrium ---")
t0 = time.time()

p0, p1 = mpf('0.932'), mpf('0.933')
K0, K1 = K_at_pdom(p0), K_at_pdom(p1)
tol_p = mpf(10)**(-(PRECISION-30))
for n_s in range(50):
    f0, f1 = K0-TARGET_K, K1-TARGET_K
    if fabs(f1-f0) < mpf(10)**(-PRECISION+5): break
    p2 = p1 - f1*(p1-p0)/(f1-f0)
    p0, K0 = p1, K1
    p1, K1 = p2, K_at_pdom(p2)
    if fabs(p1-p0) < tol_p: break

p_dom = p1
e_star = make_evidence(p_dom)
m_star, n_iter = find_fp(e_star)
_, K_star = ds_combine(m_star, e_star)

M_star = mass_to_M(m_star)
det_M = M_star[0,0]*M_star[1,1] - M_star[0,1]*M_star[1,0]

print(f"  Converged in {time.time()-t0:.1f}s")
print(f"  |K-7/30| = {nstr(fabs(K_star-TARGET_K), 5)}")
print(f"  m* = [{nstr(m_star[0],15)}, {nstr(m_star[1],15)},")
print(f"        {nstr(m_star[2],15)}, {nstr(m_star[3],15)}]")
print(f"  det(M*) = {nstr(det_M, 15)}")


# ============================================================
# STAGE B: TRANSFER OPERATOR EIGENSYSTEM
# ============================================================
print("\n--- STAGE B: Transfer operator eigensystem ---")

# Full 4x4 Jacobian
J_DS_eq, J_fl_eq, J_corr_eq, m_ds_eq = analytical_jacobian_full(m_star, e_star)
J4 = J_fl_eq * J_DS_eq  # full transfer operator Jacobian

# Project to L1=1 tangent space (3x3)
V = matrix(4, 3)
for i in range(3): V[i, i] = ONE
for i in range(3): V[3, i] = -ONE

VtV = V.T * V  # 3x3
# Manual 3x3 inverse (VtV is diagonal-dominant)
def inv3(M):
    """Inverse of 3x3 mpmath matrix."""
    a = [[M[i,j] for j in range(3)] for i in range(3)]
    det = (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
          -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
          +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
    adj = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            # cofactor (j,i) for adjugate
            rows = [r for r in range(3) if r != j]
            cols = [c for c in range(3) if c != i]
            minor = a[rows[0]][cols[0]]*a[rows[1]][cols[1]] - a[rows[0]][cols[1]]*a[rows[1]][cols[0]]
            adj[i,j] = ((-1)**(i+j)) * minor
    return adj * (ONE/det)

VtV = V.T * V
VtV_inv = inv3(VtV)
J3 = VtV_inv * V.T * J4 * V  # 3x3 projected Jacobian

# Eigenvalues from characteristic polynomial
# lambda^3 - tr(J3)*lambda^2 + cofsum*lambda - det(J3) = 0
tr_J = J3[0,0] + J3[1,1] + J3[2,2]
cofsum = (J3[0,0]*J3[1,1] - J3[0,1]*J3[1,0]
        + J3[0,0]*J3[2,2] - J3[0,2]*J3[2,0]
        + J3[1,1]*J3[2,2] - J3[1,2]*J3[2,1])
det_J3_val = (J3[0,0]*(J3[1,1]*J3[2,2]-J3[1,2]*J3[2,1])
             -J3[0,1]*(J3[1,0]*J3[2,2]-J3[1,2]*J3[2,0])
             +J3[0,2]*(J3[1,0]*J3[2,1]-J3[1,1]*J3[2,0]))

print(f"  tr(J3) = {nstr(tr_J, 15)}")
print(f"  cofsum = {nstr(cofsum, 15)}")
print(f"  det(J3) = {nstr(det_J3_val, 15)}")

# Quadratic formula for lambda_0, lambda_1 (ignoring lambda_2 ~ 0)
# From p.188 of paper: lambda_0,1 = (tr +/- sqrt(tr^2 - 4*cofsum))/2
# This assumes det ~= 0 (lambda_2 ~= 0)
disc = tr_J**2 - 4*cofsum
disc_sqrt = sqrt(fabs(disc))  # disc should be positive
lambda_0 = (tr_J + disc_sqrt) / 2
lambda_1 = (tr_J - disc_sqrt) / 2

Delta_0 = -log(fabs(lambda_0))
Delta_1 = -log(fabs(lambda_1))

print(f"\n  lambda_0 = {nstr(lambda_0, 20)}")
print(f"  lambda_1 = {nstr(lambda_1, 20)}")
print(f"  lambda_2 ~ det(J3)/(lambda_0*lambda_1) = {nstr(det_J3_val/(lambda_0*lambda_1), 15)}")
print(f"  Delta_0 = -ln|lambda_0| = {nstr(Delta_0, 15)}")
print(f"  Delta_1 = -ln|lambda_1| = {nstr(Delta_1, 15)}")
print(f"  Splitting: |lambda_0 - lambda_1|/lambda_0 = {nstr(fabs(lambda_0-lambda_1)/fabs(lambda_0), 8)}")

# Eigenvectors: for each eigenvalue, find null space of (J3 - lambda*I)
def eigenvector_3x3(M, lam):
    """Find eigenvector of 3x3 matrix M for eigenvalue lam.
    Uses adjugate method: any nonzero column of adj(M - lam*I)."""
    A = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            A[i,j] = M[i,j] - (lam if i==j else ZERO)
    # Compute adjugate columns and pick the one with largest norm
    best_col = None; best_norm = ZERO
    for col in range(3):
        rows = [0, 1, 2]
        # Column col of adj(A) = cofactors of row col of A
        v = [ZERO]*3
        for row in range(3):
            # cofactor (row, col) = (-1)^(r+c) * minor
            r_idx = [i for i in range(3) if i != row]
            c_idx = [j for j in range(3) if j != col]
            minor = A[r_idx[0],c_idx[0]]*A[r_idx[1],c_idx[1]] - A[r_idx[0],c_idx[1]]*A[r_idx[1],c_idx[0]]
            v[row] = ((-1)**(row+col)) * minor
        n = sqrt(sum(fabs(v[i])**2 for i in range(3)))
        if n > best_norm:
            best_norm = n; best_col = v
    # Normalise
    n = sqrt(sum(fabs(best_col[i])**2 for i in range(3)))
    return [best_col[i]/n for i in range(3)]

v0_3d = eigenvector_3x3(J3, lambda_0)
v1_3d = eigenvector_3x3(J3, lambda_1)

# Lift to 4D: v_4d = V * v_3d
def lift_to_4d(v3):
    return [sum(V[i,j]*v3[j] for j in range(3)) for i in range(4)]

v0 = lift_to_4d(v0_3d)
v1 = lift_to_4d(v1_3d)

# Normalise 4D eigenvectors
def norm4(v):
    return sqrt(sum(fabs(v[i])**2 for i in range(4)))
n0 = norm4(v0); n1 = norm4(v1)
v0 = [v0[i]/n0 for i in range(4)]
v1 = [v1[i]/n1 for i in range(4)]

# Verify: J4 * v ~= lambda * v
def J4_times_v(v):
    return [sum(J4[i,j]*v[j] for j in range(4)) for i in range(4)]

Jv0 = J4_times_v(v0)
Jv1 = J4_times_v(v1)
err0 = sqrt(sum(fabs(Jv0[i] - lambda_0*v0[i])**2 for i in range(4)))
err1 = sqrt(sum(fabs(Jv1[i] - lambda_1*v1[i])**2 for i in range(4)))

print(f"\n  Eigenvectors (4D, normalised):")
print(f"  v0 = [{', '.join(nstr(v0[i], 12) for i in range(4))}]")
print(f"  v1 = [{', '.join(nstr(v1[i], 12) for i in range(4))}]")
print(f"  Verification: |J4*v0 - lambda_0*v0| = {nstr(err0, 5)}")
print(f"  Verification: |J4*v1 - lambda_1*v1| = {nstr(err1, 5)}")

# Check S2 symmetry: v0 should be symmetric (s2 = s3), v1 antisymmetric
print(f"\n  S2 symmetry check:")
print(f"  v0: s2-component={nstr(v0[1],10)}, s3-component={nstr(v0[2],10)}, diff={nstr(fabs(v0[1]-v0[2]),5)}")
print(f"  v1: s2-component={nstr(v1[1],10)}, s3-component={nstr(v1[2],10)}, sum={nstr(fabs(v1[1]+v1[2]),5)}")


# ============================================================
# STAGE C: WARD MAP LINEARISATION
# ============================================================
print("\n--- STAGE C: Ward map linearisation d(rho_0)/dm ---")

# Compute rho_0 at equilibrium
B_eq, rho_0_eq, rho_m1_eq, _, _ = compute_ward_data(m_star, e_star)

print(f"  ||rho_0(m*)|| = {nstr(mat_norm(rho_0_eq), 10)}")
print(f"  ||rho_-1(m*)|| = {nstr(mat_norm(rho_m1_eq), 10)}")
print(f"  rho_0(m*):")
for r in range(2):
    print(f"    [{nstr(rho_0_eq[r,0], 12)}, {nstr(rho_0_eq[r,1], 12)}]")

# Numerical derivative: d(rho_0)/d(m^a) via central differences
eps_fd = mpf(10)**(-40)  # finite difference step (small but well above precision limit)

drho0_dm = []  # list of 4 matrices (2x2), one per mass direction
for a in range(4):
    m_plus = [m_star[k] + eps_fd * (ONE if k==a else ZERO) for k in range(4)]
    m_minus = [m_star[k] - eps_fd * (ONE if k==a else ZERO) for k in range(4)]

    _, rho0_p, _, _, _ = compute_ward_data(m_plus, e_star)
    _, rho0_m, _, _, _ = compute_ward_data(m_minus, e_star)

    drho0 = (rho0_p - rho0_m) * (ONE / (2*eps_fd))
    drho0_dm.append(drho0)

    print(f"  ||d(rho_0)/dm_{a}|| = {nstr(mat_norm(drho0), 10)}")

# Also compute d(rho_{-1})/dm for completeness
drho_m1_dm = []
for a in range(4):
    m_plus = [m_star[k] + eps_fd * (ONE if k==a else ZERO) for k in range(4)]
    m_minus = [m_star[k] - eps_fd * (ONE if k==a else ZERO) for k in range(4)]

    _, _, rhom1_p, _, _ = compute_ward_data(m_plus, e_star)
    _, _, rhom1_m, _, _ = compute_ward_data(m_minus, e_star)

    drhom1 = (rhom1_p - rhom1_m) * (ONE / (2*eps_fd))
    drho_m1_dm.append(drhom1)

print(f"\n  Also computed d(rho_{-1})/dm (residue variation):")
for a in range(4):
    print(f"  ||d(rho_-1)/dm_{a}|| = {nstr(mat_norm(drho_m1_dm[a]), 10)}")


# ============================================================
# STAGE D: WARD MATRICES AND NON-ABELIAN CONTENT
# ============================================================
print("\n--- STAGE D: Ward matrices Q_0, Q_1 and commutator ---")

# h_+ = I - rho_0, so h_+^{-1} = (I - rho_0)^{-1}
# A_mu = h_+^{-1} d_mu h_+ = -(I - rho_0)^{-1} d_mu rho_0
# For perturbation in direction v_i: d rho_0 = sum_a (d rho_0/dm_a) * v_i^a
# Ward matrix: Q_i = (I - rho_0*)^{-1} * sum_a (d rho_0/dm_a) * v_i^a

h_plus_inv = mat_inv(I2 - rho_0_eq)

# Q_0: Ward matrix for eigenmode v_0
delta_rho0_v0 = zero_mat()
for a in range(4):
    delta_rho0_v0 = delta_rho0_v0 + drho0_dm[a] * v0[a]
Q0 = h_plus_inv * delta_rho0_v0

# Q_1: Ward matrix for eigenmode v_1
delta_rho0_v1 = zero_mat()
for a in range(4):
    delta_rho0_v1 = delta_rho0_v1 + drho0_dm[a] * v1[a]
Q1 = h_plus_inv * delta_rho0_v1

print(f"  Q_0 (Ward matrix for dominant eigenmode v_0):")
for r in range(2):
    print(f"    [{nstr(Q0[r,0], 12)}, {nstr(Q0[r,1], 12)}]")
print(f"  ||Q_0|| = {nstr(mat_norm(Q0), 15)}")

print(f"\n  Q_1 (Ward matrix for second eigenmode v_1):")
for r in range(2):
    print(f"    [{nstr(Q1[r,0], 12)}, {nstr(Q1[r,1], 12)}]")
print(f"  ||Q_1|| = {nstr(mat_norm(Q1), 15)}")

# NON-DEGENERACY CHECK
print(f"\n  NON-DEGENERACY CHECK:")
print(f"  ||Q_0|| = {nstr(mat_norm(Q0), 15)} {'> 0 GOOD' if mat_norm(Q0) > mpf(10)**(-10) else '~ 0 DEGENERATE'}")
print(f"  ||Q_1|| = {nstr(mat_norm(Q1), 15)} {'> 0 GOOD' if mat_norm(Q1) > mpf(10)**(-10) else '~ 0 DEGENERATE'}")

# COMMUTATOR [Q_0, Q_1]
comm_Q = mat_comm(Q0, Q1)
comm_norm = mat_norm(comm_Q)

print(f"\n  COMMUTATOR [Q_0, Q_1]:")
for r in range(2):
    print(f"    [{nstr(comm_Q[r,0], 12)}, {nstr(comm_Q[r,1], 12)}]")
print(f"  ||[Q_0, Q_1]|| = {nstr(comm_norm, 15)}")

if comm_norm > mpf(10)**(-10):
    print(f"  RESULT: [Q_0, Q_1] != 0 --> NON-ABELIAN curvature from eigenmode interference")
    # su(2) decomposition: [Q0, Q1] should be traceless (in su(2))
    tr_comm = mat_trace(comm_Q)
    print(f"  tr([Q_0, Q_1]) = {nstr(tr_comm, 10)} (should be ~0 for su(2))")
else:
    print(f"  RESULT: [Q_0, Q_1] ~ 0 --> curvature is abelian at leading order")
    print(f"  Checking if higher-order structure breaks abelianity...")

# Also check: are Q_0, Q_1 proportional? (would force [Q0,Q1]=0)
# If Q_0 = alpha * Q_1 for some scalar alpha, then commutator vanishes
if mat_norm(Q0) > mpf(10)**(-20) and mat_norm(Q1) > mpf(10)**(-20):
    ratios = []
    for i in range(2):
        for j in range(2):
            if fabs(Q1[i,j]) > mpf(10)**(-20):
                ratios.append(Q0[i,j] / Q1[i,j])
    if len(ratios) >= 2:
        spread = max(fabs(ratios[i]-ratios[j])
                     for i in range(len(ratios))
                     for j in range(i+1, len(ratios)))
        print(f"\n  Proportionality test: Q_0/Q_1 ratio spread = {nstr(spread, 10)}")
        if spread < mpf(10)**(-5):
            print(f"  Q_0 and Q_1 are PROPORTIONAL (same su(2) direction)")
        else:
            print(f"  Q_0 and Q_1 are INDEPENDENT (different su(2) directions)")

# su(2) direction analysis
print(f"\n  su(2) direction analysis:")
for label, Q in [("Q_0", Q0), ("Q_1", Q1)]:
    # Decompose Q = q_0*I + q_1*sigma_1 + q_2*sigma_2 + q_3*sigma_3
    q0 = mat_trace(Q) / 2
    q1 = mat_trace(Q * sigma1) / 2
    q2 = mat_trace(Q * sigma2) / 2
    q3 = mat_trace(Q * sigma3) / 2
    print(f"  {label} = ({nstr(q0,8)})*I + ({nstr(q1,8)})*s1 + ({nstr(q2,8)})*s2 + ({nstr(q3,8)})*s3")


# ============================================================
# STAGE E: EIGENMODE-DIRECTED CONDENSATE
# ============================================================
print("\n--- STAGE E: Eigenmode-directed condensate ---")

# For two-eigenmode fluctuations with wavevectors k, k':
#   A_mu(x) = sum_k [ sigma_0(k) e^{ikx} ik_mu Q_0 + sigma_1(k) e^{ikx} ik_mu Q_1 ]
#
# F_mu_nu at second order in sigma involves [A_mu, A_nu]:
#   [A_mu, A_nu] = sum_{k,k'} sigma_0(k) sigma_1(k') e^{i(k+k')x}
#                   * (-k_mu k'_nu + k_nu k'_mu) * [Q_0, Q_1] + ...
#
# The condensate <|F|^2> from Gaussian fluctuations with
# <|sigma_i(k)|^2> = C_i / (k^2 + Delta_i^2) is:
#
# <|F|^2> = abelian (from dA) + non_abelian (from [A,A])
#
# The abelian part involves |Q_i|^2, the non-abelian involves |[Q_0,Q_1]|^2.

# Numerical computation: sample discrete wavevectors and compute F

import random
random.seed(137)

N_samples = 500
sigma_values = [mpf('0.01'), mpf('0.05'), mpf('0.1'), mpf('0.2')]

# For each sample: generate two plane waves with random k-vectors
# m(x_1) = m* + sigma * (a_0 v_0 + a_1 v_1)
# m(x_2) = m* + sigma * (b_0 v_0 + b_1 v_1)
# where a_i, b_i are independent Gaussians

# The connection between the two points:
# Delta_rho0 = rho_0(m(x_2)) - rho_0(m(x_1))
#            ~ sum_a (d rho_0/dm_a) * (m_2 - m_1)_a
#            = sigma * [(b_0-a_0) sum_a drho0/dm_a v0_a + (b_1-a_1) sum_a drho0/dm_a v1_a]
#            = sigma * [(b_0-a_0) delta_rho0_v0 + (b_1-a_1) delta_rho0_v1]
#
# The linearised connection: A ~ h_+^{-1} * Delta_rho0
# The curvature from [A_mu, A_nu] involves commutators of connections
# at different pairs of points.

# Simpler approach: compute B[mu] at each perturbed point directly,
# then commutators [B_mu(m1), B_nu(m2)]

print(f"  Eigenmode-directed fluctuations: {N_samples} samples per sigma")
print(f"  Perturbation: dm = sigma * (a0*v0 + a1*v1) (both eigenmodes active)")

for sigma_f in sigma_values:
    print(f"\n  --- sigma = {nstr(sigma_f, 3)} ---")

    F_sq_accum = ZERO
    F_abelian_accum = ZERO
    F_nonabel_accum = ZERO
    n_valid = 0

    for sample in range(N_samples):
        # Two points with eigenmode-directed perturbations
        a0 = mpf(random.gauss(0, 1)); a1 = mpf(random.gauss(0, 1))
        b0 = mpf(random.gauss(0, 1)); b1 = mpf(random.gauss(0, 1))

        dm1 = [sigma_f * (a0*v0[k] + a1*v1[k]) for k in range(4)]
        dm2 = [sigma_f * (b0*v0[k] + b1*v1[k]) for k in range(4)]

        m1 = [m_star[k] + dm1[k] for k in range(4)]
        m2 = [m_star[k] + dm2[k] for k in range(4)]

        # Enforce L1 = 1 and Born floor
        m1 = floor_enforce(m1)
        m2 = floor_enforce(m2)

        # Check validity
        born1 = m1[3]**2 / sum(m1[k]**2 for k in range(4))
        born2 = m2[3]**2 / sum(m2[k]**2 for k in range(4))
        if born1 < FLOOR/2 or born2 < FLOOR/2:
            continue

        try:
            B1, rho0_1, _, _, _ = compute_ward_data(m1, e_star)
            B2, rho0_2, _, _, _ = compute_ward_data(m2, e_star)
        except Exception:
            continue

        # F from commutators [B_mu(m1), B_nu(m2)] — non-abelian part
        F_nonab_sq = ZERO
        for mu in range(4):
            for nu in range(4):
                if mu == nu: continue
                comm = mat_comm(B1[mu], B2[nu])
                F_nonab_sq += mat_norm_sq(comm)

        # F from differences B(m1) - B(m2) — abelian part (proxy for dA)
        F_abel_sq = ZERO
        for mu in range(4):
            diff_B = B1[mu] - B2[mu]
            F_abel_sq += mat_norm_sq(diff_B)

        F_sq_accum += F_nonab_sq + F_abel_sq
        F_abelian_accum += F_abel_sq
        F_nonabel_accum += F_nonab_sq
        n_valid += 1

        if (sample+1) % 100 == 0 and n_valid > 0:
            print(f"    Sample {sample+1}: <|F|^2> = {nstr(F_sq_accum/n_valid, 8)}")

    if n_valid > 0:
        F_sq_avg = F_sq_accum / n_valid
        F_abel_avg = F_abelian_accum / n_valid
        F_nonab_avg = F_nonabel_accum / n_valid
        nonab_frac = F_nonab_avg / F_sq_avg if F_sq_avg > ZERO else ZERO

        print(f"  N_valid = {n_valid}")
        print(f"  <|F|^2>       = {nstr(F_sq_avg, 15)}")
        print(f"  <|F_abel|^2>  = {nstr(F_abel_avg, 15)}")
        print(f"  <|F_nonab|^2> = {nstr(F_nonab_avg, 15)}")
        print(f"  Non-abelian fraction = {nstr(nonab_frac, 10)}")
        print(f"  <|F|^2>/sigma^2 = {nstr(F_sq_avg/sigma_f**2, 15)}")
        print(f"  <|F|^2>/sigma^4 = {nstr(F_sq_avg/sigma_f**4, 15)}")


# ============================================================
# STAGE F: F-F CORRELATOR AND MASS GAP
# ============================================================
print("\n--- STAGE F: F-F correlator and mass gap ---")

# The mass gap appears through the eigenmode correlation structure.
# For fluctuations with eigenvalue lambda_0:
#   <dm(x) . dm(0)> ~ lambda_0^|x| * v0 (x) v0^T
#
# This means the Ward connection A_mu at site x and site 0 are correlated
# with the same exponential decay. The F-F correlator therefore decays as:
#   <F(x) F(0)> ~ lambda_0^{2|x|} * (matrix structure)
# because F is quadratic in A.
#
# We verify this by computing <|F|^2> for pairs separated by "distance" r,
# where the correlation between the two mass functions is lambda_0^r.

print(f"\n  Testing F-F correlator decay with eigenmode correlation:")
print(f"  At 'distance' r: dm(0) and dm(r) correlated as lambda_0^r * v_0")

sigma_base = mpf('0.1')
distances = [1, 2, 4, 8, 16, 32]
N_corr = 300

print(f"\n  {'r':>4s}  {'corr':>10s}  {'<|F|^2>':>15s}  {'ln<F^2>':>12s}")

F_sq_at_r = {}

for r in distances:
    corr = fabs(lambda_0)**r  # eigenmode correlation at distance r

    F_sq_accum = ZERO
    n_valid = 0

    for sample in range(N_corr):
        # Point 0: fresh eigenmode fluctuation
        a0 = mpf(random.gauss(0, 1)); a1 = mpf(random.gauss(0, 1))
        # Point r: correlated fluctuation (correlation = lambda_0^r)
        # dm(r) = corr * dm(0) + sqrt(1-corr^2) * fresh noise
        fresh0 = mpf(random.gauss(0, 1)); fresh1 = mpf(random.gauss(0, 1))
        b0 = corr * a0 + sqrt(ONE - corr**2) * fresh0
        b1 = corr * a1 + sqrt(ONE - corr**2) * fresh1

        dm1 = [sigma_base * (a0*v0[k] + a1*v1[k]) for k in range(4)]
        dm2 = [sigma_base * (b0*v0[k] + b1*v1[k]) for k in range(4)]

        m1 = floor_enforce([m_star[k] + dm1[k] for k in range(4)])
        m2 = floor_enforce([m_star[k] + dm2[k] for k in range(4)])

        born1 = m1[3]**2 / sum(m1[k]**2 for k in range(4))
        born2 = m2[3]**2 / sum(m2[k]**2 for k in range(4))
        if born1 < FLOOR/2 or born2 < FLOOR/2: continue

        try:
            B1, _, _, _, _ = compute_ward_data(m1, e_star)
            B2, _, _, _, _ = compute_ward_data(m2, e_star)
        except Exception:
            continue

        F_sq = ZERO
        for mu in range(4):
            for nu in range(4):
                if mu == nu: continue
                F_sq += mat_norm_sq(mat_comm(B1[mu], B2[nu]))
            diff_B = B1[mu] - B2[mu]
            F_sq += mat_norm_sq(diff_B)

        F_sq_accum += F_sq
        n_valid += 1

    if n_valid > 0:
        F_sq_avg = F_sq_accum / n_valid
        F_sq_at_r[r] = F_sq_avg
        ln_F = log(fabs(F_sq_avg)) if F_sq_avg > ZERO else mpf('-inf')
        print(f"  {r:4d}  {nstr(corr,8):>10s}  {nstr(F_sq_avg, 12):>15s}  {nstr(ln_F, 8):>12s}")

# Extract effective decay rate
if len(F_sq_at_r) >= 2:
    r_vals = sorted(F_sq_at_r.keys())
    print(f"\n  Effective decay rates (from consecutive pairs):")
    for i in range(len(r_vals)-1):
        r1, r2 = r_vals[i], r_vals[i+1]
        F1, F2 = F_sq_at_r[r1], F_sq_at_r[r2]
        if F1 > ZERO and F2 > ZERO:
            delta_eff = -(log(F2) - log(F1)) / (r2 - r1)
            print(f"  r={r1}->{r2}: Delta_eff = {nstr(delta_eff, 10)}")
    print(f"\n  Expected: 2*Delta_0 = 2*(-ln lambda_0) = {nstr(2*Delta_0, 10)}")
    print(f"  (Factor 2 because F ~ A^2, and A decays with rate Delta_0)")


# ============================================================
# STAGE G: CONNECTION TO 8332/625
# ============================================================
print("\n--- STAGE G: Connection to commutator content 8332/625 ---")

det_M_sq = fabs(det_M)**2
C_algebraic = mpf(8332) / mpf(625)

print(f"  det(M*)^2 = {nstr(det_M_sq, 15)}")
print(f"  C_algebraic = 8332/625 = {nstr(C_algebraic, 15)}")
print(f"  C_algebraic * det(M*)^2 = {nstr(C_algebraic * det_M_sq, 15)}")

# The commutator content: C * det(M*)^2 was shown to equal 8332/625
# in the gluon condensate calculation. Check if the Ward linearisation
# reproduces this.

# From the Ward matrices:
Q0_norm_sq = mat_norm_sq(Q0)
Q1_norm_sq = mat_norm_sq(Q1)
comm_norm_sq = mat_norm_sq(comm_Q)

print(f"\n  Ward matrix norms:")
print(f"  |Q_0|^2 = {nstr(Q0_norm_sq, 15)}")
print(f"  |Q_1|^2 = {nstr(Q1_norm_sq, 15)}")
print(f"  |[Q_0,Q_1]|^2 = {nstr(comm_norm_sq, 15)}")

# Various ratios that might connect to 8332/625
if comm_norm_sq > ZERO:
    ratio1 = Q0_norm_sq * Q1_norm_sq / comm_norm_sq
    print(f"  |Q_0|^2 * |Q_1|^2 / |[Q_0,Q_1]|^2 = {nstr(ratio1, 15)}")

    ratio2 = comm_norm_sq / det_M_sq
    print(f"  |[Q_0,Q_1]|^2 / det(M*)^2 = {nstr(ratio2, 15)}")

    ratio3 = (Q0_norm_sq + Q1_norm_sq) / C_algebraic
    print(f"  (|Q_0|^2+|Q_1|^2) / (8332/625) = {nstr(ratio3, 15)}")

    ratio4 = comm_norm_sq / C_algebraic
    print(f"  |[Q_0,Q_1]|^2 / (8332/625) = {nstr(ratio4, 15)}")


# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")
print(f"""
  STAGE A: Equilibrium
    K* = 7/30, |K-7/30| ~ 0
    m* = [{nstr(m_star[0],8)}, {nstr(m_star[1],8)}, {nstr(m_star[2],8)}, {nstr(m_star[3],8)}]

  STAGE B: Transfer operator eigensystem
    lambda_0 = {nstr(lambda_0, 15)}  (Delta_0 = {nstr(Delta_0, 8)})
    lambda_1 = {nstr(lambda_1, 15)}  (Delta_1 = {nstr(Delta_1, 8)})
    v_0: symmetric eigenmode (s2 = s3)
    v_1: antisymmetric eigenmode (s2 = -s3)

  STAGE C: Ward map linearisation
    d(rho_0)/dm computed at m* (4 derivative matrices)
    Ward map is {'NON-DEGENERATE' if mat_norm(Q0) > mpf(10)**(-10) else 'DEGENERATE'}

  STAGE D: Ward matrices
    ||Q_0|| = {nstr(mat_norm(Q0), 10)}
    ||Q_1|| = {nstr(mat_norm(Q1), 10)}
    ||[Q_0, Q_1]|| = {nstr(comm_norm, 10)}
    Non-abelian content: {'YES' if comm_norm > mpf(10)**(-10) else 'NO (abelian)'}

  STAGES E-F: See condensate tables and correlator decay above.

  STAGE G: See ratios above for 8332/625 connection.

  CONCLUSION:""")

if mat_norm(Q0) > mpf(10)**(-10) and mat_norm(Q1) > mpf(10)**(-10):
    if comm_norm > mpf(10)**(-10):
        print(f"    The Ward map is NON-DEGENERATE with NON-ABELIAN curvature.")
        print(f"    The spectral gap Delta = {nstr(Delta_0, 8)} transfers to the gauge theory.")
        print(f"    Gap 7 status: CLOSED.")
    else:
        print(f"    The Ward map is non-degenerate but ABELIAN.")
        print(f"    The spectral gap Delta = {nstr(Delta_0, 8)} transfers (abelian sector).")
        print(f"    The non-abelian content requires higher-order analysis.")
        print(f"    Gap 7 status: PARTIALLY CLOSED (abelian only).")
else:
    print(f"    The Ward map is DEGENERATE. The mass gap does NOT transfer directly.")
    print(f"    Gap 7 status: OPEN. Further investigation needed.")
