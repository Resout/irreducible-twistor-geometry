"""
POPOV IDENTIFICATION: IS dPhi = exp(-Hessian_Popov)?
=====================================================

The question: is the DS transfer operator the exponential of (minus) the
Popov Hessian restricted to the DS fibre?

    dPhi|_{m*} = exp(-d²S/da²|_{a*})   restricted to T_{m*}(L₁=1)

If yes: lambda_k = exp(-E_k), so Delta_DS = -ln(lambda_0) = E_0 = Delta_YM.
The DS theory IS quantised Yang-Mills.

THE COMPUTATION:
  1. Compute the equilibrium (0,1)-form a* = M*^{-1} . N . M*
     where N is the Nijenhuis data (dbar Phi, the Wirtinger Jacobian)
  2. Compute the Popov Hessian L(delta_a) = dbar_{J'} delta_a + [a*, delta_a]
     On fibre-constant forms, dbar_{J'} delta_a is the linearisation of Phi,
     i.e. the transfer operator Jacobian dPhi.
  3. The commutator [a*, .] is an additional 3x3 matrix on the tangent space.
  4. So: Hessian = dPhi + [a*, .] as a 3x3 matrix
  5. Compare eigenvalues to exp(-{1.263, 1.269, inf}) = {0.2829, 0.2813, 0}

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
# PRECISION FIRST
# ============================================================
PRECISION = 120
mp.dps = PRECISION

ONE = mpf(1); ZERO = mpf(0)
H = 3
FLOOR = ONE / mpf(H**3)
TARGET_K = mpf(7) / mpf(30)
I2 = eye(2)

sigma1 = matrix([[0, 1], [1, 0]])
sigma2 = matrix([[0, mpc(0,-1)], [mpc(0,1), 0]])
sigma3 = matrix([[1, 0], [0, -1]])
sq2 = sqrt(mpf(2))
dM_basis = [sigma1/sq2, sigma2/sq2, sigma3/sq2, I2/sq2]


# ============================================================
# 2x2 UTILITIES
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
    if total_sq < mpf(10)**(-80): return list(m)
    born = th**2 / total_sq
    if born >= FLOOR - mpf(10)**(-30): return list(m)
    ss = s1+s2+s3
    if fabs(ss) < mpf(10)**(-80): return list(m)
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
# COMPLEX DS (for Wirtinger)
# ============================================================
def complex_floor_enforce(m):
    s1, s2, s3, th = m
    s_mag_sq = fabs(s1)**2 + fabs(s2)**2 + fabs(s3)**2
    th_mag = fabs(th)
    total_mag_sq = s_mag_sq + th_mag**2
    if total_mag_sq < mpf(10)**(-80): return list(m)
    born = th_mag**2 / total_mag_sq
    if born >= FLOOR - mpf(10)**(-30): return list(m)
    S = s1 + s2 + s3
    S_mag_sq = fabs(S)**2
    if S_mag_sq < mpf(10)**(-80): return list(m)
    phase = th / th_mag if th_mag > mpf(10)**(-80) else ONE
    R = s_mag_sq / S_mag_sq
    Re_ph = phase.real if hasattr(phase, 'real') else mpf(phase)
    ac = mpf(26) - R
    bc = mpf(2) * R * Re_ph
    cc = -R
    disc = bc**2 - 4*ac*cc
    r_new = (-bc + sqrt(disc)) / (2*ac)
    if hasattr(r_new, 'real') and r_new.real <= 0:
        r_new = (-bc - sqrt(disc)) / (2*ac)
    th_new = phase * r_new
    alpha = (ONE - th_new) / S
    return [s1*alpha, s2*alpha, s3*alpha, th_new]

def complex_full_step(m, e):
    m_ds, K = ds_combine(m, e)
    return complex_floor_enforce(m_ds)


# ============================================================
# WIRTINGER ANTI-HOLOMORPHIC JACOBIAN
# ============================================================
def compute_J_anti(m_center, e, eps_w=None):
    if eps_w is None:
        eps_w = mpf(10)**(-30)
    J_anti = matrix(4, 4)
    for alpha in range(4):
        e_vec = [ZERO]*4; e_vec[alpha] = ONE
        m_pr = [m_center[k] + eps_w * e_vec[k] for k in range(4)]
        m_mr = [m_center[k] - eps_w * e_vec[k] for k in range(4)]
        phi_pr, _ = full_step(m_pr, e)
        phi_mr, _ = full_step(m_mr, e)
        df_dx = [(phi_pr[k]-phi_mr[k])/(2*eps_w) for k in range(4)]
        m_pi = [mpc(m_center[k], eps_w * e_vec[k]) for k in range(4)]
        m_mi = [mpc(m_center[k], -eps_w * e_vec[k]) for k in range(4)]
        phi_pi = complex_full_step(m_pi, e)
        phi_mi = complex_full_step(m_mi, e)
        df_dy = [(phi_pi[k]-phi_mi[k])/(2*eps_w) for k in range(4)]
        for k in range(4):
            J_anti[k, alpha] = (df_dx[k] + mpc(0,1)*df_dy[k]) / 2
    return J_anti


# ============================================================
# ANALYTICAL JACOBIAN (transfer operator)
# ============================================================
def analytical_jacobian_4x4(m, e):
    """Return J4 = J_fl * J_DS (the full transfer operator Jacobian)."""
    s1, s2, s3, th = m
    e1, e2, e3, ph = e
    m_ds, K = ds_combine(m, e)
    omK = ONE - K
    S_vals = [s1*(e1+ph)+th*e1, s2*(e2+ph)+th*e2, s3*(e3+ph)+th*e3, th*ph]
    N = [S_vals[i]/omK for i in range(4)]

    dSdm = matrix(4,4)
    for i in range(3): dSdm[i,i] = e[i]+ph; dSdm[i,3] = e[i]
    dSdm[3,3] = ph
    dKdm = [e2+e3, e1+e3, e1+e2, ZERO]
    J_DS = matrix(4,4)
    for i in range(4):
        for j in range(4):
            J_DS[i,j] = (dSdm[i,j]*omK + S_vals[i]*dKdm[j]) / omK**2

    Ss = N[0]+N[1]+N[2]
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
        for j in range(3): J_fl[i,j] = (g if i==j else ZERO) + N[i]*dgdN[j]
    J_fl[3,0]=dtdN[0]; J_fl[3,1]=dtdN[1]; J_fl[3,2]=dtdN[2]; J_fl[3,3]=ZERO

    return J_fl * J_DS


# ============================================================
# 3x3 UTILITIES
# ============================================================
def inv3(M):
    a = [[M[i,j] for j in range(3)] for i in range(3)]
    det = (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
          -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
          +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
    adj = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            rows = [r for r in range(3) if r != j]
            cols = [c for c in range(3) if c != i]
            minor = a[rows[0]][cols[0]]*a[rows[1]][cols[1]] - a[rows[0]][cols[1]]*a[rows[1]][cols[0]]
            adj[i,j] = ((-1)**(i+j)) * minor
    return adj * (ONE/det)

def eigvals3(M):
    """Eigenvalues of 3x3 real matrix via characteristic polynomial."""
    tr = M[0,0] + M[1,1] + M[2,2]
    cof = (M[0,0]*M[1,1] - M[0,1]*M[1,0]
         + M[0,0]*M[2,2] - M[0,2]*M[2,0]
         + M[1,1]*M[2,2] - M[1,2]*M[2,1])
    det = (M[0,0]*(M[1,1]*M[2,2]-M[1,2]*M[2,1])
          -M[0,1]*(M[1,0]*M[2,2]-M[1,2]*M[2,0])
          +M[0,2]*(M[1,0]*M[2,1]-M[1,1]*M[2,0]))
    return tr, cof, det


# ============================================================
# STAGE 1: EQUILIBRIUM
# ============================================================
print("=" * 70)
print("POPOV IDENTIFICATION TEST")
print(f"Working precision: {PRECISION} digits")
print("=" * 70)

print("\n--- Stage 1: Equilibrium ---")
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
m_ds_star, _ = ds_combine(m_star, e_star)

M_star = mass_to_M(m_star)
M_star_inv = mat_inv(M_star)
M_ds = mass_to_M(m_ds_star)
M_ds_inv = mat_inv(M_ds)

print(f"  Converged in {time.time()-t0:.1f}s")
print(f"  |K-7/30| = {nstr(fabs(K_star-TARGET_K), 5)}")
print(f"  m* = [{nstr(m_star[0],12)}, {nstr(m_star[1],12)}, {nstr(m_star[2],12)}, {nstr(m_star[3],12)}]")


# ============================================================
# STAGE 2: TRANSFER OPERATOR (the DS side)
# ============================================================
print("\n--- Stage 2: DS transfer operator ---")

J4 = analytical_jacobian_4x4(m_star, e_star)

# Project to L1=1 tangent space
V = matrix(4, 3)
for i in range(3): V[i, i] = ONE
for i in range(3): V[3, i] = -ONE

VtV_inv = inv3(V.T * V)
J3_DS = VtV_inv * V.T * J4 * V  # 3x3 projected transfer operator

tr_DS, cof_DS, det_DS = eigvals3(J3_DS)
disc_DS = tr_DS**2 - 4*cof_DS
lambda_0 = (tr_DS + sqrt(fabs(disc_DS))) / 2
lambda_1 = (tr_DS - sqrt(fabs(disc_DS))) / 2
lambda_2 = det_DS / (lambda_0 * lambda_1) if fabs(lambda_0 * lambda_1) > mpf(10)**(-50) else ZERO

Delta_0 = -log(fabs(lambda_0))
Delta_1 = -log(fabs(lambda_1))

print(f"  J3_DS eigenvalues:")
print(f"    lambda_0 = {nstr(lambda_0, 18)}")
print(f"    lambda_1 = {nstr(lambda_1, 18)}")
print(f"    lambda_2 = {nstr(lambda_2, 18)}")
print(f"  DS mass gaps:")
print(f"    Delta_0 = -ln(lambda_0) = {nstr(Delta_0, 15)}")
print(f"    Delta_1 = -ln(lambda_1) = {nstr(Delta_1, 15)}")

# Print the full 3x3 matrix
print(f"\n  J3_DS matrix:")
for i in range(3):
    print(f"    [{nstr(J3_DS[i,0],12)}, {nstr(J3_DS[i,1],12)}, {nstr(J3_DS[i,2],12)}]")


# ============================================================
# STAGE 3: NIJENHUIS DATA (the bridge)
# ============================================================
print("\n--- Stage 3: Nijenhuis data a* = M*^{-1} . N . M* ---")

J_anti = compute_J_anti(m_star, e_star)
print(f"  ||dbar(Phi)|| = {nstr(sqrt(sum(fabs(J_anti[i,j])**2 for i in range(4) for j in range(4))), 10)}")

# The equilibrium (0,1)-form in the bundle:
# For each mass direction alpha, the Nijenhuis data gives a 2x2 matrix:
#   N_alpha = sum_j dM_basis[j] * J_anti[j, alpha]
# The equilibrium form is a* = M_ds^{-1} * N_alpha
# (This is B[mu] from the Penrose residue script, but here we keep it
#  in mass-space indices alpha, not spacetime indices mu)

N_alpha = [zero_mat() for _ in range(4)]
a_star_alpha = [zero_mat() for _ in range(4)]
for alpha in range(4):
    for j in range(4):
        N_alpha[alpha] = N_alpha[alpha] + dM_basis[j] * J_anti[j, alpha]
    a_star_alpha[alpha] = M_ds_inv * N_alpha[alpha]

print(f"\n  Equilibrium (0,1)-form a*_alpha (4 matrices, one per mass direction):")
for alpha in range(4):
    print(f"    ||a*_{alpha}|| = {nstr(mat_norm(a_star_alpha[alpha]), 10)}")


# ============================================================
# STAGE 4: THE COMMUTATOR OPERATOR [a*, .]
# ============================================================
print("\n--- Stage 4: Commutator operator [a*, delta_a] ---")

# For a perturbation delta_m in tangent direction v (3D),
# the induced bundle perturbation is:
#   delta_a = sum_alpha (d a*/d m_alpha) * v_alpha
# But at leading order, delta_a_alpha = M_ds^{-1} * dM_basis[alpha] * v_alpha
# (the Nijenhuis data varies, but the leading term is just the Pauli embedding)

# For the commutator [a*, delta_a], we need to pick a specific representation.
# The key: a* is a 2x2 matrix for each mass direction alpha.
# A perturbation delta_m with tangent vector v = (v1, v2, v3) on L1=1
# (so v4 = -v1-v2-v3) induces:
#   delta_a = sum_{alpha=0}^{3} a_star_alpha * delta_m_alpha
#           = sum_{alpha=0}^{3} M_ds^{-1} * N_alpha * v_lifted[alpha]

# The tangent space basis: v_i has v[i]=1, v[3]=-1, others 0
tangent_vecs = []
for i in range(3):
    v = [ZERO]*4
    v[i] = ONE; v[3] = -ONE
    tangent_vecs.append(v)

# For each tangent direction, compute the induced a* perturbation (2x2 matrix)
delta_a_basis = []  # 3 matrices, one per tangent direction
for i in range(3):
    v = tangent_vecs[i]
    da = zero_mat()
    for alpha in range(4):
        da = da + a_star_alpha[alpha] * v[alpha]
    delta_a_basis.append(da)

print(f"  Induced bundle perturbations delta_a_i:")
for i in range(3):
    print(f"    ||delta_a_{i}|| = {nstr(mat_norm(delta_a_basis[i]), 10)}")

# The "average" a* that the commutator sees:
# We need a single 2x2 matrix representing a* for the commutator.
# The natural choice: a* projected onto the dominant Nijenhuis direction.
# Since J_anti has rank 1, all a*_alpha are proportional. Pick the largest:
norms = [mat_norm(a_star_alpha[alpha]) for alpha in range(4)]
dom_alpha = max(range(4), key=lambda i: norms[i])
a_star_dom = a_star_alpha[dom_alpha]
a_star_dom_norm = norms[dom_alpha]

# Actually, for the Popov Hessian, the commutator term is:
# [a*, delta_a] where a* and delta_a are both (0,1)-forms.
# In the fibre-constant case, this becomes a matrix commutator
# contracted over the form indices.
#
# The correct contraction: for fibre-constant forms,
# [a*, delta_a]_alpha = sum_beta [a*_beta, delta_a_alpha] (???)
#
# Wait -- the wedge product a* ^ delta_a involves antisymmetrisation
# over the form indices. For (0,1)-forms on CP^3 with 3 anti-holo directions:
# [a*, delta_a] = sum_{beta < alpha} [a*_beta, delta_a_alpha] - [a*_alpha, delta_a_beta]
#
# But since J_anti has rank 1 (all a*_alpha proportional to one matrix),
# all commutators [a*_alpha, a*_beta] = 0. So [a*, .] acts on perturbations
# that are NOT proportional to a*.
#
# For now, compute the full 3x3 commutator matrix:
# C_ij = "response of [a*, delta_a_j] in direction i"
# where the response is measured by projecting back onto the tangent basis.

# Method: for each tangent direction j, compute [a*_alpha, delta_a_j] for
# all alpha, form the 4-vector of results, project back to tangent space.

# But this requires choosing how to contract the form indices.
# Simplest approach: use the TOTAL a* = sum_alpha a*_alpha (sum over all directions)
a_star_total = zero_mat()
for alpha in range(4):
    a_star_total = a_star_total + a_star_alpha[alpha]

print(f"\n  Total a* = sum_alpha a*_alpha:")
print(f"    ||a*_total|| = {nstr(mat_norm(a_star_total), 10)}")
for r in range(2):
    print(f"    [{nstr(a_star_total[r,0], 12)}, {nstr(a_star_total[r,1], 12)}]")

# Commutator matrix C_ij: project [a*_total, delta_a_j] onto delta_a_i
# using Frobenius inner product: <A, B> = Re(tr(A^dag B))
C_comm = matrix(3, 3)
for j in range(3):
    comm_j = mat_comm(a_star_total, delta_a_basis[j])  # 2x2 matrix
    for i in range(3):
        # Project onto delta_a_i: Frobenius inner product
        inner = sum(delta_a_basis[i][r,c].conjugate() * comm_j[r,c]
                    for r in range(2) for c in range(2))
        norm_i = mat_norm_sq(delta_a_basis[i])
        C_comm[i, j] = inner / norm_i if norm_i > mpf(10)**(-50) else ZERO

print(f"\n  Commutator matrix C_ij = <delta_a_i, [a*_total, delta_a_j]> / ||delta_a_i||^2:")
for i in range(3):
    print(f"    [{nstr(C_comm[i,0],12)}, {nstr(C_comm[i,1],12)}, {nstr(C_comm[i,2],12)}]")

tr_C, cof_C, det_C = eigvals3(C_comm)
print(f"\n  Commutator eigenvalues:")
print(f"    tr(C) = {nstr(tr_C, 15)}")
# Eigenvalues from characteristic polynomial
disc_C = tr_C**2 - 4*cof_C
if fabs(disc_C) > mpf(10)**(-50):
    c_ev0 = (tr_C + sqrt(disc_C)) / 2
    c_ev1 = (tr_C - sqrt(disc_C)) / 2
else:
    c_ev0 = tr_C / 2
    c_ev1 = tr_C / 2
c_ev2 = det_C / (c_ev0 * c_ev1) if fabs(c_ev0 * c_ev1) > mpf(10)**(-50) else ZERO
print(f"    mu_0 = {nstr(c_ev0, 15)}")
print(f"    mu_1 = {nstr(c_ev1, 15)}")
print(f"    mu_2 = {nstr(c_ev2, 15)}")


# ============================================================
# STAGE 5: THE POPOV HESSIAN = J3_DS + C_comm
# ============================================================
print("\n--- Stage 5: Popov Hessian = J3_DS + [a*, .] ---")

H_popov = J3_DS + C_comm

print(f"  H_popov matrix:")
for i in range(3):
    print(f"    [{nstr(H_popov[i,0],12)}, {nstr(H_popov[i,1],12)}, {nstr(H_popov[i,2],12)}]")

tr_H, cof_H, det_H = eigvals3(H_popov)
disc_H = tr_H**2 - 4*cof_H
h_ev0 = (tr_H + sqrt(fabs(disc_H))) / 2
h_ev1 = (tr_H - sqrt(fabs(disc_H))) / 2
h_ev2 = det_H / (h_ev0 * h_ev1) if fabs(h_ev0 * h_ev1) > mpf(10)**(-50) else ZERO

print(f"\n  Popov Hessian eigenvalues:")
print(f"    eta_0 = {nstr(h_ev0, 18)}")
print(f"    eta_1 = {nstr(h_ev1, 18)}")
print(f"    eta_2 = {nstr(h_ev2, 18)}")

# Compare to DS eigenvalues
print(f"\n  Comparison:")
print(f"    DS:    lambda_0 = {nstr(lambda_0, 15)}, lambda_1 = {nstr(lambda_1, 15)}")
print(f"    Popov: eta_0    = {nstr(h_ev0, 15)}, eta_1    = {nstr(h_ev1, 15)}")
print(f"    Ratio: eta_0/lambda_0 = {nstr(h_ev0/lambda_0, 15)}")
print(f"    Ratio: eta_1/lambda_1 = {nstr(h_ev1/lambda_1, 15)}")

# The identification says: lambda_k = exp(-E_k)
# So if Hessian eigenvalues = E_k, then exp(-eta_k) should = lambda_k
print(f"\n  Testing lambda_k = exp(-eta_k):")
if h_ev0.real > ZERO:
    print(f"    exp(-eta_0) = {nstr(exp(-h_ev0), 15)}, lambda_0 = {nstr(lambda_0, 15)}")
if h_ev1.real > ZERO:
    print(f"    exp(-eta_1) = {nstr(exp(-h_ev1), 15)}, lambda_1 = {nstr(lambda_1, 15)}")

# Alternative: Hessian eigenvalues = -ln(lambda_k)
print(f"\n  Testing eta_k = -ln(lambda_k):")
print(f"    -ln(lambda_0) = {nstr(Delta_0, 15)}, eta_0 = {nstr(h_ev0, 15)}")
print(f"    -ln(lambda_1) = {nstr(Delta_1, 15)}, eta_1 = {nstr(h_ev1, 15)}")

# Alternative: Hessian IS the transfer operator (eta_k = lambda_k)
print(f"\n  Testing eta_k = lambda_k directly:")
print(f"    eta_0 = {nstr(h_ev0, 15)}, lambda_0 = {nstr(lambda_0, 15)}, diff = {nstr(fabs(h_ev0-lambda_0), 10)}")
print(f"    eta_1 = {nstr(h_ev1, 15)}, lambda_1 = {nstr(lambda_1, 15)}, diff = {nstr(fabs(h_ev1-lambda_1), 10)}")

# How big is the commutator correction?
print(f"\n  Size of commutator correction:")
print(f"    ||C_comm|| / ||J3_DS|| = {nstr(sqrt(sum(fabs(C_comm[i,j])**2 for i in range(3) for j in range(3))) / sqrt(sum(fabs(J3_DS[i,j])**2 for i in range(3) for j in range(3))), 10)}")


# ============================================================
# STAGE 6: ALTERNATIVE — LOGARITHM OF J3_DS
# ============================================================
print("\n--- Stage 6: Direct test — is -ln(J3_DS) the Popov Hessian? ---")

# If the identification is dPhi = exp(-H_popov), then H_popov = -ln(dPhi) = -ln(J3_DS)
# For a diagonalisable matrix with eigenvalues lambda_k:
# -ln(J3) has eigenvalues -ln(lambda_k) = {Delta_0, Delta_1, large}

# We can compute -ln(J3_DS) via eigendecomposition
# J3_DS = P * diag(lambda) * P^{-1}
# -ln(J3_DS) = P * diag(-ln(lambda)) * P^{-1}

# The eigenvalues of -ln(J3_DS):
logev_0 = -log(fabs(lambda_0))
logev_1 = -log(fabs(lambda_1))
logev_2 = -log(fabs(lambda_2)) if fabs(lambda_2) > mpf(10)**(-50) else mpf('inf')

print(f"  Eigenvalues of -ln(J3_DS):")
print(f"    E_0 = {nstr(logev_0, 15)} (= Delta_0)")
print(f"    E_1 = {nstr(logev_1, 15)} (= Delta_1)")
print(f"    E_2 = {nstr(logev_2, 15)} (~ inf if lambda_2 ~ 0)")

print(f"\n  These are the Popov energy levels IF dPhi = exp(-H_popov).")
print(f"  The mass gap is E_0 = {nstr(logev_0, 15)}.")


# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY: POPOV IDENTIFICATION TEST")
print(f"{'='*70}")
print(f"""
  DS transfer operator J3_DS:
    Eigenvalues: {nstr(lambda_0,12)}, {nstr(lambda_1,12)}, {nstr(lambda_2,8)}
    Mass gaps: {nstr(Delta_0,10)}, {nstr(Delta_1,10)}

  Commutator correction [a*, .]:
    Eigenvalues: {nstr(c_ev0,12)}, {nstr(c_ev1,12)}, {nstr(c_ev2,8)}

  Popov Hessian J3_DS + [a*, .]:
    Eigenvalues: {nstr(h_ev0,12)}, {nstr(h_ev1,12)}, {nstr(h_ev2,8)}

  KEY QUESTION: Does exp(-eta_k) = lambda_k?""")

if h_ev0.real > ZERO and h_ev1.real > ZERO:
    match_0 = fabs(exp(-h_ev0) - lambda_0)
    match_1 = fabs(exp(-h_ev1) - lambda_1)
    print(f"    |exp(-eta_0) - lambda_0| = {nstr(match_0, 10)}")
    print(f"    |exp(-eta_1) - lambda_1| = {nstr(match_1, 10)}")
    if match_0 < mpf('0.01') and match_1 < mpf('0.01'):
        print(f"\n    MATCH. The Popov Hessian eigenvalues give the DS transfer spectrum.")
        print(f"    The DS mass gap IS the Yang-Mills mass gap.")
    else:
        print(f"\n    NO MATCH at this level. The commutator correction shifts the spectrum.")
        print(f"    The identification may require a different contraction of [a*, .].")

print(f"""
  REGARDLESS: The matrix logarithm gives
    -ln(J3_DS) has eigenvalues {nstr(logev_0,10)}, {nstr(logev_1,10)}
    These ARE the mass gaps by construction.
    The question is whether THIS matrix equals the Popov Hessian.
""")
