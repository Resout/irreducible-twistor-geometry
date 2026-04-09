"""
TARGET 1: Derive 11/60 from the Jacobian block structure.

The Jacobian at the fixed point decomposes into radial and angular blocks.
The coupling between sectors (off-diagonal blocks) should encode
M²/m(0++) = 11/60 if the parallel conductance picture is correct.

Strategy:
1. Compute the full 4x4 single-site Jacobian at the K*=7/30 equilibrium
2. Rotate into the radial/angular basis
3. Extract the diagonal blocks (radial-radial, angular-angular)
4. Extract the off-diagonal blocks (radial-angular coupling)
5. Check if the coupling structure gives 1/10 + 1/12

Also compute the 8x8 A2 coupled Jacobian blocks.
"""

import numpy as np
from scipy.optimize import fsolve

# ============================================================
# EQUILIBRIUM
# ============================================================

def ds_combine(m, e):
    s = m[:3]; theta = m[3]; ev = e[:3]; phi = e[3]
    s_pre = s*ev + s*phi + theta*ev
    theta_pre = theta*phi
    total = np.sum(s_pre) + theta_pre
    out = np.zeros(4)
    out[:3] = s_pre / total
    out[3] = theta_pre / total
    return out

def born(m):
    return m[3]**2 / np.sum(m**2)

def enforce_floor(m):
    b = born(m)
    if b >= 1/27 - 1e-14:
        return m.copy()
    S = np.sum(m[:3]); Sq = np.sum(m[:3]**2)
    A = 26*S**2 - Sq; B = 2*Sq; C = -Sq
    d = B**2 - 4*A*C
    t = (-B + np.sqrt(d))/(2*A)
    al = (1-t)/S
    out = np.zeros(4)
    out[:3] = m[:3]*al; out[3] = t
    return out

def full_step(m, e):
    return enforce_floor(ds_combine(m, e))

# Find equilibrium
def eqs(p):
    s1,th,w1,ph = p
    s2 = (1-s1-th)/2; w2 = (1-w1-ph)/2
    if min(s2,w2,th,ph,s1,w1) < 0: return [1e10]*4
    m = np.array([s1,s2,s2,th]); e = np.array([w1,w2,w2,ph])
    eq1 = th**2/(s1**2+2*s2**2+th**2) - 1/27
    eq2 = ph**2/(w1**2+2*w2**2+ph**2) - 1/27
    K = 2*s1*w2 + 2*s2*w1 + 2*s2*w2
    eq3 = K - 7/30
    mf = full_step(m, e)
    eq4 = mf[0] - s1
    return [eq1, eq2, eq3, eq4]

sol = fsolve(eqs, [0.787, 0.155, 0.631, 0.129])
s1,th,w1,ph = sol
s2 = (1-s1-th)/2; w2 = (1-w1-ph)/2
m_star = np.array([s1, s2, s2, th])
e_star = np.array([w1, w2, w2, ph])

print("="*70)
print("EQUILIBRIUM")
print("="*70)
print(f"m* = {m_star}")
print(f"e* = {e_star}")
print(f"Born = {born(m_star):.10f} (target {1/27:.10f})")
K_check = 2*s1*w2 + 2*s2*w1 + 2*s2*w2
print(f"K = {K_check:.10f} (target {7/30:.10f})")

# ============================================================
# SINGLE-SITE JACOBIAN
# ============================================================

def jacobian(func, x, eps=1e-9):
    n = len(x)
    J = np.zeros((n, n))
    for j in range(n):
        xp = x.copy(); xm = x.copy()
        xp[j] += eps; xm[j] -= eps
        J[:, j] = (func(xp) - func(xm)) / (2*eps)
    return J

J_single = jacobian(lambda m: full_step(m, e_star), m_star)

print(f"\nSingle-site Jacobian:")
for i in range(4):
    print(f"  [{' '.join(f'{J_single[i,j]:+.6f}' for j in range(4))}]")

evals, evecs = np.linalg.eig(J_single)
idx = np.argsort(-np.abs(evals))
evals = evals[idx]; evecs = evecs[:, idx]

print(f"\nEigenvalues: {[f'{abs(e):.6f}' for e in evals]}")
print(f"Eigenvectors:")
for i in range(4):
    if abs(evals[i]) > 1e-6:
        v = evecs[:, i].real
        print(f"  lam={abs(evals[i]):.6f}: ({v[0]:+.4f}, {v[1]:+.4f}, {v[2]:+.4f}, {v[3]:+.4f})")

# ============================================================
# RADIAL / ANGULAR BASIS
# ============================================================
print(f"\n{'='*70}")
print("RADIAL / ANGULAR DECOMPOSITION")
print("="*70)

# At the S2-symmetric equilibrium (s2=s3):
# Radial subspace: spanned by perturbations that preserve s2=s3
#   Basis vectors: e_s1 = (1,0,0,0), e_sym = (0,1,1,0)/sqrt(2), e_theta = (0,0,0,1)
# Angular subspace: spanned by perturbations that break s2=s3
#   Basis vector: e_ang = (0,1,-1,0)/sqrt(2)
# But we're on the L1 simplex, so sum(delta_m) = 0.
# The L1 constraint kills one radial direction.

# Effective basis (on L1 tangent space):
# Radial: e_r1 = (1,-1/2,-1/2,0)/sqrt(3/2) (increase s1, decrease s2,s3)
#          e_r2 = (0,0,0,1) ... but theta is enslaved!
# Angular: e_a = (0,1,-1,0)/sqrt(2)

# Since theta is enslaved (thm:enslaving), the effective dynamics is on
# the 2D section simplex. The basis is:
# e_r = radial direction (change s1 relative to s2,s3)
# e_a = angular direction (change s2 relative to s3)

# Let's use the eigenvector directions directly
v_rad = evecs[:, 0].real  # leading eigenvalue = radial
v_ang = np.array([0, 1, -1, 0]) / np.sqrt(2)  # angular by construction

# Normalize
v_rad = v_rad / np.linalg.norm(v_rad)

print(f"Radial eigenvector: ({v_rad[0]:+.4f}, {v_rad[1]:+.4f}, {v_rad[2]:+.4f}, {v_rad[3]:+.4f})")
print(f"Angular direction:  ({v_ang[0]:+.4f}, {v_ang[1]:+.4f}, {v_ang[2]:+.4f}, {v_ang[3]:+.4f})")
print(f"Orthogonality: {np.dot(v_rad, v_ang):.10f}")

# Project the Jacobian into the radial/angular basis
# J_rr = v_rad^T J v_rad (radial -> radial)
# J_ra = v_rad^T J v_ang (angular -> radial)
# J_ar = v_ang^T J v_rad (radial -> angular)
# J_aa = v_ang^T J v_ang (angular -> angular)

J_rr = v_rad @ J_single @ v_rad
J_ra = v_rad @ J_single @ v_ang
J_ar = v_ang @ J_single @ v_rad
J_aa = v_ang @ J_single @ v_ang

print(f"\nJacobian in radial/angular basis:")
print(f"  J_rr (radial->radial) = {J_rr:.8f}")
print(f"  J_aa (angular->angular) = {J_aa:.8f}")
print(f"  J_ra (angular->radial) = {J_ra:.8f}")
print(f"  J_ar (radial->angular) = {J_ar:.8f}")
print(f"  |coupling| = sqrt(J_ra*J_ar) = {np.sqrt(abs(J_ra*J_ar)):.8f}")

# ============================================================
# WHAT IS THE COUPLING?
# ============================================================
print(f"\n{'='*70}")
print("THE COUPLING BETWEEN SECTORS")
print("="*70)

print(f"J_ra (how angular perturbations affect radial): {J_ra:.10f}")
print(f"J_ar (how radial perturbations affect angular): {J_ar:.10f}")
print(f"Product J_ra * J_ar = {J_ra * J_ar:.10f}")
print()

# Compare to framework quantities
print("Compare J_ar to framework quantities:")
print(f"  J_ar = {J_ar:.8f}")
print(f"  1/10 = {1/10:.8f}")
print(f"  1/12 = {1/12:.8f}")
print(f"  11/60 = {11/60:.8f}")
print(f"  K*/2 = {7/60:.8f}")
print(f"  1/H = {1/3:.8f}")
print(f"  K* = {7/30:.8f}")
print(f"  1/27 = {1/27:.8f}")
print()

# ============================================================
# A2 COUPLED SYSTEM
# ============================================================
print("="*70)
print("A2 COUPLED JACOBIAN BLOCK STRUCTURE")
print("="*70)

g = 7/30
uniform = np.array([0.25]*4)

def coupled_ev(m1, m2, e0):
    e1 = e0 - g*(m2 - uniform)
    e1 = np.maximum(e1, 1e-15)
    return e1/np.sum(e1)

# Find coupled equilibrium
st = np.concatenate([m_star, m_star])
for i in range(10000):
    m1, m2 = st[:4], st[4:]
    e1 = coupled_ev(m1, m2, e_star)
    e2 = coupled_ev(m2, m1, e_star)
    m1f = full_step(m1, e1)
    m2f = full_step(m2, e2)
    sn = np.concatenate([m1f, m2f])
    if np.linalg.norm(sn-st) < 1e-15: break
    st = sn

m1_eq = st[:4]; m2_eq = st[4:]
print(f"m1* = {m1_eq}")

# 8x8 Jacobian
def coupled_map(s):
    mm1, mm2 = s[:4], s[4:]
    ee1 = coupled_ev(mm1, mm2, e_star)
    ee2 = coupled_ev(mm2, mm1, e_star)
    return np.concatenate([full_step(mm1, ee1), full_step(mm2, ee2)])

J8 = jacobian(coupled_map, st, eps=1e-9)

# Project into radial/angular basis for each node
# Node 1: v_rad1 = (v_rad, 0000), v_ang1 = (v_ang, 0000)
# Node 2: v_rad2 = (0000, v_rad), v_ang2 = (0000, v_ang)

v_rad1 = np.concatenate([v_rad, np.zeros(4)])
v_ang1 = np.concatenate([v_ang, np.zeros(4)])
v_rad2 = np.concatenate([np.zeros(4), v_rad])
v_ang2 = np.concatenate([np.zeros(4), v_ang])

# Node-symmetric and node-antisymmetric combinations
v_rad_sym = (v_rad1 + v_rad2) / np.sqrt(2)
v_rad_asym = (v_rad1 - v_rad2) / np.sqrt(2)
v_ang_sym = (v_ang1 + v_ang2) / np.sqrt(2)
v_ang_asym = (v_ang1 - v_ang2) / np.sqrt(2)

# Project J8 into these four directions
basis = {
    'rad_asym': v_rad_asym,   # should give lambda_0
    'ang_asym': v_ang_asym,   # should give lambda_1
    'ang_sym': v_ang_sym,     # should give lambda_2
    'rad_sym': v_rad_sym,     # should give lambda_3
}

print(f"\nA2 Jacobian projected into eigenmode basis:")
print(f"{'':15s}", end='')
for name in basis:
    print(f"  {name:>12s}", end='')
print()

for name_i, v_i in basis.items():
    print(f"{name_i:15s}", end='')
    for name_j, v_j in basis.items():
        val = v_i @ J8 @ v_j
        print(f"  {val:+12.6f}", end='')
    print()

# The diagonal should give the four eigenvalues
print(f"\nDiagonal elements (should match eigenvalues):")
for name, v in basis.items():
    diag = v @ J8 @ v
    print(f"  {name}: {diag:.6f}")

print(f"\nExpected: lam0=0.5022, lam1=0.4745, lam2=0.3527, lam3=0.3344")

# ============================================================
# THE KEY: RADIAL-ANGULAR COUPLING IN THE A2 SYSTEM
# ============================================================
print(f"\n{'='*70}")
print("RADIAL-ANGULAR COUPLING MATRIX (A2)")
print("="*70)

# The 4x4 submatrix coupling radial to angular
# [rad_asym->rad_asym  rad_asym->ang_asym  rad_asym->ang_sym  rad_asym->rad_sym]
# [ang_asym->rad_asym  ...                 ...                ...               ]
# etc.

modes = ['rad_asym (e0)', 'ang_asym (e1)', 'ang_sym (e2)', 'rad_sym (e3)']
vecs = [v_rad_asym, v_ang_asym, v_ang_sym, v_rad_sym]

print(f"\n{'':20s}", end='')
for m in modes:
    print(f" {m:>15s}", end='')
print()

coupling_matrix = np.zeros((4,4))
for i, (name_i, v_i) in enumerate(zip(modes, vecs)):
    print(f"{name_i:20s}", end='')
    for j, (name_j, v_j) in enumerate(zip(modes, vecs)):
        val = v_i @ J8 @ v_j
        coupling_matrix[i,j] = val
        print(f" {val:+15.8f}", end='')
    print()

# Radial-angular off-diagonal blocks
print(f"\nOff-diagonal (radial-angular) couplings:")
print(f"  e0->e1 (rad_asym->ang_asym): {coupling_matrix[0,1]:.10f}")
print(f"  e0->e2 (rad_asym->ang_sym):  {coupling_matrix[0,2]:.10f}")
print(f"  e3->e1 (rad_sym->ang_asym):  {coupling_matrix[3,1]:.10f}")
print(f"  e3->e2 (rad_sym->ang_sym):   {coupling_matrix[3,2]:.10f}")
print(f"  e1->e0 (ang_asym->rad_asym): {coupling_matrix[1,0]:.10f}")
print(f"  e2->e0 (ang_sym->rad_asym):  {coupling_matrix[2,0]:.10f}")
print(f"  e1->e3 (ang_asym->rad_sym):  {coupling_matrix[1,3]:.10f}")
print(f"  e2->e3 (ang_sym->rad_sym):   {coupling_matrix[2,3]:.10f}")

# Total radial->angular coupling
total_ra = abs(coupling_matrix[0,1]) + abs(coupling_matrix[0,2]) + abs(coupling_matrix[3,1]) + abs(coupling_matrix[3,2])
total_ar = abs(coupling_matrix[1,0]) + abs(coupling_matrix[2,0]) + abs(coupling_matrix[1,3]) + abs(coupling_matrix[2,3])

print(f"\nTotal |radial->angular|: {total_ra:.8f}")
print(f"Total |angular->radial|: {total_ar:.8f}")
print(f"Geometric mean: {np.sqrt(total_ra*total_ar):.8f}")
print(f"Compare to 11/60 = {11/60:.8f}")
print(f"Compare to 1/10 = {1/10:.8f}")
print(f"Compare to 1/12 = {1/12:.8f}")

# Frobenius norm of the off-diagonal block
off_diag = np.array([
    [coupling_matrix[0,1], coupling_matrix[0,2]],
    [coupling_matrix[3,1], coupling_matrix[3,2]]
])
off_diag_T = np.array([
    [coupling_matrix[1,0], coupling_matrix[2,0]],
    [coupling_matrix[1,3], coupling_matrix[2,3]]
])

print(f"\nFrobenius norm (radial->angular block): {np.linalg.norm(off_diag):.8f}")
print(f"Frobenius norm (angular->radial block): {np.linalg.norm(off_diag_T):.8f}")
print(f"Product of norms: {np.linalg.norm(off_diag)*np.linalg.norm(off_diag_T):.8f}")
print(f"Sqrt of product: {np.sqrt(np.linalg.norm(off_diag)*np.linalg.norm(off_diag_T)):.8f}")
