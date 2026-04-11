"""
Second-order Koopman eigenfunction computation.
Verifies whether the nonlinear correction closes the 115 ppm gap
in the VP weight for 1/alpha.
"""
from mpmath import mp, mpf, matrix, sqrt, fabs, nstr, eig, chop, findroot, eye

mp.dps = 40
H = 3
FLOOR = mpf(1)/mpf(H**3)

def ds_combine(m, e):
    s = m[:3]; theta = m[3]; ev = e[:3]; phi = e[3]
    s_pre = [s[i]*ev[i] + s[i]*phi + theta*ev[i] for i in range(3)]
    theta_pre = theta * phi
    total_pre = sum(s_pre) + theta_pre
    K = mpf(1) - total_pre
    denom = mpf(1) - K
    return [sp/denom for sp in s_pre] + [theta_pre/denom], K

def born_prob(m):
    L2sq = sum(x**2 for x in m)
    if L2sq < mpf(10)**(-60): return mpf(0)
    return m[3]**2 / L2sq

def enforce_floor(m):
    b = born_prob(m)
    if b >= FLOOR - mpf(10)**(-60): return list(m)
    S = sum(m[:3])
    Sq = sum(x**2 for x in m[:3])
    A_c = mpf(26)*S**2 - Sq
    B_c = mpf(2)*Sq; C_c = -Sq
    disc = B_c**2 - 4*A_c*C_c
    t1 = (-B_c + sqrt(disc))/(2*A_c)
    t2 = (-B_c - sqrt(disc))/(2*A_c)
    cands = [t for t in [t1,t2] if mpf(0) < t < mpf(1)]
    t = min(cands, key=lambda x: fabs(x - m[3]))
    alpha_r = (mpf(1)-t)/S
    return [m[i]*alpha_r for i in range(3)] + [t]

def phi_map(m, e):
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds)

# Find equilibrium
def eq_system(*p):
    s1, theta, w1, phi = p
    s2 = (mpf(1)-s1-theta)/2
    w2 = (mpf(1)-w1-phi)/2
    m = [s1, s2, s2, theta]
    e = [w1, w2, w2, phi]
    eq1 = theta**2/(s1**2+2*s2**2+theta**2) - FLOOR
    eq2 = phi**2/(w1**2+2*w2**2+phi**2) - FLOOR
    K = 2*s1*w2 + 2*s2*(w1+w2)
    m_out = phi_map(m, e)
    return [eq1, eq2, K - mpf(7)/30, m_out[0] - s1]

print("Finding fixed point...")
sol = findroot(eq_system, [mpf('0.787'), mpf('0.155'), mpf('0.631'), mpf('0.128')])
s1s, ths, w1s, phis = sol
s2s = (mpf(1)-s1s-ths)/2
w2s = (mpf(1)-w1s-phis)/2
m_star = [s1s, s2s, s2s, ths]
e_star = [w1s, w2s, w2s, phis]
print(f"m* = ({nstr(s1s,12)}, {nstr(s2s,12)}, {nstr(s2s,12)}, {nstr(ths,12)})")

# Jacobian
print("Computing Jacobian...")
eps_j = mpf(10)**(-20)
J = matrix(4,4)
for j in range(4):
    mp_v = list(m_star); mm_v = list(m_star)
    mp_v[j] += eps_j; mm_v[j] -= eps_j
    fp = phi_map(mp_v, e_star)
    fm = phi_map(mm_v, e_star)
    for i in range(4):
        J[i,j] = (fp[i]-fm[i])/(2*eps_j)

evals_raw, evecs_raw = eig(J)
evals = [chop(e) for e in evals_raw]
idx_sorted = sorted(range(4), key=lambda k: -float(fabs(evals[k])))
lam0 = evals[idx_sorted[0]]
k0 = idx_sorted[0]

v0 = [chop(evecs_raw[i,k0]) for i in range(4)]
norm_v0 = sqrt(sum(x**2 for x in v0))
v0 = [x/norm_v0 for x in v0]
if float(v0[0].real) < 0: v0 = [-x for x in v0]

# Left eigenvector
JT = J.T
evals_T, evecs_T = eig(JT)
evals_T_clean = [chop(e) for e in evals_T]
k0_left = min(range(4), key=lambda k: float(fabs(evals_T_clean[k] - lam0)))
L0 = [chop(evecs_T[i,k0_left]) for i in range(4)]
dot_Lv = sum(L0[i]*v0[i] for i in range(4))
L0 = [x/dot_Lv for x in L0]

print(f"lam0 = {nstr(lam0, 25)}")
print(f"v0 = ({', '.join(nstr(x,10) for x in v0)})")
print(f"L0 = ({', '.join(nstr(x,10) for x in L0)})")

# Hessian
print("Computing Hessian (64 finite differences)...")
eps_h = mpf(10)**(-12)
Hess = [matrix(4,4) for _ in range(4)]
for i in range(4):
    for j in range(i, 4):
        mpp = list(m_star); mpm = list(m_star)
        mmp = list(m_star); mmm = list(m_star)
        mpp[i] += eps_h; mpp[j] += eps_h
        mpm[i] += eps_h; mpm[j] -= eps_h
        mmp[i] -= eps_h; mmp[j] += eps_h
        mmm[i] -= eps_h; mmm[j] -= eps_h
        fpp = phi_map(mpp, e_star)
        fpm = phi_map(mpm, e_star)
        fmp = phi_map(mmp, e_star)
        fmm = phi_map(mmm, e_star)
        for k in range(4):
            h_val = (fpp[k] - fpm[k] - fmp[k] + fmm[k]) / (4*eps_h**2)
            Hess[k][i,j] = h_val
            Hess[k][j,i] = h_val

# Homological equation: (J - lam0^2 I) v00 = -H(v0,v0)
print("Solving homological equation for v00...")
Hv0v0 = [sum(Hess[k][i,j]*v0[i]*v0[j] for i in range(4) for j in range(4)) for k in range(4)]
A_hom = J - lam0**2 * eye(4)
rhs_hom = matrix([[-x] for x in Hv0v0])
v00_mat = A_hom**(-1) * rhs_hom
v00 = [v00_mat[i,0] for i in range(4)]
print(f"v00 = ({', '.join(nstr(x,10) for x in v00)})")

# Sylvester equation: J^T A J - lam0 A = B
print("Solving Sylvester equation for A...")
B_mat = matrix(4,4)
for i in range(4):
    for j in range(4):
        B_mat[i,j] = -sum(L0[k]*Hess[k][i,j] for k in range(4))

M16 = matrix(16,16)
for i in range(4):
    for j in range(4):
        row = i*4 + j
        for p in range(4):
            for q in range(4):
                col = p*4 + q
                M16[row,col] = J[p,i]*J[q,j]
        M16[row,row] -= lam0

b16 = matrix(16,1)
for i in range(4):
    for j in range(4):
        b16[i*4+j,0] = B_mat[i,j]

a16 = M16**(-1) * b16
A_koop = matrix(4,4)
for i in range(4):
    for j in range(4):
        A_koop[i,j] = a16[i*4+j,0]

print("A_koopman diagonal:", [nstr(A_koop[i,i],8) for i in range(4)])

# Charge observable
c_vec = [mpf(1), mpf(-1), mpf(-1), mpf(0)]
W0 = sum(c_vec[i]*v0[i] for i in range(4))
W00 = mpf('0.5') * sum(c_vec[i]*v00[i] for i in range(4))

# What C_eff is needed?
w_exact = mpf('0.10198826550604883023')
C_eff_needed = sqrt(w_exact * 27)
Delta_needed = C_eff_needed - W0

print()
print("=" * 55)
print("RESULTS")
print("=" * 55)
print(f"W0 (linear charge amplitude) = {nstr(W0, 15)}")
print(f"W00 (quadratic charge coeff) = {nstr(W00, 15)}")
print(f"C_eff needed for exact w     = {nstr(C_eff_needed, 15)}")
print(f"Delta needed                 = {nstr(Delta_needed, 10)}")
print()
print(f"w_lin = W0^2/27              = {nstr(W0**2/27, 15)}")
print(f"w_exact                      = {nstr(w_exact, 15)}")
print()

# The correction from the curved manifold:
# On the invariant manifold m(u) = m* + u*v0 + (u^2/2)*v00
# The charge is C(u) = W0*u + W00*u^2
# The Koopman eigenfunction is psi0(u) = u (to leading order on the manifold)
# The spectral weight involves <C^2>/<psi0^2> over the invariant measure
# For u with variance sigma^2 = 1/(1-lam0^2):
sigma2 = mpf(1)/(1-lam0**2)
print(f"sigma^2 = 1/(1-lam0^2) = {nstr(sigma2, 12)}")

# <C(u)^2> = W0^2 <u^2> + 2*W0*W00 <u^3> + W00^2 <u^4>
# For symmetric distribution: <u^3> = 0, <u^4> = 3*sigma^4
# <C^2> = W0^2*sigma^2 + W00^2*3*sigma^4
# <psi0^2> = <u^2> = sigma^2
# w = <C^2>/(H^3 * <psi0^2>) = (W0^2 + 3*W00^2*sigma^2) / H^3

w_corrected = (W0**2 + 3*W00**2*sigma2) / 27
print(f"w with Gaussian u^4 correction = {nstr(w_corrected, 15)}")
print(f"gap from exact: {nstr(fabs(w_corrected - w_exact), 6)}")
print(f"gap ppm: {nstr(fabs(w_corrected/w_exact - 1)*mpf(10**6), 6)}")
print()

# Alternative: the effective amplitude approach
# C_eff = sqrt(W0^2 + 3*W00^2*sigma^2)
C_eff_gauss = sqrt(W0**2 + 3*W00**2*sigma2)
w_gauss = C_eff_gauss**2 / 27
print(f"C_eff (Gaussian) = {nstr(C_eff_gauss, 15)}")
print(f"w (Gaussian)     = {nstr(w_gauss, 15)}")
print(f"Compared to C_eff needed: {nstr(C_eff_needed, 15)}")
