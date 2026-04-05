"""
SELF-CONTAINED COMPUTATION: Second-order Birkhoff factorisation at the
DS equilibrium K*=7/30, for the Ward-reconstructed gauge connection.

PURPOSE: Compute the O(epsilon^2) correction to h+ in the perturbative
Birkhoff expansion, extract the non-abelian F+ content, and determine
whether the commutator content 8332/625 connects to the physical
condensate <Tr(F^2_Ward)> through a computable ratio.

BACKGROUND:
- The DS framework on CP^3 with Born floor produces Yang-Mills via Mason.
- The physical connection is A = h+^{-1} d h+ (Ward reconstruction), NOT
  M^{-1}dM (Maurer-Cartan, identically flat).
- At the K*=7/30 equilibrium: the floor correction epsilon = M(m*) - M(DS(m*,e*))
  satisfies ||epsilon/M_hol|| = 0.194 (small enough for perturbation theory).
- Leading-order result (COMPUTED):
    F_Ward = 0 at uniform equilibrium
    <|F+|^2> = 3.743 * sigma_eff^2/a^2 from fluctuations
    98.9% abelian, 1.1% non-abelian at O(epsilon)
    No clean algebraic ratios
- The genuinely non-abelian content enters at O(epsilon^2).

THE COMPUTATION:
The perturbative Birkhoff factorisation on each twistor line L_x ~ CP^1:
  M_tilde(zeta, zeta_bar) = M_0 * (I + delta)
  where delta = M_0^{-1} * epsilon_eff(zeta, zeta_bar) is the relative correction.

At first order:
  h+^{(1)} = M_0 * (I + P+[delta])
  h-^{(1)} = I + P-[delta]

At second order:
  h+^{(2)} = M_0 * (I + P+[delta] + P+[delta * P-[delta]] - P+[delta]*P-[delta])
  (from the multiplicative splitting M = h+^{-1} h- expanded to second order)

The physical connection at second order:
  A^{(2)} = h+^{(2),-1} * d_x h+^{(2)}

The curvature F = dA + A^2 now includes:
  - dA term: O(epsilon) [from spatial variation of first-order correction]
  - [A,A] term: O(epsilon^2) [from commutator of first-order corrections]
  The [A,A] term is where the non-abelian content lives.

WHAT TO COMPUTE:
1. The equilibrium (m*, e*) at 100-digit precision
2. The floor correction epsilon and relative correction delta = M_0^{-1} epsilon
3. The Cauchy projectors P+[delta] and P-[delta] on CP^1
4. The second-order correction P+[delta * P-[delta]]
5. The non-abelian fraction of F+ at O(epsilon^2)
6. The ratio <Tr(F^2_Ward)> / (8332/625 * sigma_eff^4 / det(M*)^2)
   (connecting the physical condensate to the commutator content)

TECHNICAL NOTES:
- The Cauchy projector on CP^1: P+[f](zeta) = (1/2pi i) oint f(w)/(w-zeta) dw
  where the contour is |w| = 1 (equator) and |zeta| < 1.
- On the equator, zeta_bar = 1/zeta. A function f(zeta, 1/zeta) has Laurent
  expansion sum_n a_n zeta^n. P+ keeps n >= 0, P- keeps n < 0.
- The incidence relation makes the floor correction's spatial derivative
  a LINEAR function of zeta on L_x. At O(epsilon): only the 1/zeta term
  contributes to F+ (via the residue). At O(epsilon^2): zeta^{-2}, zeta^{-1},
  zeta^0 terms all contribute, giving richer structure.

PRECISION: Use mpmath at 100 digits throughout.
"""

from mpmath import mp, mpf, sqrt, matrix, eye, pi, exp, nstr, fabs, mpc

mp.dps = 100
ONE = mpf(1); ZERO = mpf(0)
H = mpf(3); FLOOR = ONE / H**3; TARGET_K = mpf(7) / mpf(30)

# Pauli matrices
sigma1 = matrix([[0, 1], [1, 0]])
sigma2 = matrix([[0, -1j], [1j, 0]])
sigma3 = matrix([[1, 0], [0, -1]])
I2 = eye(2)
basis = [sigma1, sigma2, sigma3, I2]
sq2 = sqrt(mpf(2))

def to_matrix(m):
    s1, s2, s3, th = m
    return (th*I2 + s1*sigma1 + s2*sigma2 + s3*sigma3) / sq2

def enforce_floor(m):
    s1, s2, s3, th = m
    ssq = s1**2 + s2**2 + s3**2
    born = th**2 / (ssq + th**2)
    if born >= FLOOR: return [s1, s2, s3, th]
    ss = s1+s2+s3; sq = ssq; r = sq/ss**2
    disc = (2*r)**2 + 4*(26-r)*r
    t = (-2*r + sqrt(disc)) / (2*(26-r))
    alpha = (ONE-t)/ss
    return [s1*alpha, s2*alpha, s3*alpha, t]

def ds_combine(m, e):
    s1,s2,s3,th = m; e1,e2,e3,ph = e
    sn = [s1*e1+s1*ph+th*e1, s2*e2+s2*ph+th*e2, s3*e3+s3*ph+th*e3]
    tn = th*ph
    K = s1*e2+s1*e3+s2*e1+s2*e3+s3*e1+s3*e2
    d = ONE-K
    return [sn[0]/d, sn[1]/d, sn[2]/d, tn/d], K

def full_step(m, e):
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds), K

def make_evidence(p_dom):
    p_w = (ONE-p_dom)/2; sc = ONE-FLOOR
    raw = [sqrt(p_dom*sc), sqrt(p_w*sc), sqrt(p_w*sc), sqrt(FLOOR)]
    tot = sum(raw); return [r/tot for r in raw]

# ============================================================
# TASK A: Find equilibrium at 100-digit precision
# ============================================================
print("TASK A: Finding equilibrium...")

p0, p1 = mpf('0.932'), mpf('0.933')
def K_at_p(p):
    e = make_evidence(p)
    m = [mpf('0.4'), mpf('0.15'), mpf('0.15'), mpf('0.3')]
    for _ in range(10000): m, _ = full_step(m, e)
    _, K = ds_combine(m, e); return K

K0, K1 = K_at_p(p0), K_at_p(p1)
for step in range(100):
    f0, f1 = K0-TARGET_K, K1-TARGET_K
    if fabs(f1-f0) < mpf(10)**(-mp.dps): break
    p2 = p1 - f1*(p1-p0)/(f1-f0)
    p0, K0 = p1, K1; p1, K1 = p2, K_at_p(p2)
    if fabs(p1-p0) < mpf(10)**(-(mp.dps-20)): break

e_star = make_evidence(p1)
m = [mpf('0.4'), mpf('0.15'), mpf('0.15'), mpf('0.3')]
for _ in range(10000): m, _ = full_step(m, e_star)
m_star = m
m_ds, K_star = ds_combine(m_star, e_star)
print(f"  |K*-7/30| = {nstr(fabs(K_star-TARGET_K), 5)}")

# ============================================================
# TASK B: Compute all Jacobians and the floor correction
# ============================================================
print("\nTASK B: Jacobians and floor correction...")

s1,s2,s3,th = m_star; e1,e2,e3,ph = e_star
omK = ONE - K_star
S_v = [s1*(e1+ph)+th*e1, s2*(e2+ph)+th*e2, s3*(e3+ph)+th*e3, th*ph]
N = [S_v[i]/omK for i in range(4)]

# DS Jacobian
dSdm = [[ZERO]*4 for _ in range(4)]
for i in range(3): dSdm[i][i] = e_star[i]+ph; dSdm[i][3] = e_star[i]
dSdm[3][3] = ph
dKdm = [e2+e3, e1+e3, e1+e2, ZERO]
J_DS = matrix(4,4)
for i in range(4):
    for j in range(4):
        J_DS[i,j] = (dSdm[i][j]*omK + S_v[i]*dKdm[j]) / omK**2

# Floor Jacobian
Ss = N[0]+N[1]+N[2]; Sq = N[0]**2+N[1]**2+N[2]**2; R = Sq/Ss**2
u = sqrt(26*R); t_v = (u-R)/(26-R); g = (ONE-t_v)/Ss
dtdR = (338+13*R-26*u)/(u*(26-R)**2)
dRdN = [2*(N[i]*Ss-Sq)/Ss**3 for i in range(3)]
dtdN = [dtdR*dRdN[i] for i in range(3)]
dgdN = [(-dtdN[i]*Ss-(ONE-t_v))/Ss**2 for i in range(3)]
J_fl = matrix(4,4)
for i in range(3):
    for j in range(3): J_fl[i,j] = (g if i==j else ZERO) + N[i]*dgdN[j]
J_fl[3,0]=dtdN[0]; J_fl[3,1]=dtdN[1]; J_fl[3,2]=dtdN[2]; J_fl[3,3]=ZERO

J_corr = (J_fl - eye(4)) * J_DS  # Jacobian of floor correction

M_star_mat = to_matrix(m_star)
M_ds_mat = to_matrix(m_ds)
M_ds_inv = M_ds_mat**(-1)
epsilon_mat = M_star_mat - M_ds_mat
delta_mat = M_ds_inv * epsilon_mat  # relative correction

det_M_star = M_star_mat[0,0]*M_star_mat[1,1] - M_star_mat[0,1]*M_star_mat[1,0]
det_M_ds = M_ds_mat[0,0]*M_ds_mat[1,1] - M_ds_mat[0,1]*M_ds_mat[1,0]

print(f"  ||epsilon/M_0|| = {nstr(sqrt(sum(fabs(delta_mat[i,j])**2 for i in range(2) for j in range(2))), 15)}")

# ============================================================
# TASK C: Ward connection basis and incidence relation
# ============================================================
print("\nTASK C: Ward connection and incidence relation...")

B = [matrix(2,2) for _ in range(4)]
for mu in range(4):
    for j in range(4):
        B[mu] = B[mu] + M_ds_inv * basis[j] * (J_corr[j,mu] / sq2)

# Incidence: dz_bar_i/dx^mu -> A_i[mu] (constant) + C_i[mu]/zeta
A1 = [ZERO]*4; C1 = [ZERO]*4; A2 = [ZERO]*4; C2 = [ZERO]*4
isq2 = mpc(0, -1)/sq2  # -i/sqrt(2)
A1[0] = isq2;           A1[3] = isq2
C1[1] = isq2;           C1[2] = ONE/sq2
A2[1] = isq2;           A2[2] = ONE/sq2
C2[0] = isq2;           C2[3] = -isq2  # +i/sqrt(2) for mu=3

# ============================================================
# TASK D: SECOND-ORDER Birkhoff factorisation
# ============================================================
print("\nTASK D: Second-order Birkhoff factorisation via Penrose integral...")

# At first order: the integrand on L_0 is
#   rho(zeta) = sum_mu B[mu] * (A1[mu]+A2[mu] + (C1[mu]+C2[mu])/zeta)
# P+[rho] picks up the constant term; P-[rho] picks up 1/zeta term.
# F+^(1) = residue = rho_{-1} = sum_mu B[mu] * (C1[mu]+C2[mu])

rho_0 = matrix(2,2)   # constant term (in P+)
rho_m1 = matrix(2,2)  # 1/zeta term (the residue -> F+ at O(epsilon))
for mu in range(4):
    rho_0 = rho_0 + B[mu] * (A1[mu]+A2[mu])
    rho_m1 = rho_m1 + B[mu] * (C1[mu]+C2[mu])

F1_plus = rho_m1  # First-order F+
F1_sq = sum(fabs(F1_plus[i,j])**2 for i in range(2) for j in range(2))
print(f"  |F+^(1)|^2 = {nstr(F1_sq, 30)}")

# At second order: the [A,A] contribution to F.
# The Ward connection at first order has components:
#   A_mu(zeta) = B[mu] on the base (from spatial derivatives)
# The curvature F_{mu,nu} = [A_mu, A_nu] at the SAME point
# gives the commutator contribution.
#
# But through the incidence relation, each A_mu has zeta-dependent
# coefficients (from the projection of spacetime derivatives onto
# twistor-line coordinates). The FULL connection on L_0 is:
#   A(zeta) = sum_mu B[mu] * c_mu(zeta)
# where c_mu(zeta) = (A1[mu]+A2[mu]) + (C1[mu]+C2[mu])/zeta
#
# The [A,A] piece of the curvature involves:
#   [A, A] = sum_{mu<nu} [c_mu * B_mu, c_nu * B_nu]
#          = sum_{mu<nu} c_mu * c_nu * [B_mu, B_nu]
#
# The product c_mu * c_nu has Laurent terms from zeta^{-2} to zeta^0.
# The Penrose integral picks up the 1/zeta coefficient.

# Compute [B_mu, B_nu] for all pairs
print("\n  Commutator structure [B_mu, B_nu]:")
comm = [[None]*4 for _ in range(4)]
for mu in range(4):
    for nu in range(mu+1, 4):
        comm[mu][nu] = B[mu]*B[nu] - B[nu]*B[mu]
        nrm = sqrt(sum(fabs(comm[mu][nu][i,j])**2 for i in range(2) for j in range(2)))
        if nrm > mpf('1e-50'):
            print(f"    [B_{mu}, B_{nu}]: ||comm|| = {nstr(nrm, 15)}")

# The [A,A] contribution to F+:
# For each (mu, nu) pair, c_mu(zeta)*c_nu(zeta) has Laurent expansion:
#   c_mu * c_nu = (A_mu+C_mu/z)(A_nu+C_nu/z)
#              = A_mu*A_nu + (A_mu*C_nu + C_mu*A_nu)/z + C_mu*C_nu/z^2
# where A_mu = A1[mu]+A2[mu], C_mu = C1[mu]+C2[mu]
#
# The residue (1/z coefficient) of c_mu*c_nu is: A_mu*C_nu + C_mu*A_nu
# The 1/z^2 coefficient is: C_mu*C_nu (contributes to zeta * [A,A], not F+_{0'0'})

F2_AA_00 = matrix(2,2)  # [A,A] contribution to F+_{0'0'} (residue)
F2_AA_01 = matrix(2,2)  # [A,A] contribution to F+_{0'1'} (from 1/z^2 term * zeta)

for mu in range(4):
    for nu in range(mu+1, 4):
        if comm[mu][nu] is None: continue
        Am = A1[mu]+A2[mu]; Cm = C1[mu]+C2[mu]
        An = A1[nu]+A2[nu]; Cn = C1[nu]+C2[nu]

        # Residue of c_mu*c_nu at z=0: A_mu*C_nu + C_mu*A_nu
        res = Am*Cn + Cm*An
        F2_AA_00 = F2_AA_00 + comm[mu][nu] * res

        # 1/z^2 coefficient: C_mu * C_nu (contributes to F+_{0'1'} via zeta*...)
        coeff_m2 = Cm * Cn
        F2_AA_01 = F2_AA_01 + comm[mu][nu] * coeff_m2

F2_AA_sq = sum(fabs(F2_AA_00[i,j])**2 for i in range(2) for j in range(2))
print(f"\n  |[A,A]| contribution to F+_{{0'0'}}: {nstr(F2_AA_sq, 20)}")
print(f"  |[A,A]| contribution to F+_{{0'1'}}: {nstr(sum(fabs(F2_AA_01[i,j])**2 for i in range(2) for j in range(2)), 20)}")

# Total F+ at second order: F+^(1) + F+^(2)_{[A,A]}
F_total_00 = F1_plus + F2_AA_00
F_total_sq = sum(fabs(F_total_00[i,j])**2 for i in range(2) for j in range(2))

print(f"\n  |F+_total|^2 (first + second order) = {nstr(F_total_sq, 30)}")
print(f"  |F+^(1)|^2                          = {nstr(F1_sq, 30)}")
print(f"  Second-order correction:             = {nstr(F_total_sq - F1_sq, 15)}")
print(f"  Relative correction:                 = {nstr((F_total_sq - F1_sq)/F1_sq, 10)}")

# Non-abelian fraction at second order
tr_F = (F_total_00[0,0] + F_total_00[1,1]) / 2
su2_comps = [(F_total_00[0,1]+F_total_00[1,0])/2,
             (F_total_00[0,1]-F_total_00[1,0])/(2j),
             (F_total_00[0,0]-F_total_00[1,1])/2]
scalar_sq = fabs(tr_F)**2
gauge_sq = sum(fabs(c)**2 for c in su2_comps)
total_content = scalar_sq + gauge_sq

print(f"\n  su(2) decomposition (with second-order):")
print(f"    Scalar: {nstr(scalar_sq/total_content, 15)}")
print(f"    Gauge:  {nstr(gauge_sq/total_content, 15)}")

# ============================================================
# TASK E: Ratio to the commutator content 8332/625
# ============================================================
print(f"\n{'='*70}")
print("TASK E: Connection to the commutator content 8332/625")
print(f"{'='*70}")

# The commutator content from thm:condensate:
# C * det(M*)^2 = 8332/625 = 4(3*26^2 + 2*26 + 3)/25^2
# C = (8332/625) / det(M*)^2
det_M_star_sq = fabs(det_M_star)**2
C_mc = mpf(8332)/mpf(625) / det_M_star_sq
sigma_eff_sq = ONE / (ONE - mpf('0.28291')**2)
sigma_eff_4 = sigma_eff_sq**2

mc_condensate = C_mc * sigma_eff_4
ward_condensate_1 = F1_sq * sigma_eff_sq  # first order (times sigma^2/a^2)
ward_condensate_2 = F_total_sq * sigma_eff_sq

print(f"\nCommutator content C = {nstr(C_mc, 20)}")
print(f"MC <|[a,a]|^2> = C * sigma_eff^4 = {nstr(mc_condensate, 20)}")
print(f"\nWard <|F+|^2> at O(eps):   {nstr(F1_sq, 20)} * sigma_eff^2/a^2")
print(f"Ward <|F+|^2> at O(eps^2): {nstr(F_total_sq, 20)} * sigma_eff^2/a^2")

# Key ratio: Ward / MC
ratio_1 = F1_sq / C_mc
ratio_2 = F_total_sq / C_mc
print(f"\nRatio Ward(O1)/MC = {nstr(ratio_1, 20)}")
print(f"Ratio Ward(O2)/MC = {nstr(ratio_2, 20)}")

# Check for clean algebraic form of the ratio
print(f"\nAlgebraic checks on Ward(O2)/MC:")
val = ratio_2
for p in range(1, 100):
    for q in range(1, 200):
        frac = mpf(p)/mpf(q)
        if fabs(val - frac) < mpf('0.0005'):
            print(f"  {p}/{q} = {nstr(frac, 15)}, error = {nstr(fabs(val-frac), 10)}")

# Check |F+_total|^2 * det(M*)^2
ratio_det = F_total_sq * det_M_star_sq
print(f"\n|F+_total|^2 * det(M*)^2 = {nstr(ratio_det, 20)}")
print(f"Nearby: 1/3 = {nstr(ONE/3, 15)}, error = {nstr(fabs(ratio_det - ONE/3), 10)}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
