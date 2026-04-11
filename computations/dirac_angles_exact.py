"""
Dirac neutrino Koide angles: close the 3-4% gap.

From diagonal seesaw: m_Dk = sqrt(m_k * M_R), so m_Dk ∝ sqrt(m_k).
If m_k = M^2 * u_k^2 with u_k = 1 + r*tan(theta + 2*pi*k/3),
then m_Dk ∝ |u_k|.

Koide ratio Q_D = (sum |u_k|^{1/2})^2 / (3 * sum |u_k|)

This is a definite function of (r_nu, theta_nu). Compute it exactly.
"""

from mpmath import mp, mpf, pi, sqrt, tan, cos, sin, atan2, log, fabs, power
from mpmath import acos, nstr
from fractions import Fraction

mp.dps = 50

H = 3
K_star = mpf(7) / 30

print("=" * 72)
print("NEUTRINO SECTOR: Dirac Koide from diagonal seesaw")
print("=" * 72)

# Neutrino Koide parameters
Q_nu = mpf(29) / 40
theta_nu = mpf(7) / 810  # = K*/H^3
r_nu = sqrt(6 * Q_nu - 2)  # = sqrt(47/20)

print(f"\nQ_nu   = 29/40 = {Q_nu}")
print(f"theta_nu = 7/810 = {theta_nu}")
print(f"r_nu   = sqrt(47/20) = {r_nu}")

# Compute u_k = 1 + r_nu * tan(theta_nu + 2*pi*k/3)
u = []
for k in range(3):
    arg = theta_nu + 2 * pi * k / 3
    uk = 1 + r_nu * tan(arg)
    u.append(uk)
    print(f"\nu_{k}: arg = {float(arg):.6f}")
    print(f"  tan(arg) = {uk - 1}")
    print(f"  u_{k} = {uk}")
    print(f"  sign: {'positive' if uk > 0 else 'NEGATIVE'}")

print("\n" + "-" * 72)
print("Verification: Koide ratio of u_k^2 should give Q_nu = 29/40")
m_light = [uk**2 for uk in u]
sqrt_m = [fabs(uk) for uk in u]
Q_check = (sum(sqrt_m))**2 / (3 * sum(m_light))
print(f"Q(u_k^2) = {Q_check}")
print(f"29/40    = {mpf(29)/40}")
print(f"Match: {abs(Q_check - mpf(29)/40) < mpf(10)**(-40)}")

print("\n" + "=" * 72)
print("DIRAC MASSES: m_Dk ∝ |u_k|, so sqrt(m_Dk) ∝ |u_k|^{1/2}")
print("=" * 72)

# Dirac masses proportional to |u_k|
abs_u = [fabs(uk) for uk in u]
fourth_root = [sqrt(au) for au in abs_u]  # |u_k|^{1/2}

print("\n|u_k| values:")
for k in range(3):
    print(f"  |u_{k}| = {abs_u[k]}")

print("\n|u_k|^{1/2} values:")
for k in range(3):
    print(f"  |u_{k}|^{{1/2}} = {fourth_root[k]}")

# Q_D = (sum |u_k|^{1/2})^2 / (3 * sum |u_k|)
S_half = sum(fourth_root)
S_one = sum(abs_u)

Q_D = S_half**2 / (3 * S_one)
print(f"\nsum |u_k|^{{1/2}} = {S_half}")
print(f"sum |u_k|       = {S_one}")
print(f"\nQ_D = {Q_D}")
print(f"Q_D (float) = {float(Q_D):.15f}")

# Extract theta_D from the Koide decomposition of |u_k|
# If m_Dk = A*(1 + r_D*cos(theta_D + 2*pi*k/3))^2, then
# sqrt(m_Dk) = sqrt(A)*(1 + r_D*cos(theta_D + 2*pi*k/3))
# We have sqrt(m_Dk) ∝ |u_k|^{1/2}
# So we need to fit: |u_k|^{1/2} = C*(1 + r_D*cos(theta_D + 2*pi*k/3))

# From Q_D, get r_D:
r_D = sqrt(6 * Q_D - 2)
print(f"\nr_D = sqrt(6*Q_D - 2) = {r_D}")
print(f"r_D (float) = {float(r_D):.15f}")

# Extract theta_D from the mass ratios
# v_k = |u_k|^{1/2} / (sum/3) - 1 = r_D * cos(theta_D + 2*pi*k/3)
C = S_half / 3
v = [(fourth_root[k] / C - 1) / r_D for k in range(3)]

# theta_D = atan2(-sum(v_k*sin(2*pi*k/3)), sum(v_k*cos(2*pi*k/3)))
cos_part = sum(v[k] * cos(2 * pi * k / 3) for k in range(3))
sin_part = sum(v[k] * sin(2 * pi * k / 3) for k in range(3))
theta_D = atan2(-sin_part, cos_part)

print(f"\ntheta_D = {theta_D}")
print(f"theta_D (float) = {float(theta_D):.15f}")

print("\n" + "=" * 72)
print("TESTING EXACT CANDIDATES FOR Q_D")
print("=" * 72)

# Test many framework fractions
candidates_Q = [
    ("2/5", mpf(2)/5),
    ("29/70", mpf(29)/70),
    ("13/32", mpf(13)/32),
    ("3/7", mpf(3)/7),
    ("5/12", mpf(5)/12),
    ("11/27", mpf(11)/27),
    ("(H-1)/(2H-1)", mpf(H-1)/(2*H-1)),
    ("K*", K_star),
    ("29/40 - K*", mpf(29)/40 - K_star),
    ("(29/40)/sqrt(2)", mpf(29)/(40*sqrt(2))),
    ("(Q_nu+1/3)/2", (Q_nu + mpf(1)/3)/2),
    ("sqrt(Q_nu)/sqrt(3)", sqrt(Q_nu)/sqrt(3)),
    ("Q_nu/sqrt(3)", Q_nu/sqrt(3)),
    ("1-Q_nu", 1 - Q_nu),
    ("(2*Q_nu-1)", 2*Q_nu - 1),
    ("Q_nu^2", Q_nu**2),
    ("sqrt(Q_nu)", sqrt(Q_nu)),
    ("Q_nu/(1+K*)", Q_nu/(1+K_star)),
    ("(H^2-1)/(H^2+H+1)", mpf(H**2-1)/(H**2+H+1)),
    ("8/19", mpf(8)/19),
    ("8/H^3", mpf(8)/H**3),
    ("5/(3*H+3)", mpf(5)/(3*H+3)),
    ("(H^2+2)/(H^2+H+3)", mpf(H**2+2)/(H**2+H+3)),
    ("11/H^3", mpf(11)/H**3),
    ("(29-7)/(40+30)", mpf(22)/70),
    ("(Q_nu+1/H)/2", (Q_nu + mpf(1)/H)/2),
    ("7/17", mpf(7)/17),
    ("(2H-1)/(2H^2-1)", mpf(2*H-1)/(2*H**2-1)),
    ("5/(H^2+H+1)", mpf(5)/(H**2+H+1)),
    ("(Q_nu-K*)", Q_nu - K_star),
    ("29/40-7/60", mpf(29)/40 - mpf(7)/60),
    ("Q_nu*sqrt(2)/2", Q_nu*sqrt(2)/2),
    ("1/(1+sqrt(47/20))", 1/(1+r_nu)),
]

for name, val in sorted(candidates_Q, key=lambda x: abs(x[1] - Q_D)):
    gap = float((val - Q_D) / Q_D * 100)
    print(f"  {name:30s} = {float(val):.10f}  gap = {gap:+.6f}%")

print("\n" + "=" * 72)
print("TESTING EXACT CANDIDATES FOR theta_D")
print("=" * 72)

sigma = -log(mpf(23)/30)
candidates_theta = [
    ("sigma = -ln(23/30)", sigma),
    ("K*/H", K_star/H),
    ("sigma/H", sigma/H),
    ("2*theta_nu", 2*theta_nu),
    ("7/25", mpf(7)/25),
    ("7*theta_nu", 7*theta_nu),
    ("theta_nu*H^2", theta_nu*H**2),
    ("theta_nu*(H^2-2)", theta_nu*(H**2-2)),
    ("K*/sqrt(H)", K_star/sqrt(H)),
    ("1/sqrt(13)", 1/sqrt(13)),
    ("(H-1)/7", mpf(H-1)/7),
    ("2/(H^2-2)", mpf(2)/(H**2-2)),
    ("sigma*Q_nu", sigma*Q_nu),
    ("(sigma+theta_nu)/2", (sigma+theta_nu)/2),
    ("K*/(H-1)", K_star/(H-1)),
    ("sqrt(K*/H)", sqrt(K_star/H)),
    ("7/H^3*H/2", mpf(7)/(H**2*2)),
    ("theta_nu + K*/(H^3)", theta_nu + K_star/(H**3)),
    ("ln(4/3)", log(mpf(4)/3)),
    ("1/(2*pi-H)", 1/(2*pi-H)),
    ("pi/H^2 - 1/H", pi/H**2 - 1/H),
    ("sigma - theta_nu", sigma - theta_nu),
    ("K*^2*H", K_star**2*H),
    ("7/(3*H^2-2)", mpf(7)/(3*H**2-2)),
    ("7/(H*(H^2+1))", mpf(7)/(H*(H**2+1))),
    ("theta_nu*pi", theta_nu*pi),
    ("1/H - K*", mpf(1)/H - K_star),
    ("sigma*(1-K*)", sigma*(1-K_star)),
    ("sqrt(theta_nu)", sqrt(theta_nu)),
    ("2*K*/(H-1)", 2*K_star/(H-1)),
    ("pi*K*/H^2", pi*K_star/H**2),
    ("H*theta_nu + theta_nu^2", H*theta_nu + theta_nu**2),
    ("theta_nu*H*(H+1)/2", theta_nu*H*(H+1)/2),
    ("(H^2-1)*theta_nu/H", (H**2-1)*theta_nu/H),
]

for name, val in sorted(candidates_theta, key=lambda x: abs(x[1] - theta_D)):
    gap = float((val - theta_D) / theta_D * 100)
    print(f"  {name:30s} = {float(val):.10f}  gap = {gap:+.6f}%")

print("\n" + "=" * 72)
print("CHARGED LEPTON SECTOR: same analysis")
print("=" * 72)

# Charged lepton Koide: cos form
Q_e = mpf(2) / 3
theta_e = mpf(2) / 9
r_e = sqrt(2)  # sqrt(6*2/3 - 2) = sqrt(2)

print(f"\nQ_e = 2/3, theta_e = 2/9, r_e = sqrt(2)")

w = []
for k in range(3):
    arg = theta_e + 2 * pi * k / 3
    wk = 1 + r_e * cos(arg)
    w.append(wk)
    print(f"  w_{k} = {wk}")

# Dirac-type: proportional to |w_k|
abs_w = [fabs(wk) for wk in w]
sqrt_abs_w = [sqrt(aw) for aw in abs_w]

S_half_e = sum(sqrt_abs_w)
S_one_e = sum(abs_w)

Q_De = S_half_e**2 / (3 * S_one_e)
r_De = sqrt(6 * Q_De - 2)

print(f"\nQ_D(charged lepton) = {Q_De}")
print(f"Q_D(charged lepton) float = {float(Q_De):.15f}")
print(f"r_D(charged lepton) = {float(r_De):.15f}")

# Extract theta_De
C_e = S_half_e / 3
v_e = [(sqrt_abs_w[k] / C_e - 1) / r_De for k in range(3)]
cos_part_e = sum(v_e[k] * cos(2 * pi * k / 3) for k in range(3))
sin_part_e = sum(v_e[k] * sin(2 * pi * k / 3) for k in range(3))
theta_De = atan2(-sin_part_e, cos_part_e)

print(f"theta_D(charged lepton) = {theta_De}")
print(f"theta_D(charged lepton) float = {float(theta_De):.15f}")

# Test exact values for charged lepton Q_D
print("\nTesting Q_D(e) against exact values:")
candidates_Qe = [
    ("1/2", mpf(1)/2),
    ("(1+sqrt(2))/4", (1+sqrt(2))/4),
    ("sqrt(2)/2 - 1/6", sqrt(2)/2 - mpf(1)/6),
    ("sqrt(3)/3", sqrt(3)/3),
    ("(2+sqrt(2))/6", (2+sqrt(2))/6),
    ("5/9", mpf(5)/9),
    ("Q_e - 1/6", Q_e - mpf(1)/6),
    ("Q_e/sqrt(2)", Q_e/sqrt(2)),
    ("(Q_e+1/3)/2", (Q_e+mpf(1)/3)/2),
    ("sqrt(Q_e)", sqrt(Q_e)),
    ("7/12", mpf(7)/12),
    ("11/18", mpf(11)/18),
    ("(sqrt(2)+1)/(2*sqrt(2)+1)", (sqrt(2)+1)/(2*sqrt(2)+1)),
    ("4/7", mpf(4)/7),
    ("(1+sqrt(2))/(2+sqrt(2))", (1+sqrt(2))/(2+sqrt(2))),
]

for name, val in sorted(candidates_Qe, key=lambda x: abs(x[1] - Q_De)):
    gap = float((val - Q_De) / Q_De * 100)
    print(f"  {name:30s} = {float(val):.10f}  gap = {gap:+.6f}%")


print("\n" + "=" * 72)
print("ANALYTIC STRUCTURE: Koide under square root")
print("=" * 72)

# Key identity to check:
# If u_k = 1 + r*tan(phi_k) with phi_k = theta + 2*pi*k/3
# Then sum(u_k) and sum(u_k^2) have known closed forms.
# We need sum(|u_k|^{1/2}) and sum(|u_k|).
#
# Note: sum_{k=0}^{2} tan(theta+2*pi*k/3) = 3*tan(3*theta)
# So sum(u_k) = 3*(1 + r*tan(3*theta)/3*... no, that's wrong.
# sum(u_k) = 3 + r * sum(tan(phi_k)) = 3 + r * 3*tan(3*theta)

print("\nVerification of tan identity:")
tan_sum = sum(tan(theta_nu + 2*pi*k/3) for k in range(3))
print(f"sum tan(theta+2pi*k/3) = {tan_sum}")
print(f"3*tan(3*theta)         = {3*tan(3*theta_nu)}")

# Also: sum(tan^2(phi_k)) = 3 + 6/cos(6*theta) - 6/(cos(6*theta))...
# Actually: sum tan^2(phi_k) = 3*(tan^2(3*theta) + 2*sec^2(3*theta) - 2)
# Let me just verify numerically
tan2_sum = sum(tan(theta_nu + 2*pi*k/3)**2 for k in range(3))
# sum u_k^2 = sum(1 + r*tan)^2 = 3 + 2r*sum(tan) + r^2*sum(tan^2)
u2_sum = sum(uk**2 for uk in u)
print(f"\nsum(u_k^2) = {u2_sum}")
print(f"= 3 + 2r*3tan(3th) + r^2*sum(tan^2) = {3 + 2*r_nu*3*tan(3*theta_nu) + r_nu**2*tan2_sum}")

# Now the key question: is there a closed form for sum(|u_k|^{1/2})?
# This is NOT algebraic in general. But for specific (r, theta), it might simplify.

print("\n" + "=" * 72)
print("DEEPER SEARCH: varying seesaw exponent")
print("=" * 72)

# What if the seesaw isn't exactly m_D = sqrt(m_nu * M_R)?
# Type-I seesaw: m_nu = m_D^2 / M_R, so m_D = sqrt(m_nu * M_R)
# But what if there's a different power? m_D ∝ m_nu^alpha
# Then sqrt(m_D) ∝ m_nu^{alpha/2}, and the Koide becomes a function of alpha.

print("\nScan: Q_D as function of seesaw exponent alpha (m_D ∝ m_nu^alpha):")
print(f"{'alpha':>10s} {'Q_D':>15s} {'theta_D':>15s}")

target_Q_values = {
    "2/5": mpf(2)/5,
    "3/7": mpf(3)/7,
    "5/12": mpf(5)/12,
    "29/70": mpf(29)/70,
}

for alpha_num in range(1, 40):
    alpha = mpf(alpha_num) / 20  # 0.05 to 1.95
    m_D = [fabs(uk)**(2*alpha) for uk in u]  # m_k^alpha ∝ |u_k|^{2*alpha}
    sqrt_mD = [sqrt(md) for md in m_D]  # = |u_k|^alpha

    S1 = sum(sqrt_mD)
    S2 = sum(m_D)
    Q = S1**2 / (3 * S2)

    # Extract theta
    r = sqrt(6*Q - 2) if 6*Q - 2 > 0 else mpf(0)
    if r > 0:
        C_val = S1 / 3
        vv = [(sqrt_mD[k]/C_val - 1)/r for k in range(3)]
        cp = sum(vv[k] * cos(2*pi*k/3) for k in range(3))
        sp = sum(vv[k] * sin(2*pi*k/3) for k in range(3))
        th = atan2(-sp, cp)
    else:
        th = mpf(0)

    marker = ""
    for name, tgt in target_Q_values.items():
        if abs(Q - tgt) < mpf('0.005'):
            marker = f" <-- near {name}"

    print(f"  {float(alpha):8.4f}  {float(Q):15.12f}  {float(th):15.12f}{marker}")

print("\n" + "=" * 72)
print("FINE SCAN near alpha = 0.5 (standard seesaw)")
print("=" * 72)

# For standard seesaw alpha=0.5: m_D ∝ m_nu^{1/2}
# Check if Q_D matches any fraction precisely

print(f"\nPrecise Q_D at alpha=1/2: {nstr(Q_D, 40)}")
print(f"Precise theta_D at alpha=1/2: {nstr(theta_D, 40)}")

# Try PSLQ-style: is Q_D a root of a low-degree polynomial with small integer coefficients?
# Test Q_D against a^2 + b*Q_D + c = 0 type relations
print("\nSearching for minimal polynomial of Q_D:")
from mpmath import matrix, lu_solve

# Test: is Q_D algebraic of degree <= 4?
# Build [1, Q_D, Q_D^2, Q_D^3, Q_D^4] and look for integer relation
QD_powers = [Q_D**n for n in range(6)]

# Try to find integer relation using LLL/PSLQ
from mpmath import pslq
print("\nPSLQ on [1, Q_D, Q_D^2, Q_D^3, Q_D^4]:")
rel = pslq([1, Q_D, Q_D**2, Q_D**3, Q_D**4], maxcoeff=1000)
if rel:
    print(f"  Found: {rel[0]} + {rel[1]}*x + {rel[2]}*x^2 + {rel[3]}*x^3 + {rel[4]}*x^4 = 0")
else:
    print("  No relation found with coefficients <= 1000")

# Try involving sqrt(47/20), pi, and other framework constants
print("\nPSLQ on [1, Q_D, r_nu, r_nu*Q_D]:")
rel2 = pslq([1, Q_D, r_nu, r_nu*Q_D], maxcoeff=200)
if rel2:
    print(f"  Found: {rel2}")

print("\nPSLQ on [1, Q_D, theta_nu, pi]:")
rel3 = pslq([1, Q_D, theta_nu, pi], maxcoeff=200)
if rel3:
    print(f"  Found: {rel3}")

print("\nPSLQ on [1, Q_D, Q_nu]:")
rel4 = pslq([1, Q_D, Q_nu], maxcoeff=200)
if rel4:
    print(f"  Found: {rel4}")

print("\nPSLQ on [1, Q_D, K*]:")
rel5 = pslq([1, Q_D, K_star], maxcoeff=200)
if rel5:
    print(f"  Found: {rel5}")

print("\nPSLQ on [1, Q_D, Q_nu, K*, Q_nu*K*]:")
rel6 = pslq([1, Q_D, Q_nu, K_star, Q_nu*K_star], maxcoeff=200)
if rel6:
    print(f"  Found: {rel6}")

# Same for theta_D
print("\n\nPSLQ on [1, theta_D, sigma]:")
rel7 = pslq([1, theta_D, sigma], maxcoeff=200)
if rel7:
    print(f"  Found: {rel7}")

print("\nPSLQ on [1, theta_D, theta_nu, pi]:")
rel8 = pslq([1, theta_D, theta_nu, pi], maxcoeff=500)
if rel8:
    print(f"  Found: {rel8}")

print("\nPSLQ on [1, theta_D, K*, sigma]:")
rel9 = pslq([1, theta_D, K_star, sigma], maxcoeff=500)
if rel9:
    print(f"  Found: {rel9}")

print("\nPSLQ on [1, theta_D, theta_nu, K*]:")
rel10 = pslq([1, theta_D, theta_nu, K_star], maxcoeff=500)
if rel10:
    print(f"  Found: {rel10}")

print("\nPSLQ on [1, theta_D, ln(r_nu)]:")
rel11 = pslq([1, theta_D, log(r_nu)], maxcoeff=500)
if rel11:
    print(f"  Found: {rel11}")

print("\nPSLQ on [1, theta_D, sqrt(K*)]:")
rel12 = pslq([1, theta_D, sqrt(K_star)], maxcoeff=500)
if rel12:
    print(f"  Found: {rel12}")


print("\n" + "=" * 72)
print("ALTERNATIVE: Koide in TAN basis directly")
print("=" * 72)

# The neutrino Koide uses tan, not cos. What if the Dirac masses
# also follow a tan-Koide? i.e., sqrt(m_Dk) = A*(1 + r_D*tan(theta_D_tan + 2*pi*k/3))
# Then m_Dk = A^2*(1 + r_D*tan(...))^2

# We have m_Dk ∝ |u_k| = |1 + r_nu*tan(theta_nu + 2*pi*k/3)|
# So sqrt(m_Dk) ∝ |u_k|^{1/2}

# Can we write |u_k|^{1/2} = C*(1 + r'*tan(theta' + 2*pi*k/3))?
# If so, then (1 + r'*tan(theta' + 2*pi*k/3))^2 = |u_k|/C'^2
# This would require the LHS to be proportional to |u_k|, which means
# 1 + 2r'*tan(phi_k) + r'^2*tan^2(phi_k) ∝ |1 + r_nu*tan(psi_k)|
# where phi_k and psi_k have the same 2*pi*k/3 spacing but possibly different offsets.
# This is highly unlikely to hold exactly.

# Instead, let's try the COS basis for the Dirac sector
# sqrt(m_Dk) = A*(1 + r_D*cos(delta_D + 2*pi*k/3))
# This is the standard Koide parameterization.
# We already computed Q_D and theta_D in this basis above.

# But there's another possibility: the Dirac Koide might be in the TAN basis
# sqrt(m_Dk) = A*(1 + rr*tan(tt + 2*pi*k/3))
# m_Dk = A^2*(1 + rr*tan(tt + 2*pi*k/3))^2

# The Q for this is: Q_tan = (sum(1+rr*tan(tt+2pk/3)))^2 / (3*sum(1+rr*tan)^2)
# = (3 + rr*3*tan(3tt))^2 / (3*(3 + 2rr*3*tan(3tt) + rr^2*sum(tan^2)))

# For m_Dk ∝ |u_k|, we need |u_k| = A^2*(1+rr*tan(tt+2pk/3))^2
# i.e., |u_k|^{1/2} = A*(1+rr*tan(tt+2pk/3))
# So: |1 + r_nu*tan(theta_nu+2pk/3)|^{1/2} = A*(1 + rr*tan(tt+2pk/3))

# 3 equations, 3 unknowns (A, rr, tt). Let's solve:
from mpmath import findroot

def tan_koide_eqs(A, rr, tt):
    eqs = []
    for k in range(3):
        lhs = sqrt(fabs(u[k]))
        rhs = A * (1 + rr * tan(tt + 2*pi*k/3))
        eqs.append(lhs - rhs)
    return eqs

# Initial guess from cos Koide
A0 = S_half / 3
try:
    sol = findroot(tan_koide_eqs, (A0, mpf('0.1'), theta_D), tol=mpf(10)**(-40))
    A_sol, rr_sol, tt_sol = sol
    print(f"\nTan-basis fit for Dirac sector:")
    print(f"  A  = {A_sol}")
    print(f"  r' = {rr_sol}")
    print(f"  t' = {tt_sol}")

    Q_tan = (1 + rr_sol**2 / (6 - 2)) if False else None  # placeholder

    # Compute Q in tan basis: Q = (sum sqrt(m))^2 / (3 sum m)
    # with sqrt(m_k) = A*(1+r'*tan(t'+2pk/3))
    sqrt_m_tan = [A_sol * (1 + rr_sol * tan(tt_sol + 2*pi*k/3)) for k in range(3)]
    m_tan = [s**2 for s in sqrt_m_tan]
    Q_tan = sum(sqrt_m_tan)**2 / (3 * sum(m_tan))
    print(f"  Q_tan = {Q_tan}")
    print(f"  (should match Q_D = {Q_D})")

    # Now check if rr_sol and tt_sol are framework quantities
    print(f"\n  r' = {float(rr_sol):.15f}")
    print(f"  t' = {float(tt_sol):.15f}")

    # Q for tan-Koide is: 1/(1 + 2*r'^2/(9*tan^2(3t') + ...))
    # Actually just use the general formula:
    # Q_tan = (3 + r'*3*tan(3*t'))^2 / (3*(3 + 2*r'*3*tan(3*t') + r'^2*sum(tan^2)))

except Exception as e:
    print(f"  Tan-basis fit failed: {e}")


print("\n" + "=" * 72)
print("CRITICAL TEST: What if Dirac Koide uses HALF-ANGLE?")
print("=" * 72)

# theta_nu = 7/810. What if theta_D = theta_nu/2 exactly?
theta_half = theta_nu / 2  # = 7/1620
print(f"\ntheta_nu/2 = {float(theta_half):.15f}")
print(f"theta_D    = {float(theta_D):.15f}")
print(f"ratio      = {float(theta_D/theta_half):.15f}")
print(f"gap        = {float((theta_D - theta_half)/theta_D * 100):.6f}%")

# What about theta_D = theta_nu * H?
theta_H = theta_nu * H  # = 7/270
print(f"\ntheta_nu*H = {float(theta_H):.15f}")
print(f"theta_D    = {float(theta_D):.15f}")
print(f"gap        = {float((theta_D - theta_H)/theta_D * 100):.6f}%")

# What about theta_D expressed as fraction with denominator up to 1000?
print("\nBest rational approximations to theta_D:")
# Use continued fraction
from mpmath import identify
result = identify(theta_D, tol=1e-8)
if result:
    print(f"  identify: {result}")

# Manual search
best_gap = 1
best_frac = ""
for denom in range(1, 2001):
    numer = int(float(theta_D * denom) + 0.5)
    if numer > 0:
        val = mpf(numer) / denom
        gap = abs(float((val - theta_D) / theta_D))
        if gap < best_gap:
            best_gap = gap
            best_frac = f"{numer}/{denom}"
            if gap < 1e-6:
                print(f"  {numer}/{denom} = {float(val):.12f}, gap = {gap:.2e}")

print("\nBest rational approximations to Q_D:")
best_gap = 1
for denom in range(1, 2001):
    numer = int(float(Q_D * denom) + 0.5)
    if numer > 0:
        val = mpf(numer) / denom
        gap = abs(float((val - Q_D) / Q_D))
        if gap < best_gap:
            best_gap = gap
            best_frac = f"{numer}/{denom}"
            if gap < 1e-6:
                print(f"  {numer}/{denom} = {float(val):.12f}, gap = {gap:.2e}")


print("\n" + "=" * 72)
print("PSLQ ON theta_D with more framework constants")
print("=" * 72)

print("\nPSLQ on [1, theta_D, pi, theta_nu, sigma, K*, ln(2), ln(3)]:")
rel_big = pslq([1, theta_D, pi, theta_nu, sigma, K_star, log(2), log(3)], maxcoeff=100)
if rel_big:
    labels = ["1", "theta_D", "pi", "theta_nu", "sigma", "K*", "ln(2)", "ln(3)"]
    terms = [f"{c}*{l}" for c, l in zip(rel_big, labels) if c != 0]
    print(f"  Found: {' + '.join(terms)} = 0")
    print(f"  Raw: {rel_big}")

print("\nPSLQ on [1, Q_D, pi, Q_nu, K*, sqrt(2), sqrt(3)]:")
rel_big2 = pslq([1, Q_D, pi, Q_nu, K_star, sqrt(2), sqrt(3)], maxcoeff=100)
if rel_big2:
    labels2 = ["1", "Q_D", "pi", "Q_nu", "K*", "sqrt(2)", "sqrt(3)"]
    terms2 = [f"{c}*{l}" for c, l in zip(rel_big2, labels2) if c != 0]
    print(f"  Found: {' + '.join(terms2)} = 0")
    print(f"  Raw: {rel_big2}")

# Try Q_D in terms of Q_nu and tan(3*theta_nu)
t3 = tan(3 * theta_nu)
print(f"\ntan(3*theta_nu) = {t3}")
print(f"\nPSLQ on [1, Q_D, Q_nu, t3, t3^2, Q_nu*t3]:")
rel_t3 = pslq([1, Q_D, Q_nu, t3, t3**2, Q_nu*t3], maxcoeff=200)
if rel_t3:
    print(f"  Found: {rel_t3}")

# The most natural: Q_D as function of r_nu and tan(3*theta_nu)
# since sum(u_k) = 3(1 + r_nu*tan(3*theta_nu)) etc.
S1_u = 3 + r_nu * 3 * tan(3*theta_nu)
print(f"\nsum(u_k) = {sum(u)}")
print(f"3 + 3*r_nu*tan(3*theta_nu) = {S1_u}")
# These should match only if all u_k > 0
print(f"sum(|u_k|) = {sum(abs_u)}")
print(f"All u_k positive? {all(uk > 0 for uk in u)}")

if all(uk > 0 for uk in u):
    print("\nAll u_k positive! So |u_k| = u_k, and:")
    print("sum(u_k) = 3 + 3*r_nu*tan(3*theta_nu)")

    # For Q_D = (sum u_k^{1/2})^2 / (3 * sum u_k)
    # = (sum u_k^{1/2})^2 / (3*(3 + 3*r_nu*tan(3*theta_nu)))
    # = (sum u_k^{1/2})^2 / (9*(1 + r_nu*tan(3*theta_nu)))

    denom_val = 9 * (1 + r_nu * tan(3*theta_nu))
    numer_val = S_half**2
    print(f"\nQ_D = S_half^2 / denom = {numer_val} / {denom_val} = {numer_val/denom_val}")
    print(f"Matches Q_D = {Q_D}? {abs(numer_val/denom_val - Q_D) < mpf(10)**(-40)}")

    # The problem reduces to: what is sum(u_k^{1/2})?
    # u_k = 1 + r*tan(theta + 2*pi*k/3)
    # u_k^{1/2} = (1 + r*tan(theta + 2*pi*k/3))^{1/2}
    # This has no known closed form in general.

    # But numerically:
    print(f"\nsum(u_k^{{1/2}}) = {S_half}")
    print(f"sum(u_k^{{1/2}})^2 = {S_half**2}")

    # Try: is S_half^2 a nice multiple of something?
    ratio_test = S_half**2 / 9
    print(f"\nS_half^2 / 9 = {ratio_test}")
    print(f"1 + r_nu*tan(3*theta_nu) = {1 + r_nu*tan(3*theta_nu)}")
    print(f"Q_D = ratio / (1+r*t3) = {ratio_test / (1 + r_nu*tan(3*theta_nu))}")


print("\n" + "=" * 72)
print("FINAL: Is the gap physical (RG running)?")
print("=" * 72)

# The Koide relation for charged leptons holds at low energy.
# For neutrinos, the seesaw operates at M_R scale.
# The 3-4% gap might be the RG running from M_R to low energy.

# Alternatively: the seesaw might not be purely type-I diagonal.
# If there are off-diagonal elements, the Dirac Koide gets corrections.

print("""
Summary of results:
""")
print(f"Q_D  = {nstr(Q_D, 20)}")
print(f"θ_D  = {nstr(theta_D, 20)}")
print(f"r_D  = {nstr(r_D, 20)}")
print(f"")
print(f"Nearest framework values:")
print(f"  Q_D ≈ 2/5 = 0.400   (gap: {float((Q_D - mpf(2)/5)/Q_D * 100):.4f}%)")
print(f"  Q_D ≈ 29/70         (gap: {float((Q_D - mpf(29)/70)/Q_D * 100):.4f}%)")
print(f"  θ_D ≈ σ=-ln(23/30)  (gap: {float((theta_D - sigma)/theta_D * 100):.4f}%)")
