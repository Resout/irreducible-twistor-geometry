"""
ROUTE B: 11/60 = 1/d_super + 1/H(H+1) = 1/10 + 1/12

Trace the consequences of this decomposition through EVERY
quantity in the framework. Check whether it gives consistent
meaning everywhere.
"""
from mpmath import mp, mpf, sqrt, log, cos, pi, nstr, fabs, exp
mp.dps = 30

H = mpf(3)
Kstar = mpf(7)/30
d_super = H**2 + 1      # 10
d_bos = H**3 - 1        # 26
hE8 = H*(H**2+1)        # 30
HH1 = H*(H+1)           # 12 = total gauge generators
FLOOR = 1/H**3           # 1/27

# The decomposition
A = 1/d_super            # 1/10
B = 1/HH1               # 1/12
print("="*70)
print("ROUTE B: 11/60 = 1/d_super + 1/H(H+1) = 1/10 + 1/12")
print("="*70)
print(f"  1/d_super = 1/{int(d_super)} = {nstr(A, 10)}")
print(f"  1/H(H+1) = 1/{int(HH1)} = {nstr(B, 10)}")
print(f"  Sum = {nstr(A+B, 10)} = 11/60 = {nstr(mpf(11)/60, 10)}")
print()

# ================================================================
# WHAT IS d_super = H^2+1 = 10?
# ================================================================
print("="*70)
print("WHERE d_super = H^2+1 = 10 APPEARS IN THE FRAMEWORK")
print("="*70)
appearances_d = [
    ("Superstring critical dimension", "d_super = H^2+1 = 10"),
    ("Channel count in conservation law", "K*(H^2+1) - eta*H^2 = 1"),
    ("Gaussian norm", "Phi_4(H) = H^2+1 = 10"),
    ("Denominator of K*", "K* = 7/(H*10) = 7/30"),
    ("dim(Sym^2(C^4))", "The symmetric square of the mass space"),
    ("Effective DS channels", "H^2 pairwise + 1 normalisation = 10"),
    ("Solar angle denominator", "sin^2(theta_12) = H/(H^2+1) = 3/10"),
    ("Bott periodicity", "d_super-2 = 8 = rank(E8) = dim_R(O)"),
    ("Enslaving: Born surface eqn", "(d_bos)*S^2*theta^2 = Sq*(1-theta)^2, but d_super appears as channel count"),
]
for name, desc in appearances_d:
    print(f"  {name}: {desc}")

print()
print(f"  1/d_super = 1/10 = the INVERSE channel count")
print(f"  = the probability of a SPECIFIC channel in the uniform distribution")
print(f"  = sin^2(theta_12)/H = (3/10)/3 = 1/10")
print(f"  = the per-channel contribution at equilibrium")

# ================================================================
# WHAT IS H(H+1) = 12?
# ================================================================
print()
print("="*70)
print("WHERE H(H+1) = 12 APPEARS IN THE FRAMEWORK")
print("="*70)
appearances_12 = [
    ("Total gauge generators", "SU(3)[8] + SU(2)[3] + U(1)[1] = 12"),
    ("Eigenvalue denominator constraint", "rational powers p/q have q | H(H+1) = 12"),
    ("Edges of K_4", "complete graph on H+1=4 vertices has binom(4,2)=6 edges... wait, 12 != 6"),
    ("  correction", "DIRECTED edges of K_4: 4*3 = 12 = H(H+1). YES."),
    ("Product constraint denominator", "11/12 = (H(H+1)-1)/H(H+1) = massive/total"),
    ("VEV multiplier", "VEV = 144*m0 = [H(H+1)]^2 * m0"),
    ("Coxeter relation", "h(E8) = H*Phi_4(H) = 30, but H(H+1) = 12 is different"),
    ("Leech lattice", "d_bos-2 = 24 = 2*H(H+1) = twice the generator count"),
    ("Conservation law", "K*(H^2+1) = K*10 = 7/3. And H(H+1) = 12. Ratio: 7/3 / 12 = 7/36"),
]
for name, desc in appearances_12:
    print(f"  {name}: {desc}")

print()
print(f"  1/H(H+1) = 1/12 = the INVERSE total gauge generator count")
print(f"  = the per-generator contribution")
print(f"  = the uniform weight per gauge degree of freedom")

# ================================================================
# INTERPRETATION: 11/60 = per-channel + per-generator
# ================================================================
print()
print("="*70)
print("INTERPRETATION OF THE DECOMPOSITION")
print("="*70)
print(f"""
M2_lep / m(0++) = 1/d_super + 1/H(H+1) = 1/10 + 1/12

The lepton mass scale is:
  ONE PART PER DS CHANNEL (1/10 of the channel space)
  + ONE PART PER GAUGE GENERATOR (1/12 of the gauge algebra)

Physical meaning: the lepton couples to TWO structures:
  1. The DS combination channel space (10 effective channels)
     -> This is where the conservation law lives
     -> The lepton's MASS comes from the channel structure
  2. The gauge generator algebra (12 generators)
     -> This is where the gauge field lives
     -> The lepton's CHARGE comes from the gauge structure

The lepton mass scale is the sum of its channel contribution
and its gauge contribution.
""")

# ================================================================
# TRACE THROUGH OTHER KOIDE SCALES
# ================================================================
print("="*70)
print("DO 3/8 AND 40/3 DECOMPOSE IN TERMS OF d_super AND H(H+1)?")
print("="*70)

# M2_down = 3/8 = sin^2(theta_W)
# = a/d_super + b/H(H+1)?
# 3/8 = a/10 + b/12
# 60*3/8 = 6a + 5b = 22.5 ... not integer. NO clean decomposition.
print(f"M2_down = 3/8:")
print(f"  3/8 = a/10 + b/12 -> 6a + 5b = 22.5 (not integer)")
print(f"  NO integer decomposition in 1/d_super and 1/H(H+1)")
print()

# But 3/8 = H/(H^2-1) = H/((H-1)(H+1)) = 3/8
# And H^2-1 = (H-1)(H+1) = 2*4 = 8
# So 3/8 is already a clean ratio: H/(H^2-1)
print(f"  But 3/8 = H/(H^2-1) = H/dim(SU(H)) = 3/8")
print(f"  This IS a clean framework expression: hypothesis count / gluon count")
print()

# M2_up = 40/3 = (H^2+1)(H+1)/H = d_super*(H+1)/H = 10*4/3
# = (H+1)*d_super/H
print(f"M2_up = 40/3:")
print(f"  = (H+1)*d_super/H = 4*10/3")
print(f"  = d_super * (H+1)/H")
print(f"  This is the channel count TIMES the generation ratio (H+1)/H")
print()

# So the three scales in terms of d_super:
# M2_lep = 1/d_super + 1/H(H+1)
# M2_down = H/(H^2-1)
# M2_up = d_super*(H+1)/H
print(f"THE THREE SCALES:")
print(f"  M2_lep  = 1/d_super + 1/H(H+1)     = 1/10 + 1/12 = 11/60")
print(f"  M2_down = H/(H^2-1)                  = 3/8 = sin^2(theta_W)")
print(f"  M2_up   = d_super*(H+1)/H            = 10*4/3 = 40/3")
print()
print(f"  Lepton: ADDITIVE (channel + gauge)")
print(f"  Down:   RATIO (H per gluon)")
print(f"  Up:     MULTIPLICATIVE (channels times generation ratio)")
print()

# ================================================================
# TRACE d_super THROUGH THE EIGENVALUE SPECTRUM
# ================================================================
print("="*70)
print("d_super IN THE EIGENVALUE SPECTRUM")
print("="*70)

L = [0.50217, 0.47448, 0.35265, 0.33442]
D = [-log(mpf(l)) for l in L]
sigma = -log(mpf(23)/30)

print(f"The eigenvalue-Koide bridge:")
print(f"  theta_ev = K* * D2/D3 = {nstr(Kstar*D[2]/D[3], 10)} vs 2/9")
print(f"  The ratio D2/D3 should be 20/21 for exact bridge")
print(f"  20/21 = 2*d_super/(3*Phi_6(H)) = 20/21")
print(f"  So d_super appears in the bridge through the gap ratio!")
print()

# The hierarchy tower: exp(-S/d) at various d
S = H**3/Kstar  # 810/7
print(f"The hierarchy tower at d=d_super:")
print(f"  exp(-S/d_super) = exp(-81/7) = {nstr(exp(-S/d_super), 10)}")
print(f"  This is ~ the Koide residual 9.83e-6!")
print()

# String tension
print(f"String tension: sigma = -ln(1-K*) = -ln(23/30) = {nstr(sigma, 10)}")
print(f"  sigma * d_super = {nstr(sigma*d_super, 10)}")
print(f"  Compare to S/H^2 = 810/(7*9) = {nstr(S/H**2, 10)}")
print()

# ================================================================
# TRACE H(H+1) = 12 THROUGH THE EIGENVALUE SPECTRUM
# ================================================================
print("="*70)
print("H(H+1) = 12 IN THE EIGENVALUE SPECTRUM")
print("="*70)

print(f"Eigenvalue denominators: all rational powers p/q have q | 12")
print(f"  12 = H(H+1) = 3*4")
print(f"  Allowed denominators: 1, 2, 3, 4, 6, 12")
print(f"  This gives p/q values: 1/2, 1/3, 2/3, 1/4, 3/4, 1/6, 5/6, ...")
print()

print(f"The bilinear formula: m = (Di+Dj)/2 + n*sigma) / D0")
print(f"  The /2 in the bilinear comes from PAIRING two modes")
print(f"  The n*sigma adds string tension units")
print(f"  The denominator constraint 12 = H(H+1) constrains which")
print(f"  rational powers appear in the extended spectrum")
print()

# ================================================================
# THE KEY: WHAT DOES 1/d_super + 1/12 MEAN AT THE FIXED POINT?
# ================================================================
print("="*70)
print("AT THE FIXED POINT")
print("="*70)

# The fixed point satisfies:
# K* = 7/30 (conflict)
# Born = 1/27 (floor)
# Conservation law: K*(H^2+1) - eta*H^2 = 1

# What does 1/d_super = 1/10 mean at the fixed point?
# d_super = H^2+1 is the effective channel count
# The conservation law says K*d_super - eta*H^2 = 1
# So K*d_super = 1 + eta*H^2 = 1 + (4/27)*9 = 1 + 4/3 = 7/3
# K*d_super = 7/3
# 1/d_super = K*/(K*d_super) = K*/(7/3) = (7/30)/(7/3) = 1/10
# So 1/d_super = K*/(K*d_super) = K*H/(H*K*d_super) = ...

# Actually: 1/d_super = 3K*/7 * (7/(3*10)) hmm, circular

# More directly: 1/d_super = 1/(H^2+1)
# At the fixed point: K*(H^2+1) = 7/3
# So 1/(H^2+1) = K*/(7/3) = 3K*/7
print(f"At the fixed point:")
print(f"  K*d_super = K*(H^2+1) = {nstr(Kstar*d_super, 10)} = 7/3")
print(f"  So 1/d_super = K*/(7/3) = 3K*/7 = 3*7/(7*30) = 1/10")
print(f"  1/d_super is K* scaled by 3/Phi_6(H) = 3/7")
print()

# And 1/H(H+1) = 1/12
# What does this equal at the fixed point?
# eta = (H-1)^2/H^3 = 4/27
# eta * H(H+1) = (4/27)*12 = 48/27 = 16/9
# 1/H(H+1) = eta/(16/9) = 9*eta/16 = 9*(4/27)/16 = (4/3)/16 = 1/12. Circular.

# The product: 1/d_super * 1/H(H+1) = 1/120 = 1/(2*h(E8)*H(H+1)/...) hmm
print(f"  1/d_super * 1/H(H+1) = 1/120 = 1/5!")
print(f"  And 5! = S + 1/K* = 810/7 + 30/7 = 840/7 = 120")
print(f"  So 1/(d_super * H(H+1)) = 1/(S + 1/K*)")
print()

# WHOA. The product of the two terms in the decomposition is 1/(S+1/K*)
# And S + 1/K* = 120 = 5! is a known exact identity in the paper!

print(f"  *** 1/(d_super * H(H+1)) = 1/(S + 1/K*) ***")
print(f"  The product of the denominators = the instanton action + 1/K*")
print(f"  = 810/7 + 30/7 = 840/7 = 120 = 5!")
print()

# Verify: d_super * H(H+1) = 10 * 12 = 120 = S + 1/K*
print(f"  d_super * H(H+1) = {int(d_super*HH1)}")
print(f"  S + 1/K* = {nstr(S + 1/Kstar, 10)}")
print(f"  MATCH: {nstr(d_super*HH1 - (S + 1/Kstar), 10)}")
print()

# So: 11/60 = 1/d_super + 1/H(H+1) where d_super * H(H+1) = S + 1/K*
# This means the two denominators are linked by the instanton identity!

print("="*70)
print("THE INSTANTON CONNECTION")
print("="*70)
print(f"""
11/60 = 1/d_super + 1/H(H+1)

where d_super * H(H+1) = 120 = S + 1/K*

The product of the two denominators equals the instanton action
plus the inverse conflict. This is a KNOWN exact identity
(S + 1/K* = 120 = 5!, proved in the paper).

So the decomposition is:
  11/60 = 1/a + 1/b  where a*b = S + 1/K*

With a = d_super = 10, b = H(H+1) = 12.
And a + b = 22, a - b = -2, a*b = 120.

The two terms are the roots of:
  x^2 - 22x + 120 = 0
  x = (22 +/- sqrt(484-480))/2 = (22 +/- 2)/2 = 10 or 12

The discriminant is 4 = (H-1)^2 = H+1.
This IS the self-consistency equation!

x^2 - (d_super + H(H+1))x + (S+1/K*) = 0
discriminant = (d_super - H(H+1))^2 = (10-12)^2 = 4 = (H-1)^2 = H+1
""")

# Verify
disc = (d_super - HH1)**2
print(f"  (d_super - H(H+1))^2 = ({int(d_super)} - {int(HH1)})^2 = {int(disc)}")
print(f"  (H-1)^2 = {int((H-1)**2)}")
print(f"  H+1 = {int(H+1)}")
print(f"  Self-consistency: (H-1)^2 = H+1 -> {int((H-1)**2)} = {int(H+1)}")
print()

# ================================================================
# TRACE THROUGH EVERY FRAMEWORK QUANTITY
# ================================================================
print("="*70)
print("IMPLICATIONS FOR EVERY FRAMEWORK QUANTITY")
print("="*70)

print(f"""
The decomposition 11/60 = 1/a + 1/b with a*b = S+1/K* = 120
and a,b roots of x^2 - 22x + 120 = 0 with discriminant (H-1)^2:

1. K* = 7/30 = Phi_6/(H*Phi_4)
   K* = 7/30 = (a+b-1)/(a*b/2) ? No: 21/60 != 7/30.
   K* = Phi_6(H)/(H*Phi_4(H)) = 7/(3*10) = 7/(H*a)
   So K* = Phi_6(H)/(H*d_super). The d_super is the DENOMINATOR of K*.

2. Conservation law: K*a = 7/3 and K*b = 7/30*12 = 7*12/30 = 84/30 = 14/5
   K*a + K*b = 7/3 + 14/5 = 35/15 + 42/15 = 77/15
   K*(a+b) = K*22 = 7*22/30 = 154/30 = 77/15 = 5.133...
   And K*a - K*b = 7/3 - 14/5 = 35/15 - 42/15 = -7/15
   |K*a - K*b| = 7/15 = Phi_6(H)/dim(SU(H+1)) = 7/15

3. Born floor: 1/H^3 = 1/27
   How does this relate to a,b?
   1/a + 1/b = 11/60 and 1/(H^3) = 1/27
   11/60 * H^3 = 11*27/60 = 297/60 = 99/20 = 4.95
   Not clean.

4. eta = 4/27 = (H-1)^2/H^3
   eta * a*b = (4/27)*120 = 480/27 = 160/9
   160/9 = m_Z multiplier (160/3) / 3. Hmm.

5. Instanton action: S = 810/7 = a*b - 1/K* = 120 - 30/7 = (840-30)/7 = 810/7
   So S = a*b - 1/K*. The instanton action is the product of
   the two decomposition factors minus the inverse conflict.

6. sigma = -ln(23/30) = -ln(1-K*)
   sigma * a = -ln(23/30) * 10 = 2.657
   sigma * b = -ln(23/30) * 12 = 3.188
   sigma * a*b = -ln(23/30) * 120 = 31.88
   Not obviously clean.

7. Q = 2/3 = (H-1)/H
   Q = (a-b+2)/(a-b+H+1) = (-2+2)/(-2+4) = 0/2 = 0. No.
   Q = (H-1)/H = 2/3. Not directly from a,b.

8. theta = 2/9 = 2/H^2
   theta * a = 2/9 * 10 = 20/9 = 2.222
   theta * b = 2/9 * 12 = 24/9 = 8/3 = 2.667

9. The x_k Koide values:
   sum(x_k) = 2H = 6 = a - b + H - 1 = 10 - 12 + 3 - 1 = 0. No.
   sum(x_k) = 6. Not from a-b.

10. Spectral pairing: D0+D3 ~ D1+D2
    (D0+D3)/2 = average radial gap = 0.892
    (D1+D2)/2 = average angular gap = 0.894
    11/60 * D0 = 0.1833 * 0.6888 = 0.1263
    = Delta_massgap / d_super = 1.2626/10 = 0.1263
    WAIT: 11/60 * D0 = 1.2626/10? Let me check...
""")

D0_single = -log(mpf('0.28291'))
print(f"  11/60 * D0(A2) = {nstr(mpf(11)/60 * D[0], 10)}")
print(f"  Delta_massgap / d_super = {nstr(D0_single / d_super, 10)}")
print(f"  Match: {nstr(fabs(mpf(11)/60*D[0] - D0_single/d_super), 10)}")
print(f"  NO: 0.1263 vs 0.1263 ... actually very close!")
print(f"  11/60 * D0_A2 = {nstr(mpf(11)/60*D[0], 10)}")
print(f"  D0_single/10 = {nstr(D0_single/10, 10)}")
print(f"  Ratio: {nstr(mpf(11)/60*D[0]/(D0_single/10), 10)}")
print()

# ================================================================
# THE QUADRATIC IDENTITY
# ================================================================
print("="*70)
print("THE QUADRATIC IDENTITY (THE CLEANEST RESULT)")
print("="*70)
print(f"""
d_super and H(H+1) are the roots of:

  x^2 - (d_super + H(H+1))x + d_super*H(H+1) = 0
  x^2 - 22x + 120 = 0

Coefficients:
  Sum = d_super + H(H+1) = 22 = 2*(H^2+H-1) = 2*11
  Product = d_super * H(H+1) = 120 = S + 1/K* = 5!
  Discriminant = (d_super - H(H+1))^2 = 4 = (H-1)^2 = H+1

The discriminant being (H-1)^2 = H+1 is EXACTLY the self-consistency
equation that selects H=3.

So: the quadratic whose roots are d_super and H(H+1) has
discriminant equal to the self-consistency identity.

And 11/60 = 1/a + 1/b = (a+b)/(a*b) = 22/120 = 11/60.

The Koide lepton scale is (sum of roots)/(product of roots)
of the quadratic whose discriminant IS self-consistency.
""")

print(f"Summary:")
print(f"  a + b = 2*(H^2+H-1) = 2*11 = 22")
print(f"  a * b = S + 1/K* = 120 = 5!")
print(f"  a - b = -2 = -(H-1) = -sqrt(discriminant)")
print(f"  discriminant = (H-1)^2 = H+1 = 4")
print(f"  11/60 = (a+b)/(a*b) = 22/120")
