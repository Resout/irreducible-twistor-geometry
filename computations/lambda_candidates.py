"""
Lambda candidates: systematic comparison.

Three routes, three numbers. Which is the geometric Lambda?
"""
import numpy as np
from scipy.optimize import brentq
from math import gcd, log10, exp

H = 3; FLOOR = 1.0/27.0; K_STAR = 7.0/30.0; S = 810.0/7.0

# Find equilibrium
def ds_step(m, e):
    s, th = m[:3], m[3]; se, ph = e[:3], e[3]
    sn = s*se + s*ph + th*se; tn = th*ph
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
    out = np.concatenate([sn, [tn]]) / (1-K)
    born = out[3]**2/np.sum(out**2)
    if born < FLOOR:
        ss=np.sum(np.abs(out[:3])); ssq=np.sum(out[:3]**2); r=ssq/ss**2
        a,b,c = 26-r, 2*r, -r
        t = (-b+np.sqrt(b**2-4*a*c))/(2*a)
        out = np.array([out[0]*abs(out[0])*(1-t)/(ss*abs(out[0])+1e-300) if abs(out[0])>0 else 0,
                        out[1]*abs(out[1])*(1-t)/(ss*abs(out[1])+1e-300) if abs(out[1])>0 else 0,
                        out[2]*abs(out[2])*(1-t)/(ss*abs(out[2])+1e-300) if abs(out[2])>0 else 0, t])
        # simpler: just scale proportionally
        sc = (1-t)/ss if ss > 0 else 0
        out = np.array([abs(m[0])*0+abs(out[0])*0, 0, 0, t])  # wrong, redo
    return out, K

# Use the working floor from before
def ds_combine(m, e):
    s, th = m[:3], m[3]; se, ph = e[:3], e[3]
    s_new = s*se + s*ph + th*se; th_new = th*ph
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
    return np.concatenate([s_new, [th_new]]) / (1-K), K

def floor_enforce(m):
    s = m[:3]; born = m[3]**2/np.sum(m**2)
    if born < FLOOR:
        ss=np.sum(s); r=np.sum(s**2)/ss**2
        a,b,c = 26-r, 2*r, -r
        t = (-b+np.sqrt(b**2-4*a*c))/(2*a)
        return np.array([s[0]*(1-t)/ss, s[1]*(1-t)/ss, s[2]*(1-t)/ss, t])
    return m

def step(m, e): return floor_enforce(ds_combine(m,e)[0])

def make_ev(p):
    pw=(1-p)/2; sc=1-FLOOR
    raw=np.array([np.sqrt(p*sc),np.sqrt(pw*sc),np.sqrt(pw*sc),np.sqrt(FLOOR)])
    return raw/np.sum(raw)

def fp(e):
    m=np.array([0.4,0.2,0.2,0.2])
    for _ in range(5000):
        m2=step(m,e)
        if np.max(np.abs(m2-m))<1e-15: break
        m=m2
    return m2

p_dom = brentq(lambda p: ds_combine(fp(make_ev(p)),make_ev(p))[1]-7/30, 0.92, 0.94, xtol=1e-15)
e_star = make_ev(p_dom); m_star = fp(e_star)
th = m_star[3]; s = m_star[:3]

det_M = (th**2 - np.sum(s**2)) / 2

print("=" * 60)
print("LAMBDA: SYSTEMATIC CANDIDATE COMPARISON")
print("=" * 60)
print()
print(f"Equilibrium: theta* = {th:.12f}")
print(f"  det(M*) = -25*theta*^2/2 = {det_M:.12f}")
print(f"  Born(theta*) = {th**2/np.sum(m_star**2):.10f}")
print()

# Candidate 1: conformal factor from fiber area
# Each CP1 fiber restricted to 26/27 of its area
# Omega^2 = 26/27, Lambda = 3/Omega^2 = 81/26
L1 = 81.0/26.0
print(f"Candidate 1: Lambda = 3*H^3/(H^3-1) = 81/26 = {L1:.10f}")
print(f"  From: conformal factor Omega^2 = 26/27 (fiber area ratio)")
print(f"  Meaning: the Born floor reduces each fiber by 1/27,")
print(f"  the base metric compensates by 27/26, Lambda increases slightly.")
print()

# Candidate 2: hierarchy tower prefactor
# g = K*(1-K*) = 161/900
L2 = 161.0/900.0
print(f"Candidate 2: Lambda = K*(1-K*) = 161/900 = {L2:.10f}")
print(f"  From: hierarchy tower at d=1: Lambda_phys = g*exp(-S)")
print(f"  Meaning: the Bernoulli variance of the structural filter.")
print(f"  PROBLEM: g appears in ALL hierarchy entries, not just CC.")
print()

# Candidate 3: shift from round
# Lambda_shift = Lambda1 - Lambda_round = 81/26 - 3 = 3/26
L3 = 3.0/26.0
print(f"Candidate 3: Lambda = 3/(H^3-1) = 3/26 = {L3:.10f}")
print(f"  From: shift of Lambda above the round S^4 value.")
print(f"  Meaning: the fractional change 1/(H^3-1) = 1/26.")
print()

# Physical Lambda for each
print("Physical Lambda (with instanton suppression exp(-S)):")
for name, L in [("81/26", L1), ("161/900", L2), ("3/26", L3)]:
    Lp = L * exp(-S)
    print(f"  Lambda={name:>8s}: Lp = {Lp:.4e},  log10 = {log10(Lp):.2f}")
print(f"  Hierarchy tower d=1: {161/900*exp(-S):.4e},  log10 = {log10(161/900*exp(-S)):.2f}")
print()

# The DISCRIMINATOR: can any candidate be derived from the
# infinity twistor identity I_AB Z^A Z^B = (Lambda/3)*det(M)?
#
# At the equilibrium: I(m*,m*) = (Lambda/3)*det(M*)
# det(M*) = -25*theta*^2/2 (known numerically)
#
# The infinity twistor I for de Sitter is I_AB = (Lambda/3)*epsilon_AB
# Acting on Z = (omega, pi): I(Z,Z) = (Lambda/3)*(omega^0 pi_1' - omega^1 pi_0')
# = (Lambda/3)*det(M)
#
# So I(m*,m*) = (Lambda/3)*det(M*) is AUTOMATICALLY satisfied for any Lambda.
# This doesn't DETERMINE Lambda!
#
# Lambda is determined by the NORMALIZATION of I relative to the metric.
# The metric on twistor space is the Hermitian form Sigma of signature (2,2).
# The normalization condition is: I_AB I^AB = Lambda^2/9 * epsilon * epsilon = ...
# This involves contracting I with itself using the twistor metric.

print("KEY INSIGHT: I(Z,Z) = (Lambda/3)*det(M) is automatic.")
print("It holds for ANY Lambda. Lambda is determined by the")
print("NORMALIZATION of I relative to the twistor metric Sigma.")
print()
print("The twistor metric Sigma(Z,Z) = 2*Re(omega.pi) = theta^2 - |s|^2")
print(f"At equilibrium: Sigma(m*,m*) = {th**2 - np.sum(s**2):.10f}")
print(f"= 2*det(M*) = {2*det_M:.10f}")
print()

# The normalization of the infinity twistor:
# |I|^2 = I_AB I^{AB} = (Lambda/3)^2 * epsilon_AB epsilon^{AB}
# epsilon_AB epsilon^{AB} = -2 (for the standard SL(2,C) epsilon)
# Wait, for a skew 4x4 tensor: epsilon_{ABCD} epsilon^{ABCD} = 24.
# For a skew 2-form I in Lambda^2(C^4):
# |I|^2 = I_{AB} I^{AB} where indices raised by the twistor metric.
# For I = (Lambda/3)*epsilon (where epsilon is the standard volume form):
# |I|^2 = (Lambda/3)^2 * |epsilon|^2.
# |epsilon|^2 depends on the twistor metric signature.
#
# For the standard (2,2) twistor metric:
# |epsilon|^2 = Pf(Sigma)^2 = det(Sigma) = ... this is getting complicated.

# ALTERNATIVE: use the ADAMO-MASON identification directly.
# Adamo-Mason 2014 says: Lambda enters through the infinity twistor I_AB,
# with Lambda proportional to |I|^2. The SPECIFIC proportionality constant
# depends on normalization conventions.
#
# In their conventions: I_AB = (Lambda)^{1/2} * epsilon_AB (schematic).
# The action term is (Lambda/3) * I * (field).
#
# The Born floor determines I by breaking conformal symmetry.
# The SCALE of the breaking is set by Born = 1/27.
#
# In the twistor metric: Sigma(m*,m*) = 2*det(M*) = -25*theta*^2.
# The NORMALIZED infinity twistor is I such that I(m*,m*)/Sigma(m*,m*) = Lambda/3.
# But I(m*,m*) = (Lambda/3)*det(M*) and Sigma(m*,m*) = 2*det(M*).
# So I/Sigma = (Lambda/3)*det(M*)/(2*det(M*)) = Lambda/6.
# This gives Lambda/6 = Lambda/6, a tautology!

print("The ratio I/Sigma at equilibrium is tautological (Lambda/6 = Lambda/6).")
print("The infinity twistor identity alone cannot determine Lambda.")
print()
print("Lambda is determined by HOW the Born floor breaks conformal symmetry,")
print("not by the algebraic identity I(Z,Z) = (Lambda/3)*det(M).")
print()

# THE RIGHT APPROACH: Lambda is the Ricci scalar / 4 of the equilibrium metric.
# The equilibrium metric is determined by the Penrose transform.
# For a fiber-uniform section: constant conformal factor.
# Lambda = 3/Omega^2 where Omega^2 = (restricted fiber area)/(full fiber area) = 26/27.
# Lambda = 81/26.
#
# This is Candidate 1. Let me check if it's consistent with the hierarchy tower.

print("="*60)
print("CONSISTENCY CHECK: Lambda = 81/26 vs hierarchy tower")
print("="*60)
print()
print(f"If Lambda_geometric = 81/26:")
print(f"  Lambda_phys = (81/26)*exp(-S) = {81/26*exp(-S):.4e}")
print(f"  log10 = {log10(81/26*exp(-S)):.2f}")
print()
print(f"Hierarchy tower d=1: g*exp(-S) = {161/900*exp(-S):.4e}")
print(f"  log10 = {log10(161/900*exp(-S)):.2f}")
print()
print(f"Ratio: (81/26) / (161/900) = {81/26/(161/900):.6f}")
n = 81*900; d = 26*161
g_cd = gcd(n,d)
print(f"  = {n//g_cd}/{d//g_cd}")
print()

# 81*900 = 72900, 26*161 = 4186
# 72900/4186 = 36450/2093
print(f"The ratio is NOT a simple integer or fraction.")
print(f"This means either:")
print(f"  (a) Lambda = 81/26 and the tower prefactor g is NOT Lambda")
print(f"  (b) Lambda = 161/900 and the conformal factor is NOT 26/27")
print(f"  (c) Neither is exactly right")
print()
print(f"Resolution: g = K*(1-K*) is the COUPLING, not Lambda.")
print(f"Lambda = 81/26 is the CURVATURE.")
print(f"The hierarchy tower formula is: ratio = g * exp(-S/d)")
print(f"where g is the coupling and exp(-S/d) is the instanton factor.")
print(f"The CC entry (d=1) gives g*exp(-S), where g is the coupling")
print(f"and exp(-S) is the full instanton suppression.")
print(f"This is NOT 'Lambda * exp(-S)' -- it is 'coupling * tunneling'.")
print()
print(f"Lambda and g are different quantities:")
print(f"  g = K*(1-K*) = 161/900 (equilibrium coupling)")
print(f"  Lambda = 3*H^3/(H^3-1) = 81/26 (spacetime curvature)")
print(f"Both exact, both from H=3, but different physical meanings.")
