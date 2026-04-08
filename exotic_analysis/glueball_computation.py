import numpy as np
from scipy.optimize import fsolve

# --- equilibrium (same as before, abbreviated) ---
def ds_combine(m, e):
    s_pre = m[:3]*e[:3] + m[:3]*e[3] + m[3]*e[:3]
    theta_pre = m[3]*e[3]
    K = 1.0 - np.sum(s_pre) - theta_pre
    out = np.zeros(4)
    out[:3] = s_pre/(1-K); out[3] = theta_pre/(1-K)
    return out, K

def born_prob(m): return m[3]**2 / np.sum(m**2)

def enforce_floor(m):
    if born_prob(m) >= 1/27 - 1e-14: return m.copy()
    S = np.sum(m[:3]); Sq = np.sum(m[:3]**2)
    A = 26*S**2-Sq; B = 2*Sq; C = -Sq
    disc = B**2-4*A*C
    t1=(-B+np.sqrt(disc))/(2*A); t2=(-B-np.sqrt(disc))/(2*A)
    cands = [t for t in [t1,t2] if 0<t<1]
    t = min(cands, key=lambda x: abs(x-m[3]))
    out = np.zeros(4); out[:3]=m[:3]*(1-t)/S; out[3]=t
    return out

def ds_step(m,e): return enforce_floor(ds_combine(m,e)[0])

def find_eq():
    def eqs(p):
        s1,th,w1,phi = p
        s2=(1-s1-th)/2; w2=(1-w1-phi)/2
        m=np.array([s1,s2,s2,th]); e=np.array([w1,w2,w2,phi])
        eq1=th**2/(s1**2+2*s2**2+th**2)-1/27
        eq2=phi**2/(w1**2+2*w2**2+phi**2)-1/27
        eq3=2*s1*w2+2*s2*w1+2*s2*w2-7/30
        eq4=ds_step(m,e)[0]-s1
        return [eq1,eq2,eq3,eq4]
    sol=fsolve(eqs,[0.787,0.155,0.631,0.129])
    s1,th,w1,phi=sol; s2=(1-s1-th)/2; w2=(1-w1-phi)/2
    return np.array([s1,s2,s2,th]), np.array([w1,w2,w2,phi])

m_star, e_star = find_eq()

# ============================================================
# THE FLOOR KICK: F+ IS THE DISTANCE FROM LIGHT
# ============================================================

print("="*60)
print("F+ = 0  →  LIGHT  →  THE REFERENCE POINT")
print("="*60)

# The DS step (holomorphic, polynomial, integrable)
m_light, K = ds_combine(m_star, e_star)  # F+=0, this is light
m_vacuum = m_star.copy()                  # F+≠0, this is the vacuum

print(f"\nLight (post-DS, pre-floor):  {m_light}")
print(f"Vacuum (equilibrium):        {m_vacuum}")
print(f"Born(light) = {born_prob(m_light):.6f}")
print(f"Born(vacuum) = {born_prob(m_vacuum):.6f}")

# F+ IS the kick
delta_m = m_vacuum - m_light
print(f"\nδm = vacuum - light = ({delta_m[0]:.6f}, {delta_m[1]:.6f}, {delta_m[2]:.6f}, {delta_m[3]:.6f})")
print(f"  δθ  = {delta_m[3]:+.6f}  (ignorance GAINED)")
print(f"  δs₁ = {delta_m[0]:+.6f}  (dominant singleton LOST)")
print(f"  δs₂ = {delta_m[1]:+.6f}  (weak singleton lost)")
print(f"  δs₃ = {delta_m[2]:+.6f}  (weak singleton lost)")

# ============================================================
# PAULI EMBEDDING: δM = (δθ·I + δs·σ)/√2
# ============================================================

print("\n" + "="*60)
print("F+ THROUGH THE PAULI EMBEDDING")
print("="*60)

delta_theta = delta_m[3]
delta_s = delta_m[:3]

print(f"\nδM = (δθ·I + δs₁·σ₁ + δs₂·σ₂ + δs₃·σ₃)/√2")
print(f"\n  Trace (1,1):     δθ  = {delta_theta:+.6f}")
print(f"  Traceless adj:   δs  = ({delta_s[0]:+.6f}, {delta_s[1]:+.6f}, {delta_s[2]:+.6f})")
print(f"  |δθ|  = {abs(delta_theta):.6f}")
print(f"  |δs|  = {np.linalg.norm(delta_s):.6f}")

# ============================================================
# BILINEAR DECOMPOSITION: F+⊗F+ → GLUEBALL QUANTUM NUMBERS
# ============================================================

print("\n" + "="*60)
print("F+⊗F+ DECOMPOSITION INTO GLUEBALL CHANNELS")
print("="*60)

print("""
Glueballs are BILINEAR in F+.
  0++ ↔ tr(F²) = |F+|²_total = scalar contraction
  2++ ↔ symmetric traceless of s⊗s = spin-2 tensor
  0-+ ↔ tr(F∧F̃) = chirality-breaking mode (needs complex extension)

The vector δs = (δs₁, δs₂, δs₃) decomposes bilinearly:
  3⊗3 = 1 ⊕ 3 ⊕ 5
       = (spin-0) ⊕ (spin-1) ⊕ (spin-2)
""")

# Spin-0 of s⊗s: (1/3)|s|²
s_norm_sq = np.sum(delta_s**2)
T0_from_s = s_norm_sq / 3.0
T0_total = delta_theta**2 + T0_from_s  # total scalar content

print(f"|δs|² = {s_norm_sq:.8f}")
print(f"|δθ|² = {delta_theta**2:.8f}")

# Spin-2 of s⊗s: T_ij = s_i·s_j - (1/3)|s|²·δ_ij
T2 = np.outer(delta_s, delta_s) - (s_norm_sq/3)*np.eye(3)
T2_norm_sq = np.sum(T2**2)

print(f"\nSpin-0 content of s⊗s: |T₀|² = |s|⁴/3 = {s_norm_sq**2/3:.8f}")
print(f"Spin-2 content of s⊗s: |T₂|² = {T2_norm_sq:.8f}")
print(f"Ratio |T₂|²/|T₀|² = {T2_norm_sq/(s_norm_sq**2/3):.6f}")

ratio_sq = T2_norm_sq / (s_norm_sq**2/3)
ratio = np.sqrt(ratio_sq)

print(f"\n  m(2++)/m(0++) = √(|T₂|²/|T₀|²) = √{ratio_sq:.4f} = {ratio:.6f}")
print(f"  Lattice:                                          1.40")
print(f"  Deviation:                                        {abs(ratio-1.40)/1.40*100:.1f}%")

# ============================================================
# WHY IT'S √2
# ============================================================

print("\n" + "="*60)
print("WHY √2")
print("="*60)

print(f"""
At the equilibrium, s₁ = {m_star[0]:.4f} ≫ s₂ = s₃ = {m_star[1]:.4f}.
The floor kick δs is dominated by δs₁ = {delta_s[0]:.4f}.

In the limit δs ≈ (v, 0, 0):

  Spin-0 of s⊗s:  |T₀|² = |v|⁴/3

  Spin-2 of s⊗s:  T = diag(2v²/3, -v²/3, -v²/3)
                   |T₂|² = (2/3)²v⁴ + 2·(1/3)²v⁴
                         = (4/9 + 2/9)v⁴ = (2/3)v⁴

  Ratio: |T₂|²/|T₀|² = (2/3)v⁴ / (v⁴/3) = 2

  m(2++)/m(0++) = √2 = 1.41421...

This is STRUCTURAL. It follows from:
  1. H = 3 forces one dominant hypothesis at equilibrium
  2. The dominant hypothesis makes δs nearly one-dimensional  
  3. The spin-2/spin-0 norm ratio for a one-dimensional vector is exactly 2
  4. Mass ratio = √(norm ratio) = √2

Exact value at the actual equilibrium: √{ratio_sq:.6f} = {ratio:.6f}
The correction from √2 is {abs(ratio - np.sqrt(2))/np.sqrt(2)*100:.2f}%
(from the small but nonzero δs₂ = δs₃ = {delta_s[1]:.6f})
""")

# ============================================================
# NOW THE A₂ COUPLED SYSTEM
# ============================================================

print("="*60)
print("A₂ COUPLED SYSTEM: ALL FOUR STATES")
print("="*60)

g = 7/30; uniform = np.array([0.25]*4); eps = 1e-9

def coupled_ev(m1,m2,e0):
    e1=e0-g*(m2-uniform); e1=np.maximum(e1,1e-10); e1/=np.sum(e1); return e1

def coupled_step(state,e0):
    m1,m2=state[:4],state[4:]
    e1=coupled_ev(m1,m2,e0); e2=coupled_ev(m2,m1,e0)
    return np.concatenate([ds_step(m1,e1), ds_step(m2,e2)])

state = np.concatenate([m_star,m_star])
for i in range(5000):
    sn=coupled_step(state,e_star)
    if np.linalg.norm(sn-state)<1e-14: break
    state=sn

m1_eq = state[:4]; m2_eq = state[4:]

# Coupled light state (DS without floor at each node)
e1_eq = coupled_ev(m1_eq, m2_eq, e_star)
e2_eq = coupled_ev(m2_eq, m1_eq, e_star)
m1_light, K1 = ds_combine(m1_eq, e1_eq)
m2_light, K2 = ds_combine(m2_eq, e2_eq)

# Floor kick at each node
delta1 = m1_eq - m1_light
delta2 = m2_eq - m2_light

print(f"\nNode 1 light:  {m1_light}")
print(f"Node 1 vacuum: {m1_eq}")
print(f"Node 1 kick:   ({delta1[0]:+.4f}, {delta1[1]:+.4f}, {delta1[2]:+.4f}, {delta1[3]:+.4f})")

# F+ content at each node
ds1 = delta1[:3]; dth1 = delta1[3]
ds2 = delta2[:3]; dth2 = delta2[3]

# Combined s vector (both nodes)
ds_total = np.concatenate([ds1, ds2])

# Spin-2/Spin-0 ratio for node 1
s1_sq = np.sum(ds1**2)
T2_1 = np.outer(ds1, ds1) - (s1_sq/3)*np.eye(3)
T2_1_sq = np.sum(T2_1**2)
T0_1_sq = s1_sq**2 / 3

ratio_1 = np.sqrt(T2_1_sq / T0_1_sq) if T0_1_sq > 0 else 0

print(f"\nNode 1: |δs|² = {s1_sq:.6f}, |T₂|²/|T₀|² = {T2_1_sq/T0_1_sq:.4f}, √ratio = {ratio_1:.4f}")

# Same for node 2 (should be identical by symmetry)
s2_sq = np.sum(ds2**2)
T2_2 = np.outer(ds2, ds2) - (s2_sq/3)*np.eye(3)
T2_2_sq = np.sum(T2_2**2)
T0_2_sq = s2_sq**2 / 3
ratio_2 = np.sqrt(T2_2_sq / T0_2_sq) if T0_2_sq > 0 else 0

print(f"Node 2: |δs|² = {s2_sq:.6f}, |T₂|²/|T₀|² = {T2_2_sq/T0_2_sq:.4f}, √ratio = {ratio_2:.4f}")

# Jacobian eigenvalues for 0-+ and 0++*
J8 = np.zeros((8,8))
for j in range(8):
    sp=state.copy(); sp[j]+=eps; sm=state.copy(); sm[j]-=eps
    J8[:,j]=(coupled_step(sp,e_star)-coupled_step(sm,e_star))/(2*eps)

evals = np.linalg.eigvals(J8)
evals = np.sort(np.abs(evals))[::-1]
evals_sig = [e for e in evals if e > 1e-6]

Delta_0 = -np.log(evals_sig[0])

print(f"\nA₂ eigenvalues: {[f'{e:.4f}' for e in evals_sig]}")
print(f"Δ₀ = {Delta_0:.4f}")

# ============================================================
# COMPLETE SPECTRUM
# ============================================================

print("\n" + "="*60)
print("COMPLETE GLUEBALL SPECTRUM")  
print("="*60)

# 0++ from ground eigenvalue
m_0pp = Delta_0
# 2++ from √2 bilinear ratio
m_2pp = Delta_0 * ratio_1
# 0-+ from third eigenvalue (node-symmetric, s₂-s₃ antisymmetric)
m_0mp = -np.log(evals_sig[2])
# 0++* from fourth eigenvalue (node-symmetric, radial)
m_0pps = -np.log(evals_sig[3])

print(f"""
  State    Δ        m/m(0++)   Lattice   Dev

  0++      {m_0pp:.4f}    1.000      1.000     ---
  2++      {m_2pp:.4f}    {m_2pp/m_0pp:.4f}     1.40      {abs(m_2pp/m_0pp-1.40)/1.40*100:.1f}%
  0-+      {m_0mp:.4f}    {m_0mp/m_0pp:.4f}     1.50      {abs(m_0mp/m_0pp-1.50)/1.50*100:.1f}%
  0++*     {m_0pps:.4f}    {m_0pps/m_0pp:.4f}     1.56      {abs(m_0pps/m_0pp-1.56)/1.56*100:.1f}%

  The 2++ ratio is √(|T₂|²/|T₀|²) = √{T2_1_sq/T0_1_sq:.4f} = {ratio_1:.4f}
  In the one-dominant limit: √2 = {np.sqrt(2):.4f}
  Correction from equilibrium structure: {abs(ratio_1-np.sqrt(2))/np.sqrt(2)*100:.2f}%

  The 0-+ and 0++* are eigenvalues of the A₂ coupled Jacobian.
  The 2++ is the bilinear spin-2/spin-0 decomposition of F+.

  Zero free parameters.
""")

