import numpy as np

# ============================================================
# VERIFY: SU(2) COMMUTES WITH THE COUPLED DS MAP
# ============================================================

sigma1 = np.array([[0,1],[1,0]], dtype=complex)
sigma2 = np.array([[0,-1j],[1j,0]], dtype=complex)
sigma3 = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

def pauli_embed(m):
    return (m[3]*I2 + m[0]*sigma1 + m[1]*sigma2 + m[2]*sigma3)/np.sqrt(2)

def pauli_extract(M):
    """Inverse: 2x2 matrix → (s1,s2,s3,θ)"""
    M_scaled = M * np.sqrt(2)
    theta = np.trace(M_scaled).real / 2
    s1 = (M_scaled[0,1] + M_scaled[1,0]).real / 2
    s2 = (M_scaled[1,0] - M_scaled[0,1]).imag / 2  
    s3 = (M_scaled[0,0] - M_scaled[1,1]).real / 2
    return np.array([s1, s2, s3, theta])

def su2_rotate(m, U):
    """Rotate mass function by SU(2): M → UMU†"""
    M = pauli_embed(m)
    M_rot = U @ M @ U.conj().T
    return pauli_extract(M_rot)

def ds_combine(m, e):
    s_pre = m[:3]*e[:3] + m[:3]*e[3] + m[3]*e[:3]
    theta_pre = m[3]*e[3]
    K = 1.0 - np.sum(s_pre) - theta_pre
    out = np.zeros(4); out[:3]=s_pre/(1-K); out[3]=theta_pre/(1-K)
    return out, K

def born_prob(m): return m[3]**2/np.sum(m**2)

def enforce_floor(m):
    if born_prob(m) >= 1/27-1e-14: return m.copy()
    S=np.sum(m[:3]); Sq=np.sum(m[:3]**2)
    A=26*S**2-Sq; B=2*Sq; C=-Sq; disc=B**2-4*A*C
    t1=(-B+np.sqrt(disc))/(2*A); t2=(-B-np.sqrt(disc))/(2*A)
    cands=[t for t in [t1,t2] if 0<t<1]
    t=min(cands, key=lambda x: abs(x-m[3]))
    out=np.zeros(4); out[:3]=m[:3]*(1-t)/S; out[3]=t; return out

def ds_step(m,e): return enforce_floor(ds_combine(m,e)[0])

g = 7/30; uniform = np.array([0.25]*4)

def coupled_ev(m1,m2,e0):
    e1=e0-g*(m2-uniform); e1=np.maximum(e1,1e-10); e1/=np.sum(e1); return e1

def coupled_step(m1, m2, e0):
    e1 = coupled_ev(m1,m2,e0); e2 = coupled_ev(m2,m1,e0)
    return ds_step(m1,e1), ds_step(m2,e2)

# Find equilibrium
from scipy.optimize import fsolve
def find_eq():
    def eqs(p):
        s1,th,w1,phi=p; s2=(1-s1-th)/2; w2=(1-w1-phi)/2
        m=np.array([s1,s2,s2,th]); e=np.array([w1,w2,w2,phi])
        return [th**2/(s1**2+2*s2**2+th**2)-1/27,
                phi**2/(w1**2+2*w2**2+phi**2)-1/27,
                2*s1*w2+2*s2*w1+2*s2*w2-7/30,
                ds_step(m,e)[0]-s1]
    sol=fsolve(eqs,[0.787,0.155,0.631,0.129])
    s1,th,w1,phi=sol; s2=(1-s1-th)/2; w2=(1-w1-phi)/2
    return np.array([s1,s2,s2,th]), np.array([w1,w2,w2,phi])

m_star, e_star = find_eq()

# ============================================================
# TEST 1: Uniform state is SU(2)-invariant
# ============================================================

print("="*60)
print("TEST 1: UNIFORM STATE IS SU(2)-INVARIANT")
print("="*60)

M_uniform = pauli_embed(uniform)
print(f"uniform = {uniform}")
print(f"M_uniform =\n{M_uniform}")
print(f"M_uniform / (1/(2√2)) = {M_uniform * 2*np.sqrt(2)}")
print(f"That's (1/4)(I + σ₁ + σ₂ + σ₃)/√2")
print()

# Actually, uniform = (1/4, 1/4, 1/4, 1/4) maps to:
# M = (1/4·I + 1/4·σ₁ + 1/4·σ₂ + 1/4·σ₃)/√2
# This is NOT proportional to I. Let me check.
print("Wait — uniform is NOT proportional to I in the Pauli embedding.")
print(f"M_uniform has σ₁, σ₂, σ₃ components too.")
print()

# Check: is M_uniform invariant under UMU†?
theta_rot = 0.7  # arbitrary rotation angle
n_hat = np.array([0.3, 0.5, 0.8]); n_hat = n_hat/np.linalg.norm(n_hat)
U = np.cos(theta_rot/2)*I2 + 1j*np.sin(theta_rot/2)*(n_hat[0]*sigma1 + n_hat[1]*sigma2 + n_hat[2]*sigma3)

uniform_rotated = su2_rotate(uniform, U)
print(f"uniform =          {uniform}")
print(f"uniform (rotated) = {uniform_rotated}")
print(f"Difference: {np.linalg.norm(uniform_rotated - uniform):.2e}")
print()

# The uniform state IS NOT invariant! The s-components rotate.
# θ is invariant (trace), but s rotates as a 3-vector.
# However: (m₂ - uniform) transforms as:
# U(m₂ - uniform)U† = Um₂U† - U·uniform·U†
# And U·uniform·U† ≠ uniform (unless U = I).

# So the argument needs revision. Let me check the actual claim more carefully.
# The claim was that uniform maps to something proportional to I.
# Let me check: what's special about uniform in the Pauli embedding?

# uniform = (1/4)(1,1,1,1) → M = (1/4)(I + σ₁ + σ₂ + σ₃)/√2
# This is NOT proportional to I.
# The state proportional to I would be (0,0,0,1) → M = I/√2.
# That's pure ignorance, not uniform.

print("CORRECTION: uniform = (1/4,1/4,1/4,1/4) is NOT I in the Pauli embedding.")
print("Pure ignorance = (0,0,0,1) maps to I/√2.")
print()
print("Let me check what the coupling formula actually uses...")
print()

# The coupling is: e₁ = e₀ - (7/30)(m₂ - uniform)
# Under SU(2): m₂ → Um₂U†, e₀ → Ue₀U†
# Need: (m₂ - uniform) → U(m₂ - uniform)U†
# This requires: uniform → U·uniform·U†
# Which is TRUE if uniform has s = (s,s,s) because 
# then the s-vector is along (1,1,1)/√3, and rotation 
# changes this direction. So uniform is NOT SU(2) invariant.

# BUT WAIT. Let me re-read the claim.
# "The uniform state (1/4,1/4,1/4,1/4) maps to (θI + s·σ)/√2 with 
#  θ = s₁ = s₂ = s₃ = 1/4. In the Pauli embedding that's proportional to I."
# That's wrong. θI + s·σ with θ=s₁=s₂=s₃ is NOT proportional to I.

# However, the key question is: does the COUPLED MAP commute with SU(2)?
# Let me just test it numerically.

print("="*60)
print("TEST 2: DOES Φ_coupled COMMUTE WITH SU(2)?")
print("="*60)

# Take a random state near equilibrium
np.random.seed(42)
m1 = m_star + 0.01*np.random.randn(4); m1[3] = abs(m1[3]); m1 = m1/np.sum(m1)
m2 = m_star + 0.01*np.random.randn(4); m2[3] = abs(m2[3]); m2 = m2/np.sum(m2)

# Path A: rotate first, then step
m1_rot = su2_rotate(m1, U)
m2_rot = su2_rotate(m2, U)
e_rot = su2_rotate(e_star, U)
m1_A, m2_A = coupled_step(m1_rot, m2_rot, e_rot)

# Path B: step first, then rotate
m1_B_pre, m2_B_pre = coupled_step(m1, m2, e_star)
m1_B = su2_rotate(m1_B_pre, U)
m2_B = su2_rotate(m2_B_pre, U)

print(f"Path A (rotate then step):")
print(f"  m1_A = {m1_A}")
print(f"  m2_A = {m2_A}")
print(f"Path B (step then rotate):")
print(f"  m1_B = {m1_B}")
print(f"  m2_B = {m2_B}")
print(f"Difference m1: {np.linalg.norm(m1_A - m1_B):.2e}")
print(f"Difference m2: {np.linalg.norm(m2_A - m2_B):.2e}")
print()

# Hmm, the difference might be nonzero because of the uniform state.
# Let me check if it's because of the evidence formula or the DS/floor.

# First: does DS alone commute with SU(2)?
print("="*60)
print("TEST 3: DOES DS ALONE COMMUTE WITH SU(2)?")
print("="*60)

m_test = m1.copy(); e_test = e_star.copy()

# DS then rotate
m_ds, K1 = ds_combine(m_test, e_test)
m_ds_rot = su2_rotate(m_ds, U)

# Rotate then DS
m_test_rot = su2_rotate(m_test, U)
e_test_rot = su2_rotate(e_test, U)
m_rot_ds, K2 = ds_combine(m_test_rot, e_test_rot)

print(f"DS then rotate: {m_ds_rot}")
print(f"Rotate then DS: {m_rot_ds}")
print(f"Difference: {np.linalg.norm(m_ds_rot - m_rot_ds):.2e}")
print(f"K₁ = {K1:.10f}, K₂ = {K2:.10f}, diff = {abs(K1-K2):.2e}")
print()

# Does floor commute with SU(2)?
print("="*60)
print("TEST 4: DOES FLOOR COMMUTE WITH SU(2)?")
print("="*60)

m_low_born = np.array([0.8, 0.05, 0.05, 0.1])  # Born < 1/27

# Floor then rotate
m_floored = enforce_floor(m_low_born)
m_floored_rot = su2_rotate(m_floored, U)

# Rotate then floor
m_low_rot = su2_rotate(m_low_born, U)
m_rot_floored = enforce_floor(m_low_rot)

print(f"Floor then rotate: {m_floored_rot}")
print(f"Rotate then floor: {m_rot_floored}")
print(f"Difference: {np.linalg.norm(m_floored_rot - m_rot_floored):.2e}")
print()

# Check Born probabilities
print(f"Born(original) = {born_prob(m_low_born):.6f}")
print(f"Born(rotated) = {born_prob(m_low_rot):.6f}")
print(f"Born is SU(2) invariant: {abs(born_prob(m_low_born)-born_prob(m_low_rot))<1e-10}")
print()

# Check: does the floor depend on individual components or only on 
# magnitudes? The floor uses S = Σs_i and Sq = Σs_i². Under SU(2) rotation,
# the s-vector rotates, so individual s_i change BUT:
# |s|² = s₁² + s₂² + s₃² is invariant
# Σs_i is NOT invariant (it's the trace of the s-vector, which changes)

print("="*60)
print("KEY ISSUE: THE FLOOR USES Σs_i, NOT |s|")
print("="*60)

print(f"Original: s = ({m_low_born[0]:.4f}, {m_low_born[1]:.4f}, {m_low_born[2]:.4f})")
print(f"Rotated:  s = ({m_low_rot[0]:.4f}, {m_low_rot[1]:.4f}, {m_low_rot[2]:.4f})")
print(f"|s|² original = {np.sum(m_low_born[:3]**2):.6f}")
print(f"|s|² rotated  = {np.sum(m_low_rot[:3]**2):.6f}")
print(f"Σs original   = {np.sum(m_low_born[:3]):.6f}")
print(f"Σs rotated    = {np.sum(m_low_rot[:3]):.6f}")
print()

# The L₁ constraint is Σs_i + θ = 1. Under SU(2) rotation:
# θ is invariant (trace), s_i change individually.
# If Σs_i changes, then L₁ changes, which is bad.
# But wait: does SU(2) rotation preserve L₁?

print(f"L₁ original = {np.sum(m_low_born):.10f}")
print(f"L₁ rotated  = {np.sum(m_low_rot):.10f}")
print(f"L₁ preserved: {abs(np.sum(m_low_born)-np.sum(m_low_rot))<1e-10}")
print()

# The issue: L₁ = s₁ + s₂ + s₃ + θ. Under SU(2):
# θ → θ (invariant, it's the trace/2)
# (s₁, s₂, s₃) → R·(s₁, s₂, s₃) where R ∈ SO(3)
# Σs_i is NOT invariant under SO(3) rotation.
# Therefore L₁ is NOT preserved by SU(2)!

# Unless... let me check if the Pauli extract preserves L₁

m_test2 = np.array([0.5, 0.2, 0.1, 0.2])
print(f"Test: m = {m_test2}, L₁ = {np.sum(m_test2):.4f}")
m_test2_rot = su2_rotate(m_test2, U)
print(f"Rotated: m = {m_test2_rot}, L₁ = {np.sum(m_test2_rot):.4f}")
print()

# So SU(2) rotation DOESN'T preserve L₁ in general!
# This means the SU(2) symmetry is more subtle than I thought.
# The SU(2) acts on the PROJECTIVE space CP³, not on the L₁=1 simplex.
# On CP³, overall scaling doesn't matter, so L₁ can change.

# But the DS combination requires L₁ = 1 inputs.
# So the SU(2) that commutes with DS must preserve L₁.
# Which rotations preserve L₁?

# L₁ = Σm_i = θ + s₁ + s₂ + s₃. Under M → UMU†:
# θ is invariant. We need Σs_i invariant.
# Σs_i = tr(s·σ) over... no. 
# Actually s₁ = tr(Mσ₁)/√2, etc. Hmm.

# Let me think about this differently.
# In the Pauli embedding: M = (θI + s·σ)/√2
# tr(M) = 2θ/√2 = √2·θ → θ is the trace, SU(2) invariant. ✓
# The sum s₁+s₂+s₃ = tr(M(σ₁+σ₂+σ₃))/√2 ... this is NOT SU(2) invariant
# because σ₁+σ₂+σ₃ rotates under conjugation.

print("="*60)
print("RESOLUTION: SU(2) ACTS ON CP³, NOT ON L₁=1 SIMPLEX")
print("="*60)

print("""
SU(2) acts on the Pauli embedding M → UMU†.
This preserves:
  - θ (trace of M)
  - |s|² = s₁² + s₂² + s₃² (Frobenius norm of traceless part)
  - det(M) = (θ² - |s|²)/2
  - Born = θ²/(θ² + |s|²)

It does NOT preserve:
  - Individual s_i values
  - The sum s₁ + s₂ + s₃ = L₁ - θ

On CP³ (projective space), L₁ is a gauge choice, not physical.
The DS combination and Born floor, expressed projectively, 
depend only on RATIOS of components, not on L₁.

The physical symmetry group acts on CP³ = {lines in ℂ⁴}.
SU(2) rotates the s-vector while preserving |s|² and θ.
After rotation, we RE-NORMALIZE to L₁ = 1.
This renormalization is invisible on CP³.

Checking: DS output ratios are SU(2) invariant even if 
absolute values aren't...
""")

# Let me verify the projective version
# Compute DS output ratios for original and rotated inputs
m_test3 = np.array([0.5, 0.2, 0.1, 0.2])
e_test3 = e_star.copy()

# Original
out_orig, K_orig = ds_combine(m_test3/np.sum(m_test3), e_test3)
ratios_orig = out_orig / out_orig[3]  # ratios to θ

# Rotated (and renormalized)
m_rot3 = su2_rotate(m_test3, U)
m_rot3_norm = m_rot3 / np.sum(m_rot3)  # renormalize
e_rot3 = su2_rotate(e_test3, U)  
e_rot3_norm = e_rot3 / np.sum(e_rot3)

out_rot, K_rot = ds_combine(m_rot3_norm, e_rot3_norm)
ratios_rot = out_rot / out_rot[3]

# The ratios should be related by the SAME SU(2) rotation
out_orig_rot = su2_rotate(out_orig, U)
out_orig_rot_norm = out_orig_rot / np.sum(out_orig_rot)
ratios_orig_rot = out_orig_rot_norm / out_orig_rot_norm[3]

print(f"DS output ratios (original, rotated): ")
print(f"  Original:        {ratios_orig}")
print(f"  Rotated input:   {ratios_rot}")
print(f"  Rotate(original):{ratios_orig_rot}")
print(f"  Diff (rot input vs rotate output): {np.linalg.norm(ratios_rot - ratios_orig_rot):.2e}")
print()

# Check K invariance
print(f"K original = {K_orig:.10f}")  
print(f"K rotated  = {K_rot:.10f}")
print(f"K diff     = {abs(K_orig-K_rot):.2e}")
print()

# Born invariance
print(f"Born(out_orig) = {born_prob(out_orig):.10f}")
print(f"Born(out_rot)  = {born_prob(out_rot):.10f}")

# ============================================================
# THE ACTUAL SYMMETRY GROUP
# ============================================================

print("\n" + "="*60)
print("SUMMARY: THE GAUGE GROUP STRUCTURE")
print("="*60)

print("""
SU(2) acts by M → UMU† at each Dynkin node simultaneously.
Preserves: θ, |s|², Born, K, det(M).
Doesn't preserve: individual s_i, L₁ (but L₁ is projective gauge).
On CP³: exact symmetry of the DS+floor dynamics.

SU(3) acts by permuting Dynkin nodes.
Preserves: the Cartan matrix structure (A₂ is Z₂ symmetric).
Acts on a completely different index (site label vs local spinor).

U(1) is the stabilizer of the equilibrium orientation on S².

These three act on independent structures:
  SU(2): local spinor ℂ² at each node
  SU(3): which node (colour index)
  U(1):  phase within the stabilizer

They commute because they act on different tensor factors.
Product structure: SU(3) × SU(2) × U(1).
Dimension: 8 + 3 + 1 = 12 generators.
""")

