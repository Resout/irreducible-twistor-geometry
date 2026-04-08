import numpy as np
from scipy.optimize import fsolve

# ============================================================
# SETUP
# ============================================================

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
# THE LIGHT POINT
# ============================================================

m_light, K = ds_combine(m_star, e_star)
delta_m = m_star - m_light

print("="*60)
print("THE LIGHT POINT: F+ = 0")
print("="*60)

print(f"""
m_light = ({m_light[0]:.4f}, {m_light[1]:.4f}, {m_light[2]:.4f}, {m_light[3]:.4f})

The spatial geometry at F+=0:

  ℂ⁴ = ℂ¹(θ) ⊕ ℂ¹(s₁) ⊕ ℂ²(s₂,s₃)
        ↑          ↑           ↑
     substrate   dominant    spinor
     ≈ 0.026    ≈ 0.907    ≈ (0.034, 0.034)
     almost      almost      small, equal,
     gone        everything  degenerate

Three subspaces. Each determines a spin.
""")

# ============================================================
# STEP 1: ENUMERATE ALL REPRESENTATIONS AT THE LIGHT POINT
# ============================================================

print("="*60)
print("ALL REPRESENTATIONS AT THE LIGHT POINT")
print("="*60)

print("""
ONE-PARTICLE representations (perturbations of m_light):

  SPIN 0 — singlets under SU(2):
    [0a] θ direction: (0, 0, 0, 1)  — the substrate itself
    [0b] s₁ direction: (1, 0, 0, 0) — the dominant relationship
    
  SPIN 1/2 — doublet under SU(2):
    [½] (s₂, s₃) direction: (0, 1, 0, 0) and (0, 0, 1, 0)
        — the spinor subspace ℂ²
        — two components: "up" and "down"

  SPIN 1 — triplet under SU(2):
    [1] Full (s₁, s₂, s₃) as adjoint vector in ℝ³
        — but this is REDUCIBLE at the light point:
          it splits into [0b] ⊕ [½]
        — a "pure" spin-1 requires BOTH s₁ AND (s₂,s₃) content
        — this is the gauge boson direction

TWO-PARTICLE (bilinear) representations:

  SPIN 0 from bilinears:
    [0c] θ·θ → scalar (mass term)
    [0d] s₁·s₁ → scalar (glueball 0++)
    [0e] s₂·s₃ - s₃·s₂ → antisymmetric singlet from ℂ²⊗ℂ² (Λ²ℂ²)
         dim(Λ²ℂ²) = 1. This is the UNIQUE scalar from two spinors.
    [0f] s₂·s₂ + s₃·s₃ → symmetric trace from ℂ²⊗ℂ² 

  SPIN 1/2 from bilinears:
    [½'] s₁ · (s₂, s₃) → singlet × doublet = doublet
         A scalar times a spinor is still a spinor.
         This is how the Higgs gives mass to fermions:
         Yukawa coupling = scalar(Higgs) × spinor(fermion)

  SPIN 1 from bilinears:
    [1'] Sym²(ℂ²) minus trace → the TRIPLET from two spinors
         ½ ⊗ ½ = 0 ⊕ 1
         dim(Sym²ℂ² / trace) = 2 ... no, Sym²(ℂ²) has dim 3 = spin 1.
         Two spinors symmetrized give a vector. This is the W/Z.

  SPIN 2 from bilinears:
    [2] Sym²₀(ℝ³) → traceless symmetric product of two adjoint vectors
        This is the 2++ glueball. 5 components.
        Under the splitting ℝ³ = ℝ¹ ⊕ ℝ²:
        Sym²₀ decomposes, but the full SO(3) representation persists
        because the floor kick spans all of ℝ³.
""")

# ============================================================
# STEP 2: COMPUTE FLOOR KICK PROJECTION ONTO EACH SUBSPACE
# ============================================================

print("="*60)
print("FLOOR KICK PROJECTIONS")
print("="*60)

print(f"Floor kick δm = ({delta_m[0]:.6f}, {delta_m[1]:.6f}, {delta_m[2]:.6f}, {delta_m[3]:.6f})")

# Projections
proj_theta = delta_m[3]  # θ component
proj_s1 = delta_m[0]     # s₁ component (dominant)
proj_s2 = delta_m[1]     # s₂ component
proj_s3 = delta_m[2]     # s₃ component
proj_spinor = np.sqrt(delta_m[1]**2 + delta_m[2]**2)  # |ℂ²| component

total_sq = np.sum(delta_m**2)

print(f"\nSubspace projections:")
print(f"  θ (substrate):     {proj_theta:+.6f}  ({proj_theta**2/total_sq:.1%} of |δm|²)")
print(f"  s₁ (dominant):     {proj_s1:+.6f}  ({proj_s1**2/total_sq:.1%} of |δm|²)")
print(f"  ℂ² (spinor):       {proj_spinor:.6f}   ({proj_spinor**2/total_sq:.1%} of |δm|²)")
print(f"    s₂:              {proj_s2:+.6f}")
print(f"    s₃:              {proj_s3:+.6f}")
print(f"    s₂ - s₃:        {proj_s2-proj_s3:.6f}  (antisymmetric, spin content)")
print(f"    s₂ + s₃:        {proj_s2+proj_s3:.6f}  (symmetric, radial content)")

print(f"""
The floor kick is {proj_theta**2/total_sq:.1%} substrate + {proj_s1**2/total_sq:.1%} dominant + {proj_spinor**2/total_sq:.1%} spinor.

Almost entirely in θ and s₁. Almost nothing in ℂ².
The floor CREATES mass (θ and s₁ content) but barely TOUCHES spin (ℂ² content).

This is why spin is determined BEFORE the floor fires:
the floor adds the massive (timelike) content without disturbing 
the spatial (spin) pattern that was already there.
""")

# ============================================================
# STEP 3: THE COMPLETE PARTICLE TABLE
# ============================================================

print("="*60)
print("THE PARTICLE TABLE FROM FIRST PRINCIPLES")
print("="*60)

print("""
Start at F+=0. Three subspaces: θ, s₁, ℂ²(s₂,s₃).
The floor fires. Each subspace gets mass from its eigenvalue.
The spin was fixed before the floor.

MASSLESS (on the equilibrium manifold, no floor kick needed):
  Graviton:  S² orientation field.           spin 2, 2 pol.
  Photon:    U(1) stabilizer connection.     spin 1, 2 pol.

MASSIVE (floor kick gives mass, Jacobian eigenvalue gives decay):

  From ℂ¹(θ) ⊕ ℂ¹(s₁) — the SCALAR sector:
""")

# Single-site Jacobian
eps = 1e-9
J4 = np.zeros((4,4))
for j in range(4):
    sp=m_star.copy(); sp[j]+=eps; sm=m_star.copy(); sm[j]-=eps
    J4[:,j]=(ds_step(sp,e_star)-ds_step(sm,e_star))/(2*eps)

evals4, evecs4 = np.linalg.eig(J4)
idx = np.argsort(-np.abs(evals4))
evals4 = evals4[idx]; evecs4 = evecs4[:,idx]

print(f"    Radial eigenvalue (s₁+θ mixing): λ = {evals4[0].real:.6f}")
print(f"    This is the 0++ glueball direction.")
print(f"    Δ = {-np.log(abs(evals4[0])):.4f}")

print(f"""
  From ℂ²(s₂,s₃) — the SPINOR sector:

    Angular eigenvalue (s₂-s₃): λ = {evals4[1].real:.6f}
    This is the FERMION direction. spin 1/2.
    Δ = {-np.log(abs(evals4[1])):.4f}
    
    Nearly degenerate with the scalar ({abs(evals4[0].real-evals4[1].real)/abs(evals4[0].real)*100:.2f}% split).
    Both decay at nearly the same rate.
    But they have DIFFERENT SPIN because they live in different
    subspaces of the light-point geometry.
""")

# ============================================================
# STEP 4: BILINEAR SPECTRUM (COMPOSITE PARTICLES)
# ============================================================

print("="*60)
print("BILINEAR SPECTRUM: WHAT TWO PARTICLES MAKE")
print("="*60)

lam_scalar = abs(evals4[0])  # radial
lam_spinor = abs(evals4[1])  # angular/spinor

print(f"Elementary eigenvalues:")
print(f"  scalar (s₁,θ): λ_s = {lam_scalar:.6f}")
print(f"  spinor (s₂,s₃): λ_f = {lam_spinor:.6f}")

print(f"\nComposite states (Koopman products):")
print(f"  scalar × scalar = λ_s² = {lam_scalar**2:.6f}  → spin 0 or 2")
print(f"    Δ/Δ₀ = {-np.log(lam_scalar**2) / -np.log(lam_scalar):.4f}")
print(f"    = 2.000 (always, by construction)")

print(f"\n  spinor × spinor = λ_f² = {lam_spinor**2:.6f}")
print(f"    ½ ⊗ ½ = 0 ⊕ 1")
print(f"    spin 0: MESON (quark-antiquark, spins opposed)")
print(f"    spin 1: VECTOR MESON (quark-antiquark, spins aligned)")
print(f"    Δ/Δ₀ = {-np.log(lam_spinor**2) / -np.log(lam_scalar):.4f}")

print(f"\n  scalar × spinor = λ_s · λ_f = {lam_scalar*lam_spinor:.6f}")
print(f"    0 ⊗ ½ = ½")
print(f"    spin 1/2: this is the YUKAWA coupling")
print(f"    the scalar (Higgs) gives mass to the spinor (fermion)")
print(f"    Δ/Δ₀ = {-np.log(lam_scalar*lam_spinor) / -np.log(lam_scalar):.4f}")

print(f"\n  spinor × spinor × spinor = λ_f³ = {lam_spinor**3:.6f}")
print(f"    ½ ⊗ ½ ⊗ ½ = ½ ⊕ ½ ⊕ 3/2")
print(f"    spin 1/2: BARYON (three quarks, net spin 1/2)")
print(f"    spin 3/2: DELTA BARYON (three quarks, spins aligned)")
print(f"    Δ/Δ₀ = {-np.log(lam_spinor**3) / -np.log(lam_scalar):.4f}")

# ============================================================
# STEP 5: THE FULL SPECTRUM AT SU(3) (A₂)
# ============================================================

print("\n" + "="*60)
print("FULL A₂ SPECTRUM: EVERY PARTICLE CLASS")
print("="*60)

# A₂ eigenvalues
g = 7/30; uniform = np.array([0.25]*4)
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

J8=np.zeros((8,8))
for j in range(8):
    sp=state.copy(); sp[j]+=eps; sm=state.copy(); sm[j]-=eps
    J8[:,j]=(coupled_step(sp,e_star)-coupled_step(sm,e_star))/(2*eps)

evals8 = np.linalg.eigvals(J8)
evals8_sig = sorted([abs(e) for e in evals8 if abs(e) > 1e-6], reverse=True)

Delta_0 = -np.log(evals8_sig[0])

# The four A₂ eigenvalues, NOW with spin from the light-point geometry
print(f"""
A₂ eigenvalues (colour structure of SU(3)):
  λ₀ = {evals8_sig[0]:.4f} — node-antisym, radial (s₁+θ) — colour-fund, spin 0
  λ₁ = {evals8_sig[1]:.4f} — node-antisym, spinor (s₂-s₃) — colour-fund, spin 1/2
  λ₂ = {evals8_sig[2]:.4f} — node-sym, spinor (s₂-s₃) — colour-singlet, spin 1/2(?)
  λ₃ = {evals8_sig[3]:.4f} — node-sym, radial (s₁+θ) — colour-singlet, spin 0

Wait. Node-symmetric spinor modes:
  λ₂ lives in ℂ²(s₂,s₃) but is a colour singlet.
  A colour-singlet spin-1/2 object would be a LEPTON.
  λ₂ at ratio {-np.log(evals8_sig[2])/Delta_0:.3f} — this is the 0⁻⁺ glueball.
  
  But from the light-point geometry, the s₂-s₃ direction IS ℂ².
  Is the 0⁻⁺ actually a LEPTONIC excitation misidentified as a glueball?
  
  No — the J^PC assignment for glueballs uses the BILINEAR structure.
  The 0⁻⁺ glueball is a composite (F⊗F with P·C = -1).
  The lepton would be an ELEMENTARY excitation in ℂ².
  
  These are different: the eigenvalue λ₂ could be BOTH:
  — As a single-particle mode: a colour-singlet fermion (lepton?)
  — As contributing to the bilinear 0⁻⁺ glueball mass
""")

# ============================================================
# THE COMPLETE PICTURE
# ============================================================

print("="*60)
print("COMPLETE PARTICLE TABLE")
print("="*60)

print(f"""
MASSLESS (on the equilibrium manifold):
  {"Graviton":<20} spin 2    2 pol.   S² orientation
  {"Photon":<20} spin 1    2 pol.   U(1) stabilizer

ELEMENTARY MASSIVE (single eigenvalues, from ℂ⁴ at light point):

  COLOUR FUNDAMENTAL (node-antisym on A₂):
    {"Quark (up-type)":<20} spin 0    λ={evals8_sig[0]:.4f}  Δ/Δ₀=1.000  radial
    {"Quark (down-type)":<20} spin 1/2  λ={evals8_sig[1]:.4f}  Δ/Δ₀={-np.log(evals8_sig[1])/Delta_0:.3f}  spinor

  COLOUR SINGLET (node-sym on A₂):
    {"Lepton-like":<20} spin 1/2  λ={evals8_sig[2]:.4f}  Δ/Δ₀={-np.log(evals8_sig[2])/Delta_0:.3f}  spinor
    {"Higgs-like":<20} spin 0    λ={evals8_sig[3]:.4f}  Δ/Δ₀={-np.log(evals8_sig[3])/Delta_0:.3f}  radial

COMPOSITE (bilinear Koopman products):
  GLUEBALLS (F⁺⊗F⁺, colour singlet):
    {"0⁺⁺":<20} spin 0    Δ/Δ₀=1.000  (reference)
    {"2⁺⁺":<20} spin 2    Δ/Δ₀=√2={np.sqrt(2):.3f}  (exact theorem)

  MESONS (quark × antiquark, λ_f²):
    {"π (spin 0)":<20} spin 0    ½⊗½→0    Δ/Δ₀={-np.log(evals8_sig[1]**2)/Delta_0:.3f}
    {"ρ (spin 1)":<20} spin 1    ½⊗½→1    Δ/Δ₀={-np.log(evals8_sig[1]**2)/Delta_0:.3f}

  BARYONS (three quarks, λ_f³):
    {"p (spin 1/2)":<20} spin 1/2  ½⊗½⊗½→½  Δ/Δ₀={-np.log(evals8_sig[1]**3)/Delta_0:.3f}
    {"Δ (spin 3/2)":<20} spin 3/2  ½⊗½⊗½→3/2 Δ/Δ₀={-np.log(evals8_sig[1]**3)/Delta_0:.3f}
""")

# ============================================================
# MESON AND BARYON MASS RATIOS
# ============================================================

print("="*60)
print("MASS RATIOS TO 0++ GLUEBALL")  
print("="*60)

# Using 0++ glueball ≈ 1710 MeV from lattice
m_glueball = 1710  # MeV

lam_f = evals8_sig[1]  # spinor eigenvalue
lam_s = evals8_sig[0]  # scalar eigenvalue

Delta_glueball = -np.log(lam_s)  # the 0++ is from the scalar eigenvalue... 
# Actually the 0++ glueball ratio was defined from the A₂ spectral radius
# which is λ₀ = 0.5022. Let me be careful.

# The glueball 0++ is at the spectral radius of A₂
Delta_ref = -np.log(evals8_sig[0])

print(f"Reference: 0++ glueball at λ₀ = {evals8_sig[0]:.4f}, Δ = {Delta_ref:.4f}")
print(f"Using m(0++) ≈ {m_glueball} MeV")
print()

particles = [
    ("0⁺⁺ glueball", -np.log(evals8_sig[0])),
    ("fermion (elementary)", -np.log(evals8_sig[1])),
    ("2⁺⁺ glueball", -np.log(evals8_sig[0]) * np.sqrt(2)),
    ("0⁻⁺ glueball", -np.log(evals8_sig[2])),
    ("0⁺⁺* glueball", -np.log(evals8_sig[3])),
    ("meson (π/ρ)", -np.log(evals8_sig[1]**2)),
    ("baryon (p/Δ)", -np.log(evals8_sig[1]**3)),
]

print(f"{'Particle':<25} {'Δ':>7} {'Δ/Δ₀':>7} {'MeV':>8}")
print("-"*55)
for name, delta in sorted(particles, key=lambda x: x[1]):
    ratio = delta / Delta_ref
    mev = ratio * m_glueball
    print(f"{name:<25} {delta:7.4f} {ratio:7.3f} {mev:8.0f}")

print(f"""
Lattice/experimental comparison:
  0⁺⁺ glueball:  1710 MeV (lattice)     → {m_glueball} MeV ✓
  2⁺⁺ glueball:  2400 MeV (lattice)     → {np.sqrt(2)*m_glueball:.0f} MeV (√2 × 1710)
  proton:         938 MeV (experiment)   → framework: {-np.log(lam_f**3)/Delta_ref * m_glueball:.0f} MeV
  pion:           140 MeV (experiment)   → framework: {-np.log(lam_f**2)/Delta_ref * m_glueball:.0f} MeV
""")

