"""
Analytical verification of gravity claims.

For each numerical observation, determine:
  (a) Can it be proven analytically?
  (b) If so, what is the proof?
  (c) If not, is the numerical evidence sufficient?

Clay standard: every claim must be a theorem with proof,
or explicitly labelled as a conjecture with evidence.
"""
import numpy as np
from sympy import *

# ============================================================
# CLAIM 1: Raw DS is holomorphic
# ============================================================
print("=" * 60)
print("CLAIM 1: Raw DS (no normalization) is holomorphic")
print("=" * 60)

# DS combination without normalization:
#   s_new = s*se + s*θ_e + θ*se
#   θ_new = θ*θ_e
#   K = Σ_{i≠j} s_i * se_j
#   m_out = m_raw / (1 - K)
#
# Each component is a rational function of (s1, s2, s3, θ):
#   numerator and denominator are polynomials in complex variables.
# A rational function of z (not z̄) is holomorphic wherever defined.
#
# PROOF: The DS pre-normalisation products s_new_i = s_i*se_i + s_i*θ_e + θ*se_i
# and θ_new = θ*θ_e are polynomial (hence holomorphic) in the components of m.
# K = Σ_{i≠j} s_i * se_j is polynomial in s and se.
# Division by (1-K) is holomorphic where K ≠ 1.
# Therefore the raw DS map is holomorphic on {K ≠ 1}. QED.
#
# This is NOT an empirical claim. It follows directly from the definition.

print("PROVEN: Raw DS is a rational function of complex variables,")
print("        hence holomorphic on {K ≠ 1}.")
print()

# ============================================================
# CLAIM 2: L1 normalization is non-holomorphic
# ============================================================
print("=" * 60)
print("CLAIM 2: L1 normalization is non-holomorphic")
print("=" * 60)

# L1 normalization: m → m / Σ|m_i|
# |m_i| = √(m_i * m̄_i)
# This depends on m̄_i, hence is non-holomorphic.
#
# PROOF: Let f(m) = m_j / Σ_k |m_k| for some component j.
# ∂f/∂m̄_k = m_j * ∂/∂m̄_k [1/Σ|m_l|]
#           = m_j * (-1/(Σ|m_l|)²) * ∂|m_k|/∂m̄_k
#           = m_j * (-1/(Σ|m_l|)²) * m_k/(2|m_k|)
# This is generically nonzero. QED.

# Verify symbolically
m1, m2, m3, m4 = symbols('m1 m2 m3 m4')
m1b, m2b, m3b, m4b = symbols('m1b m2b m3b m4b')  # conjugates

# |m_i| = sqrt(m_i * m_ib)
abs_m1 = sqrt(m1 * m1b)
abs_m2 = sqrt(m2 * m2b)
abs_m3 = sqrt(m3 * m3b)
abs_m4 = sqrt(m4 * m4b)

L1 = abs_m1 + abs_m2 + abs_m3 + abs_m4
f1_norm = m1 / L1

# dbar w.r.t. m1b
df1_dm1b = diff(f1_norm, m1b)
print(f"∂(m1/L1)/∂m̄1 = {simplify(df1_dm1b)}")
print(f"Nonzero: {df1_dm1b != 0}")
print()
print("PROVEN: L1 normalization introduces ∂̄ ≠ 0.")
print()

# ============================================================
# CLAIM 3: Born floor is non-holomorphic
# ============================================================
print("=" * 60)
print("CLAIM 3: Born floor is non-holomorphic")
print("=" * 60)
print("ALREADY PROVEN: Theorem 19.3 in the paper.")
print("Equation (25): ∂Φ_θ/∂θ̄ = -cθ²/(2|θ|³) ≠ 0")
print()

# ============================================================
# CLAIM 4: The decomposition of ∂̄Φ into spin-0, spin-1, spin-2
# ============================================================
print("=" * 60)
print("CLAIM 4: Spin decomposition of ∂̄Φ")
print("=" * 60)

# ∂̄Φ is a 4×4 complex matrix (the anti-holomorphic Jacobian).
# Under SO(4) (the tangent space symmetry), a 4×4 matrix decomposes as:
#   4 ⊗ 4 = 1 ⊕ 6 ⊕ 9
# i.e., trace (1, scalar) + antisymmetric (6, = 3+3 under SO(3)) + symmetric traceless (9)
#
# The trace part is the conformal factor (spin-0)
# The antisymmetric part contains spin-1
# The symmetric traceless part contains spin-2 (graviton)
#
# This decomposition is representation theory, not a claim.
# The CLAIM is about the relative magnitudes.
#
# Question: can we analytically determine the spin-2 fraction?

# The ∂̄Φ matrix at the fixed point has the structure from eq (27):
# Only the θ-row and singleton rows have nonzero entries.
# The θ-row: (-cθ²/(2|θ|³), θs₁/(52|θ|c), θs₂/(52|θ|c), θs₃/(52|θ|c))
# (This is the structure for a specific parametrisation)

# At the K*=7/30 fixed point:
# m* = (0.7869, 0.0293, 0.0293, 0.1545)
# So s_dom >> s_weak, and θ is intermediate.

# The full ∂̄Φ includes L1 normalization AND floor.
# The L1 part is generic. The floor part has the structure (27).
# The spin-2 dominance depends on the specific mass function values.

# HONEST ASSESSMENT: The spin-2 fraction being ~62% is a numerical
# observation at K*=7/30 with specific evidence. It is NOT derivable
# from representation theory alone. Different mass functions give
# different fractions.

# HOWEVER: the floor-only contribution being ~77% spin-2 IS more
# structural. Let me check if this follows from the matrix structure.

print("The decomposition 4⊗4 = 1⊕6⊕9 is exact (representation theory).")
print("The spin-2 fraction at K*=7/30 (~62%) is NUMERICAL, not proven.")
print()

# Can we bound the spin-2 fraction from below?
# The floor matrix (eq 27) has the structure:
# Row 4 (θ): [θs₁/..., θs₂/..., θs₃/..., -cθ²/...]
# Other rows: contributions from si rescaling

# At the symmetric point (s1=s2=s3=s, θ):
# The θ-row cross terms are all equal: θs/(52|θ|c) for each
# The diagonal θ-term is: -cθ²/(2|θ|³)
# By symmetry, the matrix has a specific structure we can analyse.

print("Checking: at the SYMMETRIC fixed point (s1=s2=s3):")
# Use the symmetric self-evidence FP: s=0.2994, θ=0.1017
s_sym = 0.2994
th_sym = 0.1017

# The floor enforcement Jacobian structure is determined by
# the specific form of the enforcement. Let me compute it
# symbolically for the symmetric case.
print(f"At symmetric point: s={s_sym}, θ={th_sym}")
print(f"Born = {th_sym**2 / (3*s_sym**2 + th_sym**2):.6f}")
print(f"Floor active: Born < 1/27 = {th_sym**2 / (3*s_sym**2 + th_sym**2) < 1/27}")
print()

# ============================================================
# CLAIM 5: Ricci/Weyl ratio is fixed
# ============================================================
print("=" * 60)
print("CLAIM 5: Ricci/Weyl ratio at equilibrium")
print("=" * 60)

# The numerical observation: Ricci/Weyl ≈ 0.6748 across perturbations.
# Is this analytically derivable?
#
# The emergent metric g = diag(θ², -(s₁+θ)², -(s₂+θ)², -(s₃+θ)²)
# is fully determined by (s₁, s₂, s₃, θ).
# The Ricci and Weyl tensors are computable analytically from this metric.
# At the fixed point m*, all values are known.
#
# So the Ricci/Weyl ratio IS analytically computable from the metric.
# But "analytically" here means "algebra on known numbers" — not a
# clean formula.
#
# The question is: does it equal 2/3 exactly?

# Compute with higher precision
from itertools import product as iprod

def metric_components(th, s1, s2, s3):
    return np.array([th**2, -(s1+th)**2, -(s2+th)**2, -(s3+th)**2])

def full_curvature(th, s1, s2, s3, eps=1e-6):
    """Compute Ricci and Weyl tensors analytically (via finite diff)."""
    coords = np.array([th, s1, s2, s3])
    n = 4

    def g_at(c):
        return np.diag(metric_components(c[0], c[1], c[2], c[3]))

    g = g_at(coords)
    g_inv = np.diag(1.0 / np.diag(g))

    # Metric derivatives
    dg = np.zeros((n, n, n))
    for a in range(n):
        cp = coords.copy(); cp[a] += eps
        cm = coords.copy(); cm[a] -= eps
        dg[:, :, a] = (g_at(cp) - g_at(cm)) / (2 * eps)

    # Christoffel
    Gamma = np.zeros((n, n, n))
    for lam in range(n):
        for mu in range(n):
            for nu in range(n):
                for sig in range(n):
                    Gamma[lam, mu, nu] += g_inv[lam, sig] * (
                        dg[sig, nu, mu] + dg[sig, mu, nu] - dg[mu, nu, sig])
                Gamma[lam, mu, nu] *= 0.5

    # Christoffel derivatives
    def gamma_at(c):
        g_c = np.diag(metric_components(c[0], c[1], c[2], c[3]))
        gi_c = np.diag(1.0 / np.diag(g_c))
        dg_c = np.zeros((n, n, n))
        for a in range(n):
            cp2 = c.copy(); cp2[a] += eps*0.1
            cm2 = c.copy(); cm2[a] -= eps*0.1
            dg_c[:, :, a] = (np.diag(metric_components(*cp2)) - np.diag(metric_components(*cm2))) / (2*eps*0.1)
        G = np.zeros((n, n, n))
        for l in range(n):
            for m in range(n):
                for nn in range(n):
                    for ss in range(n):
                        G[l, m, nn] += gi_c[l, ss] * (dg_c[ss, nn, m] + dg_c[ss, m, nn] - dg_c[m, nn, ss])
                    G[l, m, nn] *= 0.5
        return G

    dGamma = np.zeros((n, n, n, n))
    for mu in range(n):
        cp = coords.copy(); cp[mu] += eps
        cm = coords.copy(); cm[mu] -= eps
        dGamma[:, :, :, mu] = (gamma_at(cp) - gamma_at(cm)) / (2*eps)

    # Riemann
    R = np.zeros((n, n, n, n))
    for rho in range(n):
        for sig in range(n):
            for mu in range(n):
                for nu in range(n):
                    R[rho, sig, mu, nu] = dGamma[rho, nu, sig, mu] - dGamma[rho, mu, sig, nu]
                    for lam in range(n):
                        R[rho, sig, mu, nu] += Gamma[rho, mu, lam]*Gamma[lam, nu, sig] - Gamma[rho, nu, lam]*Gamma[lam, mu, sig]

    # Ricci
    Ric = np.zeros((n, n))
    for mu in range(n):
        for nu in range(n):
            for lam in range(n):
                Ric[mu, nu] += R[lam, mu, lam, nu]

    # Scalar
    R_scalar = sum(g_inv[i, i] * Ric[i, i] for i in range(n))

    return R, Ric, R_scalar, g, g_inv

# At K*=7/30 fixed point
th_fp = 0.1545370339
s1_fp = 0.7868984462
s2_fp = 0.0292822600

R, Ric, R_s, g, g_inv = full_curvature(th_fp, s1_fp, s2_fp, s2_fp)

# Ricci squared
Ric_sq = sum(Ric[i,i]**2 * abs(g_inv[i,i]) for i in range(4))

# Weyl: compute from Kretschner - Ricci contribution
# K = |Riem|² = |Weyl|² + 2|Ric|² - R²/3  (in 4D)
# Actually: K = C² + 2(R_μν R^μν - R²/3)
# So C² = K - 2R_μν R^μν + 2R²/3

# Kretschner
K_ret = 0.0
for a, b, c, d in iprod(range(4), repeat=4):
    R_lower = np.diag(g)[a] * R[a, b, c, d]
    R_upper = R_lower * abs(g_inv[a,a]) * abs(g_inv[b,b]) * abs(g_inv[c,c]) * abs(g_inv[d,d])
    K_ret += R_lower * R_upper

# Ricci squared with both indices up
Ric_sq_2 = sum(Ric[i,i]**2 * g_inv[i,i]**2 for i in range(4))

C_sq = K_ret - 2*Ric_sq_2 + 2*R_s**2/3

print(f"Kretschner = {K_ret:.6f}")
print(f"|Ric|² = {Ric_sq_2:.6f}")
print(f"R = {R_s:.6f}")
print(f"|Weyl|² = {C_sq:.6f}")
print()

ratio = np.sqrt(abs(Ric_sq_2)) / np.sqrt(abs(C_sq)) if abs(C_sq) > 1e-10 else float('inf')
print(f"√(|Ric|²/|Weyl|²) = {ratio:.6f}")
print(f"2/3 = {2/3:.6f}")
print(f"Match: {abs(ratio - 2/3) < 0.01}")
print()

# Now check at the symmetric fixed point
print("At symmetric self-evidence FP (s=0.2994, θ=0.1017):")
R2, Ric2, R_s2, g2, g_inv2 = full_curvature(0.1017, 0.2994, 0.2994, 0.2994)

Ric_sq2 = sum(Ric2[i,i]**2 * g_inv2[i,i]**2 for i in range(4))
K_ret2 = 0.0
for a, b, c, d in iprod(range(4), repeat=4):
    R_lower = np.diag(g2)[a] * R2[a, b, c, d]
    R_upper = R_lower * abs(g_inv2[a,a]) * abs(g_inv2[b,b]) * abs(g_inv2[c,c]) * abs(g_inv2[d,d])
    K_ret2 += R_lower * R_upper

C_sq2 = K_ret2 - 2*Ric_sq2 + 2*R_s2**2/3
ratio2 = np.sqrt(abs(Ric_sq2)) / np.sqrt(abs(C_sq2)) if abs(C_sq2) > 1e-10 else float('inf')
print(f"√(|Ric|²/|Weyl|²) = {ratio2:.6f}")
print(f"2/3 = {2/3:.6f}")
print()

# Check at generic points
print("At various other mass functions:")
test_points = [
    ("Near ignorance", 0.8, 0.05, 0.05, 0.1),
    ("Democratic", 0.25, 0.25, 0.25, 0.25),
    ("Asymmetric", 0.15, 0.6, 0.1, 0.15),
    ("Extreme", 0.02, 0.9, 0.04, 0.04),
]

for name, th, s1, s2, s3 in test_points:
    try:
        Rt, Rict, R_st, gt, g_invt = full_curvature(th, s1, s2, s3)
        Ric_sqt = sum(Rict[i,i]**2 * g_invt[i,i]**2 for i in range(4))
        K_rett = 0.0
        for a, b, c, d in iprod(range(4), repeat=4):
            R_lower = np.diag(gt)[a] * Rt[a, b, c, d]
            R_upper = R_lower * abs(g_invt[a,a]) * abs(g_invt[b,b]) * abs(g_invt[c,c]) * abs(g_invt[d,d])
            K_rett += R_lower * R_upper
        C_sqt = K_rett - 2*Ric_sqt + 2*R_st**2/3
        rat = np.sqrt(abs(Ric_sqt)) / np.sqrt(abs(C_sqt)) if abs(C_sqt) > 1e-10 else float('inf')
        print(f"  {name:20s}: Ric/Weyl = {rat:.6f}")
    except:
        print(f"  {name:20s}: ERROR")

print()

# ============================================================
# CLAIM 6: Conformal invariance is broken
# ============================================================
print("=" * 60)
print("CLAIM 6: Conformal invariance breaking")
print("=" * 60)

# The Born floor introduces |θ|² = θθ̄, which is non-holomorphic.
# This is already proven (Theorem 19.3).
# But does non-holomorphicity imply conformal breaking?
#
# In Mason's framework:
# - Integrable J (holomorphic): self-dual sector only, conformal gravity
# - Non-integrable J (non-holomorphic): full theory, Yang-Mills + gravity
#
# Mason's paper treats both YM and conformal gravity.
# The jump from conformal gravity to Einstein gravity requires
# an ADDITIONAL structure: a preferred metric in each conformal class.
#
# Does the Born floor provide this?
# YES: the floor fixes Born(θ) = 1/27 at equilibrium.
# This condition picks a specific metric, not a conformal class.
# Born(Ωθ, Ωs) = Born(θ, s) (Born is scale-invariant) BUT
# the metric g depends on the actual values of θ, s (not just ratios).
# On the L1=1 simplex, Ω=1. The L1 constraint fixes the scale.
#
# So: L1 + Born floor = fixed scale = broken conformal invariance.
#
# This IS provable: it's a direct consequence of L1=1 being a
# non-projective constraint (it fixes the representative in each
# equivalence class of CP³).

print("PROVABLE: L1=1 fixes the scale on CP³.")
print("Born floor = 1/27 further constrains within the simplex.")
print("Together they pick a unique metric, breaking conformal invariance.")
print()
print("HOWEVER: this does NOT automatically give Einstein equations.")
print("It gives a preferred metric in a conformal class.")
print("Whether the dynamics (DS combination) produce Einstein equations")
print("or something else is a SEPARATE question.")
print()

# ============================================================
# SUMMARY
# ============================================================
print("=" * 60)
print("RIGOUR AUDIT")
print("=" * 60)

claims = [
    ("Raw DS is holomorphic", "PROVEN", "rational function of complex variables"),
    ("L1 normalization is non-holomorphic", "PROVEN", "∂/∂m̄ ≠ 0, symbolic verification"),
    ("Born floor is non-holomorphic", "PROVEN", "Theorem 19.3 in paper"),
    ("∂̄Φ is ~62% spin-2 (rank-2 coupled)", "NUMERICAL", "representation theory decomposition is exact, fraction is empirical"),
    ("Floor-only ∂̄Φ is ~77% spin-2", "NUMERICAL", "structural observation, not derived"),
    ("Ricci/Weyl ratio ≈ 0.6748", "COMPUTABLE", "follows from metric at FP, but not a clean formula"),
    ("Ricci/Weyl ratio = 2/3 exactly", "UNPROVEN", "numerical match to 1%, needs analytical verification"),
    ("Ricci/Weyl is FIXED (doesn't vary)", "NEEDS WORK", "confirmed numerically near FP, but only for the emergent metric"),
    ("Conformal invariance broken", "PROVABLE", "L1 + floor fixes scale, standard argument"),
    ("Conformal → Einstein specifically", "UNPROVEN", "having a preferred metric ≠ Einstein equations"),
    ("w ≈ 1 (stiff matter)", "NUMERICAL", "from emergent metric curvature, specific to FP"),
    ("Penrose transform carries floor → gravity", "NOT DONE", "formal argument only, computation needed"),
]

for claim, status, note in claims:
    marker = {"PROVEN": "✓", "PROVABLE": "~", "NUMERICAL": "?",
              "COMPUTABLE": "~", "UNPROVEN": "✗", "NEEDS WORK": "?", "NOT DONE": "✗"}
    print(f"  {marker.get(status, '?')} [{status:12s}] {claim}")
    print(f"    {note}")
print()

print("BOTTOM LINE:")
print("  3 claims proven")
print("  2 claims provable with work")
print("  4 claims numerical (need analytical derivation or honest labelling)")
print("  3 claims unproven/not done (genuine gaps)")
