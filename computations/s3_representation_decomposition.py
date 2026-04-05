"""
S₃ representation theory in the composition operator.

The eigenvalue phases split by sign representation:
  even permutations → φ₁ ≈ +0.42π
  odd permutations  → φ₁ ≈ -0.41π

This suggests the composition operator T_B(A) = compose(A, B) decomposes
under S₃ into irreducible representations. S₃ has three irreps:
  trivial (dim 1): ρ(σ) = 1 for all σ
  sign (dim 1):    ρ(σ) = sgn(σ)
  standard (dim 2): the 2D representation

Questions:
1. Does the 16-dim composition operator decompose cleanly into S₃ irreps?
2. Do the eigenvalue phase families correspond to trivial vs sign?
3. What is the standard (dim-2) component — does it carry the conjugacy info?
4. Is there a character formula connecting decay rates to conjugacy classes?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import numpy as np
from solver.algebra import H, MASS_DIM, DELTA, schmidt_number
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)


# S₃ elements as correlation matrices
S3 = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

SIGNS = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}


def build_composition_operator(B_joint):
    """Build the 16×16 linear map T_B: A → compose(A, B).

    Since compose is bilinear in components, T_B is a linear operator
    on the 16-dim space of flattened [4,4] joint masses.
    """
    dim = MASS_DIM * MASS_DIM  # 16
    T = torch.zeros(dim, dim, dtype=torch.cfloat)

    for col in range(dim):
        # Unit vector in the col-th direction
        A = torch.zeros(dim, dtype=torch.cfloat)
        A[col] = 1.0 + 0j
        A = A.reshape(MASS_DIM, MASS_DIM)

        result = compose(A, B_joint)
        T[:, col] = result.reshape(dim)

    return T


# ═══════════════════════════════════════════════════════════════
#  Build composition operators for all S₃ elements
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("S₃ REPRESENTATION DECOMPOSITION OF COMPOSITION OPERATORS")
print("=" * 80)

operators = {}
joints = {}

for name, corr in S3.items():
    ent = Entangler(corr, seed=42).build()
    joints[name] = ent.joint
    T = build_composition_operator(ent.joint)
    operators[name] = T

    # Eigenvalues of T
    eigenvalues = torch.linalg.eigvals(T)
    mags = eigenvalues.abs()
    idx = mags.argsort(descending=True)
    eigenvalues = eigenvalues[idx]

    print(f"\n  {name} (sign={SIGNS[name]:+d}):")
    print(f"    Top 6 eigenvalues (magnitude, phase/π):")
    for i in range(min(6, len(eigenvalues))):
        ev = eigenvalues[i]
        mag = ev.abs().item()
        phase = math.atan2(ev.imag.item(), ev.real.item()) / math.pi
        print(f"      λ_{i}: |λ| = {mag:.6f}, φ/π = {phase:+.4f}")


# ═══════════════════════════════════════════════════════════════
#  S₃ irrep decomposition using characters
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CHARACTER ANALYSIS")
print("="*80)

# Character table of S₃:
#           e    (12)    (123)
# trivial:  1     1       1
# sign:     1    -1       1
# standard: 2     0      -1

# For a representation ρ, the multiplicity of irrep χ is:
# m_χ = (1/|G|) Σ_g χ(g)* tr(ρ(g))

# Get traces of composition operators
traces = {}
for name in S3:
    tr = operators[name].diagonal().sum()
    traces[name] = tr
    print(f"  tr(T_{name}) = {tr.real.item():.6f} + {tr.imag.item():.6f}i")

# Characters by conjugacy class
# Conjugacy classes: {e}, {(01),(02),(12)}, {(012),(021)}
tr_e = traces["e"]
tr_trans = (traces["(01)"] + traces["(02)"] + traces["(12)"]) / 3
tr_3cyc = (traces["(012)"] + traces["(021)"]) / 2

print(f"\n  Average trace by conjugacy class:")
print(f"    e:             {tr_e.real.item():.6f}")
print(f"    transpositions: {tr_trans.real.item():.6f}")
print(f"    3-cycles:       {tr_3cyc.real.item():.6f}")

# Multiplicity formula: m_χ = (1/6)[1·χ(e)·tr(e) + 3·χ(trans)·tr(trans) + 2·χ(3cyc)·tr(3cyc)]
# Using real parts (the representation should be real)

tr_e_r = tr_e.real.item()
tr_t_r = tr_trans.real.item()
tr_c_r = tr_3cyc.real.item()

m_trivial  = (1/6) * (1*1*tr_e_r + 3*1*tr_t_r + 2*1*tr_c_r)
m_sign     = (1/6) * (1*1*tr_e_r + 3*(-1)*tr_t_r + 2*1*tr_c_r)
m_standard = (1/6) * (1*2*tr_e_r + 3*0*tr_t_r + 2*(-1)*tr_c_r)

print(f"\n  Irrep multiplicities (from character formula):")
print(f"    trivial:  {m_trivial:.3f}")
print(f"    sign:     {m_sign:.3f}")
print(f"    standard: {m_standard:.3f}")
print(f"    Total:    {m_trivial + m_sign + 2*m_standard:.3f} (should be 16)")


# ═══════════════════════════════════════════════════════════════
#  Projectors onto irreps
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("IRREP PROJECTORS")
print("="*80)

# P_χ = (dim_χ / |G|) Σ_g χ(g)* ρ(g)
dim_sq = MASS_DIM * MASS_DIM

# Need to build the GROUP representation: how S₃ acts on [4,4] joint masses
# A permutation σ acts on the joint mass M by: (σ·M)_{ij} = M_{σ⁻¹(i), σ⁻¹(j)}
# (permuting rows and columns among h0,h1,h2 — θ index stays fixed)

def permutation_matrix(perm):
    """Build 4×4 permutation matrix from a permutation of {0,1,2}.
    Index 3 (θ) is fixed."""
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0  # θ fixed
    return P

# S₃ as permutations of {0,1,2}
S3_PERMS = {
    "e":     [0, 1, 2],
    "(01)":  [1, 0, 2],
    "(02)":  [2, 1, 0],
    "(12)":  [0, 2, 1],
    "(012)": [1, 2, 0],
    "(021)": [2, 0, 1],
}

# Build 16×16 representation matrices: ρ(σ) acts on vec(M) as (P_σ ⊗ P_σ) vec(M)
def group_rep(perm):
    """16×16 matrix representing σ acting on flattened [4,4] mass."""
    P = permutation_matrix(perm)
    return torch.kron(P, P)

reps = {}
for name, perm in S3_PERMS.items():
    reps[name] = group_rep(perm)

# Verify: group multiplication
print("  Verification: (012) × (01) should give a transposition...")
product = reps["(012)"] @ reps["(01)"]
# Find which element this matches
for name, rep in reps.items():
    if torch.allclose(product, rep, atol=1e-6):
        print(f"    (012)·(01) = {name}  ✓")
        break

# Character table verification
for name in S3:
    tr = reps[name].diagonal().sum().real.item()
    print(f"  tr(ρ({name})) = {tr:.0f}")

# Projectors
chi_trivial  = {"e": 1, "(01)": 1, "(02)": 1, "(12)": 1, "(012)": 1, "(021)": 1}
chi_sign     = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}
chi_standard = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}

def build_projector(chi, dim_chi):
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for name in S3:
        P += chi[name] * reps[name]
    P *= dim_chi / 6.0
    return P

P_triv = build_projector(chi_trivial, 1)
P_sign = build_projector(chi_sign, 1)
P_std  = build_projector(chi_standard, 2)

# Verify projectors
print(f"\n  Projector ranks:")
print(f"    trivial:  rank = {torch.linalg.matrix_rank(P_triv.real.float()).item()}")
print(f"    sign:     rank = {torch.linalg.matrix_rank(P_sign.real.float()).item()}")
print(f"    standard: rank = {torch.linalg.matrix_rank(P_std.real.float()).item()}")
print(f"    Sum of ranks: {torch.linalg.matrix_rank(P_triv.real.float()).item() + torch.linalg.matrix_rank(P_sign.real.float()).item() + torch.linalg.matrix_rank(P_std.real.float()).item()}")

# Verify: P_triv + P_sign + P_std = I
total = P_triv + P_sign + P_std
identity_check = torch.eye(dim_sq, dtype=torch.cfloat)
print(f"    P_triv + P_sign + P_std = I? max error = {(total - identity_check).abs().max().item():.2e}")


# ═══════════════════════════════════════════════════════════════
#  Project crystals onto irrep subspaces
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CRYSTAL PROJECTION ONTO S₃ IRREPS")
print("="*80)

for name in S3:
    v = joints[name].reshape(dim_sq)

    v_triv = P_triv @ v
    v_sign = P_sign @ v
    v_std  = P_std @ v

    norm_total = v.abs().pow(2).sum().item()
    frac_triv = v_triv.abs().pow(2).sum().item() / norm_total
    frac_sign = v_sign.abs().pow(2).sum().item() / norm_total
    frac_std  = v_std.abs().pow(2).sum().item() / norm_total

    print(f"\n  {name} (sign={SIGNS[name]:+d}):")
    print(f"    trivial:  {frac_triv:.4f} ({frac_triv*100:.1f}%)")
    print(f"    sign:     {frac_sign:.4f} ({frac_sign*100:.1f}%)")
    print(f"    standard: {frac_std:.4f} ({frac_std*100:.1f}%)")
    print(f"    total:    {frac_triv + frac_sign + frac_std:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Composition operator restricted to each irrep subspace
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMPOSITION OPERATOR IN EACH IRREP SUBSPACE")
print("="*80)

# For the identity crystal, look at T_e restricted to each subspace
T_e = operators["e"]

# Eigenvalues of P_χ T_e P_χ restricted to im(P_χ)
for label, P_chi in [("trivial", P_triv), ("sign", P_sign), ("standard", P_std)]:
    restricted = P_chi @ T_e @ P_chi
    eigenvalues = torch.linalg.eigvals(restricted)
    # Filter out near-zero eigenvalues (from the kernel of the projector)
    nonzero = eigenvalues[eigenvalues.abs() > 1e-8]
    mags = nonzero.abs()
    idx = mags.argsort(descending=True)
    nonzero = nonzero[idx]

    print(f"\n  T_e restricted to {label} ({len(nonzero)} non-zero eigenvalues):")
    for i in range(min(6, len(nonzero))):
        ev = nonzero[i]
        mag = ev.abs().item()
        phase = math.atan2(ev.imag.item(), ev.real.item()) / math.pi
        print(f"    λ_{i}: |λ| = {mag:.6f}, φ/π = {phase:+.4f}")


# ═══════════════════════════════════════════════════════════════
#  The key test: does the sign component predict phase flipping?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SIGN COMPONENT AND PHASE STRUCTURE")
print("="*80)

# For each S₃ element, measure the sign component of its crystal
# and correlate with the eigenvalue phase under self-composition

from solver.crystals import compose as crystal_compose

for name in S3:
    joint = joints[name]
    v = joint.reshape(dim_sq)

    # Sign component
    v_sign = P_sign @ v
    sign_frac = v_sign.abs().pow(2).sum().item() / v.abs().pow(2).sum().item()

    # Phase of sign component relative to trivial
    v_triv = P_triv @ v
    triv_phase = torch.angle(v_triv[v_triv.abs().argmax()]).item() / math.pi
    sign_phase = torch.angle(v_sign[v_sign.abs().argmax()]).item() / math.pi
    relative_phase = sign_phase - triv_phase
    while relative_phase > 1: relative_phase -= 2
    while relative_phase < -1: relative_phase += 2

    # Self-composition eigenvalue phase (from previous analysis)
    current = joint.clone()
    composed = crystal_compose(current, joint)
    ev_orig = torch.linalg.eigvals(joint)
    ev_comp = torch.linalg.eigvals(composed)

    idx_orig = ev_orig.abs().argsort(descending=True)
    idx_comp = ev_comp.abs().argsort(descending=True)

    phase_orig = math.atan2(ev_orig[idx_orig[1]].imag.item(), ev_orig[idx_orig[1]].real.item()) / math.pi
    phase_comp = math.atan2(ev_comp[idx_comp[1]].imag.item(), ev_comp[idx_comp[1]].real.item()) / math.pi
    phase_step = phase_comp - phase_orig
    while phase_step > 1: phase_step -= 2
    while phase_step < -1: phase_step += 2

    print(f"  {name:>6s}: sign_frac={sign_frac:.4f}, rel_phase={relative_phase:+.3f}π, "
          f"comp_phase_step={phase_step:+.4f}π, algebraic_sign={SIGNS[name]:+d}")


# ═══════════════════════════════════════════════════════════════
#  Does the decomposition explain the rate differences?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("RATE-IRREP CONNECTION")
print("="*80)

print("""
  If the composition operator respects S₃, then the eigenvalues of
  self-composition T_B(B) = compose(B, B) should decompose as:
    - trivial sector:  eigenvalue magnitude determined by |v_triv|²
    - sign sector:     eigenvalue magnitude ∝ sign(σ) × |v_sign|²
    - standard sector: eigenvalue magnitude determined by standard rep trace

  The RATE difference between conjugacy classes would then come from
  different projections onto these sectors:
    - Identity:       symmetric, large trivial component
    - Transpositions: antisymmetric in one pair, large sign component
    - 3-cycles:       rotationally symmetric, large standard component
""")

# Compute the norm in each sector for each conjugacy class
for cls_name, members in [("identity", ["e"]), ("transpositions", ["(01)", "(02)", "(12)"]),
                           ("3-cycles", ["(012)", "(021)"])]:
    triv_fracs = []
    sign_fracs = []
    std_fracs = []

    for name in members:
        v = joints[name].reshape(dim_sq)
        norm = v.abs().pow(2).sum().item()
        triv_fracs.append((P_triv @ v).abs().pow(2).sum().item() / norm)
        sign_fracs.append((P_sign @ v).abs().pow(2).sum().item() / norm)
        std_fracs.append((P_std @ v).abs().pow(2).sum().item() / norm)

    avg_triv = sum(triv_fracs) / len(triv_fracs)
    avg_sign = sum(sign_fracs) / len(sign_fracs)
    avg_std = sum(std_fracs) / len(std_fracs)

    print(f"  {cls_name:>15s}: trivial={avg_triv:.4f}, sign={avg_sign:.4f}, standard={avg_std:.4f}")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
