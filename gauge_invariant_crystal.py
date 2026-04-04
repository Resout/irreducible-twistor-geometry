"""
The gauge-invariant crystal: S₃-symmetrization.

Physical observables must be invariant under relabeling of hypotheses.
In S₃ terms, this means they live in the trivial representation.

The S₃-symmetrized crystal is: M_sym = (1/6) Σ_{σ ∈ S₃} σ · M

Questions:
1. Does symmetrization change the Schmidt number?
2. What survives symmetrization? (Only gauge-invariant information)
3. Does the symmetrized crystal compose differently?
4. Is the symmetrized composition gap different from 2Δ?
5. What is the fixed-point manifold of symmetrized crystals?
6. Does symmetrization relate to the paper's construction?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, schmidt_number, born_probabilities
from solver.crystals import Entangler, compose, classify_relationship, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)

dim_sq = MASS_DIM * MASS_DIM

# S₃ permutations
S3_PERMS = [
    [0, 1, 2],  # e
    [1, 0, 2],  # (01)
    [2, 1, 0],  # (02)
    [0, 2, 1],  # (12)
    [1, 2, 0],  # (012)
    [2, 0, 1],  # (021)
]


def apply_perm(joint, perm):
    """Apply permutation σ to joint mass: (σ·M)_{ij} = M_{σ⁻¹(i), σ⁻¹(j)}."""
    # Build inverse permutation (for h-indices; θ index fixed)
    inv_perm = [0] * H
    for i in range(H):
        inv_perm[perm[i]] = i

    result = torch.zeros_like(joint)
    for i in range(MASS_DIM):
        for j in range(MASS_DIM):
            # Map i, j through inverse permutation (θ = index H stays fixed)
            ii = inv_perm[i] if i < H else H
            jj = inv_perm[j] if j < H else H
            result[i, j] = joint[ii, jj]
    return result


def symmetrize(joint):
    """Project onto S₃-trivial sector: (1/|G|) Σ_σ σ·M."""
    result = torch.zeros_like(joint)
    for perm in S3_PERMS:
        result += apply_perm(joint, perm)
    return result / 6.0


# ═══════════════════════════════════════════════════════════════
#  Symmetrize each crystal type
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("GAUGE-INVARIANT (S₃-SYMMETRIZED) CRYSTALS")
print("=" * 80)

for name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Entangler(corr, seed=42).build()
    joint = ent.joint
    joint_sym = symmetrize(joint)

    sc_orig = schmidt_number(joint)
    sc_sym = schmidt_number(joint_sym)

    rel_orig, conf_orig = classify_relationship(joint)
    rel_sym, conf_sym = classify_relationship(joint_sym)

    # Born probabilities
    bp_orig = born_probabilities(joint.abs().pow(2).sum(dim=1))
    bp_sym = born_probabilities(joint_sym.abs().pow(2).sum(dim=1))

    print(f"\n  {name}:")
    print(f"    Schmidt: {sc_orig:.3f} → {sc_sym:.3f}")
    print(f"    Type:    {rel_orig} → {rel_sym}")
    print(f"    Born(θ): {bp_orig[H].item():.4f} → {bp_sym[H].item():.4f}")

    # Is the symmetrized crystal symmetric in h-indices?
    # M_sym should have M[i,j] = M[σi, σj] for all σ
    # This means M[0,0] = M[1,1] = M[2,2] and M[0,1] = M[0,2] = M[1,2] etc.
    diag_h = [joint_sym[i, i] for i in range(H)]
    off_diag_h = [joint_sym[0, 1], joint_sym[0, 2], joint_sym[1, 2]]
    theta_diag = [joint_sym[i, H] for i in range(H)]

    diag_spread = max(abs(diag_h[i] - diag_h[j]) for i in range(H) for j in range(H))
    off_spread = max(abs(off_diag_h[i] - off_diag_h[j]) for i in range(3) for j in range(3))
    theta_spread = max(abs(theta_diag[i] - theta_diag[j]) for i in range(H) for j in range(H))

    print(f"    h-diagonal spread: {diag_spread.item():.2e}")
    print(f"    h-off-diagonal spread: {off_spread.item():.2e}")
    print(f"    θ-column spread: {theta_spread.item():.2e}")


# ═══════════════════════════════════════════════════════════════
#  Structure of the symmetrized crystal
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("STRUCTURE OF SYMMETRIZED CRYSTALS")
print("="*80)

# The symmetrized [4,4] crystal has the form:
#   [a  b  b  c]
#   [b  a  b  c]
#   [b  b  a  c]
#   [c  c  c  d]
# (by S₃ symmetry: all h-h diagonals equal, all h-h off-diagonals equal,
#  all h-θ equal, θ-θ is its own thing)
# That's 4 free parameters: a, b, c, d

for name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Entangler(corr, seed=42).build()
    joint_sym = symmetrize(ent.joint)

    a = joint_sym[0, 0]
    b = joint_sym[0, 1]
    c = joint_sym[0, H]
    d = joint_sym[H, H]

    print(f"\n  {name}:")
    print(f"    a (h-diag)    = {a.real.item():+.6f} {a.imag.item():+.6f}i")
    print(f"    b (h-offdiag) = {b.real.item():+.6f} {b.imag.item():+.6f}i")
    print(f"    c (h-θ cross) = {c.real.item():+.6f} {c.imag.item():+.6f}i")
    print(f"    d (θ-θ)       = {d.real.item():+.6f} {d.imag.item():+.6f}i")

    # L1 constraint: 3a + 6b + 6c + d = 1 (from symmetry)
    l1 = 3*a + 6*b + 6*c + d
    print(f"    L1 = 3a + 6b + 6c + d = {l1.real.item():.6f}")

    # Eigenvalues of the symmetrized matrix
    ev, _ = torch.linalg.eig(joint_sym)
    idx = ev.abs().argsort(descending=True)
    ev = ev[idx]
    print(f"    Eigenvalues: ", end="")
    for i in range(MASS_DIM):
        print(f"|λ_{i}|={ev[i].abs().item():.4f} ", end="")
    print()

    # For a matrix of this form, eigenvalues are:
    # λ₁ = a - b (multiplicity 2, from the standard rep of S₃ on h's)
    # λ₂ = a + 2b + 3c ... (from the trivial rep involving θ)
    # Let's check:
    pred_std = a - b  # should have multiplicity 2
    print(f"    Predicted λ_std = a - b = {pred_std.abs().item():.4f}")


# ═══════════════════════════════════════════════════════════════
#  Eigenvalue structure of the symmetrized crystal (analytical)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("ANALYTICAL EIGENVALUES OF SYMMETRIZED CRYSTAL")
print("="*80)

# The symmetrized [4,4] matrix M_sym has the block form:
#   h-h block: a*I + b*(J-I) = (a-b)*I + b*J  where J = ones(3,3)
#   h-θ cross: c * ones(3,1) and c * ones(1,3)
#   θ-θ block: d
#
# Eigenvalues of the h-h block: a-b (mult 2, standard rep) and a+2b (mult 1, trivial)
# The full 4×4 matrix couples the trivial h-component (a+2b) with θ (d) via c.
#
# The 2×2 system for the trivial+θ sector:
#   [a+2b   3c] [v_h]   = λ [v_h]
#   [c      d ] [v_θ]       [v_θ]
#
# Wait, actually: the h-trivial component is (1/√3)(|0⟩+|1⟩+|2⟩),
# and its coupling to θ is 3c (because 3 h-θ entries sum).
# Let me be more careful.

# For the symmetrized matrix:
# M_sym = [[a, b, b, c],
#          [b, a, b, c],
#          [b, b, a, c],
#          [c, c, c, d]]
#
# Standard eigenvectors (orthogonal to (1,1,1,0)):
# v_std1 = (1,-1,0,0)/√2: eigenvalue = a - b
# v_std2 = (1,1,-2,0)/√6: eigenvalue = a - b (same)
#
# Trivial + θ sector spans {(1,1,1,0)/√3, (0,0,0,1)}:
# In this 2D subspace:
# M_triv = [[a+2b, √3·c], [√3·c, d]]
#
# Eigenvalues: ((a+2b+d) ± √((a+2b-d)² + 12|c|²)) / 2

for name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Entangler(corr, seed=42).build()
    joint_sym = symmetrize(ent.joint)

    a = joint_sym[0, 0]
    b = joint_sym[0, 1]
    c = joint_sym[0, H]
    d = joint_sym[H, H]

    # Standard eigenvalue
    lambda_std = a - b

    # Trivial+θ 2×2 system
    m11 = a + 2*b
    m12 = 3**0.5 * c  # √3 c
    m22 = d

    trace = m11 + m22
    det = m11 * m22 - m12 * m12  # using M^T M symmetry... actually M is complex
    # For complex eigenvalues: use the quadratic formula
    disc = (m11 - m22)**2 + 4 * m12**2
    disc_sqrt = disc**0.5  # complex sqrt
    lambda_triv_plus = (trace + disc_sqrt) / 2
    lambda_triv_minus = (trace - disc_sqrt) / 2

    print(f"\n  {name}:")
    print(f"    λ_std = a-b = {lambda_std.abs().item():.6f} (mult 2)")
    print(f"    λ₊ = {lambda_triv_plus.abs().item():.6f}")
    print(f"    λ₋ = {lambda_triv_minus.abs().item():.6f}")

    # Compare with numerical eigenvalues
    ev, _ = torch.linalg.eig(joint_sym)
    idx = ev.abs().argsort(descending=True)
    mags = ev.abs()[idx]
    print(f"    Numerical:  {mags[0].item():.6f}, {mags[1].item():.6f}, {mags[2].item():.6f}, {mags[3].item():.6f}")

    # The Schmidt of the symmetrized crystal depends only on these 3 eigenvalues
    # (λ_std with mult 2, λ₊, λ₋)


# ═══════════════════════════════════════════════════════════════
#  Composition of symmetrized crystals
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMPOSITION OF SYMMETRIZED CRYSTALS")
print("="*80)

# Key test: is compose(M_sym, M_sym) still symmetric?
# If yes, the symmetrized crystals form a closed subalgebra.

for name in list(RELATIONSHIP_SIGNATURES.keys())[:3]:
    corr = RELATIONSHIP_SIGNATURES[name]
    ent = Entangler(corr, seed=42).build()
    joint_sym = symmetrize(ent.joint)

    composed = compose(joint_sym, joint_sym)
    composed_sym = symmetrize(composed)

    diff = (composed - composed_sym).abs().max().item()
    sc_comp = schmidt_number(composed)
    sc_sym_comp = schmidt_number(composed_sym)

    print(f"\n  {name}:")
    print(f"    max|compose(M_sym, M_sym) - symmetrize(compose(M_sym, M_sym))| = {diff:.2e}")
    print(f"    compose(M_sym, M_sym) is {'S₃-symmetric' if diff < 1e-6 else 'NOT symmetric'}!")
    print(f"    Schmidt of composition: {sc_comp:.4f}")

    # Self-composition decay
    current = joint_sym.clone()
    for step in range(8):
        sc = schmidt_number(current)
        ev, _ = torch.linalg.eig(current)
        idx = ev.abs().argsort(descending=True)
        r1 = (ev[idx[1]].abs() / ev[idx[0]].abs()).item()
        print(f"    Step {step}: Schmidt={sc:.4f}, |λ₁/λ₀|={r1:.4f}")
        current = compose(current, joint_sym)


# ═══════════════════════════════════════════════════════════════
#  The symmetrized crystal has only 4 parameters: a, b, c, d
#  Under composition, this becomes a 4-parameter dynamical system
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("4-PARAMETER DYNAMICS OF SYMMETRIZED COMPOSITION")
print("="*80)

# Track a, b, c, d through self-composition
for name in ["proportional", "inverse", "modular"]:
    corr = RELATIONSHIP_SIGNATURES[name]
    ent = Entangler(corr, seed=42).build()
    current = symmetrize(ent.joint)

    print(f"\n  {name}:")
    print(f"    {'Step':>4s} {'a':>10s} {'b':>10s} {'c':>10s} {'d':>10s} {'a-b':>10s} {'Schmidt':>8s}")
    print("    " + "-" * 65)

    for step in range(10):
        a = current[0, 0].real.item()
        b = current[0, 1].real.item()
        c = current[0, H].real.item()
        d = current[H, H].real.item()
        amb = a - b
        sc = schmidt_number(current)
        print(f"    {step:>4d} {a:>10.6f} {b:>10.6f} {c:>10.6f} {d:>10.6f} {amb:>10.6f} {sc:>8.4f}")
        current = compose(current, symmetrize(ent.joint))


# ═══════════════════════════════════════════════════════════════
#  The fixed point of symmetrized composition
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SYMMETRIZED COMPOSITION FIXED POINT")
print("="*80)

# After many compositions, a → b (standard eigenvalue → 0) and
# c, d reach equilibrium values. The fixed point should be a
# rank-1 product state that is also S₃-symmetric.

for name in ["proportional", "modular"]:
    corr = RELATIONSHIP_SIGNATURES[name]
    ent = Entangler(corr, seed=42).build()
    M = symmetrize(ent.joint)

    current = M.clone()
    for _ in range(30):
        current = compose(current, M)

    a = current[0, 0].real.item()
    b = current[0, 1].real.item()
    c = current[0, H].real.item()
    d = current[H, H].real.item()

    print(f"\n  {name} fixed point:")
    print(f"    a = {a:.8f}")
    print(f"    b = {b:.8f}")
    print(f"    c = {c:.8f}")
    print(f"    d = {d:.8f}")
    print(f"    a - b = {a-b:.8f} (should → 0)")
    print(f"    a/b = {a/b:.6f} (should → 1)")
    print(f"    c/a = {c/a:.6f}")
    print(f"    d/a = {d/a:.6f}")

    # Is the fixed point = (1/4) J (uniform)?
    print(f"    1/(H+1) = {1/(H+1):.6f}")
    print(f"    Is a ≈ 1/(H+1)? diff = {abs(a - 1/(H+1)):.6f}")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
