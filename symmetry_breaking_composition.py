"""
Does composition of S₃-symmetric crystals break S₃ symmetry?

The symmetrized crystal has a=b (all h-entries equal). Principle 32
says composition preserves S₃ symmetry. But the fixed point shows
f_std = 14-21%.

Either:
  A) Composition of real (not averaged) symmetrized crystals breaks
     symmetry through numerical effects
  B) The S₃ projector has a subtlety with complex phases
  C) The compose() function introduces asymmetry

Let me check by composing an EXACTLY symmetric [4,4] matrix with itself.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities, EPS_LOG)
from solver.crystals import compose

torch.set_grad_enabled(False)


def s3_fractions(joint):
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    signs = [1, -1, -1, -1, 1, 1]
    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M.real.float() @ P.T

    triv = sum(permute(joint, p) for p in perms) / 6
    sign = sum(s * permute(joint, p) for p, s in zip(perms, signs)) / 6
    bp = joint.abs().pow(2)
    total = bp.sum().item()
    if total < 1e-15:
        return 1.0, 0.0, 0.0
    f_triv = triv.to(joint.dtype).abs().pow(2).sum().item() / total
    f_sign = sign.to(joint.dtype).abs().pow(2).sum().item() / total
    return f_triv, f_sign, max(0, 1 - f_triv - f_sign)


def check_s3_symmetry(joint, label=""):
    """Check if joint[i,j] = joint[σ(i),σ(j)] for all σ ∈ S₃."""
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    max_diff = 0
    for p in perms[1:]:  # skip identity
        for i in range(MASS_DIM):
            for j in range(MASS_DIM):
                diff = abs(joint[i,j] - joint[p[i], p[j]])
                max_diff = max(max_diff, diff)
    return max_diff


print("=" * 80)
print("COMPOSITION AND S₃ SYMMETRY")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Build exactly S₃-symmetric crystal
# ═══════════════════════════════════════════════════════════════

# From Principle 32: M_sym = [[a,b,b,c],[b,a,b,c],[b,b,a,c],[c,c,c,d]]
a = 0.30 + 0.02j
b = 0.15 + 0.01j
c = 0.10 - 0.01j
d = 0.05 + 0.005j

exact_sym = torch.tensor([
    [a, b, b, c],
    [b, a, b, c],
    [b, b, a, c],
    [c, c, c, d],
], dtype=torch.cfloat)

# Normalize to L1 = 1
exact_sym = exact_sym / exact_sym.real.sum()

print(f"\n--- Exactly S₃-symmetric crystal ---")
print(f"  Max S₃ asymmetry: {check_s3_symmetry(exact_sym):.2e}")
f_t, f_s, f_st = s3_fractions(exact_sym)
print(f"  Irreps: triv={f_t:.6f}, sign={f_s:.6f}, std={f_st:.6f}")


# Compose with itself
result = compose(exact_sym, exact_sym)
print(f"\n  After 1 composition:")
print(f"  Max S₃ asymmetry: {check_s3_symmetry(result):.2e}")
f_t, f_s, f_st = s3_fractions(result)
print(f"  Irreps: triv={f_t:.6f}, sign={f_s:.6f}, std={f_st:.6f}")

# Compose many times
current = exact_sym.clone()
for step in range(20):
    current = compose(current, exact_sym)
    if step % 5 == 4:
        asym = check_s3_symmetry(current)
        f_t, f_s, f_st = s3_fractions(current)
        print(f"  After {step+1:2d} compositions: asym={asym:.2e}, "
              f"triv={f_t:.6f}, sign={f_s:.6f}, std={f_st:.6f}")


# ═══════════════════════════════════════════════════════════════
#  The projector on COMPLEX masses
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE PROJECTOR ON COMPLEX vs REAL MASSES")
print("="*80)

# The issue: our projector permutes M.real.float() but the crystal
# has complex entries. Complex phases break S₃ symmetry even when
# the magnitudes are symmetric.

# Check: does |M|² have S₃ symmetry?
print(f"\n  Exact symmetric crystal:")
born = exact_sym.abs().pow(2)
born_sym_check = check_s3_symmetry(born.to(torch.cfloat))
print(f"    |M|² S₃ asymmetry: {born_sym_check:.2e}")

# Check phases
print(f"\n  Phases of the symmetric crystal:")
for i in range(MASS_DIM):
    phases = [torch.angle(exact_sym[i,j]).item() / math.pi for j in range(MASS_DIM)]
    print(f"    row {i}: [{', '.join(f'{p:+.4f}π' for p in phases)}]")


# ═══════════════════════════════════════════════════════════════
#  The correct S₃ projector for complex masses
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CORRECT S₃ PROJECTOR FOR COMPLEX MASSES")
print("="*80)

def s3_fractions_complex(joint):
    """S₃ projections that properly handle complex entries."""
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    signs = [1, -1, -1, -1, 1, 1]

    def permute_complex(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T

    triv = sum(permute_complex(joint, p) for p in perms) / 6
    sign = sum(s * permute_complex(joint, p) for p, s in zip(perms, signs)) / 6
    std = joint - triv - sign

    bp = joint.abs().pow(2)
    total = bp.sum().item()
    if total < 1e-15:
        return 1.0, 0.0, 0.0

    f_triv = triv.abs().pow(2).sum().item() / total
    f_sign = sign.abs().pow(2).sum().item() / total
    f_std = std.abs().pow(2).sum().item() / total
    return f_triv, f_sign, f_std


# Compare the two projectors
print(f"\n  Exactly S₃-symmetric crystal:")
f_t1, f_s1, f_st1 = s3_fractions(exact_sym)
f_t2, f_s2, f_st2 = s3_fractions_complex(exact_sym)
print(f"    Old projector (real): triv={f_t1:.6f}, sign={f_s1:.6f}, std={f_st1:.6f}")
print(f"    New projector (cmplx): triv={f_t2:.6f}, sign={f_s2:.6f}, std={f_st2:.6f}")

# Now test on actual entangler crystals
from solver.crystals import Entangler
identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

print(f"\n  Actual identity crystals:")
for seed in [0, 42, 99]:
    ent = Entangler(identity_corr, seed=seed).build()
    crystal = ent.joint.clone()

    f_t1, f_s1, f_st1 = s3_fractions(crystal)
    f_t2, f_s2, f_st2 = s3_fractions_complex(crystal)

    print(f"    seed {seed:3d}: old std={f_st1:.5f}, new std={f_st2:.5f}, diff={f_st2-f_st1:+.5f}")

# Fixed point with correct projector
print(f"\n  Fixed points (25 compositions):")
for seed in [0, 42, 99]:
    ent = Entangler(identity_corr, seed=seed).build()
    current = ent.joint.clone()
    original = ent.joint.clone()
    for _ in range(25):
        current = compose(current, original)

    f_t1, f_s1, f_st1 = s3_fractions(current)
    f_t2, f_s2, f_st2 = s3_fractions_complex(current)

    print(f"    seed {seed:3d}: old std={f_st1:.5f}, new std={f_st2:.5f}, diff={f_st2-f_st1:+.5f}")


# ═══════════════════════════════════════════════════════════════
#  S₃-symmetrized, then compose, with CORRECT projector
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SYMMETRIZED + CORRECT PROJECTOR")
print("="*80)

perms_list = [
    [0,1,2,3], [1,0,2,3], [2,1,0,3],
    [0,2,1,3], [1,2,0,3], [2,0,1,3],
]

for seed in [0, 42, 99]:
    ent = Entangler(identity_corr, seed=seed).build()
    crystal = ent.joint.clone()

    # S₃-symmetrize (complex)
    sym = torch.zeros_like(crystal)
    for p in perms_list:
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=crystal.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        sym += P @ crystal @ P.T
    sym /= 6

    # Verify symmetry
    asym = check_s3_symmetry(sym)

    # Compose
    current = sym.clone()
    for _ in range(25):
        current = compose(current, sym)

    f_t, f_s, f_st = s3_fractions_complex(current)
    asym_fp = check_s3_symmetry(current)

    print(f"  seed {seed}: sym_asym={asym:.2e}, fp_asym={asym_fp:.2e}, "
          f"triv={f_t:.6f}, sign={f_s:.6f}, std={f_st:.6f}")


print(f"\n\n{'='*80}")
print("WHAT THIS MEANS")
print("="*80)
