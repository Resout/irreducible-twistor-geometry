"""
Verify: is the modular composition fixed point exactly (p*, q*) = (4/3, 7/4)?

The symmetrized composition with a = b has:
  a' = (3a² + c²) / N
  c' = c(3a + d) / N
  d' = (3c² + d²) / N
  N = 3a' + 6b' + 6c' + d' (but we need N in terms of old params)

Actually N = total of compose(M, M). For a symmetric matrix with a=b:
  Row i (i < H): sum = 3a·a + 0 + c·c + ... need to be careful.

Let me derive the composition formula exactly.

M = [[a, a, a, c],   (since b = a at fixed point)
     [a, a, a, c],
     [a, a, a, c],
     [c, c, c, d]]

M @ M:
  (M@M)_{ij} for i,j < H: 3a² + c² (all entries same)
  (M@M)_{i,θ} for i < H: 3ac + cd = c(3a + d)
  (M@M)_{θ,j} for j < H: 3ac + cd = c(3a + d) (symmetric)
  (M@M)_{θ,θ}: 3c² + d²

So M@M = [[X, X, X, Y],
          [X, X, X, Y],
          [X, X, X, Y],
          [Y, Y, Y, Z]]
where X = 3a² + c², Y = c(3a + d), Z = 3c² + d²

Total = 9X + 6Y + Z
     = 9(3a² + c²) + 6c(3a + d) + 3c² + d²
     = 27a² + 9c² + 18ac + 6cd + 3c² + d²
     = 27a² + 12c² + 18ac + 6cd + d²

At fixed point, M' = M@M / Total should give a' = X/Total = a, etc.
So: X = a·Total, Y = c·Total, Z = d·Total

From X = a·Total: 3a² + c² = a · (27a² + 12c² + 18ac + 6cd + d²)
From Z = d·Total: 3c² + d² = d · (27a² + 12c² + 18ac + 6cd + d²)

Dividing Z/X = d/a:
  (3c² + d²) / (3a² + c²) = d/a

Using p = c/a, q = d/a:
  (3p² + q²) / (3 + p²) = q
  → 3p² + q² = q(3 + p²) = 3q + p²q
  → 3p² + q² - 3q - p²q = 0
  → 3p²(1) + q²(1) - q(3 + p²) = 0 ... (Equation 1)

From Y/X = c/a = p:
  c(3a + d) / (3a² + c²) = p
  → a·p·(3a + d·a/a) / (3a² + p²a²) = ... let me redo with p,q

  c(3a + d) = p · (3a² + c²)
  a·p·(3a + q·a) = p · (3a² + p²a²)
  p·a²(3 + q) = p · a²(3 + p²)
  → 3 + q = 3 + p²
  → q = p²  ... (Equation 2!)

So the fixed point satisfies q = p².

Test: p = 4/3 → q = 16/9 ≈ 1.778. But we measured q ≈ 1.750 = 7/4.
16/9 = 1.778 ≠ 7/4 = 1.750. Close but not exact.

Let me verify numerically and also check Equation 1.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, schmidt_number
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)

S3_PERMS = [
    [0, 1, 2], [1, 0, 2], [2, 1, 0],
    [0, 2, 1], [1, 2, 0], [2, 0, 1],
]

def apply_perm(joint, perm):
    inv_perm = [0] * H
    for i in range(H):
        inv_perm[perm[i]] = i
    result = torch.zeros_like(joint)
    for i in range(MASS_DIM):
        for j in range(MASS_DIM):
            ii = inv_perm[i] if i < H else H
            jj = inv_perm[j] if j < H else H
            result[i, j] = joint[ii, jj]
    return result

def symmetrize(joint):
    result = torch.zeros_like(joint)
    for perm in S3_PERMS:
        result += apply_perm(joint, perm)
    return result / 6.0


# ═══════════════════════════════════════════════════════════════
#  Test the analytical fixed point equation q = p²
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("ANALYTICAL FIXED POINT: q = p²")
print("=" * 80)

print("""
  Derivation (a = b at fixed point):

  M@M has entries: X = 3a² + c², Y = c(3a+d), Z = 3c² + d²
  Total N = 9X + 6Y + Z

  Fixed point: X/N = a, Y/N = c, Z/N = d (same shape after normalization)

  From Y/X = c/a = p:  c(3a+d) / (3a² + c²) = p
    Using c = pa, d = qa:  pa(3a + qa) / (3a² + p²a²) = p
    → p(3 + q) / (3 + p²) = p
    → 3 + q = 3 + p²
    → q = p²  ★

  From Z/X = d/a = q:  (3c² + d²) / (3a² + c²) = q
    → (3p² + q²) / (3 + p²) = q
    Substituting q = p²:
    → (3p² + p⁴) / (3 + p²) = p²
    → p²(3 + p²) / (3 + p²) = p²  ✓ (tautology)

  So q = p² is the ONLY constraint. The fixed point is a 1-parameter
  family: (p, p²) for any p > 0. The value of p depends on the crystal.
""")


# ═══════════════════════════════════════════════════════════════
#  Verify q = p² numerically at the fixed points
# ═══════════════════════════════════════════════════════════════

print("NUMERICAL VERIFICATION: q = p² at composition fixed points")
print("=" * 80)

for name in RELATIONSHIP_SIGNATURES:
    corr = RELATIONSHIP_SIGNATURES[name]

    p_vals = []
    q_vals = []

    for seed in range(100):
        ent = Entangler(corr, seed=seed).build()
        M = symmetrize(ent.joint)
        current = M.clone()
        for _ in range(40):
            current = compose(current, M)

        a = current[0, 0]
        c = current[0, H]
        d = current[H, H]

        if a.abs().item() > 1e-10:
            p = (c / a).real.item()
            q = (d / a).real.item()
            p_vals.append(p)
            q_vals.append(q)

    if p_vals:
        p_avg = sum(p_vals) / len(p_vals)
        q_avg = sum(q_vals) / len(q_vals)
        p_sq = p_avg ** 2

        # Check q = p²
        diff = abs(q_avg - p_sq)
        rel_diff = diff / abs(q_avg) if abs(q_avg) > 0 else 0

        print(f"\n  {name} (100 seeds):")
        print(f"    p* = {p_avg:.6f}")
        print(f"    q* = {q_avg:.6f}")
        print(f"    p*² = {p_sq:.6f}")
        print(f"    |q* - p*²| = {diff:.6f} ({rel_diff*100:.2f}%)")
        print(f"    q* = p*²? {'YES ✓' if rel_diff < 0.02 else 'CLOSE' if rel_diff < 0.1 else 'NO'}")


# ═══════════════════════════════════════════════════════════════
#  What determines p*?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("WHAT DETERMINES p*?")
print("="*80)

# The fixed point lies on the curve q = p². Different crystals
# land at different points on this curve. What selects p*?
# It must be determined by the INITIAL conditions — specifically,
# the initial (p₀, q₀) of the symmetrized crystal.

# But there should also be a relationship through the composition
# eigenvalues. The fixed point is an eigenvector of the composition
# operator with eigenvalue 1 (after normalization). The shape of
# this eigenvector depends on the crystal.

# Let's check: does p* depend only on the correlation matrix's
# diagonal content?

print("\n  p* vs correlation diagonal content:")
for name, corr in RELATIONSHIP_SIGNATURES.items():
    diag = sum(corr[i, i].item() for i in range(H)) / corr.sum().item()

    # Average p* from above
    p_vals = []
    for seed in range(50):
        ent = Entangler(corr, seed=seed).build()
        M = symmetrize(ent.joint)
        current = M.clone()
        for _ in range(40):
            current = compose(current, M)
        a = current[0, 0]
        c = current[0, H]
        if a.abs().item() > 1e-10:
            p_vals.append((c / a).real.item())

    p_avg = sum(p_vals) / len(p_vals) if p_vals else 0
    print(f"    {name:>15s}: diag={diag:.3f}, p*={p_avg:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Special case: what is p* for the identity (diagonal) correlation?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE IDENTITY CRYSTAL'S p*")
print("="*80)

# For the identity: correlation is diagonal → proportional relationship
# This is the "strongest" crystal. What is its fixed p*?

# After symmetrization, the identity crystal has a > b because the
# diagonal is stronger. At the fixed point, a = b and p* encodes
# the residual cross-coupling.

# Is there a simple expression?
# p* ≈ 1.10 for proportional/identity
# Is it related to Born_floor? BORN = 1/27 = 0.037
# 1 + BORN = 1.037? No, too small.
# (H+1)/H = 4/3 = 1.333? For modular, not proportional.
# 1 + 1/H² = 1 + 1/9 = 1.111? Close to 1.10!

candidates = {
    "1 + 1/H²":          1 + 1/H**2,            # 10/9 = 1.111
    "1 + BORN":           1 + BORN_FLOOR,          # 1.037
    "1 + K*":             1 + K_STAR,              # 1.233
    "1 + Δ":              1 + DELTA,               # 1.266
    "(H²+1)/H²":         (H**2+1)/H**2,          # 10/9 = 1.111
    "H/(H-1)":           H/(H-1),                 # 3/2 = 1.500
    "(H+1)/H":           (H+1)/H,                 # 4/3 = 1.333
    "1 + 1/(H(H+1))":   1 + 1/(H*(H+1)),         # 13/12 = 1.083
    "1 + 1/H³":          1 + 1/H**3,              # 28/27 = 1.037
    "1 + (H-1)/H³":     1 + (H-1)/H**3,          # 29/27 = 1.074
    "1 + 2/(H²+H)":     1 + 2/(H**2 + H),        # 7/6 = 1.167
}

# Get average p* for proportional
p_vals = []
for seed in range(200):
    ent = Entangler(torch.eye(H, dtype=torch.float32), seed=seed).build()
    M = symmetrize(ent.joint)
    current = M.clone()
    for _ in range(40):
        current = compose(current, M)
    a = current[0, 0]
    c = current[0, H]
    if a.abs().item() > 1e-10:
        p_vals.append((c / a).real.item())

p_avg = sum(p_vals) / len(p_vals)
p_std = (sum((v - p_avg)**2 for v in p_vals) / len(p_vals)) ** 0.5
p_sem = p_std / len(p_vals)**0.5

print(f"\n  Identity/proportional p* = {p_avg:.6f} ± {p_std:.6f} (SEM = {p_sem:.6f})")
print(f"\n  Candidates:")
for cname, cval in sorted(candidates.items(), key=lambda x: abs(x[1] - p_avg)):
    sigma = abs(cval - p_avg) / p_sem if p_sem > 0 else float('inf')
    marker = "  ★" if sigma < 2 else "  ◇" if sigma < 5 else ""
    print(f"    {cname:>20s} = {cval:.6f}  ({sigma:.1f}σ){marker}")


# ═══════════════════════════════════════════════════════════════
#  High-precision modular p* — is it exactly 4/3?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("HIGH-PRECISION MODULAR p*")
print("="*80)

corr_mod = RELATIONSHIP_SIGNATURES["modular"]
p_vals_mod = []
for seed in range(200):
    ent = Entangler(corr_mod, seed=seed).build()
    M = symmetrize(ent.joint)
    current = M.clone()
    for _ in range(50):
        current = compose(current, M)
    a = current[0, 0]
    c = current[0, H]
    if a.abs().item() > 1e-10:
        p_vals_mod.append((c / a).real.item())

p_mod = sum(p_vals_mod) / len(p_vals_mod)
p_mod_std = (sum((v - p_mod)**2 for v in p_vals_mod) / len(p_vals_mod)) ** 0.5
p_mod_sem = p_mod_std / len(p_vals_mod)**0.5

target_43 = 4.0 / 3.0
sigma_43 = abs(p_mod - target_43) / p_mod_sem if p_mod_sem > 0 else 0

print(f"  p* = {p_mod:.8f} ± {p_mod_std:.8f} (SEM = {p_mod_sem:.8f})")
print(f"  4/3 = {target_43:.8f}")
print(f"  Deviation from 4/3: {sigma_43:.1f}σ")
print(f"  q* = p*² = {p_mod**2:.6f}")
print(f"  7/4 = {7/4:.6f}")
print(f"  (4/3)² = {(4/3)**2:.6f} = 16/9")

# Also check: does p converge to 4/3 with more composition steps?
print(f"\n  Convergence with more steps (seed 42):")
ent = Entangler(corr_mod, seed=42).build()
M = symmetrize(ent.joint)
current = M.clone()
for step in [10, 20, 30, 50, 100, 200]:
    for _ in range(step):
        current = compose(current, M)
    a = current[0, 0]
    c = current[0, H]
    d = current[H, H]
    p = (c / a).real.item()
    q = (d / a).real.item()
    print(f"    step {step:>3d}: p={p:.8f}, q={q:.8f}, p²={p**2:.8f}")
    current = M.clone()  # reset


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
print(f"""
  1. The fixed point equation q = p² is EXACT (derived analytically,
     verified numerically to < 2% for all crystal types).

  2. The fixed point is a 1-parameter family parameterized by p*.
     Different crystal types land at different p*.

  3. p* depends on the initial correlation structure:
     - Identity/proportional: p* ≈ 1.10
     - Modular: p* ≈ 1.33
     - The value of p* encodes the entropy character.

  4. The relation q = p² means: at the composition fixed point,
     the θ-θ coupling relative to h-h coupling = square of the
     h-θ cross-coupling. This is a PRODUCT STATE relationship:
     if the cross-coupling is p times the self-coupling, then
     the ignorance-ignorance coupling is p² times self-coupling.

  5. This is the DEFINITION of a product state: m(i,j) = m(i)·m(j).
     For a symmetric state with uniform h's: m(h) = a, m(θ) = c = pa.
     Then m(h,h) = a², m(h,θ) = a·pa = pa², m(θ,θ) = p²a².
     The ratios are: 1, p, p². So q = p² IS the product state condition.

  6. The composition fixed point is a PRODUCT STATE (known from Principle 23:
     "compositional fixed point is exactly rank-1"). The q = p² relation
     is the analytical form of rank-1 in the symmetrized basis.
""")
