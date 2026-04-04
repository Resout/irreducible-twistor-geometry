"""
The irreducible 2D dynamics of the trivial sector.

After symmetrization, the crystal has 4 parameters (a, b, c, d).
Under composition, a → b rapidly. Once a = b, the system is 2D:
  - p = c/a (cross-coupling ratio)
  - q = d/a (ignorance ratio)

Questions:
1. What are the dynamical equations p(n+1) = f(p(n), q(n))?
2. What is the fixed point (p*, q*)?
3. Is it the same for all crystal types?
4. What determines the convergence rate to (p*, q*)?
5. Can this be solved analytically?
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

def extract_params(joint_sym):
    """Extract (a, b, c, d) from symmetrized crystal."""
    a = joint_sym[0, 0]
    b = joint_sym[0, 1]
    c = joint_sym[0, H]
    d = joint_sym[H, H]
    return a, b, c, d


# ═══════════════════════════════════════════════════════════════
#  Derive the composition formula for (a, b, c, d)
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("COMPOSITION IN THE SYMMETRIZED BASIS")
print("=" * 80)

# compose(A, B)_{ij} = Σ_k A_{ik} B_{kj}  (then normalize)
# For A = B = M_sym:
# compose(M, M)_{00} = a*a + b*b + b*b + c*c = a² + 2b² + c²
# compose(M, M)_{01} = a*b + b*a + b*b + c*c = 2ab + b² + c²
# compose(M, M)_{03} = a*c + b*c + b*c + c*d = ac + 2bc + cd
# compose(M, M)_{33} = c*c + c*c + c*c + d*d = 3c² + d²
# Normalization: sum of all 16 entries

# Let's verify numerically
for name in ["proportional"]:
    corr = RELATIONSHIP_SIGNATURES[name]
    ent = Entangler(corr, seed=42).build()
    M = symmetrize(ent.joint)
    a, b, c, d = extract_params(M)

    # Numerical composition
    MC = compose(M, M)
    a2, b2, c2, d2 = extract_params(MC)

    # Analytical formulas (unnormalized)
    a_raw = a*a + 2*b*b + c*c
    b_raw = 2*a*b + b*b + c*c
    c_raw = a*c + 2*b*c + c*d
    d_raw = 3*c*c + d*d

    # Total for normalization
    total_raw = 3*a_raw + 6*b_raw + 6*c_raw + d_raw
    total_actual = 3*a2 + 6*b2 + 6*c2 + d2

    print(f"\n  {name}:")
    print(f"    a² + 2b² + c²:    raw={a_raw.abs().item():.6f}, actual a'={a2.abs().item():.6f}")
    print(f"    2ab + b² + c²:    raw={b_raw.abs().item():.6f}, actual b'={b2.abs().item():.6f}")
    print(f"    ac + 2bc + cd:    raw={c_raw.abs().item():.6f}, actual c'={c2.abs().item():.6f}")
    print(f"    3c² + d²:         raw={d_raw.abs().item():.6f}, actual d'={d2.abs().item():.6f}")

    # Check if compose normalizes by dividing by row sums
    # Row 0 of M@M (unnormalized): [a_raw, b_raw, b_raw, c_raw]
    row0_sum = a_raw + 2*b_raw + c_raw
    print(f"\n    Row 0 sum (unnorm): {row0_sum.abs().item():.6f}")
    print(f"    a'/a_raw = {(a2/a_raw).abs().item():.6f}")
    print(f"    b'/b_raw = {(b2/b_raw).abs().item():.6f}")

    # Is it row-normalization?
    print(f"    a_raw/row0_sum = {(a_raw/row0_sum).abs().item():.6f}")
    print(f"    Actual a' = {a2.abs().item():.6f}")

    # Or is it column-wise? Or total normalization?
    total_all = M.sum()
    print(f"\n    Total M entries: {total_all.abs().item():.6f}")
    print(f"    M@M unnorm total: {(M @ M).sum().abs().item():.6f}")
    print(f"    compose total: {MC.sum().abs().item():.6f}")


# ═══════════════════════════════════════════════════════════════
#  Actually check what compose() does to understand normalization
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMPOSE() NORMALIZATION")
print("="*80)

# Read the compose function to understand normalization
# compose(A, B) does: result = A @ B, then normalize per row (from crystals.py)

# Let's check by doing manual computation
M_manual = M @ M
# Normalize each row to sum to ... what?
for i in range(MASS_DIM):
    row_sum = M_manual[i, :].sum()
    mc_row_sum = MC[i, :].sum()
    print(f"  Row {i}: M@M sum = {row_sum.real.item():.6f}+{row_sum.imag.item():.6f}i, "
          f"compose sum = {mc_row_sum.real.item():.6f}+{mc_row_sum.imag.item():.6f}i")

# Total normalization
print(f"  Total M@M: {M_manual.sum().real.item():.6f}")
print(f"  Total compose: {MC.sum().real.item():.6f}")


# Read the actual compose function
print(f"\n  Checking compose implementation...")
from solver.crystals import compose as _compose
import inspect
src = inspect.getsource(_compose)
# Just check: does it normalize by row, total, or re-anchor?
if "sum(dim=-1" in src or "sum(dim=1" in src:
    print("    Normalizes per ROW")
elif ".sum()" in src:
    print("    Normalizes by TOTAL")
print(f"    (See crystals.py for compose implementation)")


# ═══════════════════════════════════════════════════════════════
#  Track the 2D system: p = c/a, q = d/a (once a ≈ b)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("2D DYNAMICS: p = c/a, q = d/a")
print("="*80)

for name in RELATIONSHIP_SIGNATURES:
    corr = RELATIONSHIP_SIGNATURES[name]
    ent = Entangler(corr, seed=42).build()
    M = symmetrize(ent.joint)

    current = M.clone()
    print(f"\n  {name}:")
    print(f"    {'Step':>4s} {'p=c/a':>10s} {'q=d/a':>10s} {'a-b':>12s} {'Schmidt':>8s}")
    print("    " + "-" * 50)

    for step in range(15):
        a, b, c, d = extract_params(current)
        if a.abs().item() > 1e-10:
            p = (c / a).real.item()
            q = (d / a).real.item()
        else:
            p, q = 0, 0
        amb = (a - b).abs().item()
        sc = schmidt_number(current)
        print(f"    {step:>4d} {p:>10.6f} {q:>10.6f} {amb:>12.2e} {sc:>8.4f}")
        current = compose(current, M)


# ═══════════════════════════════════════════════════════════════
#  Fixed point (p*, q*) across many seeds
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("FIXED POINT (p*, q*) ACROSS SEEDS")
print("="*80)

for name in ["proportional", "inverse", "modular"]:
    corr = RELATIONSHIP_SIGNATURES[name]
    p_stars = []
    q_stars = []

    for seed in range(50):
        ent = Entangler(corr, seed=seed).build()
        M = symmetrize(ent.joint)
        current = M.clone()
        for _ in range(30):
            current = compose(current, M)
        a, b, c, d = extract_params(current)
        if a.abs().item() > 1e-10:
            p_stars.append((c / a).real.item())
            q_stars.append((d / a).real.item())

    if p_stars:
        p_avg = sum(p_stars) / len(p_stars)
        p_std = (sum((v - p_avg)**2 for v in p_stars) / len(p_stars)) ** 0.5
        q_avg = sum(q_stars) / len(q_stars)
        q_std = (sum((v - q_avg)**2 for v in q_stars) / len(q_stars)) ** 0.5

        print(f"\n  {name} (50 seeds):")
        print(f"    p* = c*/a* = {p_avg:.6f} ± {p_std:.6f}")
        print(f"    q* = d*/a* = {q_avg:.6f} ± {q_std:.6f}")

        # Is p* = 1? (would mean c = a at fixed point)
        print(f"    p* - 1 = {p_avg - 1:.6f}")

        # Candidate expressions for q*
        candidates = {
            "1":           1.0,
            "(H-1)/H":     (H-1)/H,
            "H/(H+1)":     H/(H+1),
            "1 - BORN":    1 - BORN_FLOOR,
            "1 - K*":      1 - K_STAR,
            "(H²-1)/H²":  (H**2-1)/H**2,
        }
        print(f"    q* candidates:")
        for cname, cval in sorted(candidates.items(), key=lambda x: abs(x[1] - q_avg)):
            print(f"      {cname:>12s} = {cval:.6f} (diff = {abs(cval - q_avg):.4f})")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
