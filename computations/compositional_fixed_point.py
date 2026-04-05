"""
What survives infinite composition?

The compositional path converges. After ~5 self-compositions,
the crystal stops moving (Fisher distance plateaus). The fixed
point has Schmidt ≈ 1.04 and Fisher ≈ 0.49 from uniform.

It's NOT the uniform distribution. Something survives.
What is it? What structure does infinite composition preserve?

And: does the fixed point depend on the STARTING crystal,
or is it universal (the same for all crystals)?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (
    H, MASS_DIM, BORN_FLOOR, schmidt_number,
    born_probabilities, sym2_fingerprint,
    EPS_LOG,
)
from solver.crystals import (
    Entangler, compose, classify_relationship,
    slot_measure, RELATIONSHIP_SIGNATURES,
)

torch.set_grad_enabled(False)


def fisher_distance(j1, j2):
    p = j1.abs().pow(2).reshape(-1)
    q = j2.abs().pow(2).reshape(-1)
    p = p / p.sum().clamp(min=1e-10)
    q = q / q.sum().clamp(min=1e-10)
    bc = (p * q).sqrt().sum().item()
    return 2 * math.acos(min(bc, 1.0))


# ═══════════════════════════════════════════════════════════════
#  Compute the fixed point for each crystal type
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("THE COMPOSITIONAL FIXED POINT")
print("What survives infinite self-composition?")
print("=" * 70)

fixed_points = {}

for name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Entangler(corr).build()
    current = ent.joint.clone()

    for _ in range(20):  # 20 self-compositions (converges by ~7)
        current = compose(current, ent.joint)

    fixed_points[name] = current.clone()
    sn = schmidt_number(current)
    rel, _ = classify_relationship(current)
    sm = slot_measure(current)

    total = current.abs().pow(2).sum().item()
    hh = current[:H, :H].abs().pow(2).sum().item() / total
    theta_total = 1 - hh
    diag = sum(current[i, i].abs().pow(2).item() for i in range(MASS_DIM)) / total
    tr = sum(current[i, i] for i in range(MASS_DIM))

    print(f"\n  {name}:")
    print(f"    Schmidt = {sn:.4f}")
    print(f"    Rel = {rel}, Slot = {max(sm, key=sm.get)}")
    print(f"    Trace (real) = {tr.real.item():.6f}")
    print(f"    h×h = {hh*100:.1f}%, θ-total = {theta_total*100:.1f}%")
    print(f"    Diagonal = {diag*100:.1f}%")


# ═══════════════════════════════════════════════════════════════
#  Are all fixed points the same?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("ARE ALL FIXED POINTS THE SAME?")
print("Fisher distances between fixed points of different crystals")
print("=" * 70)

fp_names = list(fixed_points.keys())
print(f"\n  {'':>15s}", end="")
for n in fp_names:
    print(f" {n[:6]:>7s}", end="")
print()

for i, n1 in enumerate(fp_names):
    print(f"  {n1:>15s}", end="")
    for j, n2 in enumerate(fp_names):
        d = fisher_distance(fixed_points[n1], fixed_points[n2])
        print(f" {d:>7.4f}", end="")
    print()


# ═══════════════════════════════════════════════════════════════
#  The Born structure of the fixed point
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("BORN STRUCTURE OF THE FIXED POINT")
print("=" * 70)

for name in ["proportional", "inverse", "modular"]:
    fp = fixed_points[name]
    total = fp.abs().pow(2).sum().item()

    print(f"\n  {name} fixed point:")
    print(f"    Born mass matrix (% of total):")
    for i in range(MASS_DIM):
        row = []
        for j in range(MASS_DIM):
            val = fp[i, j].abs().pow(2).item() / total * 100
            row.append(f"{val:6.2f}")
        labels = ["h0", "h1", "h2", "θ"]
        print(f"      {labels[i]}: [{', '.join(row)}]")

    # The raw complex values
    print(f"    Raw joint mass:")
    for i in range(MASS_DIM):
        row = []
        for j in range(MASS_DIM):
            v = fp[i, j]
            row.append(f"{v.real.item():+.4f}{v.imag.item():+.4f}j")
        print(f"      [{', '.join(row)}]")


# ═══════════════════════════════════════════════════════════════
#  Eigenvalues of the fixed point
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("EIGENVALUES OF THE FIXED POINT")
print("=" * 70)

for name in ["proportional", "inverse", "modular"]:
    fp = fixed_points[name]
    eigenvalues, eigenvectors = torch.linalg.eig(fp)

    # Sort by magnitude
    mags = eigenvalues.abs()
    idx = mags.argsort(descending=True)
    eigenvalues = eigenvalues[idx]

    print(f"\n  {name}:")
    for i, ev in enumerate(eigenvalues):
        print(f"    λ_{i} = {ev.real.item():+.6f} {ev.imag.item():+.6f}j  "
              f"(|λ| = {abs(ev):.6f})")

    # Check: is the fixed point approximately rank-1?
    # If rank-1: only one nonzero eigenvalue
    ratio = eigenvalues[0].abs() / eigenvalues[1].abs() if eigenvalues[1].abs() > 1e-10 else float('inf')
    print(f"    |λ₀|/|λ₁| = {ratio:.2f}")
    if ratio > 100:
        print(f"    → Approximately rank-1 (product state)")


# ═══════════════════════════════════════════════════════════════
#  Cross-composition fixed points
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("CROSS-COMPOSITION: A∘B∘A∘B∘... does it converge?")
print("And if so, to WHAT?")
print("=" * 70)

pairs = [
    ("proportional", "inverse"),
    ("proportional", "modular"),
    ("inverse", "modular"),
]

for n1, n2 in pairs:
    j1 = Entangler(RELATIONSHIP_SIGNATURES[n1]).build().joint.clone()
    j2 = Entangler(RELATIONSHIP_SIGNATURES[n2]).build().joint.clone()

    current = j1.clone()
    schmidts = []
    for step in range(20):
        if step % 2 == 0:
            current = compose(current, j2)
        else:
            current = compose(current, j1)
        schmidts.append(schmidt_number(current))

    # Compare with self-composition fixed points
    d_fp1 = fisher_distance(current, fixed_points[n1])
    d_fp2 = fisher_distance(current, fixed_points[n2])

    print(f"\n  {n1} ∘ {n2} alternating:")
    print(f"    Schmidt: {schmidts[0]:.3f} → {schmidts[4]:.3f} → {schmidts[9]:.3f} → {schmidts[-1]:.3f}")
    print(f"    Fisher from {n1} fp: {d_fp1:.4f}")
    print(f"    Fisher from {n2} fp: {d_fp2:.4f}")


# ═══════════════════════════════════════════════════════════════
#  What does the fixed point KNOW?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("WHAT DOES THE FIXED POINT KNOW?")
print("The fixed point has Schmidt ~1.04, not exactly 1.")
print("What is the residual entanglement about?")
print("=" * 70)

fp = fixed_points["proportional"]

# SVD: what directions carry the residual structure?
U, S, Vh = torch.linalg.svd(fp)
print(f"\n  Singular values: [{', '.join(f'{s:.6f}' for s in S.abs().tolist())}]")

# The RATIO of σ₀ to σ₁ tells us how rank-1 the fixed point is
if S[1].abs() > 1e-10:
    ratio = S[0].abs() / S[1].abs()
    print(f"  σ₀/σ₁ = {ratio:.2f}")

# What is U[:, 0] (the dominant left singular vector)?
u0 = U[:, 0]
print(f"\n  Dominant left singular vector:")
for i in range(MASS_DIM):
    labels = ["h0", "h1", "h2", "θ"]
    print(f"    {labels[i]}: {u0[i].real.item():+.4f} {u0[i].imag.item():+.4f}j "
          f"(|u| = {u0[i].abs().item():.4f})")

# What is Vh[0, :] (the dominant right singular vector)?
v0 = Vh[0, :]
print(f"\n  Dominant right singular vector:")
for i in range(MASS_DIM):
    labels = ["h0", "h1", "h2", "θ"]
    print(f"    {labels[i]}: {v0[i].real.item():+.4f} {v0[i].imag.item():+.4f}j "
          f"(|v| = {v0[i].abs().item():.4f})")

# The rank-1 approximation σ₀ * u₀ ⊗ v₀†
rank1 = S[0] * torch.einsum("i,j->ij", u0, v0.conj())
residual = fp - rank1
residual_norm = residual.abs().pow(2).sum().sqrt().item()
total_norm = fp.abs().pow(2).sum().sqrt().item()
print(f"\n  Rank-1 approximation error: {residual_norm:.6f} ({residual_norm/total_norm*100:.1f}% of total)")

# What is the residual? It's the structure that ISN'T a product state.
# This is the genuine entanglement that survives infinite composition.
print(f"\n  Residual (what isn't product state):")
res_total = residual.abs().pow(2).sum().item()
res_hh = residual[:H, :H].abs().pow(2).sum().item() / res_total if res_total > 0 else 0
res_theta = 1 - res_hh if res_total > 0 else 0
print(f"    h×h fraction of residual: {res_hh*100:.1f}%")
print(f"    θ-total fraction of residual: {res_theta*100:.1f}%")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
