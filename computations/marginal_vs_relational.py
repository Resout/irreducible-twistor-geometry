"""
Marginal vs relational information in a crystal.

A crystal contains two types of information:
1. Marginal: what each variable looks like individually (u⊗v product)
2. Relational: how the variables correlate (entanglement)

Composition destroys relational while preserving marginal.
The Born floor constrains marginal. Nothing floors relational.

Question: for a given crystal, what fraction of its information
is marginal vs relational? Does this fraction follow any law?

The rank-1 approximation gives the best product state.
The residual (crystal minus product) is the relational part.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, BORN_FLOOR, schmidt_number
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)


def decompose_crystal(joint):
    """Decompose a crystal into marginal (product) and relational parts.

    Returns:
      product: the best rank-1 approximation (σ₀ * u₀ ⊗ v₀†)
      relational: joint - product
      marginal_fraction: ||product||² / ||joint||²
      relational_fraction: ||relational||² / ||joint||²
    """
    U, S, Vh = torch.linalg.svd(joint)

    # Best rank-1 approximation
    product = S[0] * torch.einsum("i,j->ij", U[:, 0], Vh[0, :].conj())

    # Relational residual
    relational = joint - product

    total_sq = joint.abs().pow(2).sum().item()
    product_sq = product.abs().pow(2).sum().item()
    relational_sq = relational.abs().pow(2).sum().item()

    # Note: due to SVD properties, product_sq + relational_sq ≈ total_sq
    # (Pythagorean theorem in Frobenius norm)

    return {
        "product": product,
        "relational": relational,
        "marginal_fraction": product_sq / total_sq if total_sq > 0 else 0,
        "relational_fraction": relational_sq / total_sq if total_sq > 0 else 0,
        "sigma_0": S[0].item(),
        "sigma_rest": S[1:].tolist(),
        "schmidt": schmidt_number(joint),
    }


# ═══════════════════════════════════════════════════════════════
#  Decompose all signature crystals
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("MARGINAL vs RELATIONAL INFORMATION")
print("=" * 70)

print(f"\n  {'Crystal':>15s} {'Schmidt':>8s} {'Marginal%':>10s} {'Relational%':>12s} {'σ₀':>8s} {'σ₁':>8s} {'σ₂':>8s} {'σ₃':>8s}")
print("  " + "-" * 80)

all_data = []
for name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Entangler(corr).build()
    d = decompose_crystal(ent.joint)

    sigmas = [d["sigma_0"]] + [s for s in d["sigma_rest"]]
    print(f"  {name:>15s} {d['schmidt']:>8.3f} {d['marginal_fraction']*100:>9.1f}% "
          f"{d['relational_fraction']*100:>11.1f}% "
          f"{sigmas[0]:>8.4f} {sigmas[1]:>8.4f} {sigmas[2]:>8.4f} {sigmas[3]:>8.4f}")

    all_data.append({"name": name, **d})


# ═══════════════════════════════════════════════════════════════
#  Marginal fraction vs Schmidt
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("MARGINAL FRACTION vs SCHMIDT")
print("=" * 70)

schmidts = [d["schmidt"] for d in all_data]
marg_fracs = [d["marginal_fraction"] for d in all_data]

# Is there a relationship?
def corr(xs, ys):
    n = len(xs)
    mx, my = sum(xs)/n, sum(ys)/n
    c = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx = sum((x-mx)**2 for x in xs)**0.5
    sy = sum((y-my)**2 for y in ys)**0.5
    return c/(sx*sy) if sx>1e-10 and sy>1e-10 else 0

r = corr(schmidts, marg_fracs)
print(f"\n  Correlation(Schmidt, marginal_fraction): r = {r:.3f}")

# Theoretical: for a [4,4] matrix with Schmidt S,
# the marginal fraction = σ₀² / Σσᵢ² = 1/S (by definition of Schmidt as IPR)
# So marginal_fraction should be 1/Schmidt exactly!

print(f"\n  Predicted: marginal_fraction = 1/Schmidt")
for d in all_data:
    predicted = 1.0 / d["schmidt"]
    actual = d["marginal_fraction"]
    ratio = actual / predicted if predicted > 0 else 0
    print(f"  {d['name']:>15s}: actual={actual:.4f}, predicted={predicted:.4f}, ratio={ratio:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Track the decomposition through composition
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("DECOMPOSITION THROUGH COMPOSITION")
print("How does the marginal/relational split evolve?")
print("=" * 70)

ent = Entangler(torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)).build()
current = ent.joint.clone()

print(f"\n  {'comp':>5s} {'Schmidt':>8s} {'Marginal%':>10s} {'Relational%':>12s} {'1/Schmidt':>10s}")
print("  " + "-" * 50)

for step in range(8):
    d = decompose_crystal(current)
    predicted = 1.0 / d["schmidt"]
    print(f"  {step:>5d} {d['schmidt']:>8.3f} {d['marginal_fraction']*100:>9.1f}% "
          f"{d['relational_fraction']*100:>11.1f}% {predicted*100:>9.1f}%")
    current = compose(current, ent.joint)


# ═══════════════════════════════════════════════════════════════
#  The relational content: what does it encode?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("WHAT DOES THE RELATIONAL PART ENCODE?")
print("=" * 70)

for name in ["proportional", "inverse", "modular"]:
    ent = Entangler(RELATIONSHIP_SIGNATURES[name]).build()
    d = decompose_crystal(ent.joint)

    rel = d["relational"]
    total_rel = rel.abs().pow(2).sum().item()

    if total_rel < 1e-10:
        continue

    # Block decomposition of the relational part
    rel_hh = rel[:H, :H].abs().pow(2).sum().item() / total_rel
    rel_ht = rel[:H, H].abs().pow(2).sum().item() / total_rel
    rel_th = rel[H, :H].abs().pow(2).sum().item() / total_rel
    rel_tt = rel[H, H].abs().pow(2).item() / total_rel

    print(f"\n  {name} relational part ({d['relational_fraction']*100:.1f}% of crystal):")
    print(f"    h×h: {rel_hh*100:.1f}%")
    print(f"    h×θ: {rel_ht*100:.1f}%")
    print(f"    θ×h: {rel_th*100:.1f}%")
    print(f"    θ×θ: {rel_tt*100:.1f}%")

    # Is the relational part the ENTANGLEMENT?
    # It should have Schmidt > 1 (it's what makes the whole crystal entangled)
    sn_rel = schmidt_number(rel)
    print(f"    Schmidt of relational part: {sn_rel:.3f}")


# ═══════════════════════════════════════════════════════════════
#  The information partition
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("THE INFORMATION PARTITION")
print("How is information distributed across singular values?")
print("=" * 70)

for name in ["proportional", "inverse", "modular", "exponential"]:
    ent = Entangler(RELATIONSHIP_SIGNATURES[name]).build()
    _, S, _ = torch.linalg.svd(ent.joint)

    s_sq = S.abs().pow(2)
    total = s_sq.sum().item()
    cumulative = 0

    print(f"\n  {name}:")
    for i in range(MASS_DIM):
        frac = s_sq[i].item() / total
        cumulative += frac
        print(f"    σ_{i}² / total = {frac:.4f} ({frac*100:.1f}%), cumulative = {cumulative:.4f} ({cumulative*100:.1f}%)")

    # Shannon entropy of the SV distribution
    sv_probs = s_sq / s_sq.sum()
    sv_entropy = -(sv_probs * torch.log(sv_probs.clamp(min=1e-10))).sum().item()
    max_entropy = math.log(MASS_DIM)
    print(f"    SV entropy: {sv_entropy:.4f} / {max_entropy:.4f} = {sv_entropy/max_entropy:.4f}")
    print(f"    Schmidt = 1 / Σ(σᵢ²/total)² = {1.0 / (sv_probs**2).sum().item():.4f}")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
