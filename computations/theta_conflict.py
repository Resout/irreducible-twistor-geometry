"""
What happens when the entangler receives CONFLICTING evidence?

The entangler must commit to produce high entanglement. What happens
when the correlation matrix contains contradictory signals?

Three types of contradiction:
1. Superposition: correlation is the average of two pure types
   (e.g., 0.5*identity + 0.5*inverse → [[0.5,0,0.5],[0,1,0],[0.5,0,0.5]])
2. Row conflict: some rows say proportional, others say inverse
   (e.g., [[1,0,0],[0,1,0],[1,0,0]] — first and last rows contradict)
3. Asymmetric: forward relationship differs from backward
   (correlation is NOT symmetric)

Does the entangler:
  (a) Pick one structure (winner-take-all)
  (b) Split the difference (compromise)
  (c) Collapse to noise (give up)
  (d) Something else?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
from solver.algebra import H, MASS_DIM, schmidt_number
from solver.crystals import Entangler, classify_relationship, slot_measure

torch.set_grad_enabled(False)


def theta_fingerprint(joint):
    fp = torch.stack([
        joint[0, 3], joint[1, 3], joint[2, 3],
        joint[3, 0], joint[3, 1], joint[3, 2],
    ])
    norm = fp.abs().pow(2).sum().sqrt()
    if norm > 1e-10:
        fp = fp / norm
    return fp

def theta_asymmetry(fp):
    return (fp[:3].abs() - fp[3:].abs()).abs().sum().item()

def analyze(name, corr):
    ent = Entangler(corr).build()
    sn = schmidt_number(ent.joint)
    rel, conf = classify_relationship(ent.joint)
    fp = theta_fingerprint(ent.joint)
    asymm = theta_asymmetry(fp)
    sm = slot_measure(ent.joint)
    dom = max(sm, key=sm.get)
    print(f"  {name:<40s}: Schmidt={sn:.3f}, rel={rel:13s}, "
          f"asymm={asymm:.4f}, slot={dom}")
    return ent.joint, sn, rel, asymm


# ═══════════════════════════════════════════════════════════════
#  Pure references
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("PURE REFERENCES")
print("=" * 70)

id_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
inv_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)
cyc_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)  # cyclic shift

j_id, _, _, _ = analyze("Identity (proportional)", id_corr)
j_inv, _, _, _ = analyze("Inverse", inv_corr)
j_cyc, _, _, _ = analyze("Cyclic shift", cyc_corr)


# ═══════════════════════════════════════════════════════════════
#  Type 1: Superposition of contradictory correlations
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TYPE 1: SUPERPOSITION (averaged contradictions)")
print("=" * 70)

# Identity + Inverse (contradictory: low→low vs low→high)
sup_id_inv = 0.5 * id_corr + 0.5 * inv_corr
print(f"\n  Correlation (identity+inverse)/2:")
print(f"    {sup_id_inv}")
analyze("(identity+inverse)/2", sup_id_inv)

# Identity + Cyclic (contradictory: low→low vs low→mid)
sup_id_cyc = 0.5 * id_corr + 0.5 * cyc_corr
print(f"\n  Correlation (identity+cyclic)/2:")
print(f"    {sup_id_cyc}")
analyze("(identity+cyclic)/2", sup_id_cyc)

# Inverse + Cyclic
sup_inv_cyc = 0.5 * inv_corr + 0.5 * cyc_corr
analyze("(inverse+cyclic)/2", sup_inv_cyc)

# All three bijections equally
sup_all = (id_corr + inv_corr + cyc_corr) / 3
print(f"\n  Correlation (id+inv+cyc)/3:")
print(f"    {sup_all}")
analyze("(id+inv+cyc)/3 = uniform", sup_all)

# Weighted: 80% identity + 20% inverse
for w in [0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.1]:
    mix = w * id_corr + (1-w) * inv_corr
    analyze(f"{w:.0%} identity + {1-w:.0%} inverse", mix)


# ═══════════════════════════════════════════════════════════════
#  Type 2: Row conflict (inconsistent function)
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TYPE 2: ROW CONFLICT (different rows say different things)")
print("=" * 70)

# Row 0 says proportional, row 2 says inverse
conflict1 = torch.tensor([
    [1, 0, 0],  # low → low (proportional)
    [0, 1, 0],  # mid → mid (both agree)
    [1, 0, 0],  # high → low (inverse)
], dtype=torch.float32)
analyze("row0=prop, row2=inv, row1=agree", conflict1)

# All rows contradictory
conflict2 = torch.tensor([
    [1, 0, 0],  # low → low (proportional)
    [0, 0, 1],  # mid → high (inverse-ish)
    [0, 1, 0],  # high → mid (cyclic)
], dtype=torch.float32)
analyze("all rows different direction", conflict2)

# Self-contradictory: a row that's uniform
conflict3 = torch.tensor([
    [1, 0, 0],     # low → low (clear)
    [0.33, 0.34, 0.33],  # mid → anything (uncertain)
    [0, 0, 1],     # high → high (clear)
], dtype=torch.float32)
analyze("rows 0,2 clear, row 1 uncertain", conflict3)


# ═══════════════════════════════════════════════════════════════
#  Type 3: Asymmetric correlations
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TYPE 3: ASYMMETRIC CORRELATIONS")
print("(forward relationship ≠ backward relationship)")
print("=" * 70)

# Proportional forward, uniform backward
# This represents: given X, Y is determined. But given Y, X is ambiguous.
asym1 = torch.tensor([
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
], dtype=torch.float32)
# Note: the correlation matrix IS symmetric because it's conditional P(Y|X).
# True asymmetry would need different forward/backward correlations.
# Let's construct a genuinely asymmetric one:
asym2 = torch.tensor([
    [0.8, 0.1, 0.1],  # low X → mostly low Y
    [0.1, 0.8, 0.1],  # mid X → mostly mid Y
    [0.4, 0.3, 0.3],  # high X → anything (ambiguous forward)
], dtype=torch.float32)
analyze("2 clear rows + 1 ambiguous", asym2)

asym3 = torch.tensor([
    [0.9, 0.05, 0.05],
    [0.05, 0.9, 0.05],
    [0.05, 0.05, 0.9],
], dtype=torch.float32)
analyze("Nearly identity (90%)", asym3)

asym4 = torch.tensor([
    [0.7, 0.15, 0.15],
    [0.15, 0.7, 0.15],
    [0.15, 0.15, 0.7],
], dtype=torch.float32)
analyze("Weakened identity (70%)", asym4)

asym5 = torch.tensor([
    [0.5, 0.25, 0.25],
    [0.25, 0.5, 0.25],
    [0.25, 0.25, 0.5],
], dtype=torch.float32)
analyze("Barely identity (50%)", asym5)


# ═══════════════════════════════════════════════════════════════
#  The critical experiment: entangler SEED sensitivity
#  under conflict
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SEED SENSITIVITY UNDER CONFLICT")
print("Does the entangler resolve conflict differently with different seeds?")
print("=" * 70)

# For the conflicting correlation, try many seeds
conflict_corr = 0.5 * id_corr + 0.5 * inv_corr

schmidts = []
rels = []
asymms = []
for seed in range(50):
    ent = Entangler(conflict_corr, seed=seed).build()
    sn = schmidt_number(ent.joint)
    rel, _ = classify_relationship(ent.joint)
    fp = theta_fingerprint(ent.joint)
    asymm = theta_asymmetry(fp)
    schmidts.append(sn)
    rels.append(rel)
    asymms.append(asymm)

from collections import Counter
rel_counts = Counter(rels)
avg_s = sum(schmidts) / len(schmidts)
std_s = (sum((s - avg_s)**2 for s in schmidts) / len(schmidts)) ** 0.5
avg_a = sum(asymms) / len(asymms)
std_a = (sum((a - avg_a)**2 for a in asymms) / len(asymms)) ** 0.5

print(f"\n  (identity+inverse)/2 over 50 seeds:")
print(f"    Schmidt: {avg_s:.3f} ± {std_s:.3f} [{min(schmidts):.3f}, {max(schmidts):.3f}]")
print(f"    Asymmetry: {avg_a:.4f} ± {std_a:.4f}")
print(f"    Relationship classifications: {dict(rel_counts)}")

# Compare: pure identity over 50 seeds
pure_schmidts = []
for seed in range(50):
    ent = Entangler(id_corr, seed=seed).build()
    pure_schmidts.append(schmidt_number(ent.joint))
avg_pure = sum(pure_schmidts) / len(pure_schmidts)
std_pure = (sum((s - avg_pure)**2 for s in pure_schmidts) / len(pure_schmidts)) ** 0.5
print(f"\n  Pure identity over 50 seeds:")
print(f"    Schmidt: {avg_pure:.3f} ± {std_pure:.3f} [{min(pure_schmidts):.3f}, {max(pure_schmidts):.3f}]")

# Seed sensitivity ratio
print(f"\n  Seed sensitivity: conflict σ/pure σ = {std_s/std_pure:.2f}x")
if std_s > 2 * std_pure:
    print(f"  → Conflict AMPLIFIES seed sensitivity (entangler is unstable)")
elif std_s < 0.5 * std_pure:
    print(f"  → Conflict REDUCES seed sensitivity (entangler converges)")
else:
    print(f"  → Similar sensitivity (conflict doesn't change stability)")


# ═══════════════════════════════════════════════════════════════
#  The DS resolution: does iterated combination resolve conflict?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("DS RESOLUTION OF CONFLICT")
print("What happens when we DS-combine the identity and inverse crystals?")
print("=" * 70)

from solver.algebra import ds_combine, born_probabilities, born_certainty

# Build pure crystals
j_id_flat = j_id.reshape(16)
j_inv_flat = j_inv.reshape(16)

# DS combine
combined, K = ds_combine(j_id_flat, j_inv_flat)
sn_comb = schmidt_number(combined.reshape(4, 4))
rel_comb, _ = classify_relationship(combined.reshape(4, 4))
fp_comb = theta_fingerprint(combined.reshape(4, 4))
asymm_comb = theta_asymmetry(fp_comb)

print(f"\n  DS(identity, inverse):")
print(f"    Conflict K = {K.abs().item():.4f}")
print(f"    Schmidt = {sn_comb:.3f}")
print(f"    Relationship: {rel_comb}")
print(f"    θ-asymmetry: {asymm_comb:.4f}")

# Compare with composing
composed = torch.einsum("ab,bc->ac", j_id, j_inv)
re_sum = composed.reshape(-1).real.sum()
if abs(re_sum) > 1e-8:
    composed = composed / re_sum
sn_composed = schmidt_number(composed)
rel_composed, _ = classify_relationship(composed)

print(f"\n  Composition(identity, inverse):")
print(f"    Schmidt = {sn_composed:.3f}")
print(f"    Relationship: {rel_composed}")

# DS combination PRESERVES both. Composition FILTERS.
print(f"\n  DS-combination: Schmidt {sn_comb:.3f} (preserves signal)")
print(f"  Composition:    Schmidt {sn_composed:.3f} (structural filter)")


# ═══════════════════════════════════════════════════════════════
#  Gradual conflict: how much contradiction before collapse?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("CONFLICT TOLERANCE: how much contradiction before collapse?")
print("Continuously rotate from identity toward inverse")
print("=" * 70)

# Parameterize: at angle θ, row 0 = [cos²θ, 0, sin²θ], row 2 = [sin²θ, 0, cos²θ]
import math

angles = [i * math.pi / 40 for i in range(21)]
print(f"\n  {'θ/π':>6s} {'Schmidt':>8s} {'rel':>13s} {'asymm':>8s} {'conf?':>6s}")
print("  " + "-" * 45)

for theta in angles:
    c, s = math.cos(theta), math.sin(theta)
    corr = torch.tensor([
        [c**2, 0, s**2],
        [0, 1, 0],
        [s**2, 0, c**2],
    ], dtype=torch.float32)
    ent = Entangler(corr).build()
    sn = schmidt_number(ent.joint)
    rel, _ = classify_relationship(ent.joint)
    fp = theta_fingerprint(ent.joint)
    asymm = theta_asymmetry(fp)

    # Is there genuine row conflict?
    # Row 0 and row 2 agree when θ=0 (both say proportional) or θ=π/2 (both say inverse)
    # Maximum conflict at θ=π/4
    conflict_level = min(c**2, s**2) / max(c**2, s**2) if max(c**2, s**2) > 0 else 0

    print(f"  {theta/math.pi:>6.3f} {sn:>8.3f} {rel:>13s} {asymm:>8.4f} {conflict_level:>6.2f}")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
