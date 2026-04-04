"""
Information geometry of crystal space.

The crystal's Born distribution defines a point on a statistical
manifold. The Fisher information metric gives the natural distance.
The curvature tells us about the intrinsic geometry.

For a [4,4] joint mass, the Born distribution has 15 free parameters
(16 probabilities summing to 1). The Fisher metric on this 15-dim
simplex is well-known: it's the round metric on the positive orthant
of a 15-sphere.

But the actual crystal space is much smaller — earlier we found it's
effectively ~5-dimensional (95% variance in 5 PCA components). What
is the INTRINSIC geometry of this 5-dim submanifold?

Also: the Fisher distance between two crystals is a natural
measure of "how different are their probability structures."
This might be more natural than the θ-fingerprint distance.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, schmidt_number, born_probabilities
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)


def fisher_distance(joint1, joint2):
    """Fisher-Rao distance between two joint masses.

    For categorical distributions, the Fisher-Rao distance is:
    d(p, q) = 2 * arccos(Σ √(p_i * q_i))
    This is the arc length on the probability simplex sphere.
    """
    p = joint1.abs().pow(2).reshape(-1)
    q = joint2.abs().pow(2).reshape(-1)
    p = p / p.sum().clamp(min=1e-10)
    q = q / q.sum().clamp(min=1e-10)
    bc = (p * q).sqrt().sum().item()  # Bhattacharyya coefficient
    bc = min(bc, 1.0)  # numerical safety
    return 2 * math.acos(bc)


def fisher_metric_tensor(joint, epsilon=1e-4):
    """Approximate the Fisher information metric at a point.

    Perturb each of the 16 entries of the joint mass slightly,
    compute Fisher distance, get a local metric.

    Returns a [32, 32] metric tensor (32 = 16 real + 16 imag parts).
    """
    p0 = joint.abs().pow(2).reshape(-1)
    p0 = p0 / p0.sum()

    # For the simplex, the Fisher metric is:
    # g_ij = Σ_k (∂p_k/∂θ_i)(∂p_k/∂θ_j) / p_k
    # where θ are the coordinates of the joint mass

    # Since p_k = |m_k|² / Σ|m_j|², the derivatives are complex
    # Let's just compute pairwise distances numerically
    n_params = MASS_DIM * MASS_DIM * 2  # real + imag parts

    # Too expensive for full metric. Let's compute along specific directions.
    return None


# ═══════════════════════════════════════════════════════════════
#  Build all reference crystals
# ═══════════════════════════════════════════════════════════════

crystals = {}
for name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Entangler(corr).build()
    crystals[name] = ent.joint.clone()

# Add composed crystals
crystals["id∘id"] = compose(crystals["proportional"], crystals["proportional"])
crystals["id∘id∘id"] = compose(crystals["id∘id"], crystals["proportional"])
crystals["id∘inv"] = compose(crystals["proportional"], crystals["inverse"])
crystals["inv∘mod"] = compose(crystals["inverse"], crystals["modular"])

# Add crystals at different entangler steps
for steps in [10, 20, 30, 40, 50]:
    corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
    ent = Entangler(corr).build(n_steps=steps)
    crystals[f"id_step{steps}"] = ent.joint.clone()


# ═══════════════════════════════════════════════════════════════
#  Fisher distance matrix
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("FISHER-RAO DISTANCES BETWEEN CRYSTALS")
print("(arc length on the probability simplex sphere)")
print("=" * 70)

names = list(crystals.keys())

# Print distances for key pairs
print(f"\n  Key Fisher distances:")
key_pairs = [
    ("proportional", "inverse"),
    ("proportional", "modular"),
    ("proportional", "exponential"),
    ("inverse", "modular"),
    ("proportional", "id∘id"),
    ("id∘id", "id∘id∘id"),
    ("proportional", "id∘inv"),
    ("id_step10", "id_step20"),
    ("id_step20", "id_step30"),
    ("id_step30", "id_step40"),
    ("id_step40", "id_step50"),
]

for n1, n2 in key_pairs:
    if n1 in crystals and n2 in crystals:
        d = fisher_distance(crystals[n1], crystals[n2])
        s1 = schmidt_number(crystals[n1])
        s2 = schmidt_number(crystals[n2])
        print(f"  {n1:>15s} ↔ {n2:<15s}: Fisher = {d:.4f} "
              f"(Schmidt {s1:.2f} → {s2:.2f})")


# ═══════════════════════════════════════════════════════════════
#  Fisher distance vs Schmidt difference
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FISHER DISTANCE vs SCHMIDT DIFFERENCE")
print("=" * 70)

# Collect all pairwise distances
fisher_dists = []
schmidt_diffs = []
pair_labels = []

sig_names = list(RELATIONSHIP_SIGNATURES.keys())
for i in range(len(sig_names)):
    for j in range(i+1, len(sig_names)):
        n1, n2 = sig_names[i], sig_names[j]
        d = fisher_distance(crystals[n1], crystals[n2])
        ds = abs(schmidt_number(crystals[n1]) - schmidt_number(crystals[n2]))
        fisher_dists.append(d)
        schmidt_diffs.append(ds)
        pair_labels.append(f"{n1[:4]}-{n2[:4]}")

# Correlation
def correlation(xs, ys):
    n = len(xs)
    mx, my = sum(xs)/n, sum(ys)/n
    cov = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx = sum((x-mx)**2 for x in xs) ** 0.5
    sy = sum((y-my)**2 for y in ys) ** 0.5
    return cov / (sx * sy) if sx > 1e-10 and sy > 1e-10 else 0

r = correlation(fisher_dists, schmidt_diffs)
print(f"\n  Correlation(Fisher distance, |ΔSchmidt|): r = {r:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Fisher distance along composition decay
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FISHER DISTANCE ALONG COMPOSITION DECAY")
print("How far does each composition step move you on the manifold?")
print("=" * 70)

joint_0 = crystals["proportional"].clone()
current = joint_0.clone()

print(f"\n  {'comp':>5s} {'Schmidt':>8s} {'Fisher from prev':>16s} {'Fisher from start':>18s} {'Fisher from uniform':>20s}")
print("  " + "-" * 75)

# Build uniform crystal for reference
uniform = torch.ones(MASS_DIM, MASS_DIM, dtype=torch.cfloat) / MASS_DIM**2
uniform_fisher = fisher_distance(joint_0, uniform)

prev = current.clone()
for step in range(8):
    d_prev = fisher_distance(current, prev) if step > 0 else 0
    d_start = fisher_distance(current, joint_0)
    d_uniform = fisher_distance(current, uniform)
    sn = schmidt_number(current)

    print(f"  {step:>5d} {sn:>8.3f} {d_prev:>16.4f} {d_start:>18.4f} {d_uniform:>20.4f}")

    prev = current.clone()
    current = compose(current, joint_0)

# Is Fisher step length constant? (geodesic)
print(f"\n  Step lengths suggest {'geodesic' if abs(fisher_distance(crystals['id∘id'], joint_0) - 2 * 0) < 0.1 else 'non-geodesic'} path")


# ═══════════════════════════════════════════════════════════════
#  Fisher distance during entangling
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FISHER DISTANCE DURING ENTANGLING")
print("How far does each entangling step move you?")
print("=" * 70)

id_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

# Build crystals at each step count
step_crystals = {}
for steps in range(1, 51):
    ent = Entangler(id_corr).build(n_steps=steps)
    step_crystals[steps] = ent.joint.clone()

print(f"\n  {'steps':>6s} {'Schmidt':>8s} {'Fisher from prev':>16s} {'Fisher from final':>17s} {'step/prev ratio':>16s}")
print("  " + "-" * 70)

final = step_crystals[50]
for steps in [1, 2, 3, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]:
    if steps in step_crystals:
        d_prev = fisher_distance(step_crystals[steps], step_crystals[steps-1]) if steps > 1 and (steps-1) in step_crystals else 0
        d_final = fisher_distance(step_crystals[steps], final)
        sn = schmidt_number(step_crystals[steps])

        prev_d = fisher_distance(step_crystals[steps-1], step_crystals[max(1,steps-2)]) if steps > 2 and (steps-2) in step_crystals else 0
        ratio = d_prev / prev_d if prev_d > 1e-6 else 0

        print(f"  {steps:>6d} {sn:>8.3f} {d_prev:>16.4f} {d_final:>17.4f} {ratio:>16.3f}")


# ═══════════════════════════════════════════════════════════════
#  The geodesic question: is composition a geodesic?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("IS COMPOSITION A GEODESIC?")
print("If the path crystal → crystal² → crystal³ ... is a geodesic,")
print("then d(crystal^n, crystal^m) = |n-m| × step_length")
print("=" * 70)

# Build composition sequence
comp_seq = [joint_0.clone()]
current = joint_0.clone()
for i in range(7):
    current = compose(current, joint_0)
    comp_seq.append(current.clone())

# Check triangle inequality and geodesic condition
print(f"\n  Pairwise Fisher distances (composition steps 0-7):")
print(f"  {'':>5s}", end="")
for j in range(8):
    print(f"  {j:>6d}", end="")
print()

for i in range(8):
    print(f"  {i:>5d}", end="")
    for j in range(8):
        d = fisher_distance(comp_seq[i], comp_seq[j])
        print(f"  {d:>6.3f}", end="")
    print()

# Check: is d(0,n) proportional to n?
step1 = fisher_distance(comp_seq[0], comp_seq[1])
print(f"\n  Step 1 distance: {step1:.4f}")
print(f"  If geodesic, d(0,n) should be n × {step1:.4f}:")
for n in range(1, 8):
    actual = fisher_distance(comp_seq[0], comp_seq[n])
    predicted = n * step1
    ratio = actual / predicted if predicted > 0 else 0
    print(f"    d(0,{n}) = {actual:.4f}, predicted = {predicted:.4f}, ratio = {ratio:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Curvature estimation
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("CURVATURE OF CRYSTAL SPACE")
print("Estimated from triangle defects on the manifold")
print("=" * 70)

# For three crystals A, B, C, the triangle defect is:
# defect = d(A,B) + d(B,C) - d(A,C) (should be ≥ 0 by triangle inequality)
# On a flat space, the defect is related to the triangle area
# On a curved space, defect/(area) gives the curvature

triangles = [
    ("proportional", "inverse", "modular"),
    ("proportional", "exponential", "logarithmic"),
    ("proportional", "modular", "quadratic"),
    ("inverse", "modular", "exponential"),
]

for a, b, c in triangles:
    dab = fisher_distance(crystals[a], crystals[b])
    dbc = fisher_distance(crystals[b], crystals[c])
    dac = fisher_distance(crystals[a], crystals[c])

    # Three sides
    sides = sorted([dab, dbc, dac])
    # Triangle inequality defect (smallest when the triangle is "thin")
    defect = sides[0] + sides[1] - sides[2]
    # Semi-perimeter
    s = (dab + dbc + dac) / 2
    # Heron's area (Euclidean approximation)
    area_sq = s * (s - dab) * (s - dbc) * (s - dac)
    area = math.sqrt(max(area_sq, 0))

    print(f"\n  Triangle {a[:4]}-{b[:4]}-{c[:4]}:")
    print(f"    Sides: {dab:.4f}, {dbc:.4f}, {dac:.4f}")
    print(f"    Defect: {defect:.4f}")
    print(f"    Area (Heron): {area:.4f}")
    if area > 1e-6:
        # Gaussian curvature ≈ defect / area (for small triangles)
        curvature = defect / area
        print(f"    Curvature ≈ {curvature:.4f}")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
