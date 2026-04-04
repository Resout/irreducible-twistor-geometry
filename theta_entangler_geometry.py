"""
The entangler as a map from correlation space to θ-space.

The landscape is a single hill centered on the identity/proportional
correlation. The entangler maps [3,3] correlation matrices to [4,4]
joint masses, which have θ-fingerprints. Is this map linear or curved?

If linear: the landscape shape is determined by the correlation geometry.
If curved: the entangler's nonlinearity creates structure.

Also: the identity-detector insight. The crystal framework works by
measuring how close a function is to the identity. Can we formalize
this? The correlation matrix of the identity IS the [3,3] identity.
Every other function deviates from this. The degree of deviation
determines the Schmidt number.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, schmidt_number, born_probabilities
from solver.crystals import Entangler

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


def theta_distance(fp1, fp2):
    return (fp1 - fp2).abs().pow(2).sum().sqrt().item()


def theta_asymmetry(fp):
    return (fp[:3].abs() - fp[3:].abs()).abs().sum().item()


# ═══════════════════════════════════════════════════════════════
#  Canonical correlation matrices
# ═══════════════════════════════════════════════════════════════

identity_corr = torch.tensor([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=torch.float32)
modular_corr = torch.tensor([[.4, .3, .3], [.3, .4, .3], [.3, .3, .4]], dtype=torch.float32)
inverse_corr = torch.tensor([[0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=torch.float32)
uniform_corr = torch.ones(3, 3, dtype=torch.float32) / 3
exp_corr = torch.tensor([[1, 0, 0], [.6, .4, 0], [0, 0, 1]], dtype=torch.float32)

# ═══════════════════════════════════════════════════════════════
#  Test 1: Interpolation between identity and modular
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 1: INTERPOLATION identity → modular")
print("Does the θ-path curve or stay straight?")
print("=" * 70)

# Build canonical fingerprint
ent_id = Entangler(identity_corr).build()
fp_id = theta_fingerprint(ent_id.joint)

ent_mod = Entangler(modular_corr).build()
fp_mod = theta_fingerprint(ent_mod.joint)

total_d = theta_distance(fp_id, fp_mod)
print(f"\n  Total θ-distance (identity → modular): {total_d:.4f}")

n_steps = 20
print(f"\n  {'t':>6s} {'Schmidt':>8s} {'d(id)':>8s} {'d(mod)':>8s} {'d_linear':>8s} {'deviation':>8s} {'asymm':>8s}")
print("  " + "-" * 60)

path_points = []
for i in range(n_steps + 1):
    t = i / n_steps
    corr = (1 - t) * identity_corr + t * modular_corr
    ent = Entangler(corr).build()
    fp = theta_fingerprint(ent.joint)
    sn = schmidt_number(ent.joint)

    d_id = theta_distance(fp, fp_id)
    d_mod = theta_distance(fp, fp_mod)

    # Expected distance if path were straight line in θ-space
    d_linear = t * total_d
    deviation = abs(d_id - d_linear)

    path_points.append({
        "t": t, "schmidt": sn, "d_id": d_id, "d_mod": d_mod,
        "d_linear": d_linear, "deviation": deviation,
        "fp": fp, "asymm": theta_asymmetry(fp),
    })

    if i % 2 == 0 or i == n_steps:
        print(f"  {t:>6.2f} {sn:>8.3f} {d_id:>8.4f} {d_mod:>8.4f} "
              f"{d_linear:>8.4f} {deviation:>8.4f} {path_points[-1]['asymm']:>8.4f}")

# Curvature measure: max deviation from straight line / total distance
max_dev = max(p["deviation"] for p in path_points)
print(f"\n  Max deviation from straight line: {max_dev:.4f} ({max_dev/total_d*100:.1f}% of total)")

# Step length variation
step_lengths = []
for i in range(1, len(path_points)):
    d = theta_distance(path_points[i]["fp"], path_points[i-1]["fp"])
    step_lengths.append(d)
avg_step = sum(step_lengths) / len(step_lengths)
step_cv = (sum((s - avg_step)**2 for s in step_lengths) / len(step_lengths)) ** 0.5 / avg_step if avg_step > 0 else 0
print(f"  Step length CV (0 = geodesic): {step_cv:.3f}")
min_step = min(step_lengths) if step_lengths else 0
max_step = max(step_lengths) if step_lengths else 0
ratio_step = max_step / min_step if min_step > 1e-12 else float('inf')
print(f"  Step lengths: min={min_step:.4f}, max={max_step:.4f}, ratio={ratio_step:.2f}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: Interpolation identity → inverse
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 2: INTERPOLATION identity → inverse")
print("(Both are bijections — should this be flatter?)")
print("=" * 70)

ent_inv = Entangler(inverse_corr).build()
fp_inv = theta_fingerprint(ent_inv.joint)
total_d2 = theta_distance(fp_id, fp_inv)
print(f"\n  Total θ-distance (identity → inverse): {total_d2:.4f}")

path2 = []
for i in range(n_steps + 1):
    t = i / n_steps
    corr = (1 - t) * identity_corr + t * inverse_corr
    ent = Entangler(corr).build()
    fp = theta_fingerprint(ent.joint)
    sn = schmidt_number(ent.joint)
    d_id = theta_distance(fp, fp_id)
    d_linear = t * total_d2
    path2.append({"t": t, "schmidt": sn, "d_id": d_id,
                   "deviation": abs(d_id - d_linear), "fp": fp,
                   "asymm": theta_asymmetry(fp)})
    if i % 4 == 0 or i == n_steps:
        print(f"  t={t:.2f}: Schmidt={sn:.3f}, d(id)={d_id:.4f}, "
              f"dev={abs(d_id - d_linear):.4f}, asymm={theta_asymmetry(fp):.4f}")

max_dev2 = max(p["deviation"] for p in path2)
print(f"\n  Max deviation from straight line: {max_dev2:.4f} ({max_dev2/total_d2*100:.1f}% of total)")

# Schmidt range along the path
schmidts = [p["schmidt"] for p in path2]
print(f"  Schmidt range: [{min(schmidts):.3f}, {max(schmidts):.3f}]")
print(f"  Schmidt variation: {max(schmidts) - min(schmidts):.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 3: Interpolation identity → uniform (total noise)
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 3: INTERPOLATION identity → uniform")
print("(Moving toward total noise — how does Schmidt decay?)")
print("=" * 70)

ent_uni = Entangler(uniform_corr).build()
fp_uni = theta_fingerprint(ent_uni.joint)
total_d3 = theta_distance(fp_id, fp_uni)
print(f"\n  Total θ-distance (identity → uniform): {total_d3:.4f}")

path3 = []
for i in range(n_steps + 1):
    t = i / n_steps
    corr = (1 - t) * identity_corr + t * uniform_corr
    ent = Entangler(corr).build()
    fp = theta_fingerprint(ent.joint)
    sn = schmidt_number(ent.joint)
    d_id = theta_distance(fp, fp_id)
    d_linear = t * total_d3
    path3.append({"t": t, "schmidt": sn, "d_id": d_id,
                   "deviation": abs(d_id - d_linear), "fp": fp})
    if i % 4 == 0 or i == n_steps:
        print(f"  t={t:.2f}: Schmidt={sn:.3f}, d(id)={d_id:.4f}, dev={abs(d_id - d_linear):.4f}")

max_dev3 = max(p["deviation"] for p in path3)
print(f"\n  Max deviation from straight line: {max_dev3:.4f} ({max_dev3/total_d3*100:.1f}% of total)")

# Is Schmidt a monotone function of distance from identity?
schmidts3 = [p["schmidt"] for p in path3]
dists3 = [p["d_id"] for p in path3]
# Check monotonicity
monotone = all(schmidts3[i] >= schmidts3[i+1] for i in range(len(schmidts3)-1))
print(f"  Schmidt monotonically decreasing with distance? {monotone}")

# Correlation
def correlation(xs, ys):
    n = len(xs)
    mx, my = sum(xs)/n, sum(ys)/n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = sum((x - mx)**2 for x in xs) ** 0.5
    sy = sum((y - my)**2 for y in ys) ** 0.5
    return cov / (sx * sy) if sx > 1e-10 and sy > 1e-10 else 0

r_sd = correlation(schmidts3, dists3)
print(f"  Correlation(Schmidt, d(identity)): r = {r_sd:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 4: The correlation matrix as distance from identity
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 4: CORRELATION-SPACE DISTANCE vs θ-SPACE DISTANCE")
print("=" * 70)

# For each correlation matrix, compute:
#   corr_dist = Frobenius distance from identity
#   theta_dist = θ-distance from identity crystal
# If these correlate, the entangler preserves the metric.

from solver.crystals import RELATIONSHIP_SIGNATURES

test_corrs = []
for name, sig in RELATIONSHIP_SIGNATURES.items():
    row_sums = sig.sum(dim=1, keepdim=True)
    corr = sig / row_sums.clamp(min=1e-6)
    test_corrs.append((name, corr))

# Add interpolated points
for t in [0.1, 0.2, 0.3, 0.5, 0.7, 0.9]:
    corr = (1 - t) * identity_corr + t * modular_corr
    test_corrs.append((f"id→mod t={t}", corr))
    corr2 = (1 - t) * identity_corr + t * inverse_corr
    test_corrs.append((f"id→inv t={t}", corr2))

corr_dists = []
theta_dists = []

print(f"\n  {'Name':<25s} {'corr_dist':>10s} {'θ_dist':>10s} {'Schmidt':>8s}")
print("  " + "-" * 55)

for name, corr in test_corrs:
    # Correlation-space distance from identity
    c_dist = (corr - identity_corr).pow(2).sum().sqrt().item()

    # Build crystal, get θ-distance
    ent = Entangler(corr).build()
    fp = theta_fingerprint(ent.joint)
    t_dist = theta_distance(fp, fp_id)
    sn = schmidt_number(ent.joint)

    corr_dists.append(c_dist)
    theta_dists.append(t_dist)

    print(f"  {name:<25s} {c_dist:>10.4f} {t_dist:>10.4f} {sn:>8.3f}")

r_ct = correlation(corr_dists, theta_dists)
print(f"\n  Correlation(corr_dist, θ_dist): r = {r_ct:.3f}")
if abs(r_ct) > 0.9:
    print(f"  → The entangler approximately PRESERVES the metric")
elif abs(r_ct) > 0.7:
    print(f"  → The entangler DISTORTS but roughly preserves distances")
else:
    print(f"  → The entangler strongly distorts distances")


# ═══════════════════════════════════════════════════════════════
#  Test 5: The identity as the fixed point
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 5: THE IDENTITY AS FIXED POINT OF COMPOSITION")
print("Compose any crystal with the identity crystal: what happens?")
print("=" * 70)

from solver.crystals import compose

ent_id = Entangler(identity_corr).build()
id_joint = ent_id.joint

# Compose identity with various crystals
test_joints = [
    ("inverse", Entangler(inverse_corr).build().joint),
    ("modular", Entangler(modular_corr).build().joint),
    ("exponential", Entangler(exp_corr).build().joint),
    ("identity", id_joint),
]

for name, joint in test_joints:
    composed = compose(id_joint, joint)
    sn_orig = schmidt_number(joint)
    sn_comp = schmidt_number(composed)
    fp_orig = theta_fingerprint(joint)
    fp_comp = theta_fingerprint(composed)
    d = theta_distance(fp_orig, fp_comp)

    # Also: compose the other way
    composed_rev = compose(joint, id_joint)
    sn_rev = schmidt_number(composed_rev)
    fp_rev = theta_fingerprint(composed_rev)
    d_rev = theta_distance(fp_orig, fp_rev)

    print(f"\n  identity ∘ {name}:")
    print(f"    Original Schmidt: {sn_orig:.3f}")
    print(f"    id∘X Schmidt:     {sn_comp:.3f} (θ-shift: {d:.4f})")
    print(f"    X∘id Schmidt:     {sn_rev:.3f} (θ-shift: {d_rev:.4f})")

    if abs(sn_orig - sn_comp) < 0.05 and d < 0.05:
        print(f"    → Identity preserved (approximately neutral)")
    else:
        print(f"    → Schmidt changed by {sn_comp - sn_orig:+.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 6: What makes the identity special? Eigenstructure.
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 6: EIGENSTRUCTURE OF THE IDENTITY CRYSTAL")
print("=" * 70)

# SVD of the identity crystal joint mass
U, S, V = torch.linalg.svd(id_joint)
print(f"\n  Singular values of identity crystal:")
for i, s in enumerate(S):
    print(f"    σ_{i} = {s.abs().item():.6f}")

print(f"\n  Schmidt number: {1.0 / (S.abs()**2 / S.abs().sum()**2).sum().item():.3f}")

# Compare with other crystals
for name, joint in test_joints:
    _, S2, _ = torch.linalg.svd(joint)
    sv_ratio = S2[0].abs().item() / S2[-1].abs().item() if S2[-1].abs() > 1e-10 else float('inf')
    print(f"  {name:15s}: σ ratio = {sv_ratio:.2f}, "
          f"σ = [{', '.join(f'{s.abs().item():.4f}' for s in S2)}]")


# ═══════════════════════════════════════════════════════════════
#  Test 7: Schmidt as function of correlation deviation from identity
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 7: SCHMIDT = f(||correlation - identity||)")
print("What is the functional form? Linear, quadratic, exponential?")
print("=" * 70)

# Fine-grained interpolation identity → uniform
ts = [i / 100 for i in range(101)]
data = []
for t in ts:
    corr = (1 - t) * identity_corr + t * uniform_corr
    ent = Entangler(corr).build()
    sn = schmidt_number(ent.joint)
    c_dist = t * (uniform_corr - identity_corr).pow(2).sum().sqrt().item()
    data.append({"t": t, "schmidt": sn, "corr_dist": c_dist})

# Print at key points
print(f"\n  {'t':>6s} {'Schmidt':>8s} {'corr_dist':>10s}")
for d in data[::10]:
    print(f"  {d['t']:>6.2f} {d['schmidt']:>8.3f} {d['corr_dist']:>10.4f}")

# Find where Schmidt drops below noise floor (~2.0)
for d in data:
    if d["schmidt"] < 2.0:
        print(f"\n  Schmidt drops below 2.0 at t = {d['t']:.2f} (corr_dist = {d['corr_dist']:.4f})")
        break

# Find where Schmidt drops below threshold (1.5)
for d in data:
    if d["schmidt"] < 1.5:
        print(f"  Schmidt drops below 1.5 at t = {d['t']:.2f} (corr_dist = {d['corr_dist']:.4f})")
        break

# Check: is the decay convex, linear, or concave?
# Compare midpoint Schmidt with average of endpoints
s_0 = data[0]["schmidt"]
s_100 = data[-1]["schmidt"]
s_50 = data[50]["schmidt"]
midpoint_avg = (s_0 + s_100) / 2
print(f"\n  Schmidt at t=0:    {s_0:.3f}")
print(f"  Schmidt at t=0.5:  {s_50:.3f}")
print(f"  Schmidt at t=1:    {s_100:.3f}")
print(f"  Midpoint average:  {midpoint_avg:.3f}")
if s_50 > midpoint_avg:
    print(f"  → CONCAVE decay (resists degradation in the middle)")
elif s_50 < midpoint_avg:
    print(f"  → CONVEX decay (accelerating degradation)")
else:
    print(f"  → LINEAR decay")


# ═══════════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
  1. Interpolation id→mod: max deviation {max_dev:.4f} ({max_dev/total_d*100:.1f}% of total)
     → The entangler maps {'nearly straight' if max_dev/total_d < 0.15 else 'curved'} paths in correlation space to {'nearly straight' if max_dev/total_d < 0.15 else 'curved'} paths in θ-space

  2. Interpolation id→inv: max deviation {max_dev2:.4f} ({max_dev2/total_d2*100:.1f}% of total)
     Schmidt variation along path: {max(schmidts) - min(schmidts):.3f}
     → The bijection ridge is {'flat' if max(schmidts) - min(schmidts) < 0.2 else 'not flat'}

  3. Correlation between corr_dist and θ_dist: r = {r_ct:.3f}
     → The entangler {'preserves' if abs(r_ct) > 0.8 else 'distorts'} the metric

  4. Schmidt decays {'concavely' if s_50 > midpoint_avg else 'convexly'} from identity toward uniform
""")
