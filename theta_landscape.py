"""
The θ-landscape: crystal space has attractors.

The canonical θ-fingerprint (from period-3 resonant crystals) is a
distinguished point. Every crystal's Schmidt number correlates with
its θ-distance from this point. This suggests a LANDSCAPE with peaks
at resonant configurations and valleys elsewhere.

Questions:
1. Is there a second attractor for period-2 resonance?
2. What is the shape of the landscape (spherical shells vs ridges)?
3. Does the landscape predict the computational class of a function?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
from solver.algebra import H, MASS_DIM, schmidt_number
from solver.crystals import Entangler, classify_relationship
from solver.functional import build_function_crystal, _detect_output_period
from solver import compute

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


domain = range(1, 301)


# ═══════════════════════════════════════════════════════════════
#  Build the attractor fingerprints
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("THE θ-LANDSCAPE: ATTRACTORS AND TOPOLOGY")
print("=" * 70)

# Period-3 attractor (from identity, cleanest resonance)
fc_p3 = build_function_crystal("p3", lambda x: int(x) % 3, domain)
fp_p3 = theta_fingerprint(fc_p3.joint)

# Period-2 attractor
fc_p2 = build_function_crystal("p2", lambda x: int(x) % 2, domain)
fp_p2 = theta_fingerprint(fc_p2.joint)

# Constant function (degenerate)
fc_const = build_function_crystal("const", lambda x: 0, domain)
fp_const = theta_fingerprint(fc_const.joint)

# The relationship signature crystals
from solver.crystals import RELATIONSHIP_SIGNATURES, Entangler as Ent
attractors = {
    "period-3 (resonant)": (fc_p3, fp_p3, fc_p3.schmidt),
    "period-2 (partial)": (fc_p2, fp_p2, fc_p2.schmidt),
    "constant": (fc_const, fp_const, fc_const.schmidt),
}

# Add signature crystals
for rel_name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Ent(corr).build()
    fp = theta_fingerprint(ent.joint)
    sn = schmidt_number(ent.joint)
    attractors[f"sig:{rel_name}"] = (None, fp, sn)

print("\nAttractors:")
print(f"  {'Name':<25s} {'Schmidt':>8s} {'Asymm':>8s}")
print("  " + "-" * 43)
for name, (_, fp, sn) in attractors.items():
    a = theta_asymmetry(fp)
    print(f"  {name:<25s} {sn:>8.3f} {a:>8.4f}")

# Distance matrix between attractors
att_names = list(attractors.keys())
print(f"\nθ-distances between attractors:")
print(f"  {'':>25s}", end="")
for n in att_names[:6]:
    print(f" {n[:8]:>8s}", end="")
print()
for i, n1 in enumerate(att_names[:6]):
    print(f"  {n1:<25s}", end="")
    for j, n2 in enumerate(att_names[:6]):
        d = theta_distance(attractors[n1][1], attractors[n2][1])
        print(f" {d:>8.4f}", end="")
    print()


# ═══════════════════════════════════════════════════════════════
#  Map the landscape: many functions, their position in θ-space
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("LANDSCAPE MAP: each function's position relative to attractors")
print("=" * 70)

landscape_points = []

# Modular exponentiations with various periods
for base, mod in [(2,7), (3,7), (2,13), (3,13), (5,7), (7,19),
                   (2,11), (3,11), (5,11), (2,17), (3,17), (2,19),
                   (5,13), (7,13), (2,31), (3,31)]:
    f = lambda x, b=base, m=mod: pow(b, int(x), m)
    fc = build_function_crystal(f"{base}^x mod {mod}", f, domain)
    outputs = [f(x) for x in range(1, 101)]
    # Detect true period
    for p in range(1, 50):
        if all(outputs[i] == outputs[i % p] for i in range(p, min(len(outputs), p * 5))):
            period = p
            break
    else:
        period = None

    fp = theta_fingerprint(fc.joint)
    landscape_points.append({
        "name": f"{base}^x mod {mod}",
        "family": "mod_exp",
        "period": period,
        "schmidt": fc.schmidt,
        "theta_fp": fp,
        "d_p3": theta_distance(fp, fp_p3),
        "d_p2": theta_distance(fp, fp_p2),
        "d_const": theta_distance(fp, fp_const),
        "asymm": theta_asymmetry(fp),
    })

# Pure modular (x mod n)
for n in range(2, 16):
    f = lambda x, _n=n: int(x) % _n
    fc = build_function_crystal(f"x mod {n}", f, domain)
    fp = theta_fingerprint(fc.joint)
    landscape_points.append({
        "name": f"x mod {n}",
        "family": "periodic",
        "period": n,
        "schmidt": fc.schmidt,
        "theta_fp": fp,
        "d_p3": theta_distance(fp, fp_p3),
        "d_p2": theta_distance(fp, fp_p2),
        "d_const": theta_distance(fp, fp_const),
        "asymm": theta_asymmetry(fp),
    })

# Arithmetic functions
for name, f in [
    ("φ(n)", lambda x: compute.euler_totient(max(1, int(x)))),
    ("d(n)", lambda x: compute.divisor_count(max(1, int(x)))),
    ("σ(n)", lambda x: compute.divisor_sum(max(1, int(x)))),
    ("digit_sum", lambda x: compute.digit_sum(max(1, int(x)))),
]:
    fc = build_function_crystal(name, f, domain)
    fp = theta_fingerprint(fc.joint)
    landscape_points.append({
        "name": name,
        "family": "arithmetic",
        "period": None,
        "schmidt": fc.schmidt,
        "theta_fp": fp,
        "d_p3": theta_distance(fp, fp_p3),
        "d_p2": theta_distance(fp, fp_p2),
        "d_const": theta_distance(fp, fp_const),
        "asymm": theta_asymmetry(fp),
    })

# Sort by Schmidt
landscape_points.sort(key=lambda p: -p["schmidt"])

# Report
print(f"\n{'Function':<20s} {'Period':>6s} {'Schmidt':>8s} {'d(p3)':>8s} {'d(p2)':>8s} {'d(const)':>8s} {'Asymm':>8s}")
print("-" * 80)
for pt in landscape_points:
    p_str = str(pt["period"]) if pt["period"] else "—"
    print(f"  {pt['name']:<18s} {p_str:>6s} {pt['schmidt']:>8.3f} "
          f"{pt['d_p3']:>8.4f} {pt['d_p2']:>8.4f} {pt['d_const']:>8.4f} {pt['asymm']:>8.4f}")


# ═══════════════════════════════════════════════════════════════
#  Regression: Schmidt as a function of attractor distances
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("REGRESSION: Schmidt = f(θ-distance to attractors)")
print("=" * 70)

schmidts = [p["schmidt"] for p in landscape_points]
d_p3 = [p["d_p3"] for p in landscape_points]
d_p2 = [p["d_p2"] for p in landscape_points]
d_const = [p["d_const"] for p in landscape_points]
asymms = [p["asymm"] for p in landscape_points]

def correlation(xs, ys):
    n = len(xs)
    mx, my = sum(xs)/n, sum(ys)/n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = sum((x - mx)**2 for x in xs) ** 0.5
    sy = sum((y - my)**2 for y in ys) ** 0.5
    return cov / (sx * sy) if sx > 1e-10 and sy > 1e-10 else 0

r_p3 = correlation(schmidts, d_p3)
r_p2 = correlation(schmidts, d_p2)
r_const = correlation(schmidts, d_const)
r_asymm = correlation(schmidts, asymms)

print(f"\n  Correlation(Schmidt, d(p3 attractor)):   r = {r_p3:+.3f}")
print(f"  Correlation(Schmidt, d(p2 attractor)):   r = {r_p2:+.3f}")
print(f"  Correlation(Schmidt, d(constant)):       r = {r_const:+.3f}")
print(f"  Correlation(Schmidt, θ-asymmetry):       r = {r_asymm:+.3f}")

# Check: is d_p3 the best predictor? Or is some combination better?
# Simple: which single attractor distance best predicts Schmidt?
print(f"\n  Best single predictor of Schmidt: ", end="")
rs = [("d(p3)", abs(r_p3)), ("d(p2)", abs(r_p2)),
      ("d(const)", abs(r_const)), ("asymm", abs(r_asymm))]
rs.sort(key=lambda x: -x[1])
for name, r in rs:
    print(f"{name} (|r|={r:.3f}), ", end="")
print()


# ═══════════════════════════════════════════════════════════════
#  Classify by nearest attractor
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("NEAREST ATTRACTOR CLASSIFICATION")
print("=" * 70)

# For each point, find nearest attractor
for pt in landscape_points:
    dists = {
        "period-3": pt["d_p3"],
        "period-2": pt["d_p2"],
        "constant": pt["d_const"],
    }
    nearest = min(dists, key=dists.get)
    pt["nearest"] = nearest

# Count
from collections import Counter
counts = Counter(pt["nearest"] for pt in landscape_points)
print(f"\n  Classification by nearest attractor:")
for att, count in counts.most_common():
    pts = [p for p in landscape_points if p["nearest"] == att]
    avg_schmidt = sum(p["schmidt"] for p in pts) / len(pts)
    periods = [p["period"] for p in pts if p["period"]]
    p3_frac = sum(1 for p in periods if p % 3 == 0) / len(periods) if periods else 0
    print(f"    {att:<12s}: n={count:2d}, avg Schmidt={avg_schmidt:.3f}, "
          f"3|period: {p3_frac*100:.0f}%")


# ═══════════════════════════════════════════════════════════════
#  The landscape shape: are shells spherical?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("LANDSCAPE SHAPE: CONSTANT-SCHMIDT SHELLS")
print("=" * 70)

# Group points by Schmidt level
schmidt_groups = {
    "resonant (>3.0)": [p for p in landscape_points if p["schmidt"] > 3.0],
    "partial (2.5-3.0)": [p for p in landscape_points if 2.5 <= p["schmidt"] <= 3.0],
    "weak (1.5-2.5)": [p for p in landscape_points if 1.5 <= p["schmidt"] < 2.5],
    "noise (<1.5)": [p for p in landscape_points if p["schmidt"] < 1.5],
}

for group_name, pts in schmidt_groups.items():
    if not pts:
        continue
    d_p3_vals = [p["d_p3"] for p in pts]
    d_p2_vals = [p["d_p2"] for p in pts]
    avg_d_p3 = sum(d_p3_vals) / len(d_p3_vals)
    avg_d_p2 = sum(d_p2_vals) / len(d_p2_vals)
    spread_d_p3 = (sum((d - avg_d_p3)**2 for d in d_p3_vals) / len(d_p3_vals)) ** 0.5
    spread_d_p2 = (sum((d - avg_d_p2)**2 for d in d_p2_vals) / len(d_p2_vals)) ** 0.5
    print(f"\n  {group_name} (n={len(pts)}):")
    print(f"    d(p3): {avg_d_p3:.4f} ± {spread_d_p3:.4f}")
    print(f"    d(p2): {avg_d_p2:.4f} ± {spread_d_p2:.4f}")

    # Are these shells spherical? Check if d_p3 and d_p2 are correlated
    if len(pts) > 2:
        r = correlation(d_p3_vals, d_p2_vals)
        print(f"    Correlation(d_p3, d_p2): r = {r:+.3f}", end="")
        if abs(r) > 0.5:
            print(f" → RIDGE (aligned)")
        elif abs(r) < 0.2:
            print(f" → SHELL (independent)")
        else:
            print(f" → mixed")


# ═══════════════════════════════════════════════════════════════
#  The period-3 residue: does 3|period determine attractor basin?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PERIOD STRUCTURE AND ATTRACTOR BASINS")
print("=" * 70)

# For periodic functions, classify by period mod 3
periodic_pts = [p for p in landscape_points if p["period"]]
for mod3_class in [0, 1, 2]:
    pts = [p for p in periodic_pts if p["period"] % 3 == mod3_class]
    if not pts:
        continue
    avg_d_p3 = sum(p["d_p3"] for p in pts) / len(pts)
    avg_schmidt = sum(p["schmidt"] for p in pts) / len(pts)
    periods = sorted(set(p["period"] for p in pts))
    print(f"\n  period ≡ {mod3_class} mod 3 (n={len(pts)}, periods: {periods[:8]}):")
    print(f"    avg d(p3) = {avg_d_p3:.4f}, avg Schmidt = {avg_schmidt:.3f}")
    for p in sorted(pts, key=lambda p: p["d_p3"])[:3]:
        print(f"    closest: {p['name']:20s} d(p3)={p['d_p3']:.4f} Schmidt={p['schmidt']:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Key finding test: closest to p3 ⟺ 3 | period?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("KEY TEST: does nearest-attractor = period-3 ⟺ 3|period?")
print("=" * 70)

correct = 0
total = 0
for pt in periodic_pts:
    predicted_resonant = (pt["nearest"] == "period-3")
    actually_resonant = (pt["period"] % 3 == 0)
    match = predicted_resonant == actually_resonant
    correct += match
    total += 1
    if not match:
        print(f"  MISMATCH: {pt['name']:<20s} predicted={pt['nearest']}, "
              f"period={pt['period']}, 3|p={actually_resonant}")

print(f"\n  Accuracy: {correct}/{total} ({correct/total*100:.0f}%)")
if correct == total:
    print(f"  → PERFECT: nearest attractor correctly predicts resonance capability")
