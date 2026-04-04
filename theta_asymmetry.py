"""
θ-asymmetry as invertibility detector.

The θ-fingerprint doesn't separate mathematical families. But the
PHASE STRUCTURE does: modular exponentiation has low asymmetry (0.19),
arithmetic functions have high asymmetry (0.61).

Hypothesis: asymmetry detects whether a function is invertible.
  - Bijections on Z/3Z: information flows equally both ways → low asymmetry
  - Many-to-one functions: information flows asymmetrically → high asymmetry

If this holds, the crystal detects invertibility from its ignorance
structure alone. This is a property of the FUNCTION, not the binner.

Second thread: the resonant frame is the θ-extremum. For 5^x mod 13,
r(θ-dist, Schmidt) = 0.999. Can we use this to predict the resonant
frame without trying all multipliers?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, schmidt_number
from solver.crystals import Entangler, classify_relationship
from solver.functional import (
    build_function_crystal, function_to_correlation,
    _detect_output_period, Binner,
)
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


def theta_asymmetry(joint):
    """Asymmetry between h→θ and θ→h directions."""
    fp = theta_fingerprint(joint)
    h_to_theta = fp[:3].abs()
    theta_to_h = fp[3:].abs()
    return (h_to_theta - theta_to_h).abs().sum().item()


def theta_concentration(joint):
    """How concentrated the θ-fingerprint magnitudes are."""
    fp = theta_fingerprint(joint)
    mags = fp.abs()
    mags_n = mags / mags.sum().clamp(min=1e-10)
    return (mags_n ** 2).sum().item()


domain = range(1, 301)


# ═══════════════════════════════════════════════════════════════
#  Test 1: Bijections vs many-to-one on Z/3Z
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 1: BIJECTIONS vs MANY-TO-ONE")
print("All functions Z/3Z → Z/3Z, measured through the crystal")
print("=" * 70)

# Build all 27 functions Z/3Z → Z/3Z via x mod 3 → f(x mod 3)
# using different sampling approaches

# The 6 bijections (permutations of {0,1,2})
import itertools
bijections = list(itertools.permutations([0, 1, 2]))
print(f"\nBijections (n=6):")
bij_asymmetries = []
for perm in bijections:
    table = {0: perm[0], 1: perm[1], 2: perm[2]}
    f = lambda x, t=table: t[int(x) % 3]
    fc = build_function_crystal(f"bij_{perm}", f, domain)
    a = theta_asymmetry(fc.joint)
    c = theta_concentration(fc.joint)
    bij_asymmetries.append(a)
    print(f"  {perm}: Schmidt={fc.schmidt:.3f}, asymm={a:.4f}, conc={c:.4f}")

# The 21 non-bijections
print(f"\nNon-bijections (sample of 12):")
non_bij_asymmetries = []
non_bij_examples = [
    (0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 0),
    (1, 1, 1), (2, 2, 2), (0, 1, 1), (1, 0, 0),
    (2, 0, 2), (1, 2, 1), (0, 2, 0), (2, 1, 2),
]
for vals in non_bij_examples:
    if vals in bijections:
        continue
    table = {0: vals[0], 1: vals[1], 2: vals[2]}
    f = lambda x, t=table: t[int(x) % 3]
    fc = build_function_crystal(f"nbij_{vals}", f, domain)
    a = theta_asymmetry(fc.joint)
    c = theta_concentration(fc.joint)
    non_bij_asymmetries.append(a)
    print(f"  {vals}: Schmidt={fc.schmidt:.3f}, asymm={a:.4f}, conc={c:.4f}")

avg_bij = sum(bij_asymmetries) / len(bij_asymmetries)
avg_nbij = sum(non_bij_asymmetries) / len(non_bij_asymmetries) if non_bij_asymmetries else 0
print(f"\n  Average bijection asymmetry:     {avg_bij:.4f}")
print(f"  Average non-bijection asymmetry: {avg_nbij:.4f}")
print(f"  Ratio:                           {avg_nbij / avg_bij:.2f}x" if avg_bij > 0 else "")


# ═══════════════════════════════════════════════════════════════
#  Test 2: Functions with known invertibility properties
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 2: FUNCTIONS WITH KNOWN INVERTIBILITY")
print("=" * 70)

test_cases = [
    # (name, function, invertible?)
    ("x mod 3 (identity-like)", lambda x: int(x) % 3, True),
    ("2^x mod 7 (period 3, bijective on Z/3Z)", lambda x: pow(2, int(x), 7), True),
    ("3^x mod 7 (period 6, not bijective on Z/6Z period)", lambda x: pow(3, int(x), 7), False),
    ("x² mod 7 (quadratic residues, 2-to-1)", lambda x: (int(x) ** 2) % 7, False),
    ("x³ mod 7 (cubic, bijective on Z/7Z)", lambda x: (int(x) ** 3) % 7, True),
    ("φ(n) (many-to-one)", lambda x: compute.euler_totient(max(1, int(x))), False),
    ("d(n) (many-to-one)", lambda x: compute.divisor_count(max(1, int(x))), False),
    ("x mod 2 (2-to-1 on Z/3Z)", lambda x: int(x) % 2, False),
    ("(2x+1) mod 3 (affine bijection)", lambda x: (2 * int(x) + 1) % 3, True),
    ("(x+1) mod 3 (shift, bijective)", lambda x: (int(x) + 1) % 3, True),
    ("-x mod 3 (inversion, bijective)", lambda x: (-int(x)) % 3, True),
]

print(f"\n{'Function':<45s} {'Inv?':>5s} {'Schmidt':>8s} {'Asymm':>8s} {'Conc':>8s}")
print("-" * 78)

inv_asymm = []
non_inv_asymm = []
for name, f, invertible in test_cases:
    fc = build_function_crystal(name, f, domain)
    a = theta_asymmetry(fc.joint)
    c = theta_concentration(fc.joint)
    tag = "bij" if invertible else "m→1"
    if invertible:
        inv_asymm.append(a)
    else:
        non_inv_asymm.append(a)
    print(f"  {name:<43s} {tag:>5s} {fc.schmidt:>8.3f} {a:>8.4f} {c:>8.4f}")

print(f"\n  Average invertible asymmetry:     {sum(inv_asymm)/len(inv_asymm):.4f} (n={len(inv_asymm)})")
print(f"  Average non-invertible asymmetry: {sum(non_inv_asymm)/len(non_inv_asymm):.4f} (n={len(non_inv_asymm)})")


# ═══════════════════════════════════════════════════════════════
#  Test 3: Resonant frame as θ-extremum (detailed)
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 3: RESONANT FRAME AS θ-EXTREMUM")
print("For each function, measure θ-fingerprint at every frame k.")
print("Is the resonant frame a local extremum in θ-distance from raw?")
print("=" * 70)

# Build the canonical resonant θ-fingerprint (from x mod 3)
fc_canon = build_function_crystal("x mod 3", lambda x: int(x) % 3, domain)
fp_canon = theta_fingerprint(fc_canon.joint)
print(f"\nCanonical resonant fingerprint (x mod 3): Schmidt={fc_canon.schmidt:.3f}")
print(f"  θ-fingerprint: {[f'{v.item():.4f}' for v in fp_canon.abs()]}")

# For each test function: measure θ-distance from canonical at each frame
test_functions = [
    ("3^x mod 7", lambda x: pow(3, int(x), 7), 6),
    ("2^x mod 13", lambda x: pow(2, int(x), 13), 12),
    ("5^x mod 13", lambda x: pow(5, int(x), 13), 4),
    ("2^x mod 11", lambda x: pow(2, int(x), 11), 10),
    ("x mod 6", lambda x: int(x) % 6, 6),
    ("x mod 9", lambda x: int(x) % 9, 9),
    ("x mod 12", lambda x: int(x) % 12, 12),
]

for name, f, true_period in test_functions:
    print(f"\n  ── {name} (period {true_period}) ──")

    frames = []
    for k in range(1, 19):
        fk = lambda x, _k=k, _f=f: _f(_k * int(x))
        try:
            fc_k = build_function_crystal(f"{name}_k{k}", fk, domain)
            fp_k = theta_fingerprint(fc_k.joint)
            d_canon = (fp_k - fp_canon).abs().pow(2).sum().sqrt().item()
            frames.append({
                "k": k,
                "schmidt": fc_k.schmidt,
                "theta_dist_canon": d_canon,
                "asymm": theta_asymmetry(fc_k.joint),
            })
        except Exception:
            pass

    # Sort by Schmidt
    frames.sort(key=lambda d: -d["schmidt"])

    # Report: top frames + θ-distance from canonical
    for fd in frames[:5]:
        marker = " ★" if fd["schmidt"] > 3.0 else (" ●" if fd["schmidt"] > 2.5 else "")
        print(f"    k={fd['k']:2d}: Schmidt={fd['schmidt']:.3f}, "
              f"θ-dist to canon={fd['theta_dist_canon']:.4f}, "
              f"asymm={fd['asymm']:.4f}{marker}")

    # Correlation: θ-distance from canonical vs Schmidt
    schmidts = [d["schmidt"] for d in frames]
    dists = [d["theta_dist_canon"] for d in frames]
    if len(schmidts) > 2:
        s_mean = sum(schmidts) / len(schmidts)
        d_mean = sum(dists) / len(dists)
        cov = sum((s - s_mean) * (d - d_mean) for s, d in zip(schmidts, dists))
        s_std = (sum((s - s_mean) ** 2 for s in schmidts)) ** 0.5
        d_std = (sum((d - d_mean) ** 2 for d in dists)) ** 0.5
        if s_std > 1e-10 and d_std > 1e-10:
            r = cov / (s_std * d_std)
            print(f"    Correlation(Schmidt, θ-dist to canon): r = {r:.3f}")

    # Check: is the best frame also the closest to canonical?
    closest = min(frames, key=lambda d: d["theta_dist_canon"])
    best = frames[0]
    print(f"    Best Schmidt at k={best['k']}, closest to canon at k={closest['k']}")
    if closest["k"] == best["k"]:
        print(f"    → YES: best frame = closest to canonical in θ-space")


# ═══════════════════════════════════════════════════════════════
#  Test 4: The θ-diagonal — what does joint[3,3] encode?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 4: THE θ-SELF TERM — joint[3,3]")
print("This is the ignorance-about-ignorance. What does it encode?")
print("=" * 70)

all_theta_self = []
for name, f, info in [
    ("x mod 3", lambda x: int(x) % 3, "period 3, resonant"),
    ("x mod 2", lambda x: int(x) % 2, "period 2, partial"),
    ("x mod 5", lambda x: int(x) % 5, "period 5, non-resonant"),
    ("x mod 7", lambda x: int(x) % 7, "period 7, non-resonant"),
    ("2^x mod 7", lambda x: pow(2, int(x), 7), "period 3, resonant"),
    ("2^x mod 13", lambda x: pow(2, int(x), 13), "period 12, non-resonant"),
    ("φ(n)", lambda x: compute.euler_totient(max(1, int(x))), "arithmetic"),
    ("d(n)", lambda x: compute.divisor_count(max(1, int(x))), "arithmetic"),
    ("x²", lambda x: int(x) ** 2, "quadratic"),
    ("random mod 3", lambda x: hash(str(x)) % 3, "pseudorandom"),
]:
    fc = build_function_crystal(name, f, domain)
    theta_self = fc.joint[3, 3]
    theta_self_born = abs(theta_self) ** 2
    total_born = fc.joint.abs().pow(2).sum().item()
    theta_fraction = theta_self_born / total_born

    # Also compute: how much of the joint mass is in the θ row/column?
    theta_row_mass = fc.joint[3, :].abs().pow(2).sum().item() / total_born
    theta_col_mass = fc.joint[:, 3].abs().pow(2).sum().item() / total_born

    all_theta_self.append({
        "name": name,
        "info": info,
        "theta_self": theta_self,
        "theta_fraction": theta_fraction,
        "theta_row_mass": theta_row_mass,
        "theta_col_mass": theta_col_mass,
        "schmidt": fc.schmidt,
    })
    print(f"  {name:<20s} ({info:<25s}): θ_self={theta_self:.4f}, "
          f"θ-frac={theta_fraction:.4f}, "
          f"θ-row={theta_row_mass:.3f}, θ-col={theta_col_mass:.3f}, "
          f"Schmidt={fc.schmidt:.3f}")

# Correlation: θ-fraction vs Schmidt
schmidts = [d["schmidt"] for d in all_theta_self]
fracs = [d["theta_fraction"] for d in all_theta_self]
s_mean = sum(schmidts) / len(schmidts)
f_mean = sum(fracs) / len(fracs)
cov = sum((s - s_mean) * (f - f_mean) for s, f in zip(schmidts, fracs))
s_std = (sum((s - s_mean) ** 2 for s in schmidts)) ** 0.5
f_std = (sum((f - f_mean) ** 2 for f in fracs)) ** 0.5
if s_std > 1e-10 and f_std > 1e-10:
    r = cov / (s_std * f_std)
    print(f"\n  Correlation(θ-fraction, Schmidt): r = {r:.3f}")

# θ-row vs θ-col (asymmetry of ignorance flow)
rows = [d["theta_row_mass"] for d in all_theta_self]
cols = [d["theta_col_mass"] for d in all_theta_self]
print(f"\n  θ-row mass range: [{min(rows):.3f}, {max(rows):.3f}]")
print(f"  θ-col mass range: [{min(cols):.3f}, {max(cols):.3f}]")
print(f"  Average |row - col|: {sum(abs(r-c) for r,c in zip(rows,cols))/len(rows):.4f}")


# ═══════════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
