"""
θ-fingerprint as problem classifier.

The discovery: 87% of crystal variation is in the θ↔h correlations.
The θ-fingerprint is 5-8x more discriminative than Sym² for relationship
classification. But we've only used it for CLASSIFICATION of known types.

Question: can the θ-fingerprint classify mathematical functions by their
COMPUTATIONAL TYPE — periodic vs arithmetic vs irregular? If it can,
the crystal knows what kind of problem it's looking at before we tell it.

Second question: for a raw (non-resonant) crystal, does the θ-fingerprint
predict the resonant frame? If so, the crystal can select its own
observation frame from a single construction.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import numpy as np
from collections import defaultdict
from solver.algebra import H, MASS_DIM, schmidt_number, born_probabilities
from solver.crystals import Entangler, classify_relationship, slot_measure
from solver.functional import (
    function_to_correlation, build_function_crystal,
    _detect_output_period, Binner,
)
from solver import compute

torch.set_grad_enabled(False)


# ═══════════════════════════════════════════════════════════════
#  θ-fingerprint extraction
# ═══════════════════════════════════════════════════════════════

def theta_fingerprint(joint: torch.Tensor) -> torch.Tensor:
    """Extract θ-fingerprint: the 6 off-diagonal θ↔h components.

    From the [4,4] joint mass, θ is index 3. The θ-fingerprint is:
      joint[0,3], joint[1,3], joint[2,3]  (h→θ correlations)
      joint[3,0], joint[3,1], joint[3,2]  (θ→h correlations)

    These 6 complex numbers encode how ignorance relates to each
    hypothesis. Normalized by L2 norm for comparison.
    """
    fp = torch.stack([
        joint[0, 3], joint[1, 3], joint[2, 3],
        joint[3, 0], joint[3, 1], joint[3, 2],
    ])
    norm = fp.abs().pow(2).sum().sqrt()
    if norm > 1e-10:
        fp = fp / norm
    return fp


def theta_distance(fp1: torch.Tensor, fp2: torch.Tensor) -> float:
    """L2 distance between θ-fingerprints."""
    return (fp1 - fp2).abs().pow(2).sum().sqrt().item()


def theta_phase_structure(fp: torch.Tensor) -> dict:
    """Extract phase structure from θ-fingerprint.

    The PHASES carry the non-commutative information (sign representation).
    The MAGNITUDES carry the correlation strength.
    """
    mags = fp.abs()
    phases = torch.angle(fp)

    # Asymmetry: how different are h→θ vs θ→h?
    h_to_theta = mags[:3]
    theta_to_h = mags[3:]
    asymmetry = (h_to_theta - theta_to_h).abs().sum().item()

    # Phase coherence: are all phases aligned?
    phase_spread = phases.std().item()

    # Magnitude concentration: is one hypothesis dominant?
    mag_norm = mags / mags.sum().clamp(min=1e-10)
    concentration = (mag_norm ** 2).sum().item()  # 1/6 = uniform, 1 = concentrated

    return {
        "asymmetry": asymmetry,
        "phase_spread": phase_spread,
        "concentration": concentration,
        "h_to_theta_mags": h_to_theta.tolist(),
        "theta_to_h_mags": theta_to_h.tolist(),
    }


# ═══════════════════════════════════════════════════════════════
#  Build crystals from diverse function families
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("θ-FINGERPRINT AS PROBLEM CLASSIFIER")
print("=" * 70)

domain = range(1, 301)

families = {}

# --- FAMILY 1: Modular exponentiation (periodic, period divides φ(m)) ---
print("\n── Family: Modular Exponentiation ──")
mod_exp_crystals = []
for base, mod in [(2, 7), (3, 7), (2, 13), (3, 13), (5, 7), (7, 19), (2, 11), (3, 11)]:
    f = lambda x, b=base, m=mod: pow(b, int(x), m)
    fc = build_function_crystal(f"{base}^x mod {mod}", f, domain)
    period = None
    outputs = [f(x) for x in range(1, 51)]
    period = _detect_output_period(outputs)
    fp = theta_fingerprint(fc.joint)
    mod_exp_crystals.append({
        "name": f"{base}^x mod {mod}",
        "crystal": fc,
        "theta_fp": fp,
        "schmidt": fc.schmidt,
        "period": period,
        "relationship": fc.relationship,
    })
    print(f"  {base}^x mod {mod:2d}: Schmidt={fc.schmidt:.3f}, period={period}, rel={fc.relationship}")

families["mod_exp"] = mod_exp_crystals

# --- FAMILY 2: Polynomial (non-periodic, algebraic growth) ---
print("\n── Family: Polynomials (mod small primes) ──")
poly_crystals = []
for coeffs, mod, label in [
    ([0, 1], 7, "x mod 7"),
    ([0, 0, 1], 7, "x² mod 7"),
    ([0, 0, 0, 1], 7, "x³ mod 7"),
    ([1, 2, 3], 11, "(3x²+2x+1) mod 11"),
    ([0, 1], 13, "x mod 13"),
    ([0, 0, 1], 13, "x² mod 13"),
]:
    def poly(x, c=coeffs, m=mod):
        return sum(c_i * int(x)**i for i, c_i in enumerate(c)) % m

    fc = build_function_crystal(label, poly, domain)
    outputs = [poly(x) for x in range(1, 51)]
    period = _detect_output_period(outputs)
    fp = theta_fingerprint(fc.joint)
    poly_crystals.append({
        "name": label,
        "crystal": fc,
        "theta_fp": fp,
        "schmidt": fc.schmidt,
        "period": period,
        "relationship": fc.relationship,
    })
    print(f"  {label}: Schmidt={fc.schmidt:.3f}, period={period}, rel={fc.relationship}")

families["polynomial"] = poly_crystals

# --- FAMILY 3: Arithmetic functions (irregular, number-theoretic) ---
print("\n── Family: Arithmetic Functions ──")
arith_crystals = []
for name, f in [
    ("φ(n)", lambda n: compute.euler_totient(max(1, int(n)))),
    ("d(n)", lambda n: compute.divisor_count(max(1, int(n)))),
    ("σ(n)", lambda n: compute.divisor_sum(max(1, int(n)))),
    ("digit_sum", lambda n: compute.digit_sum(max(1, int(n)))),
    ("ω(n)", lambda n: len(compute.factorize(max(2, int(n))))),
    ("Ω(n)", lambda n: sum(compute.factorize(max(2, int(n))).values())),
]:
    fc = build_function_crystal(name, f, domain)
    outputs = [f(x) for x in range(1, 51)]
    period = _detect_output_period(outputs)
    fp = theta_fingerprint(fc.joint)
    arith_crystals.append({
        "name": name,
        "crystal": fc,
        "theta_fp": fp,
        "schmidt": fc.schmidt,
        "period": period,
        "relationship": fc.relationship,
    })
    print(f"  {name}: Schmidt={fc.schmidt:.3f}, period={period}, rel={fc.relationship}")

families["arithmetic"] = arith_crystals

# --- FAMILY 4: Periodic with various periods ---
print("\n── Family: Pure Periodic (various periods) ──")
periodic_crystals = []
for p in [2, 3, 4, 5, 6, 7, 8, 9, 10, 12]:
    f = lambda x, _p=p: int(x) % _p
    fc = build_function_crystal(f"x mod {p}", f, domain)
    outputs = [f(x) for x in range(1, 51)]
    period = _detect_output_period(outputs)
    fp = theta_fingerprint(fc.joint)
    periodic_crystals.append({
        "name": f"x mod {p}",
        "crystal": fc,
        "theta_fp": fp,
        "schmidt": fc.schmidt,
        "period": period,
        "relationship": fc.relationship,
    })
    print(f"  x mod {p:2d}: Schmidt={fc.schmidt:.3f}, period={period}, rel={fc.relationship}")

families["periodic"] = periodic_crystals


# ═══════════════════════════════════════════════════════════════
#  θ-fingerprint distance matrix WITHIN and BETWEEN families
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("θ-FINGERPRINT DISTANCES: within vs between families")
print("=" * 70)

all_crystals = []
for fam_name, crystals in families.items():
    for c in crystals:
        all_crystals.append((fam_name, c))

# Within-family distances
print("\nWithin-family average θ-distances:")
for fam_name, crystals in families.items():
    if len(crystals) < 2:
        continue
    dists = []
    for i in range(len(crystals)):
        for j in range(i + 1, len(crystals)):
            d = theta_distance(crystals[i]["theta_fp"], crystals[j]["theta_fp"])
            dists.append(d)
    avg = sum(dists) / len(dists)
    std = (sum((d - avg) ** 2 for d in dists) / len(dists)) ** 0.5
    print(f"  {fam_name:15s}: avg={avg:.4f} ± {std:.4f} (n={len(dists)})")

# Between-family distances
print("\nBetween-family average θ-distances:")
fam_names = list(families.keys())
for i in range(len(fam_names)):
    for j in range(i + 1, len(fam_names)):
        dists = []
        for ci in families[fam_names[i]]:
            for cj in families[fam_names[j]]:
                d = theta_distance(ci["theta_fp"], cj["theta_fp"])
                dists.append(d)
        avg = sum(dists) / len(dists)
        print(f"  {fam_names[i]:15s} ↔ {fam_names[j]:15s}: avg={avg:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Phase structure comparison
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PHASE STRUCTURE BY FAMILY")
print("=" * 70)

for fam_name, crystals in families.items():
    asymmetries = []
    phase_spreads = []
    concentrations = []
    for c in crystals:
        ps = theta_phase_structure(c["theta_fp"])
        asymmetries.append(ps["asymmetry"])
        phase_spreads.append(ps["phase_spread"])
        concentrations.append(ps["concentration"])

    print(f"\n  {fam_name}:")
    print(f"    asymmetry:     {sum(asymmetries)/len(asymmetries):.4f} ± "
          f"{(sum((a - sum(asymmetries)/len(asymmetries))**2 for a in asymmetries)/len(asymmetries))**0.5:.4f}")
    print(f"    phase_spread:  {sum(phase_spreads)/len(phase_spreads):.4f} ± "
          f"{(sum((p - sum(phase_spreads)/len(phase_spreads))**2 for p in phase_spreads)/len(phase_spreads))**0.5:.4f}")
    print(f"    concentration: {sum(concentrations)/len(concentrations):.4f} ± "
          f"{(sum((c - sum(concentrations)/len(concentrations))**2 for c in concentrations)/len(concentrations))**0.5:.4f}")


# ═══════════════════════════════════════════════════════════════
#  The key question: does θ-fingerprint predict resonant frame?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("θ-FINGERPRINT vs RESONANT FRAME")
print("Does the raw crystal's ignorance structure predict which")
print("frame transformation will make it resonate?")
print("=" * 70)

# For each modular exponentiation, build crystals at DIFFERENT frames
# and compare their θ-fingerprints
test_functions = [
    ("3^x mod 7", lambda x: pow(3, int(x), 7), 6),    # period 6
    ("2^x mod 13", lambda x: pow(2, int(x), 13), 12),  # period 12
    ("5^x mod 13", lambda x: pow(5, int(x), 13), 4),   # period 4
    ("2^x mod 7", lambda x: pow(2, int(x), 7), 3),     # period 3 (already resonant)
    ("2^x mod 11", lambda x: pow(2, int(x), 11), 10),  # period 10
]

for name, f, true_period in test_functions:
    print(f"\n  ── {name} (period {true_period}) ──")

    # Raw crystal (k=1)
    fc_raw = build_function_crystal(f"{name}_raw", f, domain)
    fp_raw = theta_fingerprint(fc_raw.joint)
    ps_raw = theta_phase_structure(fp_raw)

    # Try frames k=1..12 and find the best
    best_k, best_schmidt = 1, fc_raw.schmidt
    frame_data = []
    for k in range(1, 13):
        fk = lambda x, _k=k, _f=f: _f(_k * int(x))
        try:
            fc_k = build_function_crystal(f"{name}_k{k}", fk, domain)
            fp_k = theta_fingerprint(fc_k.joint)
            d_from_raw = theta_distance(fp_raw, fp_k)
            frame_data.append({
                "k": k,
                "schmidt": fc_k.schmidt,
                "theta_dist": d_from_raw,
                "asymmetry": theta_phase_structure(fp_k)["asymmetry"],
                "concentration": theta_phase_structure(fp_k)["concentration"],
            })
            if fc_k.schmidt > best_schmidt:
                best_schmidt = fc_k.schmidt
                best_k = k
        except Exception:
            pass

    print(f"    Raw: Schmidt={fc_raw.schmidt:.3f}, asymm={ps_raw['asymmetry']:.4f}, "
          f"conc={ps_raw['concentration']:.4f}")
    print(f"    Best frame: k={best_k}, Schmidt={best_schmidt:.3f}")

    # Show top 3 frames by Schmidt
    frame_data.sort(key=lambda d: -d["schmidt"])
    for fd in frame_data[:4]:
        marker = " ◄ BEST" if fd["k"] == best_k else ""
        print(f"    k={fd['k']:2d}: Schmidt={fd['schmidt']:.3f}, "
              f"θ-dist from raw={fd['theta_dist']:.4f}, "
              f"asymm={fd['asymmetry']:.4f}, "
              f"conc={fd['concentration']:.4f}{marker}")

    # Check: does the resonant frame have max θ-distance from raw?
    if frame_data:
        max_dist_k = max(frame_data, key=lambda d: d["theta_dist"])
        print(f"    Max θ-distance from raw: k={max_dist_k['k']} "
              f"(dist={max_dist_k['theta_dist']:.4f})")

        # Correlation: θ-distance from raw vs Schmidt
        schmidts = [d["schmidt"] for d in frame_data]
        dists = [d["theta_dist"] for d in frame_data]
        if len(schmidts) > 2:
            s_mean = sum(schmidts) / len(schmidts)
            d_mean = sum(dists) / len(dists)
            cov = sum((s - s_mean) * (d - d_mean) for s, d in zip(schmidts, dists))
            s_std = (sum((s - s_mean) ** 2 for s in schmidts)) ** 0.5
            d_std = (sum((d - d_mean) ** 2 for d in dists)) ** 0.5
            if s_std > 1e-10 and d_std > 1e-10:
                r = cov / (s_std * d_std)
                print(f"    Correlation(θ-dist, Schmidt): r = {r:.3f}")


# ═══════════════════════════════════════════════════════════════
#  θ-fingerprint PCA: do families separate in 2D?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("θ-FINGERPRINT PCA")
print("=" * 70)

# Collect all θ-fingerprints as real vectors (real + imag parts)
all_fps = []
all_labels = []
for fam_name, crystals in families.items():
    for c in crystals:
        fp = c["theta_fp"]
        # Convert to real vector: [re0, im0, re1, im1, ...]
        real_fp = torch.stack([fp.real, fp.imag], dim=1).reshape(-1)
        all_fps.append(real_fp)
        all_labels.append(fam_name)

fps_matrix = torch.stack(all_fps)  # [N, 12]
N = fps_matrix.shape[0]

# Center
mean = fps_matrix.mean(dim=0)
centered = fps_matrix - mean

# Covariance and eigendecomposition
cov = (centered.T @ centered) / (N - 1)
eigenvalues, eigenvectors = torch.linalg.eigh(cov)

# Sort descending
idx = eigenvalues.argsort(descending=True)
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

# Variance explained
total_var = eigenvalues.sum()
print(f"\n  Total variance: {total_var:.4f}")
print(f"  Variance explained by first 4 PCs:")
for i in range(min(4, len(eigenvalues))):
    print(f"    PC{i+1}: {eigenvalues[i]:.4f} ({eigenvalues[i]/total_var*100:.1f}%)")

# Project onto first 2 PCs
projected = centered @ eigenvectors[:, :2]

# Report cluster centers
print(f"\n  Cluster centers in PC1-PC2:")
for fam_name in families:
    mask = [l == fam_name for l in all_labels]
    fam_pts = projected[mask]
    center = fam_pts.mean(dim=0)
    spread = ((fam_pts - center).pow(2).sum(dim=1).mean()).sqrt()
    print(f"    {fam_name:15s}: ({center[0]:+.3f}, {center[1]:+.3f}), spread={spread:.4f}")

# Inter-cluster distances
print(f"\n  Inter-cluster distances (PC1-PC2 center-to-center):")
centers = {}
for fam_name in families:
    mask = [l == fam_name for l in all_labels]
    centers[fam_name] = projected[mask].mean(dim=0)

fam_list = list(families.keys())
for i in range(len(fam_list)):
    for j in range(i + 1, len(fam_list)):
        d = (centers[fam_list[i]] - centers[fam_list[j]]).norm().item()
        print(f"    {fam_list[i]:15s} ↔ {fam_list[j]:15s}: {d:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

# Check: are within-family distances smaller than between?
within = {}
for fam_name, crystals in families.items():
    dists = []
    for i in range(len(crystals)):
        for j in range(i + 1, len(crystals)):
            d = theta_distance(crystals[i]["theta_fp"], crystals[j]["theta_fp"])
            dists.append(d)
    if dists:
        within[fam_name] = sum(dists) / len(dists)

between = []
for i in range(len(fam_list)):
    for j in range(i + 1, len(fam_list)):
        dists = []
        for ci in families[fam_list[i]]:
            for cj in families[fam_list[j]]:
                d = theta_distance(ci["theta_fp"], cj["theta_fp"])
                dists.append(d)
        between.append(sum(dists) / len(dists))

avg_within = sum(within.values()) / len(within) if within else 0
avg_between = sum(between) / len(between) if between else 0

print(f"\n  Average within-family θ-distance:  {avg_within:.4f}")
print(f"  Average between-family θ-distance: {avg_between:.4f}")
if avg_within > 0:
    ratio = avg_between / avg_within
    print(f"  Ratio (between/within):            {ratio:.2f}")
    if ratio > 1.5:
        print(f"  → θ-fingerprint SEPARATES families ({ratio:.1f}x)")
    elif ratio > 1.1:
        print(f"  → θ-fingerprint weakly separates ({ratio:.1f}x)")
    else:
        print(f"  → θ-fingerprint does NOT separate families ({ratio:.1f}x)")
