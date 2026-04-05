"""
The blow-up boundary: is it exactly (H+1)/H² = 4/9?

Previous session reported blow-up at discount ≈ 0.44 with 100+ steps.
Default 50 steps shows a Schmidt PEAK at d ≈ 0.40 then decline — no blow-up.
The instability may require higher step counts.

This script:
1. Coarse sweep with varying n_steps to find where instability appears
2. Fine sweep around the Schmidt peak
3. Check if the peak location is (H+1)/H² = 4/9
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, BORN_FLOOR, K_STAR, DELTA, schmidt_number, born_probabilities
from solver.crystals import Entangler, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)

target = (H + 1) / H**2  # 4/9 = 0.44444...

def safe_entangle(corr, discount, n_steps=50, seed=42):
    try:
        ent = Entangler(corr, seed=seed)
        ent.build(n_steps=n_steps, discount=discount)
        s = schmidt_number(ent.joint)
        if math.isnan(s) or math.isinf(s) or s > 1e10:
            return float('inf'), None
        return s, ent.joint
    except:
        return float('inf'), None

corr = torch.eye(H, dtype=torch.float32)

# ═══════════════════════════════════════════════════════════════
#  Step-count dependent sweep: where does blow-up appear?
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print(f"BLOW-UP vs SCHMIDT PEAK: (H+1)/H² = {target:.6f}")
print("=" * 80)

for n_steps in [50, 100, 200, 500]:
    print(f"\n--- n_steps = {n_steps} ---")
    print(f"  {'discount':>10s} {'Schmidt':>10s}")
    print("  " + "-" * 25)

    peak_d = 0
    peak_s = 0
    blowup_d = None

    for d_int in range(10, 80, 2):
        d = d_int / 100.0
        s, _ = safe_entangle(corr, d, n_steps=n_steps)

        if s == float('inf'):
            print(f"  {d:>10.3f}    BLOW-UP")
            if blowup_d is None:
                blowup_d = d
        else:
            print(f"  {d:>10.3f} {s:>10.3f}")
            if s > peak_s:
                peak_s = s
                peak_d = d

    print(f"\n  Peak: d={peak_d:.3f}, Schmidt={peak_s:.3f}")
    if blowup_d:
        print(f"  Blow-up at: d={blowup_d:.3f}")
        print(f"  Target (H+1)/H² = {target:.3f}")
    else:
        print(f"  No blow-up in range [0.10, 0.78]")


# ═══════════════════════════════════════════════════════════════
#  Fine sweep around the Schmidt peak at n_steps=50
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("FINE SWEEP AROUND SCHMIDT PEAK (n_steps=50)")
print("="*80)

print(f"  {'discount':>10s} {'Schmidt':>10s}")
print("  " + "-" * 25)

peak_d = 0
peak_s = 0
fine_data = []

for d_int in range(350, 450, 2):
    d = d_int / 1000.0
    s, _ = safe_entangle(corr, d)
    if s < float('inf'):
        print(f"  {d:>10.4f} {s:>10.5f}")
        fine_data.append((d, s))
        if s > peak_s:
            peak_s = s
            peak_d = d

print(f"\n  Peak location: d = {peak_d:.4f}")
print(f"  Peak Schmidt: {peak_s:.5f}")

# Quadratic fit around peak to find exact maximum
if len(fine_data) > 5:
    # Find the three points around the maximum
    idx_max = max(range(len(fine_data)), key=lambda i: fine_data[i][1])
    if 1 <= idx_max <= len(fine_data) - 2:
        d1, s1 = fine_data[idx_max - 1]
        d2, s2 = fine_data[idx_max]
        d3, s3 = fine_data[idx_max + 1]
        # Quadratic interpolation
        denom = (d1 - d2) * (d1 - d3) * (d2 - d3)
        if abs(denom) > 1e-15:
            a = (d3 * (s2 - s1) + d2 * (s1 - s3) + d1 * (s3 - s2)) / denom
            b = (d3**2 * (s1 - s2) + d2**2 * (s3 - s1) + d1**2 * (s2 - s3)) / denom
            d_peak = -b / (2 * a) if abs(a) > 1e-15 else d2
            print(f"  Quadratic-interpolated peak: d = {d_peak:.6f}")
            print(f"  Target (H+1)/H² = {target:.6f}")
            print(f"  Difference: {abs(d_peak - target):.6f}")

            # Other candidates for the peak
            candidates = {
                "(H+1)/H²": (H+1)/H**2,        # 4/9 = 0.4444
                "2/(H+1)":  2/(H+1),            # 2/4 = 0.5
                "H/(H²-1)": H/(H**2-1),         # 3/8 = 0.375
                "(H-1)/H":  (H-1)/H,            # 2/3 = 0.6667
                "1/ln(H)":  1/math.log(H),       # 0.9102
                "1/H":      1/H,                 # 0.3333
                "sqrt(1/H)": math.sqrt(1/H),     # 0.5774
                "ln(H)/H":  math.log(H)/H,       # 0.3662
                "2/H²":     2/H**2,              # 0.2222
            }
            print(f"\n  Candidate comparison:")
            for cname, cval in sorted(candidates.items(), key=lambda x: abs(x[1] - d_peak)):
                print(f"    {cname:>12s} = {cval:.6f}  (diff = {abs(cval - d_peak):.6f})")


# ═══════════════════════════════════════════════════════════════
#  Multi-seed peak location
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SEED DEPENDENCE OF PEAK LOCATION (20 seeds)")
print("="*80)

peak_locs = []
for seed in range(20):
    best_d = 0
    best_s = 0
    for d_int in range(35, 46):
        d = d_int / 100.0
        s, _ = safe_entangle(corr, d, seed=seed)
        if s < float('inf') and s > best_s:
            best_s = s
            best_d = d
    peak_locs.append(best_d)
    print(f"  Seed {seed:>2d}: peak at d={best_d:.2f}, Schmidt={best_s:.3f}")

avg_peak = sum(peak_locs) / len(peak_locs)
print(f"\n  Average peak location: {avg_peak:.4f}")
print(f"  Target (H+1)/H² = {target:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Type dependence
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("PEAK LOCATION BY CORRELATION TYPE")
print("="*80)

for name, corr_sig in RELATIONSHIP_SIGNATURES.items():
    best_d = 0
    best_s = 0
    for d_int in range(20, 60, 2):
        d = d_int / 100.0
        s, _ = safe_entangle(corr_sig, d)
        if s < float('inf') and s > best_s:
            best_s = s
            best_d = d
    print(f"  {name:>15s}: peak at d={best_d:.2f}, Schmidt={best_s:.3f}")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
