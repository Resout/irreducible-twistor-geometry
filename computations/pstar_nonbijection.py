"""
p* for non-bijections: does the formula generalize?

For bijections: p* = exp(Δ/2) = √(30/23), Born(θ)_fp = 10/33.
All bijections have correlation entropy ln(3) and converge to the same p*.

Non-bijections have different correlation structures:
  - Constant:    one column hot → entropy 0
  - Modular:     nearly uniform → entropy ≈ ln(9) = 2.197
  - Exponential: upper triangular → entropy ≈ 1.5
  - Quadratic:   two hot + one spread → entropy ≈ 1.2

Question: does each crystal type have its own p*, and is it
determined by its own effective Δ? Or is there a universal formula?

The composition rate 299/405 was universal across S₃ products.
But S₃ elements are all bijections. For non-bijection compositions
(e.g., modular ∘ modular), the rate might differ.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities, born_fidelity)
from solver.crystals import (Entangler, compose, RELATIONSHIP_SIGNATURES,
                             classify_relationship)

torch.set_grad_enabled(False)


# Crystal types with different correlation structures
CRYSTAL_TYPES = {
    "proportional": torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "inverse":      torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "cycle_fwd":    torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "modular":      torch.tensor([[.4,.3,.3],[.3,.4,.3],[.3,.3,.4]], dtype=torch.float32),
    "exponential":  torch.tensor([[1,0,0],[.6,.4,0],[0,0,1]], dtype=torch.float32),
    "quadratic":    torch.tensor([[1,0,0],[0,.3,.7],[0,0,1]], dtype=torch.float32),
    "logarithmic":  torch.tensor([[.8,.2,0],[0,.6,.4],[0,.1,.9]], dtype=torch.float32),
    "constant_0":   torch.tensor([[1,0,0],[1,0,0],[1,0,0]], dtype=torch.float32),
    "constant_1":   torch.tensor([[0,1,0],[0,1,0],[0,1,0]], dtype=torch.float32),
    "two_to_one":   torch.tensor([[1,0,0],[1,0,0],[0,1,0]], dtype=torch.float32),
    "uniform":      torch.tensor([[1,1,1],[1,1,1],[1,1,1]], dtype=torch.float32) / 3,
}


def correlation_entropy(corr):
    """Shannon entropy of the normalized correlation matrix."""
    c = corr / corr.sum().clamp(min=1e-10)
    flat = c.reshape(-1)
    return -sum(x.item() * math.log(max(x.item(), 1e-10)) for x in flat)


n_seeds = 30


print("=" * 80)
print("p* FOR NON-BIJECTION CRYSTALS")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Compute p* for each crystal type
# ═══════════════════════════════════════════════════════════════

print(f"\n{'type':>14s} {'H_corr':>7s} {'Schmidt':>8s} {'p*':>7s} {'Born(θ)':>8s} {'pred_p*':>8s} {'err%':>7s}")
print("-" * 70)

results = {}
for name, corr in CRYSTAL_TYPES.items():
    # Build seed-averaged crystal
    ref = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        ref += Entangler(corr, seed=seed).build().joint
    ref /= n_seeds

    sn = schmidt_number(ref)
    h_corr = correlation_entropy(corr)

    # Iterated self-composition to find fixed point
    crystal = ref.clone()
    for _ in range(40):
        crystal = compose(crystal, ref)

    bp = born_probabilities(crystal).reshape(MASS_DIM, MASS_DIM)
    marginal = bp[0]
    theta_fp = marginal[H].item()

    # Extract p*
    h_mean = sum(marginal[k].item() for k in range(H)) / H
    if h_mean > 1e-10:
        p_sq = theta_fp / h_mean
        p_star = p_sq ** 0.5
    else:
        p_star = float('inf')
        p_sq = float('inf')

    # Compute effective Δ from the Schmidt decay
    # Schmidt after 1 composition / initial Schmidt
    comp1 = compose(ref, ref)
    sn1 = schmidt_number(comp1)
    if sn > 1 and sn1 > 1:
        ratio = sn1 / sn
        eff_delta = -math.log(max(ratio, 0.01))
    else:
        eff_delta = 0
        ratio = 0

    # Predicted p* = exp(eff_delta/2)
    pred_p = math.exp(eff_delta / 2) if eff_delta > 0 else 0
    err = (p_star - pred_p) / pred_p * 100 if pred_p > 0 else 0

    results[name] = {
        'h_corr': h_corr, 'sn': sn, 'p_star': p_star,
        'theta_fp': theta_fp, 'eff_delta': eff_delta,
        'pred_p': pred_p, 'ratio': ratio,
    }

    print(f"{name:>14s} {h_corr:>7.3f} {sn:>8.3f} {p_star:>7.4f} {theta_fp:>8.5f} "
          f"{pred_p:>8.4f} {err:>+7.1f}%")


# ═══════════════════════════════════════════════════════════════
#  Correlation between h_corr and p*
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CORRELATION ENTROPY vs p*")
print("="*80)

# Sort by correlation entropy
sorted_results = sorted(results.items(), key=lambda x: x[1]['h_corr'])

h_corrs = [r['h_corr'] for _, r in sorted_results]
p_stars = [r['p_star'] for _, r in sorted_results]
thetas = [r['theta_fp'] for _, r in sorted_results]

print(f"\n  Sorted by correlation entropy:")
for name, r in sorted_results:
    print(f"    {name:>14s}: H_corr={r['h_corr']:.3f}, p*={r['p_star']:.4f}, "
          f"Born(θ)_fp={r['theta_fp']:.5f}, Schmidt_ratio={r['ratio']:.4f}")

# Check: do all bijections (proportional, inverse, cycle_fwd) give the same p*?
bij_p = [results[n]['p_star'] for n in ['proportional', 'inverse', 'cycle_fwd']]
print(f"\n  Bijection p* values: {[f'{p:.4f}' for p in bij_p]}")
print(f"  Range: {max(bij_p)-min(bij_p):.4f}")


# ═══════════════════════════════════════════════════════════════
#  Does p* follow exp(Δ_eff/2)?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("DOES p* = exp(Δ_eff/2) HOLD FOR NON-BIJECTIONS?")
print("="*80)

for name, r in sorted_results:
    if r['pred_p'] > 0 and r['p_star'] < 10:
        print(f"  {name:>14s}: p*={r['p_star']:.4f}, exp(Δ_eff/2)={r['pred_p']:.4f}, "
              f"Δ_eff={r['eff_delta']:.4f}")


# ═══════════════════════════════════════════════════════════════
#  The relationship between Schmidt and p*
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SCHMIDT NUMBER vs FIXED-POINT BORN(θ)")
print("="*80)

for name, r in sorted_results:
    if r['p_star'] < 10:
        print(f"  {name:>14s}: Schmidt={r['sn']:.3f}, Born(θ)_fp={r['theta_fp']:.5f}, "
              f"1-Born(θ)_fp={1-r['theta_fp']:.5f}")


# ═══════════════════════════════════════════════════════════════
#  Check: is Born(θ)_fp universal?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("IS Born(θ)_fp UNIVERSAL OR TYPE-DEPENDENT?")
print("="*80)

theta_vals = [r['theta_fp'] for _, r in results.items() if r['theta_fp'] < 0.9]
print(f"  Born(θ)_fp values: {[f'{t:.5f}' for t in sorted(theta_vals)]}")
print(f"  Mean: {sum(theta_vals)/len(theta_vals):.5f}")
print(f"  Range: {max(theta_vals)-min(theta_vals):.5f}")
print(f"  10/33 = {10/33:.5f}")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
