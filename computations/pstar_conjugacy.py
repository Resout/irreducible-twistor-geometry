"""
What determines p*? — The product state fixed point by conjugacy class.

Iterated self-composition converges to a product state (Schmidt → 1).
The product state has q = p² in the symmetrized basis (Principle 33).
But what is p*? Does it depend on the conjugacy class?

S₃ has 3 conjugacy classes:
  {e}:           identity
  {(01),(02),(12)}: transpositions
  {(012),(021)}:    3-cycles

If p* is conjugacy-class-specific, it carries the TYPE information
of the group element through total decorrelation. The crystal forgets
the relationship but remembers WHAT KIND of relationship it was.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities, born_fidelity)
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

S3_CORRS = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

n_seeds = 50

# Build seed-averaged crystals
refs = {}
for name, corr in S3_CORRS.items():
    s = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        s += Entangler(corr, seed=seed).build().joint
    refs[name] = s / n_seeds


print("=" * 80)
print("WHAT DETERMINES p*? — PRODUCT STATE BY CONJUGACY CLASS")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Iterated self-composition to find fixed points
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Fixed points of iterated self-composition (30 steps) ---\n")

for name in ["e", "(01)", "(02)", "(12)", "(012)", "(021)"]:
    crystal = refs[name]
    for _ in range(30):
        crystal = compose(crystal, refs[name])

    bp = born_probabilities(crystal).reshape(MASS_DIM, MASS_DIM)
    sn = schmidt_number(crystal)

    # The marginal is the row (all rows identical for product state)
    marginal = bp[0]
    h_vals = [marginal[i].item() for i in range(H)]
    theta_val = marginal[H].item()

    # p parameter: ratio of off-diagonal h to diagonal h in symmetrized basis
    # For S₃-symmetric: a = diag h, b = off-diag h, c = h-θ, d = θ-θ
    # p = c/a at the fixed point

    print(f"  {name:>6s}: marginal = [{', '.join(f'{v:.4f}' for v in h_vals)}, θ={theta_val:.4f}]"
          f"  Schmidt = {sn:.4f}")

    # Check product state: is each row the same?
    max_row_diff = 0
    for i in range(1, MASS_DIM):
        for j in range(MASS_DIM):
            diff = abs(bp[i,j].item() - bp[0,j].item())
            max_row_diff = max(max_row_diff, diff)
    print(f"         max row difference = {max_row_diff:.6f} (0 = perfect product state)")


# ═══════════════════════════════════════════════════════════════
#  Compare marginals across conjugacy classes
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Marginals by conjugacy class ---\n")

# Compute fixed-point marginals
fp_marginals = {}
for name in S3_CORRS:
    crystal = refs[name]
    for _ in range(30):
        crystal = compose(crystal, refs[name])
    bp = born_probabilities(crystal).reshape(MASS_DIM, MASS_DIM)
    fp_marginals[name] = bp[0]

# Within-class similarity
print(f"  Within-class Born fidelity:")
trans = ["(01)", "(02)", "(12)"]
for i in range(len(trans)):
    for j in range(i+1, len(trans)):
        a, b = trans[i], trans[j]
        fid = sum((fp_marginals[a][k].item() * fp_marginals[b][k].item())**0.5
                  for k in range(MASS_DIM))
        print(f"    {a} ↔ {b}: {fid:.6f}")

cycles = ["(012)", "(021)"]
fid = sum((fp_marginals[cycles[0]][k].item() * fp_marginals[cycles[1]][k].item())**0.5
          for k in range(MASS_DIM))
print(f"    (012) ↔ (021): {fid:.6f}")

# Between-class similarity
print(f"\n  Between-class Born fidelity:")
fid_et = sum((fp_marginals["e"][k].item() * fp_marginals["(01)"][k].item())**0.5
             for k in range(MASS_DIM))
print(f"    e ↔ (01): {fid_et:.6f}")
fid_ec = sum((fp_marginals["e"][k].item() * fp_marginals["(012)"][k].item())**0.5
             for k in range(MASS_DIM))
print(f"    e ↔ (012): {fid_ec:.6f}")
fid_tc = sum((fp_marginals["(01)"][k].item() * fp_marginals["(012)"][k].item())**0.5
             for k in range(MASS_DIM))
print(f"    (01) ↔ (012): {fid_tc:.6f}")


# ═══════════════════════════════════════════════════════════════
#  The p parameter from the symmetrized basis
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Extract p* from fixed-point structure ---\n")

# In the symmetrized basis:
# The product state has m(i,j) = f(i)·f(j)
# Born(i,j) = p_i · p_j where p_i = Born probability of marginal
#
# For a DIAGONAL product state: Born(i,i) = p_i²
# For an OFF-DIAGONAL:          Born(i,j) = p_i · p_j
#
# p* is defined as the ratio that parametrizes the Born probabilities:
# Born(h_k) = 1/(H + p²) for hypotheses
# Born(θ)   = p²/(H + p²) for ignorance
# where p² = Born(θ)/Born(h) · H

for name in S3_CORRS:
    m = fp_marginals[name]
    h_mean = sum(m[k].item() for k in range(H)) / H
    theta = m[H].item()

    if h_mean > 1e-10:
        p_sq = theta / h_mean
        p = p_sq ** 0.5
    else:
        p = float('inf')

    # Also compute from the Born formula: Born(θ) = p²/(H + p²)
    # → p² = H · Born(θ) / (1 - Born(θ))
    p_sq_formula = H * theta / (1 - theta) if theta < 1 else float('inf')
    p_formula = p_sq_formula ** 0.5

    # Entropy of the marginal
    entropy = -sum(m[k].item() * math.log(max(m[k].item(), 1e-10)) for k in range(MASS_DIM))

    print(f"  {name:>6s}: p* = {p:.4f} (from ratio), "
          f"p* = {p_formula:.4f} (from formula), "
          f"Born(θ) = {theta:.4f}, H_entropy = {entropy:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Is p* related to the correlation entropy?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Correlation entropy vs p* ---\n")

for name, corr in S3_CORRS.items():
    # Correlation entropy
    c = corr / corr.sum()
    h_corr = -sum(c[i,j].item() * math.log(max(c[i,j].item(), 1e-10))
                  for i in range(H) for j in range(H))

    # p* from marginal
    m = fp_marginals[name]
    theta = m[H].item()
    p_sq = H * theta / (1 - theta) if theta < 1 else 0
    p = p_sq ** 0.5

    print(f"  {name:>6s}: corr_entropy = {h_corr:.4f}, p* = {p:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Seed dependence of p*
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Does p* depend on which seed is used? ---\n")

for name in ["e", "(01)", "(012)"]:
    p_vals = []
    for seed in range(20):
        crystal = Entangler(S3_CORRS[name], seed=seed).build().joint
        for _ in range(30):
            crystal = compose(crystal, Entangler(S3_CORRS[name], seed=seed).build().joint)
        bp = born_probabilities(crystal).reshape(MASS_DIM, MASS_DIM)
        m = bp[0]
        theta = m[H].item()
        h_mean = sum(m[k].item() for k in range(H)) / H
        p = (theta/h_mean)**0.5 if h_mean > 1e-10 else 0
        p_vals.append(p)

    mean_p = sum(p_vals)/len(p_vals)
    std_p = (sum((p-mean_p)**2 for p in p_vals)/len(p_vals))**0.5
    print(f"  {name:>6s}: p* = {mean_p:.4f} ± {std_p:.4f} (over 20 seeds)")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
