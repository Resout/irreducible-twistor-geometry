"""
Gauge directions for different correlation types.

Identity correlation: evidence samples (hᵢ,hᵢ) → diagonal entries
→ these project onto A₂ ROOT directions in the standard rep.

Cycle correlation (012): evidence samples (h₀,h₁), (h₁,h₂), (h₂,h₀)
→ off-diagonal entries → do these project onto WEIGHT directions?

The A₂ root system has:
  3 positive roots at 0, 2π/3, 4π/3 (the transposition directions)
  6 weights of the fundamental rep at π/6, π/2, 5π/6, 7π/6, 3π/2, 11π/6

If cycle crystals cluster at weight positions, the correlation type
determines whether the gauge field lives in the root lattice or the
weight lattice. This would connect the crystal's evidence geometry
to the Lie algebra's representation theory.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random as _random
from solver.algebra import (H, MASS_DIM, EPS_LOG, EPS_DIV, EPS_NORM)

torch.set_grad_enabled(False)


def standard_angle(joint):
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    signs_list = [1, -1, -1, -1, 1, 1]
    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T

    triv = sum(permute(joint, p) for p in perms) / 6
    sign = sum(s * permute(joint, p) for p, s in zip(perms, signs_list)) / 6
    std = joint - triv - sign

    U, S, Vh = torch.linalg.svd(std)
    if S[0].abs() < 1e-15:
        return None
    v = U[:, 0][:H]
    e1 = torch.tensor([1, -1, 0], dtype=torch.cfloat) / math.sqrt(2)
    e2 = torch.tensor([1, 1, -2], dtype=torch.cfloat) / math.sqrt(6)
    c1 = torch.dot(v.conj(), e1).real.item()
    c2 = torch.dot(v.conj(), e2).real.item()
    return math.atan2(c2, c1)


def entangle_and_measure(correlation, seed, n_steps=50, discount=0.3):
    """Build crystal, return final gauge angle."""
    _random.seed(seed)
    corr = correlation.float()
    flat = corr.reshape(-1)
    flat = flat / flat.sum().clamp(min=EPS_LOG)

    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j

    for step in range(n_steps):
        idx = _random.choices(range(H * H), weights=flat.tolist())[0]
        hi, hj = idx // H, idx % H
        ev_a = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev_b = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev_a[hi] = 0.55 + 0.06j
        ev_b[hj] = 0.55 + 0.06j
        for k in range(H):
            noise = _random.gauss(0, 0.02)
            if k != hi: ev_a[k] = 0.05 + noise * 1j
            if k != hj: ev_b[k] = 0.05 + noise * 1j
        ev_a[H] = (1.0 - ev_a[:H].real.sum()) - ev_a[:H].imag.sum() * 1j
        ev_b[H] = (1.0 - ev_b[:H].real.sum()) - ev_b[:H].imag.sum() * 1j

        joint_dim = MASS_DIM ** 2
        jH = joint_dim - 1
        joint_ev = torch.einsum("m,n->mn", ev_a, ev_b).reshape(joint_dim)
        ev_s = joint_ev[:jH] * discount
        ev_t = 1.0 + 0j - ev_s.sum()
        disc_ev = torch.cat([ev_s, ev_t.unsqueeze(0)])

        curr = joint.reshape(joint_dim)
        cs1, ct1 = curr[:jH], curr[jH:]
        cs2, ct2 = disc_ev[:jH], disc_ev[jH:]
        comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
        comb_t = ct1 * ct2
        K = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
        norm = 1.0 - K
        safe_norm = norm if abs(norm) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
        updated = torch.cat([comb_s, comb_t]) / safe_norm
        re_sum = updated.real.sum()
        if abs(re_sum) > EPS_DIV:
            updated = updated / re_sum
        joint = updated.reshape(MASS_DIM, MASS_DIM)

    return standard_angle(joint)


# Correlation types
identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
cycle_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)
inverse_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)
modular_corr = torch.tensor([[.4,.3,.3],[.3,.4,.3],[.3,.3,.4]], dtype=torch.float32)

corrs = {
    "identity": identity_corr,
    "cycle (012)": cycle_corr,
    "inverse (02)": inverse_corr,
    "modular": modular_corr,
}

n_seeds = 400
n_bins = 12
bw = 2 * math.pi / n_bins

print("=" * 80)
print("GAUGE DIRECTIONS BY CORRELATION TYPE")
print("=" * 80)


for name, corr in corrs.items():
    angles = []
    for seed in range(n_seeds):
        angle = entangle_and_measure(corr, seed=seed)
        if angle is not None:
            angles.append(angle)

    counts = [0] * n_bins
    for a in angles:
        b = int(((a % (2*math.pi)) / bw))
        if b >= n_bins: b = n_bins - 1
        counts[b] += 1

    # Entropy
    total = sum(counts)
    probs = [c / total for c in counts] if total > 0 else [1/n_bins]*n_bins
    entropy = -sum(p * math.log(p + 1e-15) for p in probs)
    max_entropy = math.log(n_bins)

    # Circular mean
    mean_sin = sum(math.sin(a) for a in angles) / len(angles)
    mean_cos = sum(math.cos(a) for a in angles) / len(angles)
    R = math.sqrt(mean_sin**2 + mean_cos**2)

    print(f"\n  {name}: n={len(angles)}, S/S_max={entropy/max_entropy:.3f}, R={R:.3f}")
    max_c = max(counts) if counts else 1
    for i in range(n_bins):
        bar = '█' * (counts[i] * 40 // max(max_c, 1))
        print(f"    [{i*bw/math.pi:.1f}π,{(i+1)*bw/math.pi:.1f}π) {counts[i]:4d} {bar}")


# ═══════════════════════════════════════════════════════════════
#  Expected directions for each correlation type
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("EXPECTED GAUGE DIRECTIONS")
print("="*80)

# Standard rep basis: e₁ = (1,-1,0)/√2, e₂ = (1,1,-2)/√6
# The diagonal entry (i,i) in the joint mass projects to:
#   direction = eᵢ in C³ (standard basis), projected to standard rep

# For identity: samples (0,0), (1,1), (2,2)
# (0,0) → h₀ direction = (1,0,0) → standard coords (1/√2, 1/√6)
# (1,1) → h₁ direction = (0,1,0) → standard coords (-1/√2, 1/√6)
# (2,2) → h₂ direction = (0,0,1) → standard coords (0, -2/√6)

print(f"\n  Diagonal (identity) evidence directions:")
for i, label in enumerate(["h₀⊗h₀", "h₁⊗h₁", "h₂⊗h₂"]):
    v = [0, 0, 0]
    v[i] = 1
    c1 = (v[0] - v[1]) / math.sqrt(2)
    c2 = (v[0] + v[1] - 2*v[2]) / math.sqrt(6)
    angle = math.atan2(c2, c1)
    print(f"    {label}: ({c1:.4f}, {c2:.4f}), angle = {angle/math.pi:.4f}π")

# For cycle (012): samples (0,1), (1,2), (2,0)
# (0,1) → h₀ left, h₁ right. The standard projection involves both.
# In the joint mass, the entry (0,1) doesn't directly project to a single
# direction in the standard rep of the LEFT variable.
# The left variable's standard projection comes from which ROW gets mass.
# Evidence (0,1) puts mass in row 0 → same as (0,0) for the left variable.

print(f"\n  Off-diagonal (cycle) evidence directions (LEFT variable):")
for pair, label in [((0,1), "h₀⊗h₁"), ((1,2), "h₁⊗h₂"), ((2,0), "h₂⊗h₀")]:
    i = pair[0]  # left variable = row
    v = [0, 0, 0]
    v[i] = 1
    c1 = (v[0] - v[1]) / math.sqrt(2)
    c2 = (v[0] + v[1] - 2*v[2]) / math.sqrt(6)
    angle = math.atan2(c2, c1)
    print(f"    {label}: row {i} → ({c1:.4f}, {c2:.4f}), angle = {angle/math.pi:.4f}π")

print(f"\n  Key insight: the LEFT variable's standard projection depends")
print(f"  only on which ROW gets evidence, not which column.")
print(f"  Identity samples rows 0,1,2 (from (0,0),(1,1),(2,2)).")
print(f"  Cycle samples rows 0,1,2 (from (0,1),(1,2),(2,0)).")
print(f"  SAME rows → SAME root directions for the left variable!")
print(f"  The RIGHT variable might differ, but the standard projection")
print(f"  of the joint mass is dominated by the ROW structure.")


# ═══════════════════════════════════════════════════════════════
#  Compare cluster positions between identity and cycle
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CLUSTER POSITION COMPARISON")
print("="*80)

for name, corr in [("identity", identity_corr), ("cycle (012)", cycle_corr)]:
    angles = []
    for seed in range(n_seeds):
        angle = entangle_and_measure(corr, seed=seed)
        if angle is not None:
            angles.append(angle)

    # Find cluster centers using 3-fold analysis
    # Fold to [0, 2π/3) to find the single cluster center
    folded = [(a % (2*math.pi/3)) for a in angles]
    mean_sin = sum(math.sin(a) for a in folded) / len(folded)
    mean_cos = sum(math.cos(a) for a in folded) / len(folded)
    center = math.atan2(mean_sin, mean_cos)

    print(f"\n  {name}:")
    print(f"    3-fold cluster center: {center/math.pi:.4f}π")
    print(f"    Unfolded: {center/math.pi:.4f}π, {center/math.pi + 2/3:.4f}π, {center/math.pi + 4/3:.4f}π")


# ═══════════════════════════════════════════════════════════════
#  The modular crystal: uniform evidence → what gauge direction?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("MODULAR (NEAR-UNIFORM) EVIDENCE")
print("="*80)

angles_mod = []
for seed in range(n_seeds):
    angle = entangle_and_measure(modular_corr, seed=seed)
    if angle is not None:
        angles_mod.append(angle)

# Is it MORE or LESS confined than identity?
counts_mod = [0] * n_bins
for a in angles_mod:
    b = int(((a % (2*math.pi)) / bw))
    if b >= n_bins: b = n_bins - 1
    counts_mod[b] += 1

total_mod = sum(counts_mod)
probs_mod = [c / total_mod for c in counts_mod]
entropy_mod = -sum(p * math.log(p + 1e-15) for p in probs_mod)

print(f"\n  Modular (near-uniform) entropy: {entropy_mod/math.log(n_bins):.3f}")
print(f"  Identity entropy: 0.838 (from earlier)")
print(f"  If modular is MORE uniform → higher entropy (less confined)")
print(f"  If modular is SAME → evidence geometry dominates over correlation")


print(f"\n\n{'='*80}")
print("WHAT THE CORRELATION COMPARISON REVEALS")
print("="*80)
