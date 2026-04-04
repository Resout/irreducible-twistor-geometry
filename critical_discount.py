"""
critical_discount.py — Investigate entangler behavior at special discount values.

Key question: where does the discount parameter create special structure?
- K* = 7/30 ≈ 0.2333 is equilibrium conflict
- DELTA = -log(1 - K*) ≈ 0.266 is the structural filter rate
- Blow-up expected near discount = 0.5

Standalone computation. Imports from solver.algebra and solver.crystals.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch
import numpy as np
from solver.algebra import (
    H, MASS_DIM, K_STAR, DELTA, BORN_FLOOR,
    born_probabilities, schmidt_number, born_entropy,
)
from solver.crystals import Entangler

IDENTITY = torch.eye(H, dtype=torch.float32)


def run_entangler(corr, seed=42, n_steps=50, discount=0.3, track_K=False):
    """Run entangler, optionally tracking K at each step."""
    import random as _random
    from solver.algebra import EPS_LOG, EPS_NORM, EPS_DIV

    ent = Entangler(corr, seed=seed)

    if not track_K:
        ent.build(n_steps=n_steps, discount=discount)
        return ent, []

    # Reproduce build() but record K at each step
    _random.seed(ent.seed)
    flat = ent.correlation.reshape(-1)
    flat = flat / flat.sum().clamp(min=EPS_LOG)
    Ks = []

    for _ in range(n_steps):
        idx = _random.choices(range(H * H), weights=flat.tolist())[0]
        hi = idx // H
        hj = idx % H

        ev_a = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev_b = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev_a[hi] = 0.55 + 0.06j
        ev_b[hj] = 0.55 + 0.06j
        for k in range(H):
            noise = _random.gauss(0, 0.02)
            if k != hi:
                ev_a[k] = 0.05 + noise * 1j
            if k != hj:
                ev_b[k] = 0.05 + noise * 1j
        ev_a[H] = (1.0 - ev_a[:H].real.sum()) - ev_a[:H].imag.sum() * 1j
        ev_b[H] = (1.0 - ev_b[:H].real.sum()) - ev_b[:H].imag.sum() * 1j

        joint_dim = MASS_DIM ** 2
        jH = joint_dim - 1
        joint_ev = torch.einsum("m,n->mn", ev_a, ev_b).reshape(joint_dim)
        ev_s = joint_ev[:jH] * discount
        ev_t = 1.0 + 0j - ev_s.sum()
        disc_ev = torch.cat([ev_s, ev_t.unsqueeze(0)])

        curr = ent.joint.reshape(joint_dim)
        cs1, ct1 = curr[:jH], curr[jH:]
        cs2, ct2 = disc_ev[:jH], disc_ev[jH:]

        comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
        comb_t = ct1 * ct2
        K = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
        Ks.append(abs(K.item()))
        norm = 1.0 - K
        safe_norm = norm if abs(norm) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
        updated = torch.cat([comb_s, comb_t]) / safe_norm

        re_sum = updated.real.sum()
        if abs(re_sum) > EPS_DIV:
            updated = updated / re_sum

        ent.joint = updated.reshape(MASS_DIM, MASS_DIM)

    return ent, Ks


def marginal_born_theta(joint):
    """Born(theta) of the row-marginal."""
    bp = born_probabilities(joint)
    marginal = bp.sum(dim=1)  # [4]
    return marginal[-1].item()


# ═══════════════════════════════════════════════════════════════
#  1. Fine-grained discount sweep
# ═══════════════════════════════════════════════════════════════
print("=" * 72)
print("  1. DISCOUNT SWEEP  (0.05 to 0.95, identity crystal, 50 steps)")
print("=" * 72)
print(f"{'discount':>10} {'avg_K(10-40)':>14} {'Schmidt':>10} {'Born(θ)':>10} {'K - K*':>12}")
print("-" * 72)

sweep_discounts = [d / 100 for d in range(5, 96, 5)]
sweep_results = {}

for d in sweep_discounts:
    ent, Ks = run_entangler(IDENTITY, seed=42, n_steps=50, discount=d, track_K=True)
    avg_K = np.mean(Ks[10:41]) if len(Ks) > 40 else np.mean(Ks[10:])
    sc = schmidt_number(ent.joint)
    bt = marginal_born_theta(ent.joint)
    sweep_results[d] = {"avg_K": avg_K, "schmidt": sc, "born_theta": bt, "Ks": Ks}
    print(f"{d:10.3f} {avg_K:14.6f} {sc:10.4f} {bt:10.6f} {avg_K - K_STAR:12.6f}")

# Find discount where avg_K ≈ K*
best_d = min(sweep_results, key=lambda d: abs(sweep_results[d]["avg_K"] - K_STAR))
print(f"\nClosest to K* = {K_STAR:.6f}: discount = {best_d:.3f}  (avg_K = {sweep_results[best_d]['avg_K']:.6f})")


# ═══════════════════════════════════════════════════════════════
#  2. Discount = K* (0.2333): self-referential crystal
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print(f"  2. DISCOUNT = K* = {K_STAR:.6f}")
print("=" * 72)

ent_kstar, Ks_kstar = run_entangler(IDENTITY, seed=42, n_steps=50, discount=K_STAR, track_K=True)
sc_kstar = schmidt_number(ent_kstar.joint)
bt_kstar = marginal_born_theta(ent_kstar.joint)
avg_K_kstar = np.mean(Ks_kstar[10:41])

print(f"  Schmidt number:     {sc_kstar:.6f}")
print(f"  Born(θ) marginal:   {bt_kstar:.6f}")
print(f"  Avg K (steps 10-40): {avg_K_kstar:.6f}")
print(f"  K* itself:          {K_STAR:.6f}")
print(f"  discount = K* ?     discount == conflict  →  self-referential")

# Check if discount = K* means the entangler IS its own conflict
print(f"\n  K trajectory (first 15 steps):")
for i, k in enumerate(Ks_kstar[:15]):
    marker = " ← K*" if abs(k - K_STAR) < 0.01 else ""
    print(f"    step {i:2d}: K = {k:.6f}{marker}")

# Joint mass structure
bp = born_probabilities(ent_kstar.joint)
print(f"\n  Born distribution (4x4 joint):")
for i in range(MASS_DIM):
    row = "    " + " ".join(f"{bp[i,j].item():.4f}" for j in range(MASS_DIM))
    print(row)


# ═══════════════════════════════════════════════════════════════
#  3. Discount = 1 - K* (0.7667): complement
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print(f"  3. DISCOUNT = 1 - K* = {1 - K_STAR:.6f}")
print("=" * 72)

d_comp = 1.0 - K_STAR
ent_comp, Ks_comp = run_entangler(IDENTITY, seed=42, n_steps=50, discount=d_comp, track_K=True)
sc_comp = schmidt_number(ent_comp.joint)
bt_comp = marginal_born_theta(ent_comp.joint)
avg_K_comp = np.mean(Ks_comp[10:41])

print(f"  Schmidt number:     {sc_comp:.6f}")
print(f"  Born(θ) marginal:   {bt_comp:.6f}")
print(f"  Avg K (steps 10-40): {avg_K_comp:.6f}")

# Check ratio to K*
if avg_K_comp > 1e-10:
    print(f"  avg_K / K*:         {avg_K_comp / K_STAR:.6f}")
print(f"  Schmidt / H:        {sc_comp / H:.6f}")

# Is K at complement = 1 - K at K*?
print(f"\n  K(d=K*) + K(d=1-K*) = {avg_K_kstar + avg_K_comp:.6f}  (1.0 would be complement)")


# ═══════════════════════════════════════════════════════════════
#  4. Discount = DELTA (0.266): structural filter rate
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print(f"  4. DISCOUNT = DELTA = {DELTA:.6f}")
print("=" * 72)

ent_delta, Ks_delta = run_entangler(IDENTITY, seed=42, n_steps=50, discount=DELTA, track_K=True)
sc_delta = schmidt_number(ent_delta.joint)
bt_delta = marginal_born_theta(ent_delta.joint)
avg_K_delta = np.mean(Ks_delta[10:41])

print(f"  Schmidt number:     {sc_delta:.6f}")
print(f"  Born(θ) marginal:   {bt_delta:.6f}")
print(f"  Avg K (steps 10-40): {avg_K_delta:.6f}")
print(f"  -log(1 - avg_K):    {-np.log(1 - avg_K_delta):.6f}  (should be DELTA={DELTA:.6f} if self-consistent)")

# Check: does discount=DELTA produce K whose DELTA equals DELTA?
K_from_delta = avg_K_delta
delta_from_K = -np.log(1 - K_from_delta)
print(f"\n  Self-consistency check:")
print(f"    discount = DELTA = {DELTA:.6f}")
print(f"    produced avg_K   = {K_from_delta:.6f}")
print(f"    -log(1-K)        = {delta_from_K:.6f}")
print(f"    ratio to DELTA   = {delta_from_K / DELTA:.6f}")


# ═══════════════════════════════════════════════════════════════
#  5. Schmidt variability across seeds
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("  5. SCHMIDT VARIABILITY (20 seeds per discount)")
print("=" * 72)
print(f"{'discount':>10} {'mean_Sch':>10} {'std_Sch':>10} {'CV':>10} {'min':>10} {'max':>10}")
print("-" * 72)

cv_results = {}
test_discounts = [0.05, 0.10, 0.15, 0.20, K_STAR, 0.25, DELTA, 0.30, 0.35, 0.40, 0.45, 0.50]

for d in test_discounts:
    schmidts = []
    for seed in range(20):
        ent = Entangler(IDENTITY, seed=seed).build(n_steps=50, discount=d)
        schmidts.append(schmidt_number(ent.joint))
    arr = np.array(schmidts)
    mean_s = arr.mean()
    std_s = arr.std()
    cv = std_s / mean_s if mean_s > 1e-10 else float('inf')
    cv_results[d] = {"mean": mean_s, "std": std_s, "cv": cv}
    print(f"{d:10.4f} {mean_s:10.4f} {std_s:10.4f} {cv:10.6f} {arr.min():10.4f} {arr.max():10.4f}")

most_stable = min(cv_results, key=lambda d: cv_results[d]["cv"])
print(f"\nMost stable discount: {most_stable:.4f}  (CV = {cv_results[most_stable]['cv']:.6f})")
print(f"  K* = {K_STAR:.4f}, DELTA = {DELTA:.4f}")


# ═══════════════════════════════════════════════════════════════
#  6. Blow-up boundary: where does Schmidt diverge?
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("  6. BLOW-UP BOUNDARY (fine scan around 0.5)")
print("=" * 72)

fine_discounts = [d / 1000 for d in range(400, 601, 10)]
print(f"{'discount':>10} {'Schmidt':>12} {'avg_K':>12} {'max_K':>12} {'status':>10}")
print("-" * 72)

blowup_d = None
for d in fine_discounts:
    try:
        ent, Ks = run_entangler(IDENTITY, seed=42, n_steps=50, discount=d, track_K=True)
        sc = schmidt_number(ent.joint)
        avg_K = np.mean(Ks[10:41]) if len(Ks) > 40 else np.mean(Ks)
        max_K = max(Ks)

        if sc > 100 or np.isnan(sc) or np.isinf(sc):
            status = "BLOW-UP"
            if blowup_d is None:
                blowup_d = d
        elif max_K > 0.9:
            status = "DANGER"
        elif sc > 3.5:
            status = "HIGH"
        else:
            status = "ok"

        print(f"{d:10.3f} {sc:12.4f} {avg_K:12.6f} {max_K:12.6f} {status:>10}")
    except Exception as e:
        print(f"{d:10.3f} {'ERROR':>12} {str(e)[:40]}")
        if blowup_d is None:
            blowup_d = d

# Even finer scan around the boundary
if blowup_d is not None:
    print(f"\n  Approximate blow-up at discount ≈ {blowup_d:.3f}")
    print(f"  1/(H-1) = 1/2 = 0.500")
    print(f"  Refining around {blowup_d:.3f}...")

    lo, hi = blowup_d - 0.02, blowup_d + 0.02
    fine2 = [lo + i * (hi - lo) / 20 for i in range(21)]
    print(f"\n{'discount':>10} {'Schmidt':>12} {'status':>10}")
    print("-" * 40)
    for d in fine2:
        if d <= 0 or d >= 1:
            continue
        try:
            ent = Entangler(IDENTITY, seed=42).build(n_steps=50, discount=d)
            sc = schmidt_number(ent.joint)
            status = "BLOW-UP" if sc > 100 or np.isnan(sc) else "ok"
            print(f"{d:10.5f} {sc:12.4f} {status:>10}")
        except:
            print(f"{d:10.5f} {'ERROR':>12}")
else:
    print("\n  No blow-up detected in range [0.40, 0.60]")
    print("  Extending scan to 0.95...")
    for d in [0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]:
        try:
            ent, Ks = run_entangler(IDENTITY, seed=42, n_steps=50, discount=d, track_K=True)
            sc = schmidt_number(ent.joint)
            max_K = max(Ks)
            status = "BLOW-UP" if sc > 100 or np.isnan(sc) else "ok"
            print(f"  d={d:.2f}: Schmidt={sc:.4f}, max_K={max_K:.6f}  {status}")
        except Exception as e:
            print(f"  d={d:.2f}: ERROR - {e}")


# ═══════════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("  SUMMARY")
print("=" * 72)
print(f"  K* = {K_STAR:.6f} = 7/30")
print(f"  DELTA = {DELTA:.6f} = -log(1 - K*)")
print(f"  1 - K* = {1-K_STAR:.6f}")
print(f"  1/(H-1) = {1/(H-1):.6f}")
print(f"")
print(f"  Discount closest to K* equilibrium: {best_d:.3f}")
print(f"  Most stable Schmidt (lowest CV):    {most_stable:.4f}")
if blowup_d:
    print(f"  Blow-up boundary:                   {blowup_d:.3f}")
print(f"")
print(f"  d=K*:    Schmidt={sc_kstar:.4f}, avg_K={avg_K_kstar:.6f}")
print(f"  d=DELTA: Schmidt={sc_delta:.4f}, avg_K={avg_K_delta:.6f}")
print(f"  d=1-K*:  Schmidt={sc_comp:.4f},  avg_K={avg_K_comp:.6f}")
