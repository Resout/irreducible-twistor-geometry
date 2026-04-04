"""
Is the entangler's phase transition genuine or stochastic?

The previous exploration found a classification flip from proportional
to inverse at θ ≈ 0.225π when continuously rotating the correlation
matrix. But the entangler uses stochastic sampling (50 steps).

Questions:
1. Does the transition sharpen with more steps? (genuine → yes)
2. Does the transition point move with step count? (finite-size → yes)
3. Is there a deterministic critical point in the infinite-step limit?
4. Can we derive the critical angle from H=3 constants?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from collections import Counter
from solver.algebra import H, MASS_DIM, schmidt_number, K_STAR, BORN_FLOOR, DELTA
from solver.crystals import Entangler, classify_relationship

torch.set_grad_enabled(False)


def make_rotated_corr(theta):
    """Correlation matrix rotating from identity (θ=0) to inverse (θ=π/2)."""
    c, s = math.cos(theta), math.sin(theta)
    return torch.tensor([
        [c**2, 0, s**2],
        [0, 1, 0],
        [s**2, 0, c**2],
    ], dtype=torch.float32)


# ═══════════════════════════════════════════════════════════════
#  Test 1: Transition sharpness vs step count
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 1: DOES THE TRANSITION SHARPEN WITH MORE STEPS?")
print("=" * 70)

step_counts = [20, 50, 100, 200, 500]
n_angles = 41  # fine grid
n_seeds = 20   # per angle per step count

for n_steps in step_counts:
    print(f"\n  ── n_steps = {n_steps} ──")
    print(f"  {'θ/π':>7s} {'avg_S':>7s} {'std_S':>7s} {'%prop':>6s} {'%inv':>6s} {'%mod':>6s}")

    transition_angles = []

    for i in range(n_angles):
        theta = i * math.pi / (2 * (n_angles - 1))
        corr = make_rotated_corr(theta)

        schmidts = []
        rels = []
        for seed in range(n_seeds):
            ent = Entangler(corr, seed=seed).build(n_steps=n_steps, discount=0.3)
            sn = schmidt_number(ent.joint)
            rel, _ = classify_relationship(ent.joint)
            schmidts.append(sn)
            rels.append(rel)

        avg_s = sum(schmidts) / len(schmidts)
        std_s = (sum((s - avg_s)**2 for s in schmidts) / len(schmidts)) ** 0.5
        counts = Counter(rels)
        pct_prop = counts.get("proportional", 0) / n_seeds * 100
        pct_inv = counts.get("inverse", 0) / n_seeds * 100
        pct_mod = counts.get("modular", 0) / n_seeds * 100

        if i % 4 == 0 or abs(pct_prop - pct_inv) < 30:
            print(f"  {theta/math.pi:>7.4f} {avg_s:>7.3f} {std_s:>7.3f} {pct_prop:>5.0f}% {pct_inv:>5.0f}% {pct_mod:>5.0f}%")

        # Track where majority flips
        if pct_inv > pct_prop and i > 0:
            transition_angles.append(theta / math.pi)

    if transition_angles:
        print(f"  Transition at θ/π ≈ {transition_angles[0]:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: High-statistics at the critical region
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 2: HIGH-STATISTICS AT THE CRITICAL REGION")
print("200 seeds at fine angular resolution near the transition")
print("=" * 70)

n_seeds_fine = 200
n_steps_fine = 50

# Fine grid around the transition
for theta_frac in [0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28]:
    theta = theta_frac * math.pi
    corr = make_rotated_corr(theta)

    schmidts = []
    rels = []
    for seed in range(n_seeds_fine):
        ent = Entangler(corr, seed=seed).build(n_steps=n_steps_fine, discount=0.3)
        sn = schmidt_number(ent.joint)
        rel, _ = classify_relationship(ent.joint)
        schmidts.append(sn)
        rels.append(rel)

    avg_s = sum(schmidts) / len(schmidts)
    std_s = (sum((s - avg_s)**2 for s in schmidts) / len(schmidts)) ** 0.5
    counts = Counter(rels)
    pct_prop = counts.get("proportional", 0) / n_seeds_fine * 100
    pct_inv = counts.get("inverse", 0) / n_seeds_fine * 100
    pct_mod = counts.get("modular", 0) / n_seeds_fine * 100

    # Entropy of classification distribution (low entropy = sharp transition)
    probs = [pct_prop/100, pct_inv/100, pct_mod/100]
    entropy = -sum(p * math.log(p + 1e-10) for p in probs if p > 0)

    print(f"  θ/π={theta_frac:.2f}: S={avg_s:.3f}±{std_s:.3f}, "
          f"prop={pct_prop:5.1f}% inv={pct_inv:5.1f}% mod={pct_mod:5.1f}%, "
          f"H={entropy:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 3: The correlation matrix at the transition
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 3: WHAT IS SPECIAL ABOUT THE TRANSITION CORRELATION?")
print("=" * 70)

theta_c = 0.225 * math.pi
c, s = math.cos(theta_c), math.sin(theta_c)
minority_frac = s**2

print(f"\n  Critical angle θ_c ≈ 0.225π")
print(f"  cos²(θ_c) = {c**2:.6f} (majority)")
print(f"  sin²(θ_c) = {s**2:.6f} (minority)")
print(f"  Minority fraction: {minority_frac:.4f}")
print(f"  Ratio minority/majority: {s**2/c**2:.4f}")

# What H=3 constants are near this?
print(f"\n  H=3 constants for comparison:")
print(f"    K* = {K_STAR:.6f} = 7/30")
print(f"    1-K* = {1-K_STAR:.6f}")
print(f"    BORN_FLOOR = {BORN_FLOOR:.6f} = 1/27")
print(f"    1/H = {1/H:.6f}")
print(f"    1/(H+1) = {1/(H+1):.6f} = 0.25")
print(f"    1/(H²) = {1/H**2:.6f}")
print(f"    (H-1)/(H+1) = {(H-1)/(H+1):.6f}")
print(f"    (H-1)/H² = {(H-1)/H**2:.6f}")
print(f"    H/(H²+1) = {H/(H**2+1):.6f} = 3/10")
print(f"    DELTA = {DELTA:.6f}")

# The conflict K at the transition
corr_c = make_rotated_corr(theta_c)
print(f"\n  Correlation at transition:")
print(f"    {corr_c}")

# Compute theoretical K for this correlation
# For row 0: P(agree) = c⁴ + s⁴, P(disagree) = 2c²s²
p_agree_row0 = c**4 + s**4
p_disagree_row0 = 2 * c**2 * s**2
# Row 1: always (1,1) → always agrees
# Row 2: same as row 0

# The entangler samples (hi, hj) from the correlation.
# Conflict in DS comes from disagreeing singletons.
# K = Σ_i≠j m1(hi)*m2(hj) — the cross terms between different hypotheses.
# For a single DS step with evidence from this correlation:
# When (hi, hj) = (0,0): evidence agrees with identity → reinforces
# When (hi, hj) = (0,2): evidence agrees with inverse → conflicts with identity

print(f"\n  P(identity-supporting evidence | row 0): {c**2:.4f}")
print(f"  P(inverse-supporting evidence | row 0):  {s**2:.4f}")
print(f"  P(identity-supporting evidence | row 2): {s**2:.4f}")
print(f"  P(inverse-supporting evidence | row 2):  {c**2:.4f}")
print(f"  P(identity-supporting | all rows):       {(c**2 + 1 + s**2)/3:.4f} = {(1 + c**2 - s**2 + 1)/3:.4f}")

# The key: row 1 always supports identity AND inverse equally (mid→mid).
# The battle is between rows 0 and 2.
# For identity: row 0 gives c²=0.586, row 2 gives s²=0.414. Average: 0.500.
# For inverse:  row 0 gives s²=0.414, row 2 gives c²=0.586. Average: 0.500.
# It's PERFECTLY BALANCED across rows 0 and 2!
# The transition happens NOT where the minority reaches some threshold,
# but where the average across rows reaches 50%.

identity_support = (c**2 + s**2) / 2  # average of row 0 and row 2 identity support
inverse_support = (s**2 + c**2) / 2   # average of row 0 and row 2 inverse support
print(f"\n  Average identity support (rows 0,2): {identity_support:.4f}")
print(f"  Average inverse support (rows 0,2):  {inverse_support:.4f}")
print(f"  → Always 50/50! The transition is NOT at average balance.")
print(f"  → The transition is a SYMMETRY BREAKING by the entangler.")


# ═══════════════════════════════════════════════════════════════
#  Test 4: Track the entangler's internal dynamics
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 4: ENTANGLER DYNAMICS AT THE TRANSITION")
print("Track Schmidt after each step to see how the commitment builds")
print("=" * 70)

import random as _random

def entangler_trajectory(corr, n_steps=100, discount=0.3, seed=42):
    """Run entangler step by step, recording Schmidt after each."""
    _random.seed(seed)
    from solver.algebra import EPS_LOG, EPS_DIV, EPS_NORM

    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j

    flat = corr.reshape(-1)
    flat = flat / flat.sum().clamp(min=EPS_LOG)

    trajectory = []

    for step in range(n_steps):
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

        sn = schmidt_number(joint)
        rel, _ = classify_relationship(joint)
        trajectory.append({"step": step, "schmidt": sn, "rel": rel,
                          "hi": hi, "hj": hj, "K": abs(K.item()) if isinstance(K, torch.Tensor) else abs(K)})

    return trajectory


# Run at transition point with two seeds that give different outcomes
print(f"\n  At θ/π = 0.225 (transition point):")

corr_transition = make_rotated_corr(0.225 * math.pi)

# Find seeds that give different classifications
prop_seed = None
inv_seed = None
for seed in range(100):
    ent = Entangler(corr_transition, seed=seed).build(n_steps=50)
    rel, _ = classify_relationship(ent.joint)
    if rel == "proportional" and prop_seed is None:
        prop_seed = seed
    if rel == "inverse" and inv_seed is None:
        inv_seed = seed
    if prop_seed is not None and inv_seed is not None:
        break

print(f"  Seed {prop_seed} → proportional, Seed {inv_seed} → inverse")

for seed, label in [(prop_seed, "→ proportional"), (inv_seed, "→ inverse")]:
    traj = entangler_trajectory(corr_transition, n_steps=100, seed=seed)
    print(f"\n  Seed {seed} ({label}):")
    print(f"  {'step':>5s} {'Schmidt':>8s} {'rel':>13s} {'hi,hj':>6s} {'K':>8s}")
    for t in traj:
        if t["step"] % 10 == 0 or t["step"] < 5 or t["step"] == 49:
            print(f"  {t['step']:>5d} {t['schmidt']:>8.3f} {t['rel']:>13s} "
                  f"({t['hi']},{t['hj']}) {t['K']:>8.4f}")

    # When does the classification LOCK IN?
    final_rel = traj[-1]["rel"]
    lock_step = 0
    for i in range(len(traj) - 1, -1, -1):
        if traj[i]["rel"] != final_rel:
            lock_step = i + 1
            break
    print(f"  Classification locks in at step {lock_step}")

    # Count identity vs inverse evidence
    id_evidence = sum(1 for t in traj if t["hi"] == t["hj"])
    inv_evidence = sum(1 for t in traj if (t["hi"], t["hj"]) in [(0,2),(2,0)])
    mid_evidence = sum(1 for t in traj if t["hi"] == 1 and t["hj"] == 1)
    print(f"  Evidence: {id_evidence} identity-supporting, {inv_evidence} inverse-supporting, "
          f"{mid_evidence} mid-mid (neutral)")


# ═══════════════════════════════════════════════════════════════
#  Test 5: Is the transition at a universal point?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 5: IS THE TRANSITION ANGLE UNIVERSAL?")
print("Try different parameterizations of the id→inv interpolation")
print("=" * 70)

# Parameterization 1: linear interpolation (already tested)
# Parameterization 2: different row structure

# Rotate only row 0, keep row 2 fixed as identity
print("\n  Rotating ONLY row 0 (row 2 stays identity):")
for theta_frac in [0.00, 0.10, 0.20, 0.30, 0.40, 0.50]:
    theta = theta_frac * math.pi
    c, s = math.cos(theta), math.sin(theta)
    corr = torch.tensor([
        [c**2, 0, s**2],
        [0, 1, 0],
        [0, 0, 1],  # row 2 stays identity
    ], dtype=torch.float32)

    rels = []
    schmidts = []
    for seed in range(50):
        ent = Entangler(corr, seed=seed).build()
        rel, _ = classify_relationship(ent.joint)
        rels.append(rel)
        schmidts.append(schmidt_number(ent.joint))
    counts = Counter(rels)
    avg_s = sum(schmidts) / len(schmidts)
    print(f"  θ/π={theta_frac:.2f}: S={avg_s:.3f}, "
          f"prop={counts.get('proportional',0):2d}/50, "
          f"inv={counts.get('inverse',0):2d}/50, "
          f"mod={counts.get('modular',0):2d}/50")

# Rotate all three rows independently
print("\n  All rows rotate together (full matrix rotation):")
for theta_frac in [0.00, 0.10, 0.15, 0.20, 0.22, 0.24, 0.25, 0.30, 0.40, 0.50]:
    theta = theta_frac * math.pi
    c, s = math.cos(theta), math.sin(theta)
    corr = torch.tensor([
        [c**2, 0, s**2],
        [s**2, 0, c**2],  # row 1 ALSO rotates (not fixed at [0,1,0])
        [s**2, 0, c**2],
    ], dtype=torch.float32)

    rels = []
    schmidts = []
    for seed in range(50):
        ent = Entangler(corr, seed=seed).build()
        rel, _ = classify_relationship(ent.joint)
        rels.append(rel)
        schmidts.append(schmidt_number(ent.joint))
    counts = Counter(rels)
    avg_s = sum(schmidts) / len(schmidts)
    print(f"  θ/π={theta_frac:.2f}: S={avg_s:.3f}, "
          f"prop={counts.get('proportional',0):2d}/50, "
          f"inv={counts.get('inverse',0):2d}/50, "
          f"mod={counts.get('modular',0):2d}/50")


# ═══════════════════════════════════════════════════════════════
#  Test 6: The critical K
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 6: CONFLICT K AT THE TRANSITION")
print("What is K when the entangler switches?")
print("=" * 70)

# For each angle, compute the average K across the trajectory
for theta_frac in [0.00, 0.10, 0.20, 0.225, 0.25, 0.30, 0.40, 0.50]:
    theta = theta_frac * math.pi
    corr = make_rotated_corr(theta)

    all_Ks = []
    for seed in range(10):
        traj = entangler_trajectory(corr, n_steps=50, seed=seed)
        Ks = [t["K"] for t in traj]
        all_Ks.extend(Ks)

    avg_K = sum(all_Ks) / len(all_Ks)
    std_K = (sum((k - avg_K)**2 for k in all_Ks) / len(all_Ks)) ** 0.5

    c = math.cos(theta)
    s = math.sin(theta)
    # Theoretical single-step K for this correlation
    # When evidence is (0,0): singleton is at position 0. K = m(other singletons) * evidence(0)
    # This is complex — just report empirical
    print(f"  θ/π={theta_frac:.3f}: avg K = {avg_K:.4f} ± {std_K:.4f}, "
          f"minority={s**2:.3f}")


print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
