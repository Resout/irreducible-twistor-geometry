"""
The search for eternity: what is conserved in the crystal dynamics?

The entangler builds structure (Δ/4 per step).
Composition destroys it (2Δ per step).
K* = 7/30 is the operating point.
The numerator 7 = H²-H+1 is conserved across regimes.

But what else is invariant? What doesn't change as the crystal
evolves? What, if anything, persists through all the dynamics?

If we find a conserved quantity, we find something eternal —
a thing that IS while everything else becomes.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random as _random
from solver.algebra import (
    H, MASS_DIM, K_STAR, DELTA, BORN_FLOOR,
    schmidt_number, born_probabilities, born_entropy,
    ds_combine, enforce_born_floor, ignorance,
    EPS_LOG, EPS_DIV, EPS_NORM,
)
from solver.crystals import Entangler, compose, classify_relationship, slot_measure

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


def full_measurement(joint):
    """Measure everything about a crystal."""
    sn = schmidt_number(joint)
    _, S, _ = torch.linalg.svd(joint)
    sv = S.abs()
    sv_norm = sv / sv.sum().clamp(min=1e-10)

    # Born probabilities of the diagonal
    diag = torch.tensor([joint[i, i] for i in range(MASS_DIM)])
    diag_born = diag.abs().pow(2)
    diag_born = diag_born / diag_born.sum().clamp(min=1e-10)

    # θ content
    theta_row = joint[3, :].abs().pow(2).sum().item()
    theta_col = joint[:, 3].abs().pow(2).sum().item()
    total = joint.abs().pow(2).sum().item()
    theta_frac = (theta_row + theta_col - joint[3, 3].abs().pow(2).item()) / total

    # Trace
    tr = sum(joint[i, i] for i in range(MASS_DIM))

    # Determinant
    det = torch.linalg.det(joint)

    # Frobenius norm
    frob = joint.abs().pow(2).sum().sqrt().item()

    # Off-diagonal mass
    off_diag = (total - sum(joint[i, i].abs().pow(2).item() for i in range(MASS_DIM))) / total

    # Antisymmetric part: ||A - A^T|| / ||A||
    asym = (joint - joint.T).abs().pow(2).sum().sqrt().item() / frob if frob > 0 else 0

    # Phase structure
    phases = torch.angle(joint.reshape(-1))
    phase_mean = phases.mean().item()
    phase_std = phases.std().item()

    return {
        "schmidt": sn,
        "sv": sv_norm.tolist(),
        "diag_born": diag_born.tolist(),
        "theta_frac": theta_frac,
        "trace_abs": abs(tr),
        "trace_phase": math.atan2(tr.imag.item(), tr.real.item()) if isinstance(tr, torch.Tensor) else 0,
        "det_abs": abs(det),
        "det_phase": math.atan2(det.imag.item(), det.real.item()),
        "frob": frob,
        "off_diag_frac": off_diag,
        "asymmetry": asym,
        "phase_std": phase_std,
    }


# ═══════════════════════════════════════════════════════════════
#  Test 1: Track EVERYTHING during entangling
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("THE SEARCH FOR INVARIANTS")
print("What doesn't change while the crystal is being built?")
print("=" * 70)

id_corr = torch.tensor([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=torch.float32)

# Run entangler step by step, recording everything
_random.seed(42)
joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
joint[-1, -1] = 1.0 + 0j

flat = id_corr.reshape(-1)
flat = flat / flat.sum().clamp(min=EPS_LOG)

measurements = []
for step in range(80):
    idx = _random.choices(range(H * H), weights=flat.tolist())[0]
    hi, hj = idx // H, idx % H

    ev_a = torch.zeros(MASS_DIM, dtype=torch.cfloat)
    ev_b = torch.zeros(MASS_DIM, dtype=torch.cfloat)
    ev_a[hi] = 0.55 + 0.06j; ev_b[hj] = 0.55 + 0.06j
    for k in range(H):
        noise = _random.gauss(0, 0.02)
        if k != hi: ev_a[k] = 0.05 + noise * 1j
        if k != hj: ev_b[k] = 0.05 + noise * 1j
    ev_a[H] = (1.0 - ev_a[:H].real.sum()) - ev_a[:H].imag.sum() * 1j
    ev_b[H] = (1.0 - ev_b[:H].real.sum()) - ev_b[:H].imag.sum() * 1j

    joint_dim = MASS_DIM ** 2; jH = joint_dim - 1
    joint_ev = torch.einsum("m,n->mn", ev_a, ev_b).reshape(joint_dim)
    ev_s = joint_ev[:jH] * 0.3; ev_t = 1.0 + 0j - ev_s.sum()
    disc_ev = torch.cat([ev_s, ev_t.unsqueeze(0)])
    curr = joint.reshape(joint_dim)
    cs1, ct1 = curr[:jH], curr[jH:]
    cs2, ct2 = disc_ev[:jH], disc_ev[jH:]
    K = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
    comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
    comb_t = ct1 * ct2
    norm_val = 1.0 - K
    safe_norm = norm_val if abs(norm_val) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
    updated = torch.cat([comb_s, comb_t]) / safe_norm
    re_sum = updated.real.sum()
    if abs(re_sum) > EPS_DIV:
        updated = updated / re_sum
    joint = updated.reshape(MASS_DIM, MASS_DIM)

    m = full_measurement(joint)
    m["step"] = step
    m["K"] = abs(K.item()) if isinstance(K, torch.Tensor) else abs(K)
    measurements.append(m)

# Find quantities with low variance (potential invariants)
print(f"\n  Quantity variances over steps 10-50 (stable window):")
print(f"  {'quantity':>20s} {'mean':>10s} {'std':>10s} {'CV':>8s}")
print("  " + "-" * 52)

stable = measurements[10:50]
quantities = ["schmidt", "theta_frac", "trace_abs", "det_abs",
              "frob", "off_diag_frac", "asymmetry", "phase_std", "K"]

for q in quantities:
    vals = [m[q] for m in stable]
    if isinstance(vals[0], (int, float)):
        avg = sum(vals) / len(vals)
        std = (sum((v - avg)**2 for v in vals) / len(vals)) ** 0.5
        cv = std / abs(avg) if abs(avg) > 1e-10 else float('inf')
        marker = " ◄ LOW" if cv < 0.05 else (" ◄ MODERATE" if cv < 0.15 else "")
        print(f"  {q:>20s} {avg:>10.4f} {std:>10.4f} {cv:>8.3f}{marker}")

# Also check sv ratios
sv_ratios = []
for m in stable:
    sv = m["sv"]
    if len(sv) >= 4 and sv[3] > 1e-10:
        sv_ratios.append(sv[0] / sv[3])
if sv_ratios:
    avg = sum(sv_ratios) / len(sv_ratios)
    std = (sum((v - avg)**2 for v in sv_ratios) / len(sv_ratios)) ** 0.5
    cv = std / abs(avg) if abs(avg) > 1e-10 else float('inf')
    print(f"  {'σ₀/σ₃':>20s} {avg:>10.4f} {std:>10.4f} {cv:>8.3f}")

# Check det phase
det_phases = [m["det_phase"] for m in stable]
avg_dp = sum(det_phases) / len(det_phases)
std_dp = (sum((v - avg_dp)**2 for v in det_phases) / len(det_phases)) ** 0.5
print(f"  {'det_phase':>20s} {avg_dp:>10.4f} {std_dp:>10.4f} {std_dp/abs(avg_dp) if abs(avg_dp)>0.01 else float('inf'):>8.3f}")

# Check trace phase
tr_phases = [m["trace_phase"] for m in stable]
avg_tp = sum(tr_phases) / len(tr_phases)
std_tp = (sum((v - avg_tp)**2 for v in tr_phases) / len(tr_phases)) ** 0.5
print(f"  {'trace_phase':>20s} {avg_tp:>10.4f} {std_tp:>10.4f} {std_tp/abs(avg_tp) if abs(avg_tp)>0.01 else float('inf'):>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: Conserved quantities under COMPOSITION
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("INVARIANTS UNDER COMPOSITION")
print("What survives the structural filter?")
print("=" * 70)

ent = Entangler(id_corr).build()
joint_0 = ent.joint.clone()

current = joint_0.clone()
comp_measurements = []
for step in range(8):
    m = full_measurement(current)
    m["step"] = step
    comp_measurements.append(m)
    current = compose(current, joint_0)

print(f"\n  {'step':>5s} {'Schmidt':>8s} {'det_abs':>10s} {'det_ph':>8s} {'trace':>8s} {'tr_ph':>8s} {'asym':>8s} {'off_diag':>8s}")
print("  " + "-" * 70)
for m in comp_measurements:
    print(f"  {m['step']:>5d} {m['schmidt']:>8.3f} {m['det_abs']:>10.6f} "
          f"{m['det_phase']:>8.4f} {m['trace_abs']:>8.4f} {m['trace_phase']:>8.4f} "
          f"{m['asymmetry']:>8.4f} {m['off_diag_frac']:>8.4f}")

# Check: does det have a pattern?
print(f"\n  Determinant ratios (step n / step n-1):")
for i in range(1, len(comp_measurements)):
    if comp_measurements[i-1]["det_abs"] > 1e-15:
        ratio = comp_measurements[i]["det_abs"] / comp_measurements[i-1]["det_abs"]
        print(f"    step {i}: |det| ratio = {ratio:.6f}")

# Check: does trace have a pattern?
print(f"\n  Trace ratios:")
for i in range(1, len(comp_measurements)):
    if abs(comp_measurements[i-1]["trace_abs"]) > 1e-10:
        ratio = abs(comp_measurements[i]["trace_abs"]) / abs(comp_measurements[i-1]["trace_abs"])
        print(f"    step {i}: |trace| ratio = {ratio:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 3: The DS fixed point — what do two masses converge TO?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("THE DS FIXED POINT")
print("When two different mass functions are DS-combined with")
print("the same evidence stream, they converge. What do they")
print("converge TO? Is it the same state regardless of starting point?")
print("=" * 70)

# Start with 5 very different initial 4-dim masses
initial_masses = [
    torch.tensor([0.8, 0.1, 0.05, 0.05], dtype=torch.cfloat),  # concentrated on h0
    torch.tensor([0.05, 0.05, 0.8, 0.1], dtype=torch.cfloat),  # concentrated on h2
    torch.tensor([0.0, 0.0, 0.0, 1.0], dtype=torch.cfloat),    # pure ignorance
    torch.tensor([0.33, 0.33, 0.34, 0.0], dtype=torch.cfloat),  # uniform singletons, no θ
    torch.tensor([0.25, 0.25, 0.25, 0.25], dtype=torch.cfloat), # everything equal
]

# Generate a shared evidence stream
_random.seed(42)
evidence_stream = []
for _ in range(100):
    # Evidence concentrated on one hypothesis
    hi = _random.randint(0, H-1)
    ev = torch.zeros(MASS_DIM, dtype=torch.cfloat)
    ev[hi] = 0.6 + 0j
    for k in range(H):
        if k != hi:
            ev[k] = 0.05 + 0j
    ev[H] = 1.0 - ev[:H].sum()
    evidence_stream.append(ev)

# Apply the same evidence stream to all initial masses
print(f"\n  Convergence of 5 initial masses under shared evidence:")
print(f"  {'step':>5s}", end="")
for i in range(len(initial_masses)):
    print(f"  {'mass'+str(i):>10s}", end="")
print(f"  {'max_dist':>10s}")

masses = [m.clone() for m in initial_masses]
for step in range(100):
    ev = evidence_stream[step]

    for i in range(len(masses)):
        # DS combine with discounted evidence
        disc_ev = torch.zeros_like(ev)
        disc_ev[:H] = ev[:H] * 0.3
        disc_ev[H] = 1.0 - disc_ev[:H].sum()

        combined, K = ds_combine(masses[i], disc_ev)
        masses[i] = combined

    if step % 10 == 0 or step < 5:
        # Report Born(h0) for each mass as a proxy
        borns = [born_probabilities(m)[0].item() for m in masses]
        max_dist = max(borns) - min(borns)
        print(f"  {step:>5d}", end="")
        for b in borns:
            print(f"  {b:>10.4f}", end="")
        print(f"  {max_dist:>10.6f}")

# Final state comparison
print(f"\n  Final states (step 99):")
for i, m in enumerate(masses):
    bp = born_probabilities(m)
    print(f"    mass{i}: Born = [{bp[0]:.4f}, {bp[1]:.4f}, {bp[2]:.4f}, {bp[3]:.4f}]")

# Do they converge to the SAME state?
final_borns = [born_probabilities(m) for m in masses]
max_diff = 0
for i in range(len(final_borns)):
    for j in range(i+1, len(final_borns)):
        diff = (final_borns[i] - final_borns[j]).abs().max().item()
        max_diff = max(max_diff, diff)
print(f"\n  Maximum Born difference between any two final states: {max_diff:.6f}")
if max_diff < 0.001:
    print(f"  → CONVERGED to the same state (universal attractor)")
else:
    print(f"  → NOT fully converged (max diff = {max_diff:.4f})")


# ═══════════════════════════════════════════════════════════════
#  Test 4: What IS the fixed point?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("WHAT IS THE FIXED POINT?")
print("The converged state — what are its properties?")
print("=" * 70)

# Use mass0 as representative (they should all be the same)
fixed_mass = masses[0]
bp = born_probabilities(fixed_mass)
print(f"\n  Fixed point Born probabilities:")
print(f"    Born(h0) = {bp[0].item():.6f}")
print(f"    Born(h1) = {bp[1].item():.6f}")
print(f"    Born(h2) = {bp[2].item():.6f}")
print(f"    Born(θ)  = {bp[3].item():.6f}")
print(f"    Born floor = {BORN_FLOOR:.6f}")

# Entropy
ent_val = born_entropy(fixed_mass).item()
max_ent = math.log(MASS_DIM)
print(f"\n  Born entropy: {ent_val:.4f} (max = {max_ent:.4f})")
print(f"  Entropy ratio: {ent_val/max_ent:.4f}")

# Is the fixed point the evidence?
# If we take the "average evidence" = uniform over h0, h1, h2
avg_evidence = torch.tensor([0.6/3, 0.6/3, 0.6/3, 0.2], dtype=torch.cfloat)
bp_ev = born_probabilities(avg_evidence)
print(f"\n  Average evidence Born: [{bp_ev[0]:.4f}, {bp_ev[1]:.4f}, {bp_ev[2]:.4f}, {bp_ev[3]:.4f}]")

# Is it determined by the specific evidence sequence?
# Run with different seeds
print(f"\n  Dependence on evidence sequence (different seeds):")
for seed in [0, 1, 2, 3, 4]:
    _random.seed(seed)
    evidence_2 = []
    for _ in range(100):
        hi = _random.randint(0, H - 1)
        ev = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev[hi] = 0.6; ev[H] = 0.4
        for k in range(H):
            if k != hi: ev[k] = 0.05
        ev[H] = 1.0 - ev[:H].sum()
        evidence_2.append(ev)

    m = torch.tensor([0.25, 0.25, 0.25, 0.25], dtype=torch.cfloat)
    for ev in evidence_2:
        disc_ev = torch.zeros_like(ev); disc_ev[:H] = ev[:H] * 0.3; disc_ev[H] = 1.0 - disc_ev[:H].sum()
        m, _ = ds_combine(m, disc_ev)

    bp_s = born_probabilities(m)
    # Which hypothesis has highest Born?
    dominant = bp_s[:H].argmax().item()
    # Count evidence: how many times was each hypothesis the dominant in evidence?
    ev_counts = [0, 0, 0]
    for ev in evidence_2:
        ev_counts[ev[:H].abs().argmax().item()] += 1

    print(f"    seed {seed}: Born=[{bp_s[0]:.3f},{bp_s[1]:.3f},{bp_s[2]:.3f},{bp_s[3]:.3f}], "
          f"dominant=h{dominant}, evidence counts={ev_counts}")


# ═══════════════════════════════════════════════════════════════
#  Test 5: Equilibrium between building and destroying
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("THE EQUILIBRIUM: building vs destroying")
print("Alternate: one step of entangling, one step of composition.")
print("Does Schmidt reach a steady state?")
print("=" * 70)

id_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

# Build one crystal as the "composition partner"
ent_partner = Entangler(id_corr, seed=99).build()
partner_joint = ent_partner.joint.clone()

# Start from ignorance
joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
joint[-1, -1] = 1.0 + 0j

_random.seed(42)
flat = id_corr.reshape(-1)
flat = flat / flat.sum().clamp(min=EPS_LOG)

print(f"\n  {'cycle':>6s} {'after_build':>12s} {'after_compose':>14s} {'net_change':>11s}")
print("  " + "-" * 50)

for cycle in range(30):
    # BUILD: one entangling step
    idx = _random.choices(range(H * H), weights=flat.tolist())[0]
    hi, hj = idx // H, idx % H
    ev_a = torch.zeros(MASS_DIM, dtype=torch.cfloat)
    ev_b = torch.zeros(MASS_DIM, dtype=torch.cfloat)
    ev_a[hi] = 0.55 + 0.06j; ev_b[hj] = 0.55 + 0.06j
    for k in range(H):
        noise = _random.gauss(0, 0.02)
        if k != hi: ev_a[k] = 0.05 + noise * 1j
        if k != hj: ev_b[k] = 0.05 + noise * 1j
    ev_a[H] = (1.0 - ev_a[:H].real.sum()) - ev_a[:H].imag.sum() * 1j
    ev_b[H] = (1.0 - ev_b[:H].real.sum()) - ev_b[:H].imag.sum() * 1j

    joint_dim = MASS_DIM ** 2; jH = joint_dim - 1
    joint_ev = torch.einsum("m,n->mn", ev_a, ev_b).reshape(joint_dim)
    ev_s = joint_ev[:jH] * 0.3; ev_t = 1.0 + 0j - ev_s.sum()
    disc_ev = torch.cat([ev_s, ev_t.unsqueeze(0)])
    curr = joint.reshape(joint_dim)
    cs1, ct1 = curr[:jH], curr[jH:]
    cs2, ct2 = disc_ev[:jH], disc_ev[jH:]
    K = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
    comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
    comb_t = ct1 * ct2
    norm_val = 1.0 - K
    safe_norm = norm_val if abs(norm_val) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
    updated = torch.cat([comb_s, comb_t]) / safe_norm
    re_sum = updated.real.sum()
    if abs(re_sum) > EPS_DIV: updated = updated / re_sum
    joint = updated.reshape(MASS_DIM, MASS_DIM)

    sn_after_build = schmidt_number(joint)

    # DESTROY: one composition step
    joint = compose(joint, partner_joint)
    sn_after_compose = schmidt_number(joint)

    net = sn_after_compose - sn_after_build

    if cycle % 3 == 0 or cycle < 5:
        print(f"  {cycle:>6d} {sn_after_build:>12.4f} {sn_after_compose:>14.4f} {net:>+11.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 6: The deeper equilibrium — build N, compose 1
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("BUILD N, COMPOSE 1: finding the balance ratio")
print("How many build steps per compose step gives equilibrium?")
print("=" * 70)

for build_ratio in [1, 2, 4, 8, 16, 32]:
    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j

    _random.seed(42)
    schmidts = []

    for cycle in range(20):
        # Build N steps
        for _ in range(build_ratio):
            idx = _random.choices(range(H * H), weights=flat.tolist())[0]
            hi, hj = idx // H, idx % H
            ev_a = torch.zeros(MASS_DIM, dtype=torch.cfloat)
            ev_b = torch.zeros(MASS_DIM, dtype=torch.cfloat)
            ev_a[hi] = 0.55 + 0.06j; ev_b[hj] = 0.55 + 0.06j
            for k in range(H):
                noise = _random.gauss(0, 0.02)
                if k != hi: ev_a[k] = 0.05 + noise * 1j
                if k != hj: ev_b[k] = 0.05 + noise * 1j
            ev_a[H] = (1.0 - ev_a[:H].real.sum()) - ev_a[:H].imag.sum() * 1j
            ev_b[H] = (1.0 - ev_b[:H].real.sum()) - ev_b[:H].imag.sum() * 1j
            joint_dim = MASS_DIM ** 2; jH = joint_dim - 1
            joint_ev = torch.einsum("m,n->mn", ev_a, ev_b).reshape(joint_dim)
            ev_s = joint_ev[:jH] * 0.3; ev_t = 1.0 + 0j - ev_s.sum()
            disc_ev = torch.cat([ev_s, ev_t.unsqueeze(0)])
            curr = joint.reshape(joint_dim)
            cs1, ct1 = curr[:jH], curr[jH:]
            cs2, ct2 = disc_ev[:jH], disc_ev[jH:]
            K_val = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
            comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
            comb_t = ct1 * ct2
            nv = 1.0 - K_val
            sn_v = nv if abs(nv) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
            upd = torch.cat([comb_s, comb_t]) / sn_v
            rs = upd.real.sum()
            if abs(rs) > EPS_DIV: upd = upd / rs
            joint = upd.reshape(MASS_DIM, MASS_DIM)

        # Compose 1 step
        joint = compose(joint, partner_joint)
        schmidts.append(schmidt_number(joint))

    # Check if Schmidt has stabilized in the last 5 cycles
    late = schmidts[-5:]
    avg_late = sum(late) / len(late)
    std_late = (sum((s - avg_late)**2 for s in late) / len(late)) ** 0.5

    print(f"  build:compose = {build_ratio:2d}:1 → avg Schmidt (last 5) = {avg_late:.4f} ± {std_late:.4f}")

    # The growth/decay ratio is (Δ/4 * N) vs (2Δ * 1)
    # Equilibrium when Δ/4 * N = 2Δ → N = 8
    # So build:compose = 8:1 should be the equilibrium


print("\n  Predicted equilibrium: build:compose = 8:1")
print(f"  (from growth rate Δ/4 per step, decay rate 2Δ per composition)")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
