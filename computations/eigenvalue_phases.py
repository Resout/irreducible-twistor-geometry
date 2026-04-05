"""
Eigenvalue phases: conservation law and representation theory.

The sub-dominant eigenvalue phases sum to ≈ π per composition step.
Is this exact? What determines how the π distributes among the three
sub-dominant eigenvalues?

And: the two phase families (0.12π vs 0.85π) — are these connected
to S₃ representations? S₃ has three irreps: trivial (dim 1),
sign (dim 1), standard (dim 2). The crystals might decompose
into these irreps, with eigenvalue phases as characters.

Also: the eigenvalue as natural coordinate. If eigenvalues evolve
independently under composition, they're the FUNDAMENTAL description
of the crystal. Everything else (Schmidt, θ-fingerprint, slot)
should be derivable from the eigenvalue spectrum.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, schmidt_number
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)


def get_sorted_eigenvalues(joint):
    eigenvalues, eigenvectors = torch.linalg.eig(joint)
    mags = eigenvalues.abs()
    idx = mags.argsort(descending=True)
    return eigenvalues[idx], eigenvectors[:, idx]


# ═══════════════════════════════════════════════════════════════
#  Test 1: Phase sum conservation
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 1: DO SUB-DOMINANT PHASES SUM TO π?")
print("=" * 70)

for name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Entangler(corr).build()
    ev, _ = get_sorted_eigenvalues(ent.joint)

    # Phases
    phases = [math.atan2(e.imag.item(), e.real.item()) for e in ev]

    # Phase sum (all 4)
    total_phase = sum(phases)
    # Sub-dominant sum (phases 1-3)
    sub_sum = sum(phases[1:])

    # Phase of determinant
    det_val = torch.linalg.det(ent.joint)
    det_phase = math.atan2(det_val.imag.item(), det_val.real.item())

    print(f"\n  {name}:")
    print(f"    Phases/π: [{', '.join(f'{p/math.pi:+.4f}' for p in phases)}]")
    print(f"    Sum(all)/π = {total_phase/math.pi:+.6f}")
    print(f"    Sum(sub)/π = {sub_sum/math.pi:+.6f}")
    print(f"    det phase/π = {det_phase/math.pi:+.6f}")
    print(f"    Sum(all) = det phase? diff = {abs(total_phase - det_phase):.6f}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: Phase steps per composition
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 2: PHASE STEP SUMS UNDER COMPOSITION")
print("Does the sum of phase steps per composition = π?")
print("=" * 70)

for name in ["proportional", "inverse", "modular"]:
    corr = RELATIONSHIP_SIGNATURES[name]
    ent = Entangler(corr).build()
    joint = ent.joint.clone()

    current = joint.clone()
    prev_phases = None

    print(f"\n  {name}:")
    print(f"  {'step':>5s} {'Δφ₀/π':>8s} {'Δφ₁/π':>8s} {'Δφ₂/π':>8s} {'Δφ₃/π':>8s} {'Σ(Δφ)/π':>8s} {'Σ(sub)/π':>8s}")

    for step in range(10):
        ev, _ = get_sorted_eigenvalues(current)
        phases = [math.atan2(e.imag.item(), e.real.item()) for e in ev]

        if prev_phases is not None:
            deltas = []
            for i in range(MASS_DIM):
                dp = phases[i] - prev_phases[i]
                # Unwrap to [-π, π]
                while dp > math.pi: dp -= 2 * math.pi
                while dp < -math.pi: dp += 2 * math.pi
                deltas.append(dp)

            total_step = sum(deltas)
            sub_step = sum(deltas[1:])
            print(f"  {step:>5d} {deltas[0]/math.pi:>+8.4f} {deltas[1]/math.pi:>+8.4f} "
                  f"{deltas[2]/math.pi:>+8.4f} {deltas[3]/math.pi:>+8.4f} "
                  f"{total_step/math.pi:>+8.4f} {sub_step/math.pi:>+8.4f}")

        prev_phases = phases
        current = compose(current, joint)


# ═══════════════════════════════════════════════════════════════
#  Test 3: Eigenvalue basis as natural coordinates
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 3: EIGENVALUE COORDINATES OF THE CRYSTAL SPACE")
print("Map each crystal to (|λ₁/λ₀|, |λ₂/λ₀|, |λ₃/λ₀|, φ₁, φ₂, φ₃)")
print("=" * 70)

print(f"\n  {'Crystal':>15s} {'|λ₁/λ₀|':>8s} {'|λ₂/λ₀|':>8s} {'|λ₃/λ₀|':>8s} "
      f"{'φ₁/π':>8s} {'φ₂/π':>8s} {'φ₃/π':>8s} {'Schmidt':>8s}")
print("  " + "-" * 75)

ev_coords = {}
for name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Entangler(corr).build()
    ev, _ = get_sorted_eigenvalues(ent.joint)

    ratios = [ev[i].abs().item() / ev[0].abs().item() for i in range(1, MASS_DIM)]
    phases = [math.atan2(ev[i].imag.item(), ev[i].real.item()) / math.pi for i in range(1, MASS_DIM)]
    # Phase relative to λ₀
    phi0 = math.atan2(ev[0].imag.item(), ev[0].real.item())
    rel_phases = [(math.atan2(ev[i].imag.item(), ev[i].real.item()) - phi0) / math.pi
                  for i in range(1, MASS_DIM)]
    # Wrap to [-1, 1]
    rel_phases = [p - 2*round(p/2) for p in rel_phases]

    sn = schmidt_number(ent.joint)

    ev_coords[name] = {
        "ratios": ratios,
        "phases": rel_phases,
        "schmidt": sn,
    }

    print(f"  {name:>15s} {ratios[0]:>8.4f} {ratios[1]:>8.4f} {ratios[2]:>8.4f} "
          f"{rel_phases[0]:>+8.4f} {rel_phases[1]:>+8.4f} {rel_phases[2]:>+8.4f} {sn:>8.3f}")

# Can Schmidt be predicted from eigenvalue ratios alone?
# Schmidt = 1 / Σ(σᵢ/Σσⱼ)² = 1/IPR of singular values
# For normal matrices, |λᵢ| = σᵢ. For non-normal, they differ.
# But the eigenvalue ratios should still correlate with Schmidt.

all_ratios_sum = []
all_schmidts = []
for name, ec in ev_coords.items():
    # Sum of squared ratios ≈ 1/Schmidt (if normal)
    ipr = 1.0 + sum(r**2 for r in ec["ratios"])
    total = 1.0 + sum(ec["ratios"])
    normalized_ipr = ipr / (total**2)
    predicted_schmidt = 1.0 / normalized_ipr if normalized_ipr > 0 else 1

    all_ratios_sum.append(sum(ec["ratios"]))
    all_schmidts.append(ec["schmidt"])

    print(f"    {name}: predicted Schmidt = {predicted_schmidt:.3f}, actual = {ec['schmidt']:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 4: Seed variation in eigenvalue coordinates
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 4: ARE EIGENVALUE COORDINATES SEED-INVARIANT?")
print("If eigenvalues are fundamental, they should be stable across seeds")
print("=" * 70)

id_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

all_r1 = []; all_r2 = []; all_r3 = []
all_p1 = []; all_p2 = []; all_p3 = []

for seed in range(50):
    ent = Entangler(id_corr, seed=seed).build()
    ev, _ = get_sorted_eigenvalues(ent.joint)
    phi0 = math.atan2(ev[0].imag.item(), ev[0].real.item())

    for i, (rl, pl) in enumerate([(all_r1, all_p1), (all_r2, all_p2), (all_r3, all_p3)]):
        rl.append(ev[i+1].abs().item() / ev[0].abs().item())
        p = (math.atan2(ev[i+1].imag.item(), ev[i+1].real.item()) - phi0) / math.pi
        p = p - 2*round(p/2)
        pl.append(p)

def stats(vals):
    avg = sum(vals) / len(vals)
    std = (sum((v-avg)**2 for v in vals) / len(vals))**0.5
    return avg, std

for label, rl, pl in [("λ₁", all_r1, all_p1), ("λ₂", all_r2, all_p2), ("λ₃", all_r3, all_p3)]:
    ra, rs = stats(rl)
    pa, ps = stats(pl)
    print(f"  {label}: |λ/λ₀| = {ra:.4f} ± {rs:.4f} (CV={rs/ra*100:.1f}%), "
          f"φ/π = {pa:+.4f} ± {ps:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 5: The S₃ decomposition
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 5: S₃ PERMUTATION EIGENVALUES")
print("Do different S₃ elements produce different eigenvalue spectra?")
print("=" * 70)

s3_elements = {
    "identity": torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "swap(01)": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "swap(02)": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "swap(12)": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "cyc(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "cyc(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

print(f"\n  {'Element':>12s} {'order':>5s} {'|λ₁/λ₀|':>8s} {'|λ₂/λ₀|':>8s} {'|λ₃/λ₀|':>8s} "
      f"{'φ₁/π':>8s} {'φ₂/π':>8s} {'φ₃/π':>8s} {'Schmidt':>8s}")
print("  " + "-" * 75)

for name, corr in s3_elements.items():
    ent = Entangler(corr).build()
    ev, _ = get_sorted_eigenvalues(ent.joint)
    sn = schmidt_number(ent.joint)

    ratios = [ev[i].abs().item() / ev[0].abs().item() for i in range(1, MASS_DIM)]
    phi0 = math.atan2(ev[0].imag.item(), ev[0].real.item())
    rel_phases = [(math.atan2(ev[i].imag.item(), ev[i].real.item()) - phi0) / math.pi
                  for i in range(1, MASS_DIM)]
    rel_phases = [p - 2*round(p/2) for p in rel_phases]

    # Order of the permutation
    order = 1
    perm = list(range(3))
    current_perm = [corr[i].argmax().item() for i in range(3)]
    temp = list(current_perm)
    for _ in range(6):
        order += 1
        temp = [current_perm[int(t)] for t in temp]
        if temp == list(range(3)):
            break

    print(f"  {name:>12s} {order:>5d} {ratios[0]:>8.4f} {ratios[1]:>8.4f} {ratios[2]:>8.4f} "
          f"{rel_phases[0]:>+8.4f} {rel_phases[1]:>+8.4f} {rel_phases[2]:>+8.4f} {sn:>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 6: Eigenvalue evolution visualized
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 6: EIGENVALUE TRAJECTORIES IN THE COMPLEX PLANE")
print("Under composition, eigenvalues spiral inward. Show the spiral.")
print("=" * 70)

for name in ["proportional", "inverse"]:
    corr = RELATIONSHIP_SIGNATURES[name]
    ent = Entangler(corr).build()
    joint = ent.joint.clone()
    current = joint.clone()

    print(f"\n  {name} — eigenvalue trajectories (re, im):")
    print(f"  {'step':>5s}  {'λ₀':>20s}  {'λ₁':>20s}  {'λ₂':>20s}  {'λ₃':>20s}")

    for step in range(8):
        ev, _ = get_sorted_eigenvalues(current)
        parts = []
        for e in ev:
            parts.append(f"({e.real.item():+.4f},{e.imag.item():+.4f})")
        print(f"  {step:>5d}  {'  '.join(parts)}")
        current = compose(current, joint)


# ═══════════════════════════════════════════════════════════════
#  Test 7: The det phase is exactly the initial det phase per step
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 7: DET PHASE ACCUMULATION")
print("det(A^n) = det(A)^n (up to normalization)")
print("=" * 70)

for name in ["proportional", "inverse", "modular"]:
    corr = RELATIONSHIP_SIGNATURES[name]
    ent = Entangler(corr).build()
    joint = ent.joint.clone()

    det_0 = torch.linalg.det(joint)
    det_phase_0 = math.atan2(det_0.imag.item(), det_0.real.item())

    current = joint.clone()
    print(f"\n  {name} (det phase₀ = {det_phase_0/math.pi:.4f}π):")

    for step in range(8):
        current = compose(current, joint)
        det_n = torch.linalg.det(current)
        det_phase_n = math.atan2(det_n.imag.item(), det_n.real.item())
        # Predicted: (step+1) * det_phase_0, mod 2π
        predicted = ((step + 2) * det_phase_0) % (2 * math.pi)
        if predicted > math.pi: predicted -= 2 * math.pi
        diff = abs(det_phase_n - predicted)
        if diff > math.pi: diff = 2 * math.pi - diff

        print(f"    comp {step+1}: det phase = {det_phase_n/math.pi:+.4f}π, "
              f"predicted = {predicted/math.pi:+.4f}π, diff = {diff/math.pi:.4f}π")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
