"""
S₃ crystal geometry: fermion masses from phases, CKM from distances.

Thread 1 (handover): "The masses might be encoded in the PHASES, not just
the magnitudes." Eigenvalue phases: transpositions vs 3-cycles differ.

Thread 7 (handover): "The Cabibbo angle might be a GEOMETRIC angle in
crystal space (FS distance between quark crystals?)."
V_us ≈ 1/(S/26) to 0.14%, V_cb ≈ 1/S^(2/3) to 0.21%.

Hypothesis: the S₃ conjugacy classes {e}, {transpositions}, {3-cycles}
define three "flavor" states. The FS distances between them give mixing
angles. The eigenvalue phases give mass ratios.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import numpy as np
from solver.algebra import (
    H, MASS_DIM, K_STAR, BORN_FLOOR, DELTA,
    born_probabilities, born_fidelity, schmidt_number,
    sym2_fingerprint, sym2_distance,
)
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)

S = H**3 / K_STAR  # 810/7
S_gauge = S / H     # 270/7


# ═══════════════════════════════════════════════════════════════════
#  S₃ elements as correlation matrices
# ═══════════════════════════════════════════════════════════════════

S3_ELEMENTS = {
    "e":       torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":    torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":    torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":    torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)":   torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)":   torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

CONJUGACY_CLASSES = {
    "identity":      ["e"],
    "transpositions": ["(01)", "(02)", "(12)"],
    "3-cycles":      ["(012)", "(021)"],
}


def get_eigenvalues(joint):
    eigenvalues, eigenvectors = torch.linalg.eig(joint)
    mags = eigenvalues.abs()
    idx = mags.argsort(descending=True)
    return eigenvalues[idx], eigenvectors[:, idx]


def fs_distance(j1, j2):
    """Fubini-Study distance between two crystals (via Born fidelity)."""
    fid = born_fidelity(j1, j2)
    fid = min(fid, 1.0)  # numerical safety
    return math.acos(max(fid, -1.0))


# ═══════════════════════════════════════════════════════════════════
print("=" * 80)
print("PART 1: SEED-AVERAGED EIGENVALUE SPECTRA BY CONJUGACY CLASS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════════

N_SEEDS = 200

class_spectra = {}  # class_name -> {ratios: [...], phases: [...]}

for class_name, elements in CONJUGACY_CLASSES.items():
    all_ratios = [[] for _ in range(3)]
    all_phases = [[] for _ in range(3)]
    all_det_phases = []

    for elem_name in elements:
        corr = S3_ELEMENTS[elem_name]
        for seed in range(N_SEEDS):
            ent = Entangler(corr, seed=seed).build()
            ev, _ = get_eigenvalues(ent.joint)

            phi0 = math.atan2(ev[0].imag.item(), ev[0].real.item())
            for i in range(3):
                r = ev[i+1].abs().item() / max(ev[0].abs().item(), 1e-15)
                all_ratios[i].append(r)
                p = (math.atan2(ev[i+1].imag.item(), ev[i+1].real.item()) - phi0) / math.pi
                p = p - 2*round(p/2)
                all_phases[i].append(p)

            det_val = torch.linalg.det(ent.joint)
            dp = math.atan2(det_val.imag.item(), det_val.real.item())
            all_det_phases.append(dp)

    # Statistics
    avg_r = [np.mean(r) for r in all_ratios]
    std_r = [np.std(r) for r in all_ratios]
    avg_p = [np.mean(p) for p in all_phases]
    std_p = [np.std(p) for p in all_phases]
    avg_det = np.mean(all_det_phases)
    std_det = np.std(all_det_phases)

    class_spectra[class_name] = {
        "avg_ratios": avg_r,
        "std_ratios": std_r,
        "avg_phases": avg_p,
        "std_phases": std_p,
        "avg_det_phase": avg_det,
        "std_det_phase": std_det,
    }

    n_samples = len(all_ratios[0])
    print(f"\n  {class_name} ({n_samples} samples):")
    for i in range(3):
        print(f"    λ_{i+1}/λ₀: {avg_r[i]:.4f} ± {std_r[i]:.4f}  "
              f"  φ_{i+1}/π: {avg_p[i]:+.4f} ± {std_p[i]:.4f}")
    print(f"    det phase/π: {avg_det/math.pi:+.4f} ± {std_det/math.pi:.4f}")


# ═══════════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 2: CONJUGACY CLASS AVERAGES — THE THREE FLAVOR CRYSTALS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════════

# Build seed-averaged crystals for each conjugacy class
class_crystals = {}

for class_name, elements in CONJUGACY_CLASSES.items():
    avg_joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    count = 0
    for elem_name in elements:
        corr = S3_ELEMENTS[elem_name]
        for seed in range(N_SEEDS):
            ent = Entangler(corr, seed=seed).build()
            avg_joint += ent.joint
            count += 1
    avg_joint /= count
    # Re-normalize
    re_sum = avg_joint.reshape(-1).real.sum()
    if abs(re_sum) > 1e-10:
        avg_joint = avg_joint / re_sum
    class_crystals[class_name] = avg_joint

    sn = schmidt_number(avg_joint)
    bp = born_probabilities(avg_joint.reshape(-1))

    print(f"\n  {class_name}:")
    print(f"    Schmidt = {sn:.4f}")
    ev, _ = get_eigenvalues(avg_joint)
    print(f"    |λ₀| = {ev[0].abs().item():.6f}")
    for i in range(1, MASS_DIM):
        r = ev[i].abs().item() / max(ev[0].abs().item(), 1e-15)
        phi0 = math.atan2(ev[0].imag.item(), ev[0].real.item())
        p = (math.atan2(ev[i].imag.item(), ev[i].real.item()) - phi0) / math.pi
        p = p - 2*round(p/2)
        print(f"    λ_{i}/λ₀: |{r:.6f}| ∠ {p:+.6f}π")


# ═══════════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 3: FUBINI-STUDY DISTANCE MATRIX")
print("Between conjugacy class crystals (the geometry of flavor space)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════════

names = list(class_crystals.keys())
print(f"\n  {'':>20s}", end="")
for n in names:
    print(f"  {n:>15s}", end="")
print()

dist_matrix = {}
for i, n1 in enumerate(names):
    print(f"  {n1:>20s}", end="")
    for j, n2 in enumerate(names):
        d = fs_distance(class_crystals[n1], class_crystals[n2])
        dist_matrix[(n1, n2)] = d
        print(f"  {d:>15.6f}", end="")
    print()

# Also compute element-level distances
print(f"\n  Full 6×6 element distance matrix:")
elem_names = list(S3_ELEMENTS.keys())
elem_crystals = {}
for name, corr in S3_ELEMENTS.items():
    # Seed-average over 200 seeds
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(N_SEEDS):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= N_SEEDS
    avg /= avg.reshape(-1).real.sum()
    elem_crystals[name] = avg

print(f"\n  {'':>10s}", end="")
for n in elem_names:
    print(f"  {n:>10s}", end="")
print()
for n1 in elem_names:
    print(f"  {n1:>10s}", end="")
    for n2 in elem_names:
        d = fs_distance(elem_crystals[n1], elem_crystals[n2])
        print(f"  {d:>10.6f}", end="")
    print()


# ═══════════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 4: FS DISTANCES vs CKM MATRIX ELEMENTS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════════

# CKM matrix elements (PDG 2024):
CKM = {
    "V_ud": 0.97435,  "V_us": 0.22500, "V_ub": 0.00369,
    "V_cd": 0.22486,  "V_cs": 0.97349, "V_cb": 0.04182,
    "V_td": 0.00857,  "V_ts": 0.04110, "V_tb": 0.999118,
}

# Wolfenstein parameters:
lam_W = 0.22500   # = sin(θ_C) ≈ V_us
A_W   = 0.826
rho_W = 0.159
eta_W = 0.348

# CKM angles (radians):
theta_12 = math.asin(lam_W)          # Cabibbo angle = 0.2271 rad = 13.01°
theta_23 = math.asin(A_W * lam_W**2) # = 0.04182 rad = 2.40°
theta_13 = math.asin(lam_W**3 * A_W * math.sqrt(rho_W**2 + eta_W**2))

print(f"\n  CKM angles:")
print(f"    θ₁₂ (Cabibbo): {theta_12:.6f} rad = {math.degrees(theta_12):.4f}°")
print(f"    θ₂₃:           {theta_23:.6f} rad = {math.degrees(theta_23):.4f}°")
print(f"    θ₁₃:           {theta_13:.6f} rad = {math.degrees(theta_13):.4f}°")

print(f"\n  FS distances between conjugacy classes:")
d_et = dist_matrix[("identity", "transpositions")]
d_e3 = dist_matrix[("identity", "3-cycles")]
d_t3 = dist_matrix[("transpositions", "3-cycles")]
print(f"    d(identity, transpositions) = {d_et:.6f} rad = {math.degrees(d_et):.4f}°")
print(f"    d(identity, 3-cycles)       = {d_e3:.6f} rad = {math.degrees(d_e3):.4f}°")
print(f"    d(transpositions, 3-cycles) = {d_t3:.6f} rad = {math.degrees(d_t3):.4f}°")

# Compare various ratios/combinations to CKM
print(f"\n  Testing ratios against CKM:")
for name, d in [("d(e,t)", d_et), ("d(e,3)", d_e3), ("d(t,3)", d_t3)]:
    for ckm_name, angle in [("θ_12", theta_12), ("θ_23", theta_23)]:
        ratio = d / angle
        print(f"    {name} / {ckm_name} = {ratio:.4f}")
    # Also sin of the distance
    print(f"    sin({name}) = {math.sin(d):.6f}  (V_us = {CKM['V_us']:.6f})")

# Now test: are the FS distances PROPORTIONAL to CKM angles
# through some crystal constant?
print(f"\n  Direct value comparisons:")
for name, d in [("d(e,t)", d_et), ("d(e,3)", d_e3), ("d(t,3)", d_t3)]:
    print(f"    {name} = {d:.6f}")
    print(f"      /{S:.2f} = {d/S:.6f}  (V_ub = {CKM['V_ub']:.6f})")
    print(f"      ×H = {d*H:.6f}")
    print(f"      /K* = {d/K_STAR:.6f}")
    print(f"      /π = {d/math.pi:.6f}")


# ═══════════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 5: EIGENVALUE PHASE DIFFERENCES BETWEEN CONJUGACY CLASSES")
print("(the mass ratios?)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════════

# The eigenvalue phases differ between conjugacy classes.
# Phase differences might encode mass RATIOS.

for i in range(3):
    p_e = class_spectra["identity"]["avg_phases"][i]
    p_t = class_spectra["transpositions"]["avg_phases"][i]
    p_3 = class_spectra["3-cycles"]["avg_phases"][i]

    dp_et = abs(p_e - p_t)
    dp_e3 = abs(p_e - p_3)
    dp_t3 = abs(p_t - p_3)

    print(f"\n  φ_{i+1} differences:")
    print(f"    identity:      {p_e:+.4f}π")
    print(f"    transpositions:{p_t:+.4f}π")
    print(f"    3-cycles:      {p_3:+.4f}π")
    print(f"    Δφ(e,t) = {dp_et:.4f}π = {dp_et*180:.2f}°")
    print(f"    Δφ(e,3) = {dp_e3:.4f}π = {dp_e3*180:.2f}°")
    print(f"    Δφ(t,3) = {dp_t3:.4f}π = {dp_t3*180:.2f}°")

# Ratio magnitudes between classes
print(f"\n  Eigenvalue ratio differences (|λ/λ₀| averaged):")
for i in range(3):
    r_e = class_spectra["identity"]["avg_ratios"][i]
    r_t = class_spectra["transpositions"]["avg_ratios"][i]
    r_3 = class_spectra["3-cycles"]["avg_ratios"][i]
    print(f"    λ_{i+1}/λ₀: identity={r_e:.4f}, trans={r_t:.4f}, 3cyc={r_3:.4f}")
    if r_e > 0.01 and r_t > 0.01 and r_3 > 0.01:
        print(f"      ratios: trans/identity = {r_t/r_e:.4f}, "
              f"3cyc/identity = {r_3/r_e:.4f}, "
              f"3cyc/trans = {r_3/r_t:.4f}")


# ═══════════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 6: FERMION MASS RATIOS FROM CRYSTAL INVARIANTS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════════

# Known fermion mass ratios (at MZ scale, PDG):
# Down-type quarks:
m_d = 4.67e-3    # GeV
m_s = 93.4e-3    # GeV
m_b = 4.18       # GeV

# Up-type quarks:
m_u = 2.16e-3    # GeV
m_c = 1.27       # GeV
m_t = 172.69     # GeV

# Charged leptons:
m_e = 0.511e-3   # GeV
m_mu = 0.1057    # GeV
m_tau = 1.777    # GeV

print(f"\n  Fermion mass ratios within generations:")
print(f"    Down-type: m_d/m_b = {m_d/m_b:.6f}, m_s/m_b = {m_s/m_b:.6f}")
print(f"    Up-type:   m_u/m_t = {m_u/m_t:.2e}, m_c/m_t = {m_c/m_t:.6f}")
print(f"    Leptons:   m_e/m_τ = {m_e/m_tau:.6f}, m_μ/m_τ = {m_mu/m_tau:.6f}")

print(f"\n  Cross-generation ratios:")
print(f"    m_e/m_μ = {m_e/m_mu:.6f}")
print(f"    m_μ/m_τ = {m_mu/m_tau:.6f}")
print(f"    m_d/m_s = {m_d/m_s:.6f}")
print(f"    m_s/m_b = {m_s/m_b:.6f}")
print(f"    m_u/m_c = {m_u/m_c:.6f}")
print(f"    m_c/m_t = {m_c/m_t:.6f}")

# The eigenvalue ratio |λ₁/λ₀| is the dominant crystal parameter.
# Under n compositions, the ratio becomes (|λ₁/λ₀|)^n.
# Could fermion masses be exp(-n × crystal_gap) for different n?

print(f"\n  If mass = exp(-n × S_gauge/H^k):")
print(f"    S_gauge = {S_gauge:.4f}")
print(f"    S_gauge/H = {S_gauge/H:.4f}")
print(f"    S_gauge/H² = {S_gauge/H**2:.4f}")

# The Koide formula: (m_e + m_μ + m_τ) / (√m_e + √m_μ + √m_τ)² = 2/3
koide = (m_e + m_mu + m_tau) / (math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(m_tau))**2
print(f"\n  Koide formula: (Σm)/(Σ√m)² = {koide:.6f}  (theory: 2/3 = {2/3:.6f})")

# Crystal version: 2/3 = (H-1)/H = 1 - 1/H
print(f"    2/3 = (H-1)/H at H=3: {(H-1)/H:.6f} ✓")


# ═══════════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 7: THE COMPOSITION OPERATOR SPECTRUM BY CONJUGACY CLASS")
print("How does each class evolve under self-composition?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════════

for class_name in CONJUGACY_CLASSES:
    crystal = class_crystals[class_name]
    current = crystal.clone()

    print(f"\n  {class_name}:")
    print(f"  {'step':>5s} {'Schmidt':>8s} {'|λ₁/λ₀|':>10s} {'|λ₂/λ₀|':>10s} "
          f"{'|λ₃/λ₀|':>10s} {'φ₁/π':>10s}")

    for step in range(12):
        sn = schmidt_number(current)
        ev, _ = get_eigenvalues(current)
        ratios = [ev[i].abs().item() / max(ev[0].abs().item(), 1e-15)
                  for i in range(1, MASS_DIM)]
        phi0 = math.atan2(ev[0].imag.item(), ev[0].real.item())
        p1 = (math.atan2(ev[1].imag.item(), ev[1].real.item()) - phi0) / math.pi
        p1 = p1 - 2*round(p1/2)

        print(f"  {step:>5d} {sn:>8.4f} {ratios[0]:>10.6f} {ratios[1]:>10.6f} "
              f"{ratios[2]:>10.6f} {p1:>+10.4f}")

        current = compose(current, crystal)


# ═══════════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 8: DS CONFLICT BETWEEN CONJUGACY CLASSES")
print("K(class_i, class_j) — the conflict matrix of flavor space")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════════

from solver.algebra import ds_combine

print(f"\n  Conflict matrix K(i,j):")
print(f"  {'':>20s}", end="")
for n in names:
    print(f"  {n:>15s}", end="")
print()

for n1 in names:
    print(f"  {n1:>20s}", end="")
    for n2 in names:
        # Flatten to 1D mass vectors and combine
        m1 = class_crystals[n1].reshape(-1)
        m2 = class_crystals[n2].reshape(-1)
        # Re-normalize for DS combine
        re1 = m1.real.sum()
        re2 = m2.real.sum()
        if abs(re1) > 1e-10: m1 = m1 / re1
        if abs(re2) > 1e-10: m2 = m2 / re2
        _, K = ds_combine(m1.unsqueeze(0), m2.unsqueeze(0))
        print(f"  {K.item():>15.6f}", end="")
    print()


# ═══════════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 9: CROSS-COMPOSITION — THE MIXING MATRIX")
print("compose(class_i crystal, class_j crystal) → what relationship?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════════

# If S₃ classes are flavors, cross-composition IS the CKM matrix
# (rotation from one flavor basis to another)

print(f"\n  Cross-composition eigenvalue spectra:")
for n1 in names:
    for n2 in names:
        result = compose(class_crystals[n1], class_crystals[n2])
        sn = schmidt_number(result)
        ev, _ = get_eigenvalues(result)
        ratios = [ev[i].abs().item() / max(ev[0].abs().item(), 1e-15)
                  for i in range(1, MASS_DIM)]
        # Fidelity with each class crystal
        fids = {n: born_fidelity(result, class_crystals[n]) for n in names}
        best = max(fids, key=fids.get)
        print(f"  {n1:>15s} ○ {n2:<15s}: Schmidt={sn:.4f}  "
              f"|λ₁/λ₀|={ratios[0]:.4f}  closest={best} (fid={fids[best]:.4f})")


# ═══════════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 10: THE NUMBER (H-1)/H = 2/3 AND THE KOIDE FORMULA")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════════

# Koide: (m_e + m_μ + m_τ) / (√m_e + √m_μ + √m_τ)² = 2/3
# This is (H-1)/H. What if the three masses come from the crystal's
# three sub-dominant eigenvalue magnitudes?

# For each conjugacy class, the three sub-dominant eigenvalue ratios
# (relative to dominant) define three "mass coordinates"

for class_name in CONJUGACY_CLASSES:
    crystal = class_crystals[class_name]
    ev, _ = get_eigenvalues(crystal)
    r = [ev[i].abs().item() / max(ev[0].abs().item(), 1e-15) for i in range(1, MASS_DIM)]

    # Koide-like quantity
    sum_r = sum(r)
    sum_sqrt_r = sum(math.sqrt(ri) for ri in r)
    koide_r = sum_r / sum_sqrt_r**2 if sum_sqrt_r > 0 else 0

    # Using r² (Born probabilities of eigenvalues)
    r2 = [ri**2 for ri in r]
    sum_r2 = sum(r2)
    sum_sqrt_r2 = sum(math.sqrt(ri) for ri in r2)
    koide_r2 = sum_r2 / sum_sqrt_r2**2 if sum_sqrt_r2 > 0 else 0

    print(f"\n  {class_name}:")
    print(f"    Eigenvalue ratios: {r[0]:.6f}, {r[1]:.6f}, {r[2]:.6f}")
    print(f"    Koide(ratios):     {koide_r:.6f}  (target: {2/3:.6f})")
    print(f"    Koide(ratios²):    {koide_r2:.6f}  (target: {2/3:.6f})")

    # Physical mass comparison:
    # If r[i] ∝ √(m_i/m_3), then Koide holds iff it holds for the masses
    # The actual lepton masses give Koide = 0.6667
    # The crystal eigenvalues give...?


# ═══════════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("FINDINGS")
print("=" * 80)
