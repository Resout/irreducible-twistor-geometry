"""
Thread 97b: CKM mixing angles, quark masses, and the √2 mechanism.

Four sub-threads:
A. Verify the √2 phase mechanism with per-seed averaging (fix transposition)
B. 3-cycle exponent a₁ with proper per-seed measurement
C. Standard sector doublets → quark mass connection
D. CKM mixing angles from individual element FS distances
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import numpy as np
from solver.algebra import (H, MASS_DIM, K_STAR, BORN_FLOOR, DELTA,
                            born_fidelity, sym2_fingerprint, sym2_distance)
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

dim_sq = MASS_DIM ** 2

# S₃ structure
S3_PERMS = {
    "e":     [0, 1, 2],
    "(01)":  [1, 0, 2],
    "(02)":  [2, 1, 0],
    "(12)":  [0, 2, 1],
    "(012)": [1, 2, 0],
    "(021)": [2, 0, 1],
}

S3_CORR = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

CONJUGACY_CLASSES = {
    "identity": ["e"],
    "3-cycles": ["(012)", "(021)"],
    "transpositions": ["(01)", "(02)", "(12)"],
}

def perm_mat_4x4(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return P

reps = {n: torch.kron(perm_mat_4x4(p), perm_mat_4x4(p)) for n, p in S3_PERMS.items()}

chi_s = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}
chi_t = {"e": 1, "(01)": 1, "(02)": 1, "(12)": 1, "(012)": 1, "(021)": 1}
chi_d = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}

def proj(chi, dim_chi):
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for n in S3_PERMS:
        P += chi[n] * reps[n]
    return P * dim_chi / 6.0

P_sign = proj(chi_s, 1)
P_triv = proj(chi_t, 1)
P_std = proj(chi_d, 2)


def build_T_B(B_joint):
    T = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for col in range(dim_sq):
        A = torch.zeros(dim_sq, dtype=torch.cfloat)
        A[col] = 1.0
        result = compose(A.reshape(MASS_DIM, MASS_DIM), B_joint)
        T[:, col] = result.reshape(dim_sq)
    return T


def sign_coupling(T):
    restricted = P_sign @ T @ P_sign
    evals = torch.linalg.eigvals(restricted)
    idx = evals.abs().argmax()
    return evals[idx]


def std_spectrum(T):
    """Leading standard-sector eigenvalues."""
    restricted = P_std @ T @ P_std
    evals = torch.linalg.eigvals(restricted)
    mags = evals.abs().sort(descending=True).values
    return mags[mags > 1e-8]


def triv_coupling(T):
    restricted = P_triv @ T @ P_triv
    evals = torch.linalg.eigvals(restricted)
    idx = evals.abs().argmax()
    return evals[idx]


N_SEEDS = 200

# Physical constants
m_e = 0.51100; m_mu = 105.66; m_tau = 1776.86
m_u = 2.16; m_d = 4.67; m_s = 93.4
m_c = 1270; m_b = 4180; m_t = 172500

# CKM matrix elements (magnitudes)
V_CKM = {
    "ud": 0.97370, "us": 0.2245, "ub": 0.00382,
    "cd": 0.221, "cs": 0.987, "cb": 0.0410,
    "td": 0.0080, "ts": 0.0388, "tb": 0.999,
}

# Cabibbo angle
theta_C = math.asin(V_CKM["us"])


# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("PART A: THE √2 PHASE MECHANISM (per-seed measurement)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print("\n  Per-seed sign eigenvalue measurement:")

all_data = {}
for class_name, members in CONJUGACY_CLASSES.items():
    mags = []
    phases = []
    for name in members:
        for seed in range(N_SEEDS):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T = build_T_B(ent.joint)
            lam = sign_coupling(T)
            mags.append(lam.abs().item())
            phases.append(math.atan2(lam.imag.item(), lam.real.item()))

    avg_mag = np.mean(mags)
    std_mag = np.std(mags) / np.sqrt(len(mags))
    avg_phase = np.mean(phases)
    std_phase = np.std(phases) / np.sqrt(len(phases))

    all_data[class_name] = {
        "mag": avg_mag, "mag_std": std_mag,
        "phase": avg_phase, "phase_std": std_phase,
    }

    a = -math.log(avg_mag) / math.log(6) if avg_mag > 1e-20 else float('inf')
    print(f"\n  {class_name:>15s}:")
    print(f"    |λ_sign| = {avg_mag:.6f} ± {std_mag:.6f}")
    print(f"    a = -log₆|λ| = {a:.6f}")
    print(f"    phase = {avg_phase:.6f} rad = {avg_phase/math.pi:.6f}π")

# The √2 connection
l_id = all_data["identity"]
theta_crystal = (1 - K_STAR) / H  # 23/90
sqrt2_theta = math.sqrt(2) * theta_crystal

print(f"\n  √2 × θ_crystal = √2 × (1-K*)/H = {sqrt2_theta:.6f} rad = {sqrt2_theta/math.pi:.6f}π")
print(f"  Identity phase  = {l_id['phase']:.6f} rad = {l_id['phase']/math.pi:.6f}π")
print(f"  Ratio: phase/(√2·θ_c) = {l_id['phase']/sqrt2_theta:.6f}")
print(f"  Match: {abs(1 - l_id['phase']/sqrt2_theta)*100:.4f}%")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART B: PRECISE EXPONENTS (per-seed magnitude averaging)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

a0 = -math.log(all_data["identity"]["mag"]) / math.log(6)
a1 = -math.log(all_data["3-cycles"]["mag"]) / math.log(6)
a2 = -math.log(all_data["transpositions"]["mag"]) / math.log(6)

print(f"\n  Exponents |λ_sign| = 6^{{-a}}:")
print(f"    a₀(identity)       = {a0:.8f}  (1/3 = {1/3:.8f}, dev = {abs(a0-1/3)/(1/3)*100:.4f}%)")
print(f"    a₁(3-cycles)       = {a1:.8f}")
print(f"    a₂(transpositions) = {a2:.8f}  (59/30 = {59/30:.8f}, dev = {abs(a2-59/30)/(59/30)*100:.4f}%)")

print(f"\n  Relationships:")
print(f"    a₂ - a₀ = {a2-a0:.8f}  (59/30 - 1/3 = {59/30-1/3:.8f} = {(59-10)/30:.8f} = 49/30)")
print(f"    a₂ + a₀ = {a2+a0:.8f}  (59/30 + 1/3 = {59/30+1/3:.8f} = {(59+10)/30:.8f} = 69/30)")
print(f"    a₂ / a₀ = {a2/a0:.6f}  (59/10 = {59/10:.6f})")
print(f"    a₁ - a₀ = {a1-a0:.8f}")
print(f"    a₂ - a₁ = {a2-a1:.8f}")
print(f"    a₁/a₀ = {a1/a0:.6f}")
print(f"    (a₂-a₁)/(a₁-a₀) = {(a2-a1)/(a1-a0):.6f}")
print(f"    a₀ + a₁ + a₂ = {a0+a1+a2:.6f}")

# If a₀ = 1/3 and a₂ = 59/30, check what a₁ would need to be
# for various relations
print(f"\n  What determines a₁?")
print(f"    Geometric mean √(a₀·a₂) = {math.sqrt(1/3 * 59/30):.6f}  (a₁ = {a1:.6f}, dev = {abs(math.sqrt(1/3*59/30)-a1)/a1*100:.2f}%)")
print(f"    Arithmetic mean (a₀+a₂)/2 = {(1/3+59/30)/2:.6f}  (dev = {abs((1/3+59/30)/2-a1)/a1*100:.2f}%)")
print(f"    Harmonic mean = {2/((1/(1/3))+(1/(59/30))):.6f}  (dev = {abs(2/((1/(1/3))+(1/(59/30)))-a1)/a1*100:.2f}%)")

# Check: a₁ = (a₀ × h + a₂ × h') / h'' for Coxeter numbers?
# a₁ = (a₀ × 30 + something) / something_else?
# a₁ × 30 ≈ 21.5, a₁ × 90 ≈ 64.6

# Let me look for a₁ as a fraction with denominator 30 or 90
for d in [6, 10, 15, 18, 20, 30, 42, 60, 90, 120]:
    n = round(a1 * d)
    dev = abs(n/d - a1)
    if dev / a1 < 0.005:
        print(f"    {n}/{d} = {n/d:.8f}  (dev = {dev/a1*100:.4f}%)")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART C: STANDARD SECTOR — 5 DOUBLETS AND QUARKS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# The standard sector has dimension 10 = 5 doublets of the 2D standard rep.
# Question: do the 5 doublets correspond to 5 quarks (u,d,s,c,b)?
# The top quark is special (the only quark heavier than W) — it might be
# the remaining sign-sector direction.

# Measure standard sector spectrum per conjugacy class
print("\n  Standard sector eigenvalue spectra per class:")

for class_name, members in CONJUGACY_CLASSES.items():
    all_spectra = []
    for name in members:
        for seed in range(min(50, N_SEEDS)):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T = build_T_B(ent.joint)
            spec = std_spectrum(T)
            all_spectra.append(spec[:5].tolist() if len(spec) >= 5 else spec.tolist())

    # Average spectra
    max_len = max(len(s) for s in all_spectra)
    avg_spec = []
    for i in range(min(max_len, 5)):
        vals = [s[i] for s in all_spectra if len(s) > i]
        avg_spec.append(np.mean(vals))

    print(f"\n  {class_name:>15s}: {', '.join(f'{v:.4f}' for v in avg_spec)}")

    # Do the leading eigenvalues form Koide triplets?
    if len(avg_spec) >= 3:
        from itertools import combinations
        for i, j, k in combinations(range(min(len(avg_spec), 5)), 3):
            s = math.sqrt(avg_spec[i]) + math.sqrt(avg_spec[j]) + math.sqrt(avg_spec[k])
            if s > 1e-10:
                Q = (avg_spec[i] + avg_spec[j] + avg_spec[k]) / s**2
                if abs(Q - 2/3) < 0.05:
                    print(f"    Koide triplet ({i},{j},{k}): Q = {Q:.6f}  ← near 2/3")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART D: CKM MIXING FROM INDIVIDUAL ELEMENT FS DISTANCES")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Build individual element crystals
print("\n  Building individual S₃ element crystals (50 seeds each)...")

element_joints = {}
for name in S3_PERMS:
    joints = []
    for seed in range(50):
        ent = Entangler(S3_CORR[name], seed=seed).build()
        joints.append(ent.joint)
    # Average joint
    element_joints[name] = torch.stack(joints).mean(dim=0)

# FS distance matrix (born_fidelity is the Bhattacharyya coefficient)
print("\n  Born fidelity (Bhattacharyya) between individual elements:")
elements = list(S3_PERMS.keys())
print(f"  {'':>8s}", end="")
for e in elements:
    print(f"  {e:>8s}", end="")
print()

fid_matrix = {}
for e1 in elements:
    print(f"  {e1:>8s}", end="")
    for e2 in elements:
        f = born_fidelity(element_joints[e1], element_joints[e2])
        fid_matrix[(e1, e2)] = f
        print(f"  {f:>8.4f}", end="")
    print()

# FS distance = arccos(fidelity) — Bhattacharyya angle
print(f"\n  FS (Bhattacharyya) angle between elements:")
print(f"  {'':>8s}", end="")
for e in elements:
    print(f"  {e:>8s}", end="")
print()

for e1 in elements:
    print(f"  {e1:>8s}", end="")
    for e2 in elements:
        f = fid_matrix[(e1, e2)]
        angle = math.acos(min(f, 1.0))
        print(f"  {angle:>8.4f}", end="")
    print()

# Sym² distances
print(f"\n  Sym² (fingerprint) distances between elements:")
fps = {n: sym2_fingerprint(element_joints[n]) for n in elements}

print(f"  {'':>8s}", end="")
for e in elements:
    print(f"  {e:>8s}", end="")
print()

for e1 in elements:
    print(f"  {e1:>8s}", end="")
    for e2 in elements:
        d = sym2_distance(fps[e1], fps[e2])
        print(f"  {d:>8.4f}", end="")
    print()

# Now: do the inter-class distances correspond to CKM angles?
print(f"\n\n  CKM mixing angle comparison:")
print(f"  Physical Cabibbo angle: θ_C = {theta_C:.6f} rad = {math.degrees(theta_C):.4f}°")
print(f"  Physical |V_us| = {V_CKM['us']:.4f}, |V_cb| = {V_CKM['cb']:.4f}, |V_ub| = {V_CKM['ub']:.5f}")

# Average intra-class and inter-class distances
intra_distances = {}
inter_distances = {}

for c1_name, c1_members in CONJUGACY_CLASSES.items():
    for c2_name, c2_members in CONJUGACY_CLASSES.items():
        key = f"{c1_name}-{c2_name}"
        dists = []
        for e1 in c1_members:
            for e2 in c2_members:
                if e1 != e2:
                    f = fid_matrix[(e1, e2)]
                    dists.append(math.acos(min(f, 1.0)))
        if dists:
            avg = np.mean(dists)
            if c1_name == c2_name:
                intra_distances[c1_name] = avg
            else:
                inter_distances[key] = avg

print(f"\n  Intra-class FS angles:")
for name, d in intra_distances.items():
    print(f"    {name:>20s}: {d:.6f} rad = {math.degrees(d):.4f}°")

print(f"\n  Inter-class FS angles:")
for name, d in sorted(inter_distances.items(), key=lambda x: x[1]):
    print(f"    {name:>35s}: {d:.6f} rad = {math.degrees(d):.4f}°")

# The key comparison: inter-class FS distances as mixing angles
# Hypothesis: id-transposition distance ∝ Cabibbo angle
# id-3cycle distance ∝ V_cb-like angle
print(f"\n  Inter-class distances normalized by Cabibbo angle:")
for name, d in sorted(inter_distances.items(), key=lambda x: x[1]):
    ratio = d / theta_C
    print(f"    {name:>35s}: {d:.6f} / θ_C = {ratio:.4f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART E: INDIVIDUAL TRANSPOSITION FS DISTANCES → CKM ELEMENTS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# The three transpositions (01), (02), (12) correspond to three
# "directions" of flavor mixing. Their mutual FS distances might
# encode the CKM matrix structure.

trans = ["(01)", "(02)", "(12)"]
print(f"\n  Transposition mutual FS angles:")
for t1 in trans:
    for t2 in trans:
        if t1 < t2:
            f = fid_matrix[(t1, t2)]
            angle = math.acos(min(f, 1.0))
            print(f"    {t1} ↔ {t2}: {angle:.6f} rad = {math.degrees(angle):.4f}°")

# 3-cycle mutual distances
cyc = ["(012)", "(021)"]
f12 = fid_matrix[(cyc[0], cyc[1])]
angle12 = math.acos(min(f12, 1.0))
print(f"\n  3-cycle mutual FS angle:")
print(f"    (012) ↔ (021): {angle12:.6f} rad = {math.degrees(angle12):.4f}°")

# Cross: each transposition to each 3-cycle
print(f"\n  Transposition ↔ 3-cycle FS angles:")
for t in trans:
    for c in cyc:
        f = fid_matrix[(t, c)]
        angle = math.acos(min(f, 1.0))
        print(f"    {t} ↔ {c}: {angle:.6f} rad = {math.degrees(angle):.4f}°")

# Each transposition to identity
print(f"\n  Identity ↔ transposition FS angles:")
for t in trans:
    f = fid_matrix[("e", t)]
    angle = math.acos(min(f, 1.0))
    print(f"    e ↔ {t}: {angle:.6f} rad = {math.degrees(angle):.4f}°")

# Compare ratios with CKM hierarchy
print(f"\n  Wolfenstein hierarchy check:")
print(f"  λ_Wolfenstein = sin(θ_C) = {V_CKM['us']:.4f}")
print(f"  λ² = {V_CKM['us']**2:.4f}")
print(f"  λ³ = {V_CKM['us']**3:.5f}")
print(f"  |V_cb| = {V_CKM['cb']:.4f} ≈ λ² = {V_CKM['us']**2:.4f}")
print(f"  |V_ub| = {V_CKM['ub']:.5f} ≈ λ³ = {V_CKM['us']**3:.5f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART F: QUARK KOIDE — STANDARD SECTOR COUPLINGS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Physical quark Koide ratios
def koide_Q(m1, m2, m3):
    s = math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    return (m1 + m2 + m3) / s**2 if s > 1e-30 else float('nan')

print(f"\n  Physical quark Koide ratios:")
print(f"    Charged leptons (e,μ,τ):  Q = {koide_Q(m_e, m_mu, m_tau):.6f}")
print(f"    Heavy quarks (c,b,t):     Q = {koide_Q(m_c, m_b, m_t):.6f}")
print(f"    Up-type (u,c,t):          Q = {koide_Q(m_u, m_c, m_t):.6f}")
print(f"    Down-type (d,s,b):        Q = {koide_Q(m_d, m_s, m_b):.6f}")
print(f"    Light quarks (u,d,s):     Q = {koide_Q(m_u, m_d, m_s):.6f}")

# The crystal's standard sector has 5 doublets.
# Extract the 5 leading standard eigenvalues for each conjugacy class
print(f"\n  Standard sector eigenvalue extraction (per-seed):")

for class_name, members in CONJUGACY_CLASSES.items():
    all_mags = [[] for _ in range(5)]
    for name in members:
        for seed in range(100):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T = build_T_B(ent.joint)
            spec = std_spectrum(T)
            for i in range(min(5, len(spec))):
                all_mags[i].append(spec[i].item())

    print(f"\n  {class_name:>15s} standard eigenvalues:")
    avgs = []
    for i in range(5):
        if all_mags[i]:
            avg = np.mean(all_mags[i])
            avgs.append(avg)
            print(f"    σ_{i}: {avg:.6f} (6^{{-{-math.log(avg)/math.log(6):.4f}}})")

    # Koide triplets within the standard sector
    if len(avgs) >= 3:
        from itertools import combinations
        print(f"    Koide scan:")
        for i, j, k in combinations(range(len(avgs)), 3):
            Q = koide_Q(avgs[i], avgs[j], avgs[k])
            if abs(Q - 2/3) < 0.05:
                print(f"      ({i},{j},{k}): Q = {Q:.6f} ←")
            elif abs(Q - 2/3) < 0.1:
                print(f"      ({i},{j},{k}): Q = {Q:.6f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART G: SECTOR COUPLING CROSS-MATRIX")
print("How does each conjugacy class couple sign→standard→trivial?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print(f"\n  Sector coupling matrix (per-seed averaged):")
print(f"  {'class':>15s}  {'|λ_triv|':>10s}  {'|λ_sign|':>10s}  {'|λ_std₀|':>10s}  {'ratio s/t':>10s}")
print(f"  {'-'*60}")

for class_name, members in CONJUGACY_CLASSES.items():
    triv_mags = []
    sign_mags = []
    std_mags = []
    for name in members:
        for seed in range(100):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T = build_T_B(ent.joint)
            triv_mags.append(triv_coupling(T).abs().item())
            sign_mags.append(sign_coupling(T).abs().item())
            spec = std_spectrum(T)
            if len(spec) > 0:
                std_mags.append(spec[0].item())

    at = np.mean(triv_mags)
    as_ = np.mean(sign_mags)
    astd = np.mean(std_mags) if std_mags else 0
    ratio = as_ / at if at > 1e-10 else float('inf')
    print(f"  {class_name:>15s}  {at:>10.6f}  {as_:>10.6f}  {astd:>10.6f}  {ratio:>10.6f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("SUMMARY")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════
print("""
  Results collected. See output above.
""")
