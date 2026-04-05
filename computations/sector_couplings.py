"""
Sector couplings: how the composition operator acts WITHIN and BETWEEN sectors.

The sign sector is 1-dimensional. T_B restricted to it is a scalar λ_sign(B).
This scalar IS the chiral coupling of crystal B: how strongly B interacts
with the fermion sector.

The standard sector is 10-dimensional. T_B restricted to it has 10 eigenvalues.
These determine how B interacts with gauge fields.

The CROSS-SECTOR couplings: how T_B maps one sector into another. These are
the "forbidden" transitions (sector-changing interactions).

Key questions:
  1. What are the chiral couplings λ_sign(B) for each S₃ element?
  2. Do the ratios of chiral couplings give mass ratios?
  3. What are the sector-changing amplitudes (off-diagonal blocks)?
  4. Does standard × standard → sign give pair production rates?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import numpy as np
from solver.algebra import H, MASS_DIM, K_STAR, BORN_FLOOR, DELTA
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

dim_sq = MASS_DIM ** 2  # 16

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

# Build projectors
def perm_mat_4x4(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return P

reps = {n: torch.kron(perm_mat_4x4(p), perm_mat_4x4(p)) for n, p in S3_PERMS.items()}

chi_t = {"e": 1, "(01)": 1, "(02)": 1, "(12)": 1, "(012)": 1, "(021)": 1}
chi_s = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}
chi_d = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}

def proj(chi, dim_chi):
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for n in S3_PERMS:
        P += chi[n] * reps[n]
    return P * dim_chi / 6.0

P_triv = proj(chi_t, 1)
P_sign = proj(chi_s, 1)
P_std  = proj(chi_d, 2)


def build_T_B(B_joint):
    """16×16 composition operator T_B: A → compose(A, B)."""
    T = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for col in range(dim_sq):
        A = torch.zeros(dim_sq, dtype=torch.cfloat)
        A[col] = 1.0
        result = compose(A.reshape(MASS_DIM, MASS_DIM), B_joint)
        T[:, col] = result.reshape(dim_sq)
    return T


# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("PART 1: CHIRAL COUPLING λ_sign(B) FOR EACH S₃ ELEMENT")
print("(sign sector is 1D → T_B restricted to sign is a scalar)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

N_SEEDS = 50

for name in S3_PERMS:
    corr = S3_CORR[name]
    sign_couplings = []

    for seed in range(N_SEEDS):
        ent = Entangler(corr, seed=seed).build()
        T = build_T_B(ent.joint)

        # Restrict to sign sector: λ_sign = (P_sign T P_sign) restricted to im(P_sign)
        # Since sign is 1D, this is a scalar
        restricted = P_sign @ T @ P_sign
        # The non-zero eigenvalue
        evals = torch.linalg.eigvals(restricted)
        nonzero = evals[evals.abs() > 1e-10]
        if len(nonzero) > 0:
            sign_couplings.append(nonzero[0])
        else:
            sign_couplings.append(torch.tensor(0.0 + 0j))

    mags = [c.abs().item() for c in sign_couplings]
    phases = [math.atan2(c.imag.item(), c.real.item()) / math.pi for c in sign_couplings]

    print(f"\n  {name:>6s}: |λ_sign| = {np.mean(mags):.6f} ± {np.std(mags):.6f}"
          f"  φ/π = {np.mean(phases):+.4f} ± {np.std(phases):.4f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 2: SECTOR-RESTRICTED SPECTRA FOR ALL S₃ ELEMENTS")
print("(eigenvalues of P_sector T_B P_sector for each B)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Use seed=42 for each
for name in S3_PERMS:
    ent = Entangler(S3_CORR[name], seed=42).build()
    T = build_T_B(ent.joint)

    print(f"\n  {name}:")
    for label, P_sec in [("trivial", P_triv), ("sign", P_sign), ("standard", P_std)]:
        restricted = P_sec @ T @ P_sec
        evals = torch.linalg.eigvals(restricted)
        nonzero = evals[evals.abs() > 1e-8]
        mags = nonzero.abs()
        idx = mags.argsort(descending=True)
        nonzero = nonzero[idx]

        if len(nonzero) == 0:
            print(f"    {label:>10s}: (empty)")
            continue

        top = min(4, len(nonzero))
        vals = "  ".join(f"|{nonzero[i].abs().item():.4f}|∠{math.atan2(nonzero[i].imag.item(), nonzero[i].real.item())/math.pi:+.3f}π"
                        for i in range(top))
        print(f"    {label:>10s}: {vals}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 3: CROSS-SECTOR AMPLITUDES")
print("How much does T_B map between sectors?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# The cross-sector amplitude from sector i to sector j under T_B:
# A_{i→j}(B) = ||P_j T_B P_i|| / ||P_i||
# This measures how much composing with B changes the sector content

for name in S3_PERMS:
    ent = Entangler(S3_CORR[name], seed=42).build()
    T = build_T_B(ent.joint)

    print(f"\n  {name} — sector transfer amplitudes:")
    header = "from\\to"
    print(f"  {header:>12s} {'trivial':>10s} {'sign':>10s} {'standard':>10s}")

    for from_label, P_from in [("trivial", P_triv), ("sign", P_sign), ("standard", P_std)]:
        row = f"  {from_label:>12s}"
        for to_label, P_to in [("trivial", P_triv), ("sign", P_sign), ("standard", P_std)]:
            cross = P_to @ T @ P_from
            amp = cross.abs().pow(2).sum().sqrt().item()
            # Normalize by source sector dimension
            src_dim = torch.linalg.matrix_rank(P_from.real.float()).item()
            row += f" {amp/src_dim:>10.6f}"
        print(row)


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 4: THE SIGN COUPLING AS MASS GENERATION")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

S = H**3 / K_STAR

# If the chiral coupling λ_sign determines mass, then:
# mass ∝ |λ_sign|^n for n compositions
# Different conjugacy classes have different |λ_sign|

# Seed-averaged chiral couplings by class
class_couplings = {}
for class_name, members in [("identity", ["e"]),
                             ("transpositions", ["(01)", "(02)", "(12)"]),
                             ("3-cycles", ["(012)", "(021)"])]:
    all_mags = []
    for name in members:
        for seed in range(N_SEEDS):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T = build_T_B(ent.joint)
            restricted = P_sign @ T @ P_sign
            evals = torch.linalg.eigvals(restricted)
            nonzero = evals[evals.abs() > 1e-10]
            if len(nonzero) > 0:
                all_mags.append(nonzero[0].abs().item())

    avg_mag = np.mean(all_mags) if all_mags else 0
    std_mag = np.std(all_mags) if all_mags else 0
    class_couplings[class_name] = avg_mag

    print(f"  {class_name:>15s}: |λ_sign| = {avg_mag:.6f} ± {std_mag:.6f}"
          f"  -ln|λ| = {-math.log(avg_mag) if avg_mag > 0 else float('inf'):.4f}")

# Ratios
if all(v > 0 for v in class_couplings.values()):
    le = class_couplings["identity"]
    lt = class_couplings["transpositions"]
    l3 = class_couplings["3-cycles"]

    print(f"\n  Ratios:")
    print(f"    |λ_trans/λ_identity| = {lt/le:.6f}")
    print(f"    |λ_3cyc/λ_identity| = {l3/le:.6f}")
    print(f"    |λ_3cyc/λ_trans|    = {l3/lt:.6f}")

    # If mass ∝ |λ_sign|^n:
    # The chiral coupling DECREASES mass (smaller coupling = lighter)
    # Heaviest class has largest |λ_sign|
    print(f"\n  Mass ordering (larger |λ_sign| = heavier):")
    ordered = sorted(class_couplings.items(), key=lambda x: x[1], reverse=True)
    for cn, cv in ordered:
        print(f"    {cn:>15s}: |λ_sign| = {cv:.6f}")

    # Compare with physical mass ratios
    print(f"\n  If three generations map to three classes:")
    print(f"    m_τ/m_e ∝ (|λ_heavy|/|λ_light|)^n")
    ratio_heavy_light = ordered[0][1] / ordered[-1][1]
    m_tau_m_e = 1777/0.511
    n_needed = math.log(m_tau_m_e) / math.log(ratio_heavy_light) if ratio_heavy_light > 1 else float('inf')
    print(f"    Ratio of extremes: {ratio_heavy_light:.4f}")
    print(f"    m_τ/m_e = {m_tau_m_e:.1f}")
    print(f"    n needed: {n_needed:.1f} compositions")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 5: FUSION RULES — SECTOR × SECTOR → SECTOR")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# S₃ tensor product (fusion) rules:
# trivial ⊗ X = X for any X
# sign ⊗ sign = trivial
# sign ⊗ standard = standard
# standard ⊗ standard = trivial ⊕ sign ⊕ standard

print(f"""
  S₃ fusion rules (from representation ring):

  ⊗           | trivial | sign    | standard
  ------------+---------+---------+----------
  trivial     | trivial | sign    | standard
  sign        | sign    | trivial | standard
  standard    | standard| standard| triv ⊕ sign ⊕ std

  Physical interpretation:
    gravity × gravity    = gravity       (gravity is self-coupling)
    chirality × chirality = gravity      (fermion-antifermion → graviton)
    chirality × gauge     = gauge        (charged fermion = gauge + chiral)
    gauge × gauge         = ALL SECTORS  (gluon-gluon → everything)

  The crucial fusion: standard ⊗ standard ⊃ sign
  → gauge × gauge produces chirality (pair production!)

  Verify numerically: compose two transposition crystals, measure sign content.
""")

# Compose two transposition crystals and check sign content
for n1 in ["(01)", "(02)", "(12)"]:
    for n2 in ["(01)", "(02)", "(12)"]:
        j1 = Entangler(S3_CORR[n1], seed=42).build().joint
        j2 = Entangler(S3_CORR[n2], seed=42).build().joint
        result = compose(j1, j2)

        v = result.reshape(dim_sq)
        norm = v.abs().pow(2).sum().item()
        fs = (P_sign @ v).abs().pow(2).sum().item() / norm if norm > 0 else 0

        if n1 == n2:
            # Same transposition: should give identity (pure trivial)
            ft = (P_triv @ v).abs().pow(2).sum().item() / norm if norm > 0 else 0
            print(f"  {n1} ○ {n2}: sign = {fs:.4f}, trivial = {ft:.4f}")
        elif fs > 0.01:
            print(f"  {n1} ○ {n2}: sign = {fs:.4f} ← PAIR PRODUCTION")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 6: THE STANDARD SECTOR EIGENVALUE SPLITTING")
print("(10 eigenvalues in the standard sector — do they split into generations?)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# The standard sector has 10 dimensions = 5 copies of the 2D standard rep.
# Under composition, these might split into distinct groups.

ent_id = Entangler(S3_CORR["e"], seed=42).build()
T_e = build_T_B(ent_id.joint)
restricted_std = P_std @ T_e @ P_std
evals_std = torch.linalg.eigvals(restricted_std)
nonzero_std = evals_std[evals_std.abs() > 1e-8]
mags_std = nonzero_std.abs()
idx = mags_std.argsort(descending=True)
nonzero_std = nonzero_std[idx]

print(f"\n  T_e restricted to standard sector ({len(nonzero_std)} eigenvalues):")
print(f"  {'idx':>5s} {'|λ|':>10s} {'φ/π':>10s} {'|λ|/|λ₀|':>10s} {'gap':>10s}")

prev_mag = None
for i in range(len(nonzero_std)):
    ev = nonzero_std[i]
    mag = ev.abs().item()
    phase = math.atan2(ev.imag.item(), ev.real.item()) / math.pi
    ratio = mag / nonzero_std[0].abs().item()
    gap = f"{-math.log(mag/prev_mag):.4f}" if prev_mag is not None and mag > 0 else ""
    prev_mag = mag
    print(f"  {i:>5d} {mag:>10.6f} {phase:>+10.4f} {ratio:>10.6f} {gap:>10s}")

# Look for groupings
print(f"\n  Eigenvalue magnitudes cluster at:")
unique_mags = []
for ev in nonzero_std:
    mag = ev.abs().item()
    found = False
    for um in unique_mags:
        if abs(mag - um[0]) / mag < 0.01:
            um[1] += 1
            found = True
            break
    if not found:
        unique_mags.append([mag, 1])

for mag, count in sorted(unique_mags, key=lambda x: -x[0]):
    print(f"    |λ| ≈ {mag:.6f} (×{count})")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("FINDINGS")
print("=" * 80)
