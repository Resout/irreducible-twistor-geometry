"""
Representation-theoretic correction to the 2Δ spectral gap.

The composition gap is approximately 2Δ for bijections, but varies
by conjugacy class:
  identity:      2.15Δ
  transpositions: 2.03Δ
  3-cycles:       2.21Δ

Hypothesis: each S₃ irrep sector has its own spectral gap, and the
effective gap for a crystal is a weighted average:
  gap_eff = Σ_χ (fraction_χ × gap_χ)

This would make the conjugacy-class-specific rates a PREDICTION of
the irrep decomposition, not just a correlation.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, DELTA, schmidt_number
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)


# ═══════════════════════════════════════════════════════════════
#  Build the S₃ group representation on [4,4] joint masses
# ═══════════════════════════════════════════════════════════════

S3_PERMS = {
    "e":     [0, 1, 2],
    "(01)":  [1, 0, 2],
    "(02)":  [2, 1, 0],
    "(12)":  [0, 2, 1],
    "(012)": [1, 2, 0],
    "(021)": [2, 0, 1],
}

S3_CORRS = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

SIGNS = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}
CONJ = {"e": "e", "(01)": "trans", "(02)": "trans", "(12)": "trans",
        "(012)": "3-cyc", "(021)": "3-cyc"}

dim_sq = MASS_DIM * MASS_DIM  # 16

def permutation_matrix(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return P

def group_rep_16(perm):
    P = permutation_matrix(perm)
    return torch.kron(P, P)

# Build projectors
reps = {name: group_rep_16(perm) for name, perm in S3_PERMS.items()}

chi = {
    "trivial":  {"e": 1, "(01)": 1, "(02)": 1, "(12)": 1, "(012)": 1, "(021)": 1},
    "sign":     {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1},
    "standard": {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1},
}
dims = {"trivial": 1, "sign": 1, "standard": 2}

projectors = {}
for irrep_name in chi:
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for g_name in S3_PERMS:
        P += chi[irrep_name][g_name] * reps[g_name]
    P *= dims[irrep_name] / 6.0
    projectors[irrep_name] = P


# ═══════════════════════════════════════════════════════════════
#  For each S₃ element: compute composition operator restricted
#  to each irrep sector and measure its spectral gap
# ═══════════════════════════════════════════════════════════════

def build_composition_operator(B_joint):
    dim = MASS_DIM * MASS_DIM
    T = torch.zeros(dim, dim, dtype=torch.cfloat)
    for col in range(dim):
        A = torch.zeros(dim, dtype=torch.cfloat)
        A[col] = 1.0 + 0j
        A = A.reshape(MASS_DIM, MASS_DIM)
        result = compose(A, B_joint)
        T[:, col] = result.reshape(dim)
    return T

def spectral_gap_in_sector(T, P):
    """Spectral gap of T restricted to image(P)."""
    # T restricted to sector: P T P (acts on image of P)
    restricted = P @ T @ P
    ev = torch.linalg.eigvals(restricted)
    mags = ev.abs()
    idx = mags.argsort(descending=True)
    mags = mags[idx]

    # Filter near-zero
    nonzero = mags[mags > 1e-6]
    if len(nonzero) < 2:
        return 0, nonzero[0].item() if len(nonzero) > 0 else 0, 0

    # Find first significant drop
    l0 = nonzero[0].item()
    for i in range(1, len(nonzero)):
        if nonzero[i].item() < 0.95 * l0:
            l1 = nonzero[i].item()
            return math.log(l0 / l1), l0, l1

    return 0, l0, nonzero[-1].item()


print("=" * 80)
print("IRREP-RESOLVED SPECTRAL GAPS")
print("=" * 80)

sector_gaps = {}

for name, corr in S3_CORRS.items():
    ent = Entangler(corr, seed=42).build()
    T = build_composition_operator(ent.joint)

    # Overall gap
    overall_gap, l0_all, l1_all = spectral_gap_in_sector(T, torch.eye(dim_sq, dtype=torch.cfloat))

    # Per-sector gaps
    gaps = {}
    for irrep_name, P in projectors.items():
        gap, l0, l1 = spectral_gap_in_sector(T, P)
        gaps[irrep_name] = {"gap": gap, "l0": l0, "l1": l1}

    sector_gaps[name] = gaps

    print(f"\n  {name} (class={CONJ[name]}, sign={SIGNS[name]:+d}):")
    print(f"    Overall: gap = {overall_gap:.4f} = {overall_gap/DELTA:.3f}Δ")
    for irrep_name in ["trivial", "sign", "standard"]:
        g = gaps[irrep_name]
        if g["gap"] > 0:
            print(f"    {irrep_name:>10s}: gap = {g['gap']:.4f} = {g['gap']/DELTA:.3f}Δ  (λ₀={g['l0']:.4f}, λ₁={g['l1']:.4f})")
        else:
            print(f"    {irrep_name:>10s}: no gap (dim too small or degenerate)")


# ═══════════════════════════════════════════════════════════════
#  Test: weighted average of sector gaps = effective gap?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("WEIGHTED AVERAGE TEST")
print("="*80)

for name, corr in S3_CORRS.items():
    ent = Entangler(corr, seed=42).build()
    v = ent.joint.reshape(dim_sq)
    norm_total = v.abs().pow(2).sum().item()

    fracs = {}
    for irrep_name, P in projectors.items():
        v_proj = P @ v
        fracs[irrep_name] = v_proj.abs().pow(2).sum().item() / norm_total

    # Weighted gap
    gaps = sector_gaps[name]
    weighted = sum(fracs[ir] * gaps[ir]["gap"] for ir in gaps if gaps[ir]["gap"] > 0)

    # Actual eigenvalue decay rate (from composition)
    ev0 = torch.linalg.eigvals(ent.joint)
    idx0 = ev0.abs().argsort(descending=True)
    ev0 = ev0[idx0]

    composed = compose(ent.joint, ent.joint)
    ev1 = torch.linalg.eigvals(composed)
    idx1 = ev1.abs().argsort(descending=True)
    ev1 = ev1[idx1]

    r0 = (ev0[1].abs() / ev0[0].abs()).item()
    r1 = (ev1[1].abs() / ev1[0].abs()).item()
    actual_rate = math.log(r0 / r1) if r0 > 1e-12 and r1 > 1e-12 else 0

    print(f"\n  {name}:")
    print(f"    Fractions: triv={fracs['trivial']:.3f}, sign={fracs['sign']:.3f}, std={fracs['standard']:.3f}")
    print(f"    Weighted gap: {weighted:.4f} = {weighted/DELTA:.3f}Δ")
    print(f"    Actual rate:  {actual_rate:.4f} = {actual_rate/DELTA:.3f}Δ")
    print(f"    Match: {'YES' if abs(weighted - actual_rate) < 0.05 else 'NO'} (diff = {abs(weighted - actual_rate):.4f})")


# ═══════════════════════════════════════════════════════════════
#  Can we predict the rate from the irrep fractions alone?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("RATE PREDICTION FROM IRREP FRACTIONS")
print("="*80)

# Use the AVERAGE sector gap across all S₃ elements, then predict
# individual rates from individual fractions

avg_sector_gap = {}
for irrep_name in ["trivial", "sign", "standard"]:
    vals = [sector_gaps[n][irrep_name]["gap"] for n in S3_CORRS if sector_gaps[n][irrep_name]["gap"] > 0]
    if vals:
        avg_sector_gap[irrep_name] = sum(vals) / len(vals)
        print(f"  Average {irrep_name} gap: {avg_sector_gap[irrep_name]:.4f} = {avg_sector_gap[irrep_name]/DELTA:.3f}Δ")

print(f"\n  Predictions using average sector gaps:")
for name in S3_CORRS:
    ent = Entangler(S3_CORRS[name], seed=42).build()
    v = ent.joint.reshape(dim_sq)
    norm_total = v.abs().pow(2).sum().item()

    fracs = {}
    for irrep_name, P in projectors.items():
        v_proj = P @ v
        fracs[irrep_name] = v_proj.abs().pow(2).sum().item() / norm_total

    predicted = sum(fracs[ir] * avg_sector_gap.get(ir, 0) for ir in fracs)

    # Actual
    ev0 = torch.linalg.eigvals(ent.joint)
    idx0 = ev0.abs().argsort(descending=True)
    composed = compose(ent.joint, ent.joint)
    ev1 = torch.linalg.eigvals(composed)
    idx1 = ev1.abs().argsort(descending=True)
    r0 = (ev0[idx0[1]].abs() / ev0[idx0[0]].abs()).item()
    r1 = (ev1[idx1[1]].abs() / ev1[idx1[0]].abs()).item()
    actual = math.log(r0 / r1) if r0 > 1e-12 and r1 > 1e-12 else 0

    print(f"  {name:>6s}: predicted={predicted/DELTA:.3f}Δ, actual={actual/DELTA:.3f}Δ, diff={abs(predicted-actual)/DELTA:.3f}Δ")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
