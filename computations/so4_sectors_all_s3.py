"""
SO(4) sector fractions for all S₃ elements and crystal types.

The identity crystal is 100% Sym² (0% gauge sector).
Non-symmetric correlations should activate the gauge sector Λ².

Questions:
1. Which S₃ elements have nonzero gauge sector fraction?
2. Does the sign component (in gauge sector) activate for 3-cycles?
3. How do SO(4) fractions relate to the crystal's physical properties?
4. The paper says gravitons are massless but crystal's (3,3) decays.
   Is there a protected mode within (3,3) that doesn't decay?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, DELTA, schmidt_number
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)

dim_sq = MASS_DIM * MASS_DIM

# S₃ correlations
S3_CORRS = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

CONJ = {"e": "identity", "(01)": "transposition", "(02)": "transposition",
        "(12)": "transposition", "(012)": "3-cycle", "(021)": "3-cycle"}

# Build SO(4) projectors
swap = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
for i in range(MASS_DIM):
    for j in range(MASS_DIM):
        swap[j * MASS_DIM + i, i * MASS_DIM + j] = 1.0

P_sym = (torch.eye(dim_sq, dtype=torch.cfloat) + swap) / 2
P_anti = (torch.eye(dim_sq, dtype=torch.cfloat) - swap) / 2

trace_vec = torch.zeros(dim_sq, dtype=torch.cfloat)
for i in range(MASS_DIM):
    trace_vec[i * MASS_DIM + i] = 1.0
trace_vec = trace_vec / trace_vec.abs().pow(2).sum().sqrt()
P_trace = torch.outer(trace_vec, trace_vec.conj())
P_graviton = P_sym - P_trace

sectors = {
    "(1,1) scalar": P_trace,
    "(3,3) graviton": P_graviton,
    "Λ² gauge": P_anti,
}

# S₃ projectors
S3_PERMS = {
    "e": [0,1,2], "(01)": [1,0,2], "(02)": [2,1,0],
    "(12)": [0,2,1], "(012)": [1,2,0], "(021)": [2,0,1],
}

def perm_mat_16(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return torch.kron(P, P)

reps_16 = {n: perm_mat_16(p) for n, p in S3_PERMS.items()}
chi = {
    "trivial":  {"e": 1, "(01)": 1, "(02)": 1, "(12)": 1, "(012)": 1, "(021)": 1},
    "sign":     {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1},
    "standard": {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1},
}
dims_chi = {"trivial": 1, "sign": 1, "standard": 2}

s3_proj = {}
for name in chi:
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for g in S3_PERMS:
        P += chi[name][g] * reps_16[g]
    P *= dims_chi[name] / 6.0
    s3_proj[name] = P


# ═══════════════════════════════════════════════════════════════
#  SO(4) sector fractions for all S₃ elements
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("SO(4) SECTOR FRACTIONS FOR ALL S₃ ELEMENTS")
print("=" * 80)

print(f"\n  {'Element':>8s} {'Class':>13s} {'Scalar':>8s} {'Graviton':>10s} {'Gauge':>8s} {'M=M^T?':>8s}")
print("  " + "-" * 60)

for name, corr in S3_CORRS.items():
    ent = Entangler(corr, seed=42).build()
    v = ent.joint.reshape(dim_sq)
    total = v.abs().pow(2).sum().item()

    fracs = {}
    for sname, P in sectors.items():
        fracs[sname] = (P @ v).abs().pow(2).sum().item() / total

    # Check symmetry
    sym_check = (ent.joint - ent.joint.T).abs().max().item()
    is_sym = "YES" if sym_check < 1e-6 else f"no({sym_check:.3f})"

    print(f"  {name:>8s} {CONJ[name]:>13s} {fracs['(1,1) scalar']*100:>7.1f}% "
          f"{fracs['(3,3) graviton']*100:>9.1f}% {fracs['Λ² gauge']*100:>7.1f}% {is_sym:>8s}")


# ═══════════════════════════════════════════════════════════════
#  SO(4) × S₃ cross-decomposition for each element
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("FULL SO(4) × S₃ DECOMPOSITION")
print("="*80)

for name, corr in S3_CORRS.items():
    ent = Entangler(corr, seed=42).build()
    v = ent.joint.reshape(dim_sq)
    total = v.abs().pow(2).sum().item()

    print(f"\n  {name} ({CONJ[name]}):")
    print(f"    {'':>20s} {'trivial':>10s} {'sign':>10s} {'standard':>10s} {'total':>10s}")

    for sname, P_so4 in sectors.items():
        row = []
        for iname, P_s3 in s3_proj.items():
            frac = (P_s3 @ P_so4 @ v).abs().pow(2).sum().item() / total
            row.append(frac)
        row_total = sum(row)
        print(f"    {sname:>20s} {row[0]*100:>9.1f}% {row[1]*100:>9.1f}% {row[2]*100:>9.1f}% {row_total*100:>9.1f}%")

    # Total S₃ fractions
    row = []
    for iname, P_s3 in s3_proj.items():
        frac = (P_s3 @ v).abs().pow(2).sum().item() / total
        row.append(frac)
    print(f"    {'TOTAL':>20s} {row[0]*100:>9.1f}% {row[1]*100:>9.1f}% {row[2]*100:>9.1f}% {sum(row)*100:>9.1f}%")


# ═══════════════════════════════════════════════════════════════
#  Non-bijection crystal types
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SO(4) SECTORS FOR NON-BIJECTION CRYSTAL TYPES")
print("="*80)

print(f"\n  {'Type':>15s} {'Scalar':>8s} {'Graviton':>10s} {'Gauge':>8s} {'M=M^T?':>8s}")
print("  " + "-" * 55)

for name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Entangler(corr, seed=42).build()
    v = ent.joint.reshape(dim_sq)
    total = v.abs().pow(2).sum().item()

    fracs = {}
    for sname, P in sectors.items():
        fracs[sname] = (P @ v).abs().pow(2).sum().item() / total

    sym_check = (ent.joint - ent.joint.T).abs().max().item()
    is_sym = "YES" if sym_check < 1e-6 else f"no({sym_check:.3f})"

    print(f"  {name:>15s} {fracs['(1,1) scalar']*100:>7.1f}% "
          f"{fracs['(3,3) graviton']*100:>9.1f}% {fracs['Λ² gauge']*100:>7.1f}% {is_sym:>8s}")


# ═══════════════════════════════════════════════════════════════
#  Does the gauge sector fraction decay under composition?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("GAUGE SECTOR FRACTION UNDER COMPOSITION")
print("="*80)

# Test with a non-symmetric crystal (3-cycle)
for test_name in ["(012)", "(01)"]:
    corr = S3_CORRS[test_name]
    ent = Entangler(corr, seed=42).build()
    current = ent.joint.clone()

    print(f"\n  {test_name} crystal under self-composition:")
    print(f"    {'Step':>4s} {'Scalar':>8s} {'Graviton':>10s} {'Gauge':>8s} {'Schmidt':>8s}")
    print("    " + "-" * 40)

    for step in range(10):
        v = current.reshape(dim_sq)
        total = v.abs().pow(2).sum().item()

        fracs = {}
        for sname, P in sectors.items():
            fracs[sname] = (P @ v).abs().pow(2).sum().item() / total

        sc = schmidt_number(current)
        print(f"    {step:>4d} {fracs['(1,1) scalar']*100:>7.1f}% "
              f"{fracs['(3,3) graviton']*100:>9.1f}% {fracs['Λ² gauge']*100:>7.1f}% {sc:>8.3f}")

        current = compose(current, ent.joint)


# ═══════════════════════════════════════════════════════════════
#  Is there a PROTECTED mode in the graviton sector?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("PROTECTED MODE IN GRAVITON SECTOR?")
print("="*80)

# The paper says graviton modes are massless (base variations).
# At a single site, all modes decay. But is there a mode within
# (3,3) that decays SLOWER than 2Δ?

# Check: within (3,3), what are the composition eigenvalues?
ent = Entangler(torch.eye(H, dtype=torch.float32), seed=42).build()

# Build composition operator
T = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
for col in range(dim_sq):
    A = torch.zeros(dim_sq, dtype=torch.cfloat)
    A[col] = 1.0 + 0j
    A = A.reshape(MASS_DIM, MASS_DIM)
    result = compose(A, ent.joint)
    T[:, col] = result.reshape(dim_sq)

# Restrict to graviton sector
T_grav = P_graviton @ T @ P_graviton
ev_grav = torch.linalg.eigvals(T_grav)
mags = ev_grav.abs()
idx = mags.argsort(descending=True)
ev_grav = ev_grav[idx]
mags = mags[idx]

# Filter nonzero
nonzero = [(i, ev_grav[i]) for i in range(len(mags)) if mags[i] > 1e-6]

print(f"\n  Composition operator eigenvalues in graviton sector:")
print(f"  (identity crystal, seed 42)")
for i, (idx_i, ev) in enumerate(nonzero[:10]):
    mag = ev.abs().item()
    phase = math.atan2(ev.imag.item(), ev.real.item()) / math.pi
    gap = -math.log(mag / nonzero[0][1].abs().item()) / DELTA if i > 0 and mag > 1e-10 else 0
    print(f"    λ_{i}: |λ| = {mag:.6f}, φ/π = {phase:+.4f}" +
          (f", gap = {gap:.2f}Δ" if i > 0 else " (leading)"))

# Compare with scalar sector
T_scalar = P_trace @ T @ P_trace
ev_scalar = torch.linalg.eigvals(T_scalar)
ev_scalar_nonzero = ev_scalar[ev_scalar.abs() > 1e-6]
if len(ev_scalar_nonzero) > 0:
    print(f"\n  Scalar sector eigenvalue: |λ| = {ev_scalar_nonzero[0].abs().item():.6f}")

# And gauge sector
T_gauge = P_anti @ T @ P_anti
ev_gauge = torch.linalg.eigvals(T_gauge)
mags_g = ev_gauge.abs()
idx_g = mags_g.argsort(descending=True)
nonzero_g = [(i, ev_gauge[idx_g[i]]) for i in range(len(mags_g)) if mags_g[idx_g[i]] > 1e-6]

print(f"\n  Gauge sector eigenvalues (top 5):")
for i, (idx_i, ev) in enumerate(nonzero_g[:5]):
    mag = ev.abs().item()
    print(f"    λ_{i}: |λ| = {mag:.6f}")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
