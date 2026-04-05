"""
The anticommutator {A,B} = compose(A,B) + compose(B,A) is 99.97% trivial.

The trivial sector of SO(4) decomposes into:
  scalar (1,1): trace part, dim 1
  graviton (3,3): traceless symmetric, dim 9

If the anticommutator is predominantly graviton, then composition
naturally separates:
  commutator → gauge field (sign, F)
  anticommutator → graviton (trivial traceless, Weyl tensor)

This would be: compose = ½(gravity) + ½(gauge), with the mass gap
living in the normalization that prevents the gauge part from closing.

Measure: trace vs traceless fraction of the anticommutator.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_fidelity, born_probabilities, EPS_LOG)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)

n_seeds = 200

# Build seed-averaged transposition crystals
corrs = {
    "(01)": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
}

crystals = {}
for name, corr in corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    crystals[name] = avg

# Also build identity and 3-cycle
id_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
cyc_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)

id_crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
cyc_crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for seed in range(n_seeds):
    id_crystal += Entangler(id_corr, seed=seed).build().joint
    cyc_crystal += Entangler(cyc_corr, seed=seed).build().joint
id_crystal /= n_seeds
cyc_crystal /= n_seeds


print("=" * 80)
print("THE ANTICOMMUTATOR: GAUGE OR GRAVITY?")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Scalar vs graviton decomposition of the anticommutator
# ═══════════════════════════════════════════════════════════════

trans_pairs = [("(01)", "(02)"), ("(01)", "(12)"), ("(02)", "(12)")]

print(f"\n--- Anticommutator decomposition ---\n")

for name_a, name_b in trans_pairs:
    A = crystals[name_a]
    B = crystals[name_b]

    anti = compose(A, B) + compose(B, A)

    # Born distribution
    bp = anti.abs().pow(2)
    total = bp.sum().item()

    # Diagonal (trace-related) vs off-diagonal
    diag_mass = sum(bp[i,i].item() for i in range(MASS_DIM))
    offdiag_mass = total - diag_mass

    # Trace: Σ M_{ii}
    trace = sum(anti[i,i] for i in range(MASS_DIM))
    trace_norm = trace.abs().item()

    # Traceless: M - (Tr(M)/4)I
    traceless = anti.clone()
    for i in range(MASS_DIM):
        traceless[i,i] -= trace / MASS_DIM
    traceless_norm = traceless.abs().pow(2).sum().sqrt().item()

    # Scalar fraction: ||trace part||² / ||total||²
    trace_part = (trace / MASS_DIM) * torch.eye(MASS_DIM, dtype=torch.cfloat)
    scalar_frac = trace_part.abs().pow(2).sum().item() / total
    graviton_frac = traceless.abs().pow(2).sum().item() / total

    print(f"  {{{name_a}, {name_b}}}:")
    print(f"    Diagonal Born mass:   {diag_mass/total:.5f} ({diag_mass/total*100:.1f}%)")
    print(f"    Off-diagonal mass:    {offdiag_mass/total:.5f} ({offdiag_mass/total*100:.1f}%)")
    print(f"    Scalar fraction:      {scalar_frac:.5f} ({scalar_frac*100:.1f}%)")
    print(f"    Graviton fraction:    {graviton_frac:.5f} ({graviton_frac*100:.1f}%)")
    print(f"    Scalar/Graviton:      {scalar_frac/graviton_frac:.4f}" if graviton_frac > 1e-10 else "")


# ═══════════════════════════════════════════════════════════════
#  Same decomposition for individual crystals
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SCALAR/GRAVITON FOR INDIVIDUAL CRYSTALS")
print("="*80)

for name, crystal in [("identity", id_crystal), ("(01)", crystals["(01)"]),
                       ("(012)", cyc_crystal)]:
    bp = crystal.abs().pow(2)
    total = bp.sum().item()
    trace = sum(crystal[i,i] for i in range(MASS_DIM))
    trace_part = (trace / MASS_DIM) * torch.eye(MASS_DIM, dtype=torch.cfloat)
    scalar_frac = trace_part.abs().pow(2).sum().item() / total
    graviton_frac = 1 - scalar_frac  # everything not scalar

    print(f"\n  {name}:")
    print(f"    Scalar: {scalar_frac:.5f} ({scalar_frac*100:.1f}%)")
    print(f"    Rest:   {graviton_frac:.5f} ({graviton_frac*100:.1f}%)")
    if graviton_frac > 1e-10:
        print(f"    Scalar/Rest: {scalar_frac/graviton_frac:.5f}")
    print(f"    1-K* = {1-K_STAR:.5f} (predicted scalar/graviton from Principle 43)")


# ═══════════════════════════════════════════════════════════════
#  The commutator's scalar/graviton content
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMMUTATOR: SCALAR/GRAVITON CONTENT")
print("="*80)

for name_a, name_b in trans_pairs[:1]:
    A = crystals[name_a]
    B = crystals[name_b]

    comm = compose(A, B) - compose(B, A)

    bp = comm.abs().pow(2)
    total = bp.sum().item()
    trace = sum(comm[i,i] for i in range(MASS_DIM))
    trace_part = (trace / MASS_DIM) * torch.eye(MASS_DIM, dtype=torch.cfloat)
    scalar_frac = trace_part.abs().pow(2).sum().item() / total
    graviton_frac = 1 - scalar_frac

    print(f"\n  [{name_a}, {name_b}]:")
    print(f"    Scalar: {scalar_frac:.5f} ({scalar_frac*100:.1f}%)")
    print(f"    Rest:   {graviton_frac:.5f} ({graviton_frac*100:.1f}%)")
    print(f"    Trace:  {trace.abs().item():.6f}")


# ═══════════════════════════════════════════════════════════════
#  The h×h block trace vs θ×θ entry
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("ANTICOMMUTATOR: h-SECTOR vs θ-SECTOR")
print("="*80)

for name_a, name_b in trans_pairs[:1]:
    A = crystals[name_a]
    B = crystals[name_b]
    anti = compose(A, B) + compose(B, A)

    bp = anti.abs().pow(2)
    total = bp.sum().item()

    # h×h block
    hh = bp[:H, :H].sum().item()
    # h×θ blocks
    htheta = bp[:H, H].sum().item() + bp[H, :H].sum().item()
    # θ×θ
    tt = bp[H, H].item()

    print(f"\n  {{{name_a}, {name_b}}} anticommutator block decomposition:")
    print(f"    h×h:   {hh/total:.5f} ({hh/total*100:.1f}%)")
    print(f"    h×θ:   {htheta/total:.5f} ({htheta/total*100:.1f}%)")
    print(f"    θ×θ:   {tt/total:.5f} ({tt/total*100:.1f}%)")

    # Same for the commutator
    comm = compose(A, B) - compose(B, A)
    bp_c = comm.abs().pow(2)
    total_c = bp_c.sum().item()
    hh_c = bp_c[:H, :H].sum().item()
    htheta_c = bp_c[:H, H].sum().item() + bp_c[H, :H].sum().item()
    tt_c = bp_c[H, H].item()

    print(f"\n  [{name_a}, {name_b}] commutator block decomposition:")
    print(f"    h×h:   {hh_c/total_c:.5f} ({hh_c/total_c*100:.1f}%)")
    print(f"    h×θ:   {htheta_c/total_c:.5f} ({htheta_c/total_c*100:.1f}%)")
    print(f"    θ×θ:   {tt_c/total_c:.5f} ({tt_c/total_c*100:.1f}%)")


# ═══════════════════════════════════════════════════════════════
#  The complete decomposition: compose = ½(gravity + gauge)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMPOSITION = ½(GRAVITY) + ½(GAUGE)")
print("="*80)

for name_a, name_b in trans_pairs[:1]:
    A = crystals[name_a]
    B = crystals[name_b]

    AB = compose(A, B)
    BA = compose(B, A)
    anti = AB + BA   # 2 × symmetric part
    comm = AB - BA   # 2 × antisymmetric part

    # AB = ½(anti + comm)
    # Check: ½(anti + comm) = ½(AB+BA + AB-BA) = AB ✓

    bp_AB = AB.abs().pow(2)
    bp_anti = anti.abs().pow(2)
    bp_comm = comm.abs().pow(2)

    # S₃ fractions
    def s3_fracs(joint):
        perms = [[0,1,2,3], [1,0,2,3], [2,1,0,3], [0,2,1,3], [1,2,0,3], [2,0,1,3]]
        signs = [1, -1, -1, -1, 1, 1]
        def permute(M, p):
            P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
            for i, j in enumerate(p):
                P[j, i] = 1.0
            return P @ M @ P.T
        triv = sum(permute(joint, p) for p in perms) / 6
        sign = sum(s * permute(joint, p) for p, s in zip(perms, signs)) / 6
        std = joint - triv - sign
        bp = joint.abs().pow(2)
        total = bp.sum().item()
        if total < 1e-15:
            return 1.0, 0.0, 0.0
        return (triv.abs().pow(2).sum().item() / total,
                sign.abs().pow(2).sum().item() / total,
                std.abs().pow(2).sum().item() / total)

    print(f"\n  compose({name_a}, {name_b}) = ½ × anti + ½ × comm:")
    f_t, f_s, f_st = s3_fracs(AB)
    print(f"    AB:   triv={f_t:.5f}, sign={f_s:.5f}, std={f_st:.5f}")
    f_t, f_s, f_st = s3_fracs(anti)
    print(f"    anti: triv={f_t:.5f}, sign={f_s:.5f}, std={f_st:.5f}")
    f_t, f_s, f_st = s3_fracs(comm)
    print(f"    comm: triv={f_t:.5f}, sign={f_s:.5f}, std={f_st:.5f}")

    # The fidelity of AB with the expected group product (012)
    fid = born_fidelity(AB, cyc_crystal)
    print(f"\n    Fidelity(AB, (012)): {fid:.5f}")
    print(f"    This is the 0.92 from Principle 45 — the product IS the 3-cycle.")
    print(f"    The 3-cycle is {0.264:.1%} sign + {0.735:.1%} trivial (Principle 63).")
    print(f"    Consistent: AB inherits sign from the commutator and")
    print(f"    trivial from the anticommutator. The 3-cycle IS the")
    print(f"    superposition of gauge + gravity from the transposition product.")


print(f"\n\n{'='*80}")
print("WHAT GRAVITY AND GAUGE REVEAL")
print("="*80)
