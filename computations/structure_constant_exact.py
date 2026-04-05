"""
Is the structure constant N = √(4/3) exact?

From structure_constants.py: N ≈ 1.161 across all transposition pairs.
Candidates:
  √(4/3) ≈ 1.15470
  7/6    ≈ 1.16667
  √(27/20) ≈ 1.16190
  √(Δ_comp/Δ_DS) where Δ_comp = 2Δ

Approach:
  1. High-precision computation (many seeds)
  2. Test N² against simple fractions
  3. Check if N depends on entangler parameters (discount, iterations)
  4. Derive from L1 norms and bracket structure
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, EPS_LOG,
                            schmidt_number, born_fidelity)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


def s3_project(joint):
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    signs = [1, -1, -1, -1, 1, 1]
    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T
    triv = sum(permute(joint, p) for p in perms) / 6
    sign = sum(s * permute(joint, p) for p, s in zip(perms, signs)) / 6
    std = joint - triv - sign
    return triv, sign, std


# ═══════════════════════════════════════════════════════════════
#  High-precision N measurement
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("STRUCTURE CONSTANT N: HIGH-PRECISION MEASUREMENT")
print("=" * 80)

corrs = {
    "(01)": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

# Test with increasing seed counts
for n_seeds in [200, 500, 1000, 2000]:
    crystals = {}
    for name, corr in corrs.items():
        avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
        for seed in range(n_seeds):
            avg += Entangler(corr, seed=seed).build().joint
        avg /= n_seeds
        crystals[name] = avg

    # Measure N for each transposition pair
    Ns = []
    for name_a, name_b in [("(01)", "(02)"), ("(01)", "(12)"), ("(02)", "(12)")]:
        comm = compose(crystals[name_a], crystals[name_b]) - compose(crystals[name_b], crystals[name_a])
        _, sign_comm, _ = s3_project(comm)
        comm_sign_norm = sign_comm.abs().pow(2).sum().sqrt().item()

        _, sign_012, _ = s3_project(crystals["(012)"])
        _, sign_021, _ = s3_project(crystals["(021)"])
        avg_cycle_sign = (sign_012.abs().pow(2).sum().sqrt().item() +
                         sign_021.abs().pow(2).sum().sqrt().item()) / 2

        N = comm_sign_norm / avg_cycle_sign
        Ns.append(N)

    avg_N = sum(Ns) / len(Ns)
    N_sq = avg_N ** 2
    print(f"\n  Seeds={n_seeds:4d}: N = {avg_N:.8f}, N² = {N_sq:.8f}")
    print(f"    spread: {min(Ns):.8f} to {max(Ns):.8f}")

    # Test candidates for N²
    candidates = {
        "4/3": 4/3,
        "27/20": 27/20,
        "7/6": (7/6)**2,  # for N, not N²
        "N=7/6": 7/6,
        "1+K*": 1 + K_STAR,
        "1/(1-K*)": 1/(1-K_STAR),
        "30/23": 30/23,
        "H/(H-1)": H/(H-1),
        "(H²+1)/H²": (H**2+1)/H**2,
        "exp(2Δ)-1": math.exp(2*DELTA)-1,
        "1+2Δ": 1 + 2*DELTA,
    }

    if n_seeds == 2000:
        print(f"\n  --- Testing N² = {N_sq:.8f} against candidates ---")
        for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]-N_sq)):
            diff = N_sq - val
            print(f"    {name:>15s} = {val:.8f}  diff = {diff:+.6f}")

        print(f"\n  --- Testing N = {avg_N:.8f} against candidates ---")
        n_candidates = {
            "√(4/3)": math.sqrt(4/3),
            "7/6": 7/6,
            "√(27/20)": math.sqrt(27/20),
            "√(30/23)": math.sqrt(30/23),
            "1+Δ/3": 1 + DELTA/3,
            "1+K*/2": 1 + K_STAR/2,
            "(H²+H+1)/(H²+1)": (H**2+H+1)/(H**2+1),
            "13/11": 13/11,
            "√(H/(H-1))": math.sqrt(H/(H-1)),
            "exp(Δ/2)": math.exp(DELTA/2),
            "1/(1-Δ/2)": 1/(1-DELTA/2),
            "√((H²+H+1)/H²)": math.sqrt((H**2+H+1)/H**2),
        }
        for name, val in sorted(n_candidates.items(), key=lambda x: abs(x[1]-avg_N)):
            diff = avg_N - val
            print(f"    {name:>25s} = {val:.8f}  diff = {diff:+.8f}")


# ═══════════════════════════════════════════════════════════════
#  Does N depend on entangler parameters?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("PARAMETER DEPENDENCE OF N")
print("="*80)

n_seeds = 500

for n_iters in [20, 50, 100]:
    for discount in [0.1, 0.3, 0.5]:
        if n_iters == 100 and discount == 0.5:
            continue  # known to blow up

        crystals_param = {}
        for name, corr in list(corrs.items())[:3]:  # just transpositions
            avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
            for seed in range(n_seeds):
                avg += Entangler(corr, seed=seed).build(n_steps=n_iters, discount=discount).joint
            avg /= n_seeds
            crystals_param[name] = avg

        # 3-cycles
        for name in ["(012)", "(021)"]:
            avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
            for seed in range(n_seeds):
                avg += Entangler(corrs[name], seed=seed).build(n_steps=n_iters, discount=discount).joint
            avg /= n_seeds
            crystals_param[name] = avg

        comm = compose(crystals_param["(01)"], crystals_param["(02)"]) - \
               compose(crystals_param["(02)"], crystals_param["(01)"])
        _, sign_comm, _ = s3_project(comm)

        _, sign_012, _ = s3_project(crystals_param["(012)"])
        _, sign_021, _ = s3_project(crystals_param["(021)"])
        avg_csn = (sign_012.abs().pow(2).sum().sqrt().item() +
                   sign_021.abs().pow(2).sum().sqrt().item()) / 2

        N = sign_comm.abs().pow(2).sum().sqrt().item() / avg_csn if avg_csn > 1e-10 else 0
        print(f"  iters={n_iters:3d}, discount={discount:.1f}: N = {N:.6f}")


# ═══════════════════════════════════════════════════════════════
#  Alternative: N from raw (un-normalized) composition
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("N FROM RAW (UN-NORMALIZED) COMPOSITION")
print("="*80)

# Build with default params, high seed count
n_seeds = 1000
crystals_raw = {}
for name, corr in corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    crystals_raw[name] = avg

# Raw commutator
A, B = crystals_raw["(01)"], crystals_raw["(02)"]
raw_AB = torch.einsum("ab,bc->ac", A, B)
raw_BA = torch.einsum("ab,bc->ac", B, A)
raw_comm = raw_AB - raw_BA

_, sign_raw_comm, _ = s3_project(raw_comm)
_, sign_012, _ = s3_project(crystals_raw["(012)"])
_, sign_021, _ = s3_project(crystals_raw["(021)"])

avg_csn = (sign_012.abs().pow(2).sum().sqrt().item() +
           sign_021.abs().pow(2).sum().sqrt().item()) / 2
N_raw = sign_raw_comm.abs().pow(2).sum().sqrt().item() / avg_csn

print(f"\n  N (raw composition) = {N_raw:.8f}")
print(f"  N (normalized)      = {Ns[-1] if Ns else 0:.8f}")
print(f"  Ratio N_norm/N_raw  = {Ns[-1]/N_raw if N_raw > 1e-10 else 0:.8f}")


# ═══════════════════════════════════════════════════════════════
#  The L1 factor
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("L1 NORMALIZATION FACTOR")
print("="*80)

L1_AB = raw_AB.reshape(-1).real.sum().item()
L1_BA = raw_BA.reshape(-1).real.sum().item()

print(f"  L1(A·B) = {L1_AB:.8f}")
print(f"  L1(B·A) = {L1_BA:.8f}")
print(f"  L1_AB / L1_BA = {L1_AB/L1_BA:.8f}")

# The normalized commutator = AB/L1_AB - BA/L1_BA
# The raw commutator = AB - BA
# So N_norm/N_raw should relate to the L1 ratio

# Actually: compose(A,B) = raw(A,B)/L1(A,B)
# So [A,B]_norm = raw(AB)/L1_AB - raw(BA)/L1_BA
# And [A,B]_raw = raw(AB) - raw(BA)

# If L1_AB ≈ L1_BA ≈ L, then [A,B]_norm ≈ [A,B]_raw / L
# So N_norm / N_raw ≈ 1/L × (correction from L1 difference)

print(f"\n  1/L1_AB = {1/L1_AB:.8f}")
print(f"  Expected N_norm = N_raw / L1 = {N_raw/L1_AB:.8f}")
print(f"  Actual N_norm               = {Ns[-1] if Ns else 0:.8f}")


print(f"\n\n{'='*80}")
print("WHAT THE STRUCTURE CONSTANT REVEALS")
print("="*80)
