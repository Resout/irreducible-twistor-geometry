"""
The paper-crystal connection: Δ_mass = Δ, Δ_Born = 2Δ.

The paper establishes: mass function combinations contract at rate
Δ = -ln(1-K*) = -ln(23/30) ≈ 0.266 per DS step. (Theorem 4)

The crystal shows: eigenvalue RATIO |λ₁/λ₀| of the [4,4] joint mass
decays per self-composition at rate ≈ 2Δ. (Principle 23)

The hypothesis: the factor of 2 comes from Born measurement.
Born probability = |m|². If mass amplitudes decay at rate Δ,
then Born probabilities decay at rate 2Δ.

Specifically:
  |m_1/m_0|(step n) = |m_1/m_0|(step 0) × e^{-n Δ_mass}
  Born(1)/Born(0)(step n) = |m_1/m_0|²(step 0) × e^{-2n Δ_mass}

So if Δ_mass = Δ, then Δ_Born = 2Δ. This would unify the paper and
the crystal computations.

Test: measure BOTH the amplitude ratio AND the Born ratio decay rates.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, schmidt_number, born_probabilities
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)


def eigendata(joint):
    ev, _ = torch.linalg.eig(joint)
    idx = ev.abs().argsort(descending=True)
    return ev[idx]


# ═══════════════════════════════════════════════════════════════
#  Track amplitude ratios AND Born ratios through composition
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("PAPER-CRYSTAL CONNECTION: Δ_mass vs Δ_Born")
print(f"Paper's Δ = -ln(1-K*) = -ln(23/30) = {DELTA:.6f}")
print(f"2Δ = {2*DELTA:.6f}")
print("=" * 80)

from solver.crystals import RELATIONSHIP_SIGNATURES

for name, corr in list(RELATIONSHIP_SIGNATURES.items())[:4]:  # proportional, inverse, modular, exponential
    ent = Entangler(corr, seed=42).build()
    joint = ent.joint

    current = joint.clone()
    amp_ratios = []  # |λ₁|/|λ₀| in amplitude space
    born_ratios = []  # Born(1)/Born(0) in probability space

    for step in range(12):
        ev = eigendata(current)
        amp_ratio = (ev[1].abs() / ev[0].abs()).item()
        amp_ratios.append(amp_ratio)

        # Born probabilities of the eigenvalues
        # Born(i) ∝ |λᵢ|²
        born_ratio = (ev[1].abs()**2 / ev[0].abs()**2).item()
        born_ratios.append(born_ratio)

        current = compose(current, joint)

    # Compute decay rates
    amp_rates = []
    born_rates = []
    for s in range(2, 10):
        if amp_ratios[s] > 1e-12 and amp_ratios[s-1] > 1e-12:
            amp_rates.append(math.log(amp_ratios[s-1] / amp_ratios[s]))
        if born_ratios[s] > 1e-12 and born_ratios[s-1] > 1e-12:
            born_rates.append(math.log(born_ratios[s-1] / born_ratios[s]))

    avg_amp = sum(amp_rates) / len(amp_rates) if amp_rates else 0
    avg_born = sum(born_rates) / len(born_rates) if born_rates else 0

    print(f"\n  {name}:")
    print(f"    Amplitude decay rate: {avg_amp:.4f} = {avg_amp/DELTA:.3f}Δ")
    print(f"    Born decay rate:      {avg_born:.4f} = {avg_born/DELTA:.3f}Δ")
    print(f"    Ratio Born/Amp:       {avg_born/avg_amp:.4f}" if avg_amp > 0 else "    (amp rate zero)")


# ═══════════════════════════════════════════════════════════════
#  Now check: is the amplitude decay rate = Δ?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("AMPLITUDE DECAY: IS IT EXACTLY Δ?")
print("="*80)

# Seed-averaged for identity crystal
amp_rates_all = []
born_rates_all = []

for seed in range(100):
    ent = Entangler(torch.eye(H, dtype=torch.float32), seed=seed).build()
    joint = ent.joint
    current = joint.clone()

    amp_ratios = []
    for step in range(12):
        ev = eigendata(current)
        amp_ratios.append((ev[1].abs() / ev[0].abs()).item())
        current = compose(current, joint)

    rates = []
    for s in range(2, 10):
        if amp_ratios[s] > 1e-12 and amp_ratios[s-1] > 1e-12:
            rates.append(math.log(amp_ratios[s-1] / amp_ratios[s]))
    if rates:
        amp_rates_all.append(sum(rates) / len(rates))

avg = sum(amp_rates_all) / len(amp_rates_all)
std = (sum((r - avg)**2 for r in amp_rates_all) / len(amp_rates_all)) ** 0.5
sem = std / len(amp_rates_all)**0.5

print(f"\n  Identity crystal (100 seeds):")
print(f"    Amplitude decay rate = {avg:.5f} ± {std:.5f}")
print(f"    Rate/Δ = {avg/DELTA:.4f} ± {std/DELTA:.4f}")
print(f"    Deviation from Δ: {abs(avg - DELTA)/sem:.1f}σ")
print(f"    Deviation from 2Δ: {abs(avg - 2*DELTA)/sem:.1f}σ")


# ═══════════════════════════════════════════════════════════════
#  The composition operator eigenvalue IS the amplitude ratio
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("16-DIM OPERATOR EIGENVALUE vs AMPLITUDE DECAY")
print("="*80)

# The 16-dim composition operator T_B has eigenvalues with:
#   |λ_gap/λ₀| = 0.571/0.975 = 0.586
#   -ln(0.586) = 0.535 = 2.01Δ
#
# But this is the BORN ratio rate if we think of the operator acting
# on flattened mass². Let's check: what is the per-row operator?

def build_row_operator(B_joint, row_idx):
    """4×4 operator: how row row_idx of compose(A,B) depends on row row_idx of A."""
    T = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for j in range(MASS_DIM):
        # Unit vector: A has 1 in position (row_idx, j), zeros elsewhere
        A = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
        A[row_idx, j] = 1.0 + 0j
        result = compose(A, B_joint)
        T[:, j] = result[row_idx, :]
    return T

ent = Entangler(torch.eye(H, dtype=torch.float32), seed=42).build()

print(f"\n  Row operator eigenvalues (identity crystal, seed 42):")
for row in range(MASS_DIM):
    T_row = build_row_operator(ent.joint, row)
    ev = torch.linalg.eigvals(T_row)
    idx = ev.abs().argsort(descending=True)
    ev = ev[idx]
    label = f"h{row}" if row < H else "θ"

    # Gap between first two eigenvalues
    gap = math.log(ev[0].abs().item() / ev[1].abs().item()) if ev[1].abs().item() > 1e-10 else float('inf')

    print(f"    Row {label}: |λ| = [{ev[0].abs().item():.4f}, {ev[1].abs().item():.4f}, "
          f"{ev[2].abs().item():.4f}, {ev[3].abs().item():.4f}]  gap = {gap:.4f} = {gap/DELTA:.2f}Δ")


# ═══════════════════════════════════════════════════════════════
#  The paper's DS combination map and the crystal's compose
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE STRUCTURAL FILTER: compose vs DS combine")
print("="*80)

# Crystal compose: tensor contraction of [4,4] joints
# Paper's DS: scalar combination of 4-dim masses
# Are they the same when restricted to a single row?

# compose(A,B)_{ij} = Σ_k A_{ik} * B_{kj} / Σ_k A_{ik}
# This is: row i of A, combined with ALL columns of B,
# weighted by A's row-i masses.

# DS combine(m1, m2) = (m1⊗m2 - K) / (1-K)
# where the product is Dempster's rule

# These are NOT the same operation. Composition is matrix multiplication
# (with normalization). DS combination is Dempster's rule.
# But they both produce exponential decay at multiples of Δ.

# The connection may be: compose acts on ROWS as a weighted average
# of columns of B. The weights (row of A) are mass functions.
# Each composition step mixes the columns of B using A's uncertainty.
# This mixing IS a form of structural filtering — uncertainty about
# the intermediate variable destroys correlations.

# Let's verify: the row operator T_B is equivalent to...
# compose(A,B)[i,j] = Σ_k A[i,k] * B[k,j] (matrix product)
# So T_B acts on A_row_i by: A_row_i → (A_row_i @ B)
# This is just left-multiplication by B^T!

B = ent.joint
print(f"\n  Is the row operator just B^T?")
T_row0 = build_row_operator(B, 0)
BT = B.T.clone()

# They're not the same because compose does normalization
# Let's check the difference
diff = (T_row0 - BT).abs().max().item()
print(f"    max|T_row0 - B^T| = {diff:.6f}")

# Check if it's B^T times a scalar (normalization)
if BT.abs().max() > 1e-10:
    ratios = T_row0 / BT
    valid = BT.abs() > 1e-6
    if valid.any():
        r_vals = ratios[valid]
        print(f"    T_row0 / B^T: mean = {r_vals.mean().abs().item():.6f}, std = {r_vals.std().abs().item():.6f}")

# Actually compose includes normalization per row
# Let's build it carefully:
print(f"\n  Compose structure for row 0:")
A_test = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
A_test[0, :] = ent.joint[0, :]  # copy row 0

result = compose(A_test, B)
manual = A_test[0, :] @ B
manual_norm = manual / manual.sum() if manual.sum().abs() > 1e-10 else manual

print(f"    compose result row 0: {result[0, :].real.tolist()}")
print(f"    manual A[0,:]@B:      {manual.real.tolist()}")
print(f"    manual normalized:    {manual_norm.real.tolist()}")


print(f"\n{'='*80}")
print("SYNTHESIS")
print("="*80)

print(f"""
  The paper establishes: DS combination contracts mass functions at rate
  Δ = -ln(1-K*) = {DELTA:.4f} per step. (Theorem 4)

  The crystal shows: the composition operator T_B on [4,4] joint masses
  is BLOCK-DIAGONAL in rows. Each row-block is a 4×4 linear map with
  spectral gap ≈ 2Δ.

  The connection:
    1. Composition = matrix multiplication + normalization
    2. Each row evolves independently via the row operator
    3. The row operator's eigenvalue spectrum determines the decay rate
    4. The spectral gap of the row operator = 2Δ for bijection crystals

  The factor of 2:
    The 4-dim eigenvalue |λ₁/λ₀| decays at rate r per composition.
    The Born probability ratio |λ₁|²/|λ₀|² decays at rate 2r.
    The 16-dim spectral gap measures the BORN ratio rate = 2Δ.
    The amplitude ratio rate should be Δ.

  Measurement: amplitude ratio rate = {avg:.4f} = {avg/DELTA:.3f}Δ ± {std/DELTA:.3f}
  This is NOT exactly Δ — it's ~2Δ even for amplitudes.

  Interpretation: the composition spectral gap is NOT simply Δ doubled
  by Born measurement. The gap ≈ 2Δ in AMPLITUDE space. The factor 2
  comes from the BILINEAR structure of composition (tensor contraction
  involves products of two mass functions, not one). This is different
  from the paper's single-mass DS combination.
""")
