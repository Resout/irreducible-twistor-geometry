"""
The normalization drives Born(θ) to 10/33.

Composition: C^{n+1} = C^n @ C / ||C^n @ C||_{L1}

Without normalization: C^n → σ₀^n u₀ v₀†, Born(θ) → 0.264 (SVD limit)
With normalization: Born(θ) → 0.303 = 10/33

The normalization step (dividing by real L1) redistributes mass.
Every time we normalize, we compress the singletons and inflate θ
relative to the SVD limit. Over many steps, this converges to 10/33.

The mechanism: composition produces a raw result with L1 < 1
(because the conflict K removes mass). Normalizing back to L1=1
uniformly scales all entries, but Born(θ) = |θ|²/Σ|m|² is NOT
linearly affected — it depends on the RATIO of θ to singletons.

Actually, normalization scales all entries equally, so ratios
are preserved... unless the Born floor enforcement kicks in.
That's the nonlinear step! enforce_born_floor() scales singletons
DOWN when Born(θ) drops below 1/27, which pushes Born(θ) UP.

But in the compose() function there's NO born floor enforcement —
it's just L1 normalization. So where does the 10/33 come from?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities)
from solver.crystals import Entangler

torch.set_grad_enabled(False)

EPS_DIV = 1e-8

# Build a clean identity crystal
n_seeds = 50
crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for seed in range(n_seeds):
    crystal += Entangler(
        torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
        seed=seed
    ).build().joint
crystal /= n_seeds


print("=" * 80)
print("HOW NORMALIZATION DRIVES Born(θ) TO 10/33")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Track Born(θ) step by step: with vs without normalization
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Step-by-step Born(θ) of C^n ---\n")

def compose_raw(c1, c2):
    """Compose WITHOUT normalization."""
    return torch.einsum("ab,bc->ac", c1, c2)

def compose_norm(c1, c2):
    """Compose WITH L1 normalization (standard)."""
    result = torch.einsum("ab,bc->ac", c1, c2)
    re_sum = result.reshape(-1).real.sum()
    if abs(re_sum) > EPS_DIV:
        result = result / re_sum
    return result

# Track with normalization
cn_norm = crystal.clone()
# Track without normalization (but measure Born from magnitude ratios)
cn_raw = crystal.clone()

print(f"  {'step':>4s} {'Born(θ)_norm':>13s} {'Born(θ)_raw':>12s} {'L1_raw':>10s} "
      f"{'Schmidt_n':>10s} {'Schmidt_r':>10s}")
print("-" * 70)

for step in range(20):
    # Born probabilities
    bp_norm = born_probabilities(cn_norm).reshape(MASS_DIM, MASS_DIM)
    marg_norm = bp_norm.sum(dim=1)

    bp_raw = cn_raw.abs().pow(2).reshape(MASS_DIM, MASS_DIM)
    bp_raw = bp_raw / bp_raw.sum()
    marg_raw = bp_raw.sum(dim=1)

    l1_raw = cn_raw.reshape(-1).real.sum().item()

    sn_n = schmidt_number(cn_norm) if step < 15 else 0  # skip expensive for late steps
    sn_r = schmidt_number(cn_raw / cn_raw.reshape(-1).real.sum()) if step < 15 else 0

    print(f"  {step:4d} {marg_norm[H].item():>13.6f} {marg_raw[H].item():>12.6f} "
          f"{l1_raw:>10.4e} {sn_n:>10.4f} {sn_r:>10.4f}")

    # Next step
    cn_norm = compose_norm(cn_norm, crystal)
    cn_raw = compose_raw(cn_raw, crystal)


# ═══════════════════════════════════════════════════════════════
#  What happens to L1 at each step?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- L1 norm at each step (before normalization) ---\n")

cn = crystal.clone()
l1_vals = []
for step in range(12):
    raw = torch.einsum("ab,bc->ac", cn, crystal)
    l1 = raw.reshape(-1).real.sum().item()
    l1_vals.append(l1)
    cn = raw / l1  # normalize

    born_theta = born_probabilities(cn).reshape(MASS_DIM, MASS_DIM).sum(dim=1)[H].item()
    print(f"  Step {step:2d}: L1_before = {l1:.6f}, Born(θ) = {born_theta:.6f}")

# The L1 values should converge to σ₀ (the dominant singular value)
print(f"\n  Mean L1 (last 5): {sum(l1_vals[-5:])/5:.6f}")
print(f"  σ₀ (from SVD) ≈ 0.264")


# ═══════════════════════════════════════════════════════════════
#  Direct check: is Born(θ) ratio-invariant under normalization?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Is Born(θ) invariant under uniform scaling? ---\n")

# Born(θ) = |m_θ|² / Σ|m_i|² — this IS scale-invariant!
# So normalization (uniform scaling) CANNOT change Born(θ).
# The 0.264 → 0.303 shift must come from somewhere else.

m_test = torch.tensor([0.5+0.1j, 0.3+0.05j, 0.1+0.02j, 0.1+0.03j], dtype=torch.cfloat)
bp1 = born_probabilities(m_test)
bp2 = born_probabilities(m_test * 2.5)
bp3 = born_probabilities(m_test * 0.001)

print(f"  Original:   Born(θ) = {bp1[H].item():.6f}")
print(f"  × 2.5:      Born(θ) = {bp2[H].item():.6f}")
print(f"  × 0.001:    Born(θ) = {bp3[H].item():.6f}")
print(f"  (All identical — Born is scale-invariant)")

# So WHERE does the Born(θ) shift come from?
# It must come from the COMPOSITION ITSELF, not the normalization!
# The einsum("ab,bc->ac") operation mixes entries nonlinearly
# (it's matrix multiplication of complex numbers).

# Let's check: does composition change Born(θ) even without normalization?
print(f"\n  Checking: does raw composition change Born(θ)?")
raw1 = torch.einsum("ab,bc->ac", crystal, crystal)
bp_raw1 = raw1.abs().pow(2).reshape(-1)
bp_raw1 = bp_raw1 / bp_raw1.sum()
born_theta_raw = bp_raw1.reshape(MASS_DIM, MASS_DIM).sum(dim=1)[H].item()

bp_orig = born_probabilities(crystal).reshape(MASS_DIM, MASS_DIM).sum(dim=1)
born_theta_orig = bp_orig[H].item()

print(f"  Born(θ) of C:   {born_theta_orig:.6f}")
print(f"  Born(θ) of C²:  {born_theta_raw:.6f}")
print(f"  Change: {born_theta_raw - born_theta_orig:+.6f}")


# ═══════════════════════════════════════════════════════════════
#  Track Born(θ) of the RAW C^n (no normalization at any step)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Born(θ) of C^n WITHOUT any normalization ---\n")

cn = crystal.clone()
for step in range(12):
    bp = cn.abs().pow(2).reshape(MASS_DIM, MASS_DIM)
    bp = bp / bp.sum()
    marg = bp.sum(dim=1)
    sn = schmidt_number(cn / cn.reshape(-1).real.sum()) if step < 10 else 0
    print(f"  C^{step+1:2d}: Born(θ) = {marg[H].item():.6f}, Schmidt = {sn:.4f}")
    cn = torch.einsum("ab,bc->ac", cn, crystal)  # NO normalization


# ═══════════════════════════════════════════════════════════════
#  The answer: matrix multiplication mixes θ entries
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE MECHANISM: MATRIX MULTIPLICATION AMPLIFIES θ")
print("="*80)

# In C² = Σ_b C_{ab} C_{bc}:
# The (θ,θ) entry gets contributions from ALL intermediate values:
# C²_{θ,θ} = Σ_b C_{θ,b} C_{b,θ}
# This includes C_{θ,h} × C_{h,θ} cross-terms.
# Since θ is always populated (Born floor), these cross-terms
# systematically increase the θ-sector relative to h-h entries.

# Show the θ-row and θ-column of C
print(f"\n  Crystal θ-row (C[θ,:]):")
theta_row = crystal[H, :]
bp_tr = theta_row.abs().pow(2) / theta_row.abs().pow(2).sum()
print(f"    Born = [{', '.join(f'{b.item():.4f}' for b in bp_tr)}]")

print(f"\n  Crystal θ-column (C[:,θ]):")
theta_col = crystal[:, H]
bp_tc = theta_col.abs().pow(2) / theta_col.abs().pow(2).sum()
print(f"    Born = [{', '.join(f'{b.item():.4f}' for b in bp_tc)}]")

# The θ-θ element of C²
c2 = torch.einsum("ab,bc->ac", crystal, crystal)
c2_tt = c2[H, H]
c_tt = crystal[H, H]
print(f"\n  C[θ,θ] = {c_tt.abs().item():.6f}")
print(f"  C²[θ,θ] = {c2_tt.abs().item():.6f}")
print(f"  Ratio |C²[θ,θ]|/|C[θ,θ]| = {c2_tt.abs().item() / c_tt.abs().item():.4f}")

# Compare to a diagonal h-h element
c2_00 = c2[0, 0]
c_00 = crystal[0, 0]
print(f"\n  C[0,0] = {c_00.abs().item():.6f}")
print(f"  C²[0,0] = {c2_00.abs().item():.6f}")
print(f"  Ratio |C²[0,0]|/|C[0,0]| = {c2_00.abs().item() / c_00.abs().item():.4f}")

# The θ-θ entry grows FASTER than the h-h entries under composition!
# This is because the θ-row and θ-column have more "spread" weight.
# More spread → more cross-term contributions → faster growth of θ.

# Show that the relative θ weight increases per step
print(f"\n  Relative θ-weight per composition step:")
cn = crystal.clone()
for step in range(8):
    total_hh = cn[:H, :H].abs().pow(2).sum().item()
    total_ht = cn[:H, H].abs().pow(2).sum().item()
    total_th = cn[H, :H].abs().pow(2).sum().item()
    total_tt = cn[H, H].abs().pow(2).item()
    total = cn.abs().pow(2).sum().item()

    frac_tt = total_tt / total
    frac_theta = (total_ht + total_th + total_tt) / total

    print(f"  C^{step+1}: frac(θ,θ)={frac_tt:.5f}, frac(any θ)={frac_theta:.5f}")
    cn = compose_norm(cn, crystal)


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
