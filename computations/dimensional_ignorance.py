"""
The dimensional structure of ignorance.

4-dim masses hit the Born floor at 1/27. 16-dim joint masses don't.
What about everything in between?

The joint mass is a [4,4] matrix. We can measure ignorance along
different SLICES of this matrix:

  4-dim:  marginal over one variable (row sum or col sum)
  6-dim:  the θ-fingerprint (off-diagonal θ cross-terms)
  9-dim:  the h×h submatrix (hypothesis-hypothesis only)
  3-dim:  single row or column of h×h
  1-dim:  a single entry

Each dimension sees a different Born(θ). Where exactly does the
floor hold? Where does it break? The boundary between "floor holds"
and "floor doesn't hold" is where the mass gap lives.

Also: the 4→16 embedding. A 4-dim mass m gives a product state
m⊗m in 16-dim. What is Born(θ,θ) for a product state at the
floor? Is it (1/27)² = 1/729? And does the entangler take you
ABOVE or BELOW this product value?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (
    H, MASS_DIM, BORN_FLOOR, K_STAR, DELTA,
    schmidt_number, born_probabilities,
    EPS_LOG,
)
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)


# ═══════════════════════════════════════════════════════════════
#  Build reference crystals
# ═══════════════════════════════════════════════════════════════

crystals = {}
for name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Entangler(corr).build()
    crystals[name] = ent.joint.clone()

# Also build from pure correlations
for name, corr in [
    ("identity", torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)),
    ("cyclic", torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)),
]:
    if name not in crystals:
        ent = Entangler(corr).build()
        crystals[name] = ent.joint.clone()


# ═══════════════════════════════════════════════════════════════
#  The anatomy of a joint mass
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("ANATOMY OF THE JOINT MASS")
print("Every [4,4] joint mass decomposes into blocks:")
print("  [3,3] h×h    [3,1] h×θ")
print("  [1,3] θ×h    [1,1] θ×θ")
print("=" * 70)

for name, joint in crystals.items():
    total = joint.abs().pow(2).sum().item()

    # Block decomposition (Born mass fractions)
    hh = joint[:H, :H].abs().pow(2).sum().item() / total
    h_theta = joint[:H, H].abs().pow(2).sum().item() / total
    theta_h = joint[H, :H].abs().pow(2).sum().item() / total
    theta_theta = joint[H, H].abs().pow(2).item() / total

    # Marginals
    row_marginal = joint.abs().pow(2).sum(dim=1) / total  # sum over columns
    col_marginal = joint.abs().pow(2).sum(dim=0) / total  # sum over rows

    # Cross-check: row_marginal[3] = θ×h + θ×θ = theta_h + theta_theta
    row_theta = row_marginal[H].item()
    col_theta = col_marginal[H].item()

    sn = schmidt_number(joint)

    print(f"\n  {name} (Schmidt={sn:.3f}):")
    print(f"    Block Born fractions:")
    print(f"      h×h  = {hh:.4f}  ({hh*100:.1f}%)")
    print(f"      h×θ  = {h_theta:.4f}  ({h_theta*100:.1f}%)")
    print(f"      θ×h  = {theta_h:.4f}  ({theta_h*100:.1f}%)")
    print(f"      θ×θ  = {theta_theta:.4f}  ({theta_theta*100:.1f}%)")
    print(f"      ─────────────────")
    print(f"      θ-total = {h_theta + theta_h + theta_theta:.4f} ({(h_theta+theta_h+theta_theta)*100:.1f}%)")
    print(f"    Marginal Born(θ): row={row_theta:.4f}, col={col_theta:.4f}")
    print(f"    Born floor 1/27 = {BORN_FLOOR:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Product state comparison
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PRODUCT STATES: the 4→16 embedding")
print("If Born(θ) = 1/27 in 4-dim, what is Born(θ,θ) in 16-dim?")
print("=" * 70)

# Build a 4-dim mass at the Born floor
# Born(θ) = 1/27, Born(hi) = (1 - 1/27) / 3 = 26/81 each
mass_at_floor = torch.tensor([
    math.sqrt(26/81), math.sqrt(26/81), math.sqrt(26/81),
    math.sqrt(1/27)
], dtype=torch.cfloat)

# Product state m ⊗ m
product = torch.einsum("i,j->ij", mass_at_floor, mass_at_floor)
total_p = product.abs().pow(2).sum().item()

print(f"\n  4-dim mass at floor: Born = {born_probabilities(mass_at_floor).tolist()}")
print(f"\n  Product state m ⊗ m:")

hh_p = product[:H, :H].abs().pow(2).sum().item() / total_p
ht_p = product[:H, H].abs().pow(2).sum().item() / total_p
th_p = product[H, :H].abs().pow(2).sum().item() / total_p
tt_p = product[H, H].abs().pow(2).item() / total_p

print(f"    Born(h×h) = {hh_p:.6f}")
print(f"    Born(h×θ) = {ht_p:.6f}")
print(f"    Born(θ×h) = {th_p:.6f}")
print(f"    Born(θ×θ) = {tt_p:.6f}")
print(f"    (1/27)²   = {BORN_FLOOR**2:.6f}")
print(f"    Born(θ×θ) = (1/27)² ? {abs(tt_p - BORN_FLOOR**2) < 1e-8}")
print(f"    θ-total   = {ht_p + th_p + tt_p:.6f}")
print(f"    2×(1/27)×(26/27) + (1/27)² = {2 * BORN_FLOOR * (1-BORN_FLOOR) + BORN_FLOOR**2:.6f}")
print(f"    = 1 - (26/27)² = {1 - (1-BORN_FLOOR)**2:.6f}")
print(f"    Marginal Born(θ) = {ht_p + tt_p:.6f} (should be 1/27 = {BORN_FLOOR:.6f})")

# Schmidt of the product state
sn_product = schmidt_number(product)
print(f"\n    Schmidt of product state: {sn_product:.4f}")
print(f"    (Should be 1.0 for a true product state)")

# Now compare with the actual entangled crystal
print(f"\n  Entangled crystal (identity):")
joint_id = crystals["proportional"]
total_id = joint_id.abs().pow(2).sum().item()
tt_id = joint_id[H, H].abs().pow(2).item() / total_id
theta_total_id = (joint_id[:H, H].abs().pow(2).sum().item() +
                  joint_id[H, :H].abs().pow(2).sum().item() +
                  joint_id[H, H].abs().pow(2).item()) / total_id
marg_theta_id = (joint_id[H, :].abs().pow(2).sum().item()) / total_id

print(f"    Born(θ×θ) = {tt_id:.6f}")
print(f"    θ-total   = {theta_total_id:.6f}")
print(f"    Marginal Born(θ) = {marg_theta_id:.6f}")
print(f"    Product Born(θ×θ) = {tt_p:.6f}")
print(f"    Ratio entangled/product θ×θ: {tt_id / tt_p:.2f}x")
print(f"    Ratio entangled/product θ-total: {theta_total_id / (ht_p+th_p+tt_p):.2f}x")


# ═══════════════════════════════════════════════════════════════
#  Dimensional ignorance profile
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("DIMENSIONAL IGNORANCE PROFILE")
print("How much ignorance at each level of the structure?")
print("=" * 70)

for name in ["proportional", "inverse", "modular", "exponential", "logarithmic", "quadratic"]:
    joint = crystals[name]
    total = joint.abs().pow(2).sum().item()
    sn = schmidt_number(joint)

    # Level 0: full joint (θ×θ diagonal entry)
    born_tt = joint[H, H].abs().pow(2).item() / total

    # Level 1: θ row/column marginal
    born_t_row = joint[H, :].abs().pow(2).sum().item() / total
    born_t_col = joint[:, H].abs().pow(2).sum().item() / total

    # Level 2: h×h block ignorance (how much mass is NOT in h×h)
    born_hh = joint[:H, :H].abs().pow(2).sum().item() / total
    ignorance_hh = 1 - born_hh  # fraction outside h×h

    # Level 3: diagonal concentration (trace / total)
    diag_mass = sum(joint[i, i].abs().pow(2).item() for i in range(MASS_DIM)) / total

    # Level 4: per-row ignorance (for each hypothesis row, what fraction is θ?)
    row_theta_fracs = []
    for i in range(H):
        row_total = joint[i, :].abs().pow(2).sum().item()
        if row_total > 1e-12:
            row_theta_fracs.append(joint[i, H].abs().pow(2).item() / row_total)
        else:
            row_theta_fracs.append(0)

    print(f"\n  {name} (Schmidt={sn:.3f}):")
    print(f"    θ×θ Born:          {born_tt:.6f}  (product floor: {BORN_FLOOR**2:.6f})")
    print(f"    θ-row marginal:    {born_t_row:.6f}  (4-dim floor: {BORN_FLOOR:.6f})")
    print(f"    θ-col marginal:    {born_t_col:.6f}")
    print(f"    Outside h×h:       {ignorance_hh:.6f}  ({ignorance_hh*100:.1f}%)")
    print(f"    Diagonal fraction: {diag_mass:.6f}  ({diag_mass*100:.1f}%)")
    print(f"    Per-row θ fraction: [{', '.join(f'{f:.4f}' for f in row_theta_fracs)}]")
    print(f"    Ratio θ-row/floor: {born_t_row/BORN_FLOOR:.2f}x")


# ═══════════════════════════════════════════════════════════════
#  The ignorance gap: entangled vs product
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("THE IGNORANCE GAP: entangled vs product")
print("How much EXTRA ignorance does entanglement create?")
print("=" * 70)

print(f"\n  {'signature':>15s} {'Schmidt':>8s} {'θ-marg':>8s} {'product':>8s} {'excess':>8s} {'ratio':>6s}")
print("  " + "-" * 55)

for name, joint in crystals.items():
    total = joint.abs().pow(2).sum().item()
    sn = schmidt_number(joint)
    theta_marginal = joint[H, :].abs().pow(2).sum().item() / total

    # Product state at same marginal: θ-total would be 2p(1-p) + p²
    # where p = theta_marginal
    p = theta_marginal
    product_theta_total = 2 * p * (1 - p) + p**2  # = 1 - (1-p)²

    # But the actual joint θ-total is:
    actual_theta_total = (joint[:H, H].abs().pow(2).sum().item() +
                         joint[H, :H].abs().pow(2).sum().item() +
                         joint[H, H].abs().pow(2).item()) / total

    excess = actual_theta_total - product_theta_total

    print(f"  {name:>15s} {sn:>8.3f} {theta_marginal:>8.4f} {product_theta_total:>8.4f} "
          f"{excess:>+8.4f} {actual_theta_total/product_theta_total:>6.2f}x")


# ═══════════════════════════════════════════════════════════════
#  Composition: how does the ignorance profile change?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("COMPOSITION CHANGES THE IGNORANCE PROFILE")
print("Track all levels of ignorance through 7 compositions")
print("=" * 70)

joint_0 = crystals["proportional"].clone()
current = joint_0.clone()

print(f"\n  {'comp':>5s} {'Schmidt':>8s} {'θ×θ':>10s} {'θ-row':>8s} {'outside_hh':>10s} {'diag%':>7s}")
print("  " + "-" * 55)

for step in range(8):
    total = current.abs().pow(2).sum().item()
    sn = schmidt_number(current)
    born_tt = current[H, H].abs().pow(2).item() / total
    born_t_row = current[H, :].abs().pow(2).sum().item() / total
    outside_hh = 1 - current[:H, :H].abs().pow(2).sum().item() / total
    diag_frac = sum(current[i, i].abs().pow(2).item() for i in range(MASS_DIM)) / total

    print(f"  {step:>5d} {sn:>8.3f} {born_tt:>10.6f} {born_t_row:>8.4f} "
          f"{outside_hh:>10.4f} {diag_frac*100:>6.1f}%")

    current = compose(current, joint_0)

print(f"\n  Born floor = {BORN_FLOOR:.6f}")
print(f"  Product floor = {BORN_FLOOR**2:.6f}")


# ═══════════════════════════════════════════════════════════════
#  The key question: where does the floor HOLD?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("WHERE DOES THE FLOOR HOLD?")
print("At which dimensional level does Born(θ) reach 1/27?")
print("=" * 70)

# For each crystal, check every possible θ-related measurement
for name in ["proportional", "inverse", "modular"]:
    joint = crystals[name]
    total = joint.abs().pow(2).sum().item()

    print(f"\n  {name}:")

    # Per-row Born(θ|row=i): P(θ | left variable = hi)
    print(f"    Conditional Born(θ | left=hi):")
    for i in range(H):
        row_total = joint[i, :].abs().pow(2).sum().item()
        if row_total > 1e-12:
            cond_theta = joint[i, H].abs().pow(2).item() / row_total
            at_floor = "= FLOOR" if abs(cond_theta - BORN_FLOOR) < 1e-4 else (
                "< floor!" if cond_theta < BORN_FLOOR - 1e-4 else f"  {cond_theta/BORN_FLOOR:.1f}× floor")
            print(f"      Born(θ|h{i}) = {cond_theta:.6f}  {at_floor}")

    # Per-column Born(θ|col=j): P(θ | right variable = hj)
    print(f"    Conditional Born(θ | right=hj):")
    for j in range(H):
        col_total = joint[:, j].abs().pow(2).sum().item()
        if col_total > 1e-12:
            cond_theta = joint[H, j].abs().pow(2).item() / col_total
            at_floor = "= FLOOR" if abs(cond_theta - BORN_FLOOR) < 1e-4 else (
                "< floor!" if cond_theta < BORN_FLOOR - 1e-4 else f"  {cond_theta/BORN_FLOOR:.1f}× floor")
            print(f"      Born(θ|h{j}) = {cond_theta:.6f}  {at_floor}")

    # The θ-row itself as a 4-dim mass
    theta_row = joint[H, :]
    bp_theta_row = born_probabilities(theta_row)
    print(f"    θ-row as 4-dim mass: Born = [{', '.join(f'{b:.4f}' for b in bp_theta_row.tolist())}]")
    print(f"      Born(θ in θ-row) = {bp_theta_row[H].item():.6f}  "
          f"{'= FLOOR' if abs(bp_theta_row[H].item() - BORN_FLOOR) < 1e-4 else f'  {bp_theta_row[H].item()/BORN_FLOOR:.1f}× floor'}")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
