"""
Does the graviton gap close in the multi-site limit?

Single-site graviton sector has internal gap ≈ 0.34Δ (smallest in framework).
Paper says gravitons are massless (base variations preserving K=K*).

Test: build a 2-site system. A "site" is a crystal. Coupling is through
DS propagation (shared variable). The graviton mode is a coherent
orientation change across sites — both crystals rotating together.

If the graviton gap closes from 0.34Δ (1 site) to < 0.34Δ (2 sites),
the massless graviton is emerging.

Approach:
1. Build two identical crystals at two "sites"
2. Form the 2-site state as a tensor product: M₁ ⊗ M₂ (256-dim)
3. Define the 2-site composition operator
4. Measure the graviton sector gap
5. Compare with single-site gap

Actually, 256 dims is too large for direct operator construction.
Simpler approach: use the crystal graph's DS propagation between
two sites connected by a shared variable, and measure how the
graviton mode propagates.

Even simpler: the paper's argument is that base variations are
GLOBAL rotations. At a single site, a rotation is a perturbation
that decays. At multiple sites, a UNIFORM rotation costs no energy
(it's a symmetry). The gap should close because the uniform mode
becomes a zero-cost direction.

Test this with composition: compose the SAME perturbation at
multiple sites (uniform rotation) vs DIFFERENT perturbations
(non-uniform rotation). Uniform should decay slower.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, schmidt_number, ds_combine, born_probabilities
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)

dim_sq = MASS_DIM * MASS_DIM


# ═══════════════════════════════════════════════════════════════
#  Single-site: graviton sector perturbation and decay
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("GRAVITON MODE: SINGLE SITE vs MULTI-SITE")
print("=" * 80)

# Build the base crystal (identity, seed 42)
ent = Entangler(torch.eye(H, dtype=torch.float32), seed=42).build()
B = ent.joint  # the "equilibrium" crystal

# A "graviton perturbation" is a smooth rotation of the hypothesis labels.
# In S₃ terms, it's a small interpolation from identity toward a permutation.
# In the (3,3) sector, it's a symmetric traceless perturbation.

# Build a small perturbation in the graviton sector:
# Take a symmetric traceless matrix and add it to B.
# The simplest: δM_{ij} = ε(δ_{i0}δ_{j0} - δ_{i1}δ_{j1}) for h-indices
# This is symmetric, traceless within h-block, and lives in (3,3).

def graviton_perturbation(epsilon):
    """Small symmetric traceless perturbation in the h-h block."""
    delta = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    delta[0, 0] = epsilon
    delta[1, 1] = -epsilon
    # traceless: 1 + (-1) = 0
    # symmetric: yes
    return delta

def gauge_perturbation(epsilon):
    """Small antisymmetric perturbation in the h-h block."""
    delta = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    delta[0, 1] = epsilon
    delta[1, 0] = -epsilon
    return delta


# ═══════════════════════════════════════════════════════════════
#  Decay of graviton vs gauge perturbations at single site
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Single-site perturbation decay ---")

for pert_name, pert_fn in [("graviton (sym traceless)", graviton_perturbation),
                             ("gauge (antisymmetric)", gauge_perturbation)]:
    epsilon = 0.01
    B_perturbed = B + pert_fn(epsilon)
    # Re-normalize
    re_sum = B_perturbed.reshape(-1).real.sum()
    B_perturbed = B_perturbed / re_sum

    print(f"\n  {pert_name} (ε={epsilon}):")
    print(f"    {'Step':>4s} {'dist from B':>12s} {'decay rate':>12s} {'Schmidt':>8s}")
    print("    " + "-" * 40)

    current = B_perturbed.clone()
    prev_dist = None

    for step in range(12):
        # Distance from equilibrium
        diff = current - B
        dist = diff.abs().pow(2).sum().sqrt().item()

        rate = ""
        if prev_dist is not None and prev_dist > 1e-15 and dist > 1e-15:
            r = math.log(prev_dist / dist)
            rate = f"{r/DELTA:.2f}Δ"

        sc = schmidt_number(current)
        print(f"    {step:>4d} {dist:>12.2e} {rate:>12s} {sc:>8.3f}")

        prev_dist = dist
        current = compose(current, B)


# ═══════════════════════════════════════════════════════════════
#  Multi-site: DS propagation of graviton mode
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("MULTI-SITE: DS PROPAGATION OF PERTURBATIONS")
print("="*80)

# Model: N sites, each with crystal B.
# Site i has state Bᵢ = B + δᵢ.
# "Propagation" = compose adjacent sites.
# UNIFORM perturbation: all δᵢ = δ (graviton mode)
# NON-UNIFORM perturbation: δ₁ = δ, δ₂ = -δ (glueball mode)

# With 2 sites:
# Site 1: B + δ
# Site 2: B + δ (uniform) or B - δ (non-uniform)
# Composed: compose(B₁, B₂)

# After composition, measure the residual perturbation.

epsilon = 0.01
delta = graviton_perturbation(epsilon)

print(f"\n  --- 2 sites, graviton perturbation, ε={epsilon} ---")

# Uniform: both sites get same perturbation
B1_uniform = B + delta
B2_uniform = B + delta
B1_uniform = B1_uniform / B1_uniform.reshape(-1).real.sum()
B2_uniform = B2_uniform / B2_uniform.reshape(-1).real.sum()

# Non-uniform: opposite perturbations
B1_nonunif = B + delta
B2_nonunif = B - delta
B1_nonunif = B1_nonunif / B1_nonunif.reshape(-1).real.sum()
B2_nonunif = B2_nonunif / B2_nonunif.reshape(-1).real.sum()

# Compose
C_uniform = compose(B1_uniform, B2_uniform)
C_nonunif = compose(B1_nonunif, B2_nonunif)
C_unperturbed = compose(B, B)

dist_uniform = (C_uniform - C_unperturbed).abs().pow(2).sum().sqrt().item()
dist_nonunif = (C_nonunif - C_unperturbed).abs().pow(2).sum().sqrt().item()

# Initial perturbation magnitude
dist_initial = delta.abs().pow(2).sum().sqrt().item()

print(f"    Initial perturbation |δ|: {dist_initial:.6f}")
print(f"    After compose (uniform):     |C - C₀| = {dist_uniform:.6f} ({dist_uniform/dist_initial:.4f}× initial)")
print(f"    After compose (non-uniform): |C - C₀| = {dist_nonunif:.6f} ({dist_nonunif/dist_initial:.4f}× initial)")
print(f"    Ratio non-uniform/uniform: {dist_nonunif/dist_uniform:.4f}")

if dist_nonunif > dist_uniform:
    print(f"    ★ Non-uniform LARGER — non-uniform modes are AMPLIFIED by composition")
else:
    print(f"    ★ Uniform LARGER — uniform modes persist MORE through composition")


# ═══════════════════════════════════════════════════════════════
#  Chain of N sites: how does graviton mode scale?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("N-SITE CHAIN: GRAVITON vs GLUEBALL SURVIVAL")
print("="*80)

epsilon = 0.01

for N in [2, 3, 4, 5, 8]:
    # Uniform chain: all sites perturbed the same way
    current_uniform = B + graviton_perturbation(epsilon)
    current_uniform = current_uniform / current_uniform.reshape(-1).real.sum()

    # Non-uniform chain: alternating perturbations
    current_nonunif = B + graviton_perturbation(epsilon)
    current_nonunif = current_nonunif / current_nonunif.reshape(-1).real.sum()

    unperturbed = B.clone()

    for site in range(1, N):
        # Next site
        if True:  # uniform
            next_site = B + graviton_perturbation(epsilon)
        next_site_u = next_site / next_site.reshape(-1).real.sum()
        current_uniform = compose(current_uniform, next_site_u)

        if site % 2 == 1:  # alternating sign
            next_site_n = B - graviton_perturbation(epsilon)
        else:
            next_site_n = B + graviton_perturbation(epsilon)
        next_site_n = next_site_n / next_site_n.reshape(-1).real.sum()
        current_nonunif = compose(current_nonunif, next_site_n)

        unperturbed = compose(unperturbed, B)

    dist_u = (current_uniform - unperturbed).abs().pow(2).sum().sqrt().item()
    dist_n = (current_nonunif - unperturbed).abs().pow(2).sum().sqrt().item()

    print(f"  N={N}: uniform residual = {dist_u:.6f}, non-uniform = {dist_n:.6f}, ratio = {dist_n/dist_u:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Compare: how much does COMPOSITION decay uniform vs non-uniform?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("ITERATED COMPOSITION: UNIFORM vs NON-UNIFORM DECAY RATES")
print("="*80)

# Instead of a chain, do iterated self-composition with the perturbed crystal.
# This isolates the composition dynamics from multi-site topology.

for pert_type, sign2 in [("uniform (both +δ)", 1), ("non-uniform (+δ, -δ)", -1)]:
    B1 = B + graviton_perturbation(epsilon)
    B2 = B + sign2 * graviton_perturbation(epsilon)
    B1 = B1 / B1.reshape(-1).real.sum()
    B2 = B2 / B2.reshape(-1).real.sum()

    # Build the "effective 2-site crystal" via composition
    C = compose(B1, B2)
    C_ref = compose(B, B)

    # Now iterate composition of C with itself
    current = C.clone()
    ref = C_ref.clone()

    print(f"\n  {pert_type}:")
    print(f"    {'Step':>4s} {'|C^n - ref^n|':>14s} {'decay':>8s}")
    print("    " + "-" * 30)

    prev_dist = None
    for step in range(8):
        dist = (current - ref).abs().pow(2).sum().sqrt().item()
        rate = ""
        if prev_dist and prev_dist > 1e-15 and dist > 1e-15:
            rate = f"{math.log(prev_dist/dist)/DELTA:.2f}Δ"
        print(f"    {step:>4d} {dist:>14.2e} {rate:>8s}")
        prev_dist = dist
        current = compose(current, C)
        ref = compose(ref, C_ref)


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
