"""
D₄ triality from crystal graph topology.

D₄ has the Dynkin diagram:
        α₂
        |
  α₁ — α₀ — α₃
        |
        α₄

The central node α₀ connects to 3 branches. The outer symmetry
group of this diagram is S₃ — permuting the three legs.

This is the crystal's OWN symmetry group. If we build a crystal
graph with D₄ topology, the crystal framework should see its own
S₃ acting as triality — the crystal examining its own structure.

Questions:
1. Build a 4-node crystal graph with D₄ topology
2. Does the composition operator preserve the S₃ triality?
3. What is the spectral gap structure? (D₄ has 3 degenerate legs)
4. Does the triality mix the three representations (vector, spinor+, spinor-)?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities, born_fidelity,
                            ds_combine, enforce_born_floor)
from solver.crystals import Entangler, compose, slot_measure, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)

dim_sq = MASS_DIM * MASS_DIM

# S₃ elements and their correlation matrices
S3_CORRS = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}


def make_crystal(corr, seed=42):
    return Entangler(corr, seed=seed).build().joint


# ═══════════════════════════════════════════════════════════════
#  Build D₄ crystal graph
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("D₄ TRIALITY FROM CRYSTAL GRAPH TOPOLOGY")
print("=" * 80)

# D₄ has 4 nodes. The central node (α₀) connects to three outer nodes.
# In the crystal framework, each edge is a crystal.
# The three legs are related by S₃ (triality).
#
# For D₄ to emerge, the three outer edges should be the three
# transpositions of S₃: (01), (02), (12).
# The central node mediates their interaction.
#
# Crystal graph:
#   leg₁: (01) crystal connecting node₀—node₁
#   leg₂: (02) crystal connecting node₀—node₂
#   leg₃: (12) crystal connecting node₀—node₃

print(f"\n--- Building D₄ crystal graph ---")
print(f"  Central node: α₀ (identity)")
print(f"  Leg 1: α₁ via (01) transposition")
print(f"  Leg 2: α₂ via (02) transposition")
print(f"  Leg 3: α₃ via (12) transposition")

n_seeds = 50

# Build seed-averaged crystals for each transposition
legs = {}
for trans_name in ["(01)", "(02)", "(12)"]:
    joint_sum = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        joint_sum += make_crystal(S3_CORRS[trans_name], seed=seed)
    legs[trans_name] = joint_sum / n_seeds
    sn = schmidt_number(legs[trans_name])
    print(f"  {trans_name}: Schmidt = {sn:.4f}")

# Also build the identity crystal (central node's self-loop)
id_sum = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for seed in range(n_seeds):
    id_sum += make_crystal(S3_CORRS["e"], seed=seed)
identity = id_sum / n_seeds
print(f"  identity: Schmidt = {schmidt_number(identity):.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 1: Triality symmetry — are the three legs equivalent?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("TEST 1: TRIALITY SYMMETRY — ARE THE THREE LEGS EQUIVALENT?")
print("="*80)

# Born fidelity between each pair of legs
print(f"\n  Born fidelity between legs:")
leg_names = ["(01)", "(02)", "(12)"]
for i in range(len(leg_names)):
    for j in range(i+1, len(leg_names)):
        bf = born_fidelity(legs[leg_names[i]], legs[leg_names[j]])
        print(f"    {leg_names[i]} ↔ {leg_names[j]}: {bf:.6f}")

# Schmidt numbers should be identical
print(f"\n  Schmidt numbers:")
for name in leg_names:
    print(f"    {name}: {schmidt_number(legs[name]):.6f}")

# Slot measures should be identical
print(f"\n  Slot measures:")
for name in leg_names:
    sm = slot_measure(legs[name])
    print(f"    {name}: gen={sm['generator']:.4f} spec={sm['spectral']:.4f} orb={sm['orbit']:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: Composition through the hub — D₄ path composition
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("TEST 2: COMPOSITION THROUGH THE HUB")
print("="*80)

# In D₄, a path from node₁ to node₂ goes through the central node:
#   node₁ —(01)→ α₀ —(02)→ node₂
# This composition is compose(C₁, C₂) where C₁ = (01), C₂ = (02)
#
# Group theory: (02) ∘ (01) = (012) (a 3-cycle!)
# The composition of two transpositions gives a 3-cycle.

print(f"\n  Path compositions (leg_i through hub to leg_j):")
for i in range(len(leg_names)):
    for j in range(len(leg_names)):
        if i == j:
            continue
        comp = compose(legs[leg_names[i]], legs[leg_names[j]])
        sn = schmidt_number(comp)
        rel_type, rel_conf = slot_measure(comp), born_fidelity(comp, legs[leg_names[0]])
        print(f"    {leg_names[i]} ∘ {leg_names[j]}: Schmidt = {sn:.4f}", end="")

        # Compare to 3-cycle crystals
        c012 = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
        c021 = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
        for seed in range(n_seeds):
            c012 += make_crystal(S3_CORRS["(012)"], seed=seed)
            c021 += make_crystal(S3_CORRS["(021)"], seed=seed)
        c012 /= n_seeds
        c021 /= n_seeds

        bf_012 = born_fidelity(comp, c012)
        bf_021 = born_fidelity(comp, c021)
        print(f"  fidelity to (012)={bf_012:.4f}, (021)={bf_021:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 3: The composition algebra of D₄ legs
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("TEST 3: COMPOSITION ALGEBRA — DOES IT CLOSE TO S₃?")
print("="*80)

# If the three legs are the three transpositions,
# their compositions should generate all of S₃:
#   (01)∘(02) = (012) or (021)
#   (02)∘(01) = (021) or (012)
#   (01)∘(12) = ... etc.
#
# The full multiplication table of S₃ should emerge.

# Build reference crystals for all S₃ elements
refs = {}
for name in S3_CORRS:
    ref_sum = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        ref_sum += make_crystal(S3_CORRS[name], seed=seed)
    refs[name] = ref_sum / n_seeds

# For each product of legs, find the closest S₃ element
print(f"\n  Multiplication table (best match by Born fidelity):")
print(f"  {'':>6s}", end="")
for j in leg_names:
    print(f"  {j:>8s}", end="")
print()

for i in leg_names:
    print(f"  {i:>6s}", end="")
    for j in leg_names:
        comp = compose(legs[i], legs[j])
        best_match, best_fid = "?", 0.0
        for ref_name, ref_joint in refs.items():
            fid = born_fidelity(comp, ref_joint)
            if fid > best_fid:
                best_fid = fid
                best_match = ref_name
        print(f"  {best_match:>5s}({best_fid:.2f})", end="")
    print()

# Expected S₃ multiplication table for transpositions:
# (01)∘(01) = e,   (01)∘(02) = (012), (01)∘(12) = (021)
# (02)∘(01) = (021), (02)∘(02) = e,   (02)∘(12) = (012)
# (12)∘(01) = (012), (12)∘(02) = (021), (12)∘(12) = e
print(f"\n  Expected S₃ table:")
print(f"  {'':>6s}   {'(01)':>8s}  {'(02)':>8s}  {'(12)':>8s}")
print(f"  {'(01)':>6s}       e     (012)    (021)")
print(f"  {'(02)':>6s}    (021)       e     (012)")
print(f"  {'(12)':>6s}    (012)    (021)       e")


# ═══════════════════════════════════════════════════════════════
#  Test 4: D₄ root system — 12 roots from 3 simple roots
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("TEST 4: D₄ ROOT STRUCTURE")
print("="*80)

# D₄ has rank 4, with 24 roots (12 positive + 12 negative).
# Simple roots: α₁, α₂, α₃, α₄ with Cartan matrix:
#   A = [[ 2,  0, -1,  0],
#        [ 0,  2, -1,  0],
#        [-1, -1,  2, -1],
#        [ 0,  0, -1,  2]]
#
# But we only have 3 legs (3 transpositions).
# The 4th dimension comes from the CENTRAL NODE.
# The identity crystal IS the 4th root (the central node).
#
# Actually, in our hub-and-spoke topology:
#   α₀ = central node (identity / hub)
#   α₁, α₂, α₃ = three legs (transpositions)
#
# The Cartan matrix should have:
#   ⟨αᵢ, αⱼ⟩ = 0 for i,j both legs (they don't share an edge)
#   ⟨α₀, αᵢ⟩ < 0 for legs (they share an edge with the hub)

# Build S₃ projectors for the standard sector
def perm_mat_16(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return torch.kron(P, P)

S3_PERMS = {
    "e": [0,1,2], "(01)": [1,0,2], "(02)": [2,1,0],
    "(12)": [0,2,1], "(012)": [1,2,0], "(021)": [2,0,1],
}

reps = {n: perm_mat_16(p) for n, p in S3_PERMS.items()}
chi_std = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}

P_std = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
for g in S3_PERMS:
    P_std += chi_std[g] * reps[g]
P_std *= 2 / 6.0

# Project all crystals onto standard sector
def project_std(joint):
    v = joint.reshape(dim_sq)
    v_std = P_std @ v
    if v_std.abs().max() > 1e-10:
        phase = torch.angle(v_std[v_std.abs().argmax()])
        v_std = v_std * torch.exp(-1j * phase)
    norm = v_std.abs().pow(2).sum().sqrt().item()
    if norm > 1e-10:
        v_std = v_std / norm
    return v_std, norm

print(f"\n  Standard-sector projections:")
all_vecs = {}
for name in ["e", "(01)", "(02)", "(12)", "(012)", "(021)"]:
    v, norm = project_std(refs[name])
    all_vecs[name] = v
    print(f"    {name:>6s}: ||v_std|| = {norm:.6f}")

# Inner products between ALL pairs
print(f"\n  Inner product matrix in standard sector:")
all_names = ["e", "(01)", "(02)", "(12)", "(012)", "(021)"]
print(f"  {'':>8s}", end="")
for j in all_names:
    print(f" {j:>7s}", end="")
print()

for i in all_names:
    print(f"  {i:>8s}", end="")
    for j in all_names:
        ip = (all_vecs[i].conj() @ all_vecs[j]).real.item()
        print(f" {ip:>7.3f}", end="")
    print()


# ═══════════════════════════════════════════════════════════════
#  Test 5: DS combination of the three legs (multi-path)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("TEST 5: DS COMBINATION OF THREE LEGS (TRIALITY FUSION)")
print("="*80)

# D₄ triality says the three 8-dim representations (vector, spinor+, spinor-)
# are permuted by the outer automorphism S₃.
#
# In the crystal: if we DS-combine the three leg crystals,
# the result should be S₃-invariant (triality-averaged).
# This is the "multi-path through the hub" computation from Principle 40.

leg_joints = [legs[n] for n in leg_names]

# Flatten to 4D masses for DS combination
leg_masses_4d = []
for j in leg_joints:
    # Marginal: sum over rows → 4D mass
    marginal = born_probabilities(j).sum(dim=0)
    # Convert to complex mass with L1=1
    m = marginal.to(torch.cfloat)
    m = m / m.real.sum()
    leg_masses_4d.append(m)

# DS-combine the three legs
combined = leg_masses_4d[0].unsqueeze(0)
for m in leg_masses_4d[1:]:
    combined, K = ds_combine(combined, m.unsqueeze(0))
    print(f"  After combining: K = {K.item():.6f}")

bp = born_probabilities(combined.squeeze(0))
print(f"\n  Born probabilities of triality fusion:")
print(f"    h₀: {bp[0].item():.6f}")
print(f"    h₁: {bp[1].item():.6f}")
print(f"    h₂: {bp[2].item():.6f}")
print(f"    θ:  {bp[3].item():.6f}")

# Is it S₃-symmetric? (all h_i should be equal)
h_vals = [bp[i].item() for i in range(H)]
h_mean = sum(h_vals) / H
h_spread = max(h_vals) - min(h_vals)
print(f"\n  S₃ symmetry check:")
print(f"    h mean: {h_mean:.6f}")
print(f"    h spread: {h_spread:.6f} (0 = perfectly symmetric)")
print(f"    Born(θ): {bp[H].item():.6f} (floor = {BORN_FLOOR:.6f})")


# ═══════════════════════════════════════════════════════════════
#  Test 6: Composition depth — 2-hop through D₄ hub
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("TEST 6: DEEP COMPOSITION — REPEATED TRAVERSAL OF D₄")
print("="*80)

# Start at leg 1, go through hub, out to leg 2, back through hub, out to leg 3...
# This is iterated composition: (12) ∘ (02) ∘ (01) = ?

chain = compose(legs["(01)"], legs["(02)"])
chain = compose(chain, legs["(12)"])
sn_chain = schmidt_number(chain)
print(f"  (12) ∘ (02) ∘ (01): Schmidt = {sn_chain:.4f}")

# Compare to all S₃ elements
for ref_name, ref_joint in refs.items():
    fid = born_fidelity(chain, ref_joint)
    if fid > 0.8:
        print(f"    matches {ref_name} with fidelity {fid:.4f}")

# Full loop: back through all three
loop = compose(chain, legs["(01)"])
loop = compose(loop, legs["(02)"])
loop = compose(loop, legs["(12)"])
sn_loop = schmidt_number(loop)
print(f"\n  Full D₄ loop (6 compositions): Schmidt = {sn_loop:.4f}")

for ref_name, ref_joint in refs.items():
    fid = born_fidelity(loop, ref_joint)
    if fid > 0.8:
        print(f"    matches {ref_name} with fidelity {fid:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 7: The hub crystal under triality action
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("TEST 7: HOW DOES THE IDENTITY CRYSTAL TRANSFORM UNDER TRIALITY?")
print("="*80)

# Triality acts by permuting the three legs.
# The central node should be FIXED by triality (it's the fixed point).
# Check: compose identity with each transposition.

for name in leg_names:
    comp = compose(identity, legs[name])
    fid_id = born_fidelity(comp, identity)
    fid_trans = born_fidelity(comp, legs[name])
    sn = schmidt_number(comp)
    print(f"  e ∘ {name}: Schmidt = {sn:.4f}, fidelity(e)={fid_id:.4f}, fidelity({name})={fid_trans:.4f}")


# ═══════════════════════════════════════════════════════════════
#  FINDINGS
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
print("""
The D₄ crystal graph with three transposition legs meeting at an
identity hub. Key questions:

1. Are the three legs equivalent under triality? (Born fidelity,
   Schmidt, slot measures should be identical up to seed noise)

2. Does the multiplication table close to S₃? (Composition of
   legs should give 3-cycles, self-composition gives identity)

3. Is the triality fusion S₃-symmetric? (DS combination of three
   legs should give equal h_i probabilities)

4. What is the root structure? (Inner products in standard sector
   reveal whether D₄ root system angles emerge)
""")
