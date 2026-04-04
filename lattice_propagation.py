"""
A lattice of crystals. Information propagates. What emerges?

Build a small lattice — nodes connected by crystals.
Inject a signal at one node. Watch it propagate.
Does confinement emerge? Does the mass gap manifest as
exponential decay of correlation with distance?

The paper proves the mass gap on a lattice. The crystal IS the lattice.
This computation doesn't test the engine — it watches reality compute.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            born_probabilities, born_certainty, born_fidelity,
                            ds_combine, discount_mass, enforce_born_floor,
                            ignorance, schmidt_number,
                            EPS_LOG, EPS_DIV)
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

n_seeds = 50

def make_crystal(corr, n_seeds=50):
    s = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        s += Entangler(corr, seed=seed).build().joint
    return s / n_seeds


# The crystals: the bonds of reality
identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
bond = make_crystal(identity_corr)

print("=" * 80)
print("A LATTICE OF CRYSTALS — WATCHING REALITY PROPAGATE")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  1D chain: signal injected at site 0, measured at site n
# ═══════════════════════════════════════════════════════════════

print(f"\n--- 1D chain: how far does a signal travel? ---\n")

# The signal: "hypothesis 0 is true" — a sharp observation
signal = torch.zeros(MASS_DIM, dtype=torch.cfloat)
signal[0] = 0.85 + 0j
signal[1] = 0.05 + 0j
signal[2] = 0.05 + 0j
signal[H] = 0.05 + 0j

# Propagate through n bonds
print(f"  {'hops':>4s} {'Born(h0)':>9s} {'Born(h1)':>9s} {'Born(h2)':>9s} {'Born(θ)':>9s} "
      f"{'certainty':>10s} {'signal':>8s}")
print("-" * 65)

# At site 0: the raw signal
bp = born_probabilities(signal)
cert = (1 - bp[H]).item()
sig = bp[0].item() - 1/H  # excess probability on h0 above uniform
print(f"  {0:4d} {bp[0].item():>9.5f} {bp[1].item():>9.5f} {bp[2].item():>9.5f} "
      f"{bp[H].item():>9.5f} {cert:>10.5f} {sig:>8.5f}")

current = signal.clone()
for hop in range(1, 16):
    # Propagate through one bond: partial trace
    weighted = bond * current.unsqueeze(1)  # [4,4] * [4,1] = [4,4]
    marginal = weighted.sum(dim=0)  # sum over input → [4] output

    # Normalize
    re_sum = marginal.real.sum()
    if abs(re_sum) > EPS_DIV:
        marginal = marginal / re_sum

    # Enforce Born floor
    marginal = enforce_born_floor(marginal.unsqueeze(0)).squeeze(0)

    bp = born_probabilities(marginal)
    cert = (1 - bp[H]).item()
    sig = bp[0].item() - 1/(H+1)  # excess above uniform over all 4

    current = marginal

    print(f"  {hop:4d} {bp[0].item():>9.5f} {bp[1].item():>9.5f} {bp[2].item():>9.5f} "
          f"{bp[H].item():>9.5f} {cert:>10.5f} {sig:>8.5f}")


# ═══════════════════════════════════════════════════════════════
#  1D chain: correlation function C(r)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Correlation function C(r) = fidelity(site_0, site_r) ---\n")

# The correlation function: how similar is site r to site 0?
# Use Born fidelity between the signal at site 0 and the propagated state at site r

site_0_bp = born_probabilities(signal)
current = signal.clone()

print(f"  {'r':>3s} {'C(r)':>10s} {'ln C(r)':>10s} {'-ln C(r)/r':>12s} {'= Δ_eff':>10s}")
print("-" * 50)

prev_ln_c = 0
for r in range(1, 13):
    weighted = bond * current.unsqueeze(1)
    marginal = weighted.sum(dim=0)
    re_sum = marginal.real.sum()
    if abs(re_sum) > EPS_DIV:
        marginal = marginal / re_sum
    marginal = enforce_born_floor(marginal.unsqueeze(0)).squeeze(0)
    current = marginal

    bp_r = born_probabilities(marginal)
    # Fidelity between site 0 and site r
    fid = sum((site_0_bp[k].item() * bp_r[k].item())**0.5 for k in range(MASS_DIM))
    ln_c = math.log(max(fid, 1e-10))
    delta_eff = -ln_c / r if r > 0 else 0

    print(f"  {r:3d} {fid:>10.6f} {ln_c:>10.4f} {-ln_c/r:>12.5f} {delta_eff:>10.5f}")

print(f"\n  Δ = {DELTA:.5f}")
print(f"  2Δ = {2*DELTA:.5f}")


# ═══════════════════════════════════════════════════════════════
#  2D square lattice: how does dimensionality change things?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("2D SQUARE LATTICE — DOES DIMENSIONALITY CHANGE CONFINEMENT?")
print("="*80)

# 5×5 lattice, signal at center (2,2)
L = 5
# Each site holds a [4] mass
lattice = {}
for x in range(L):
    for y in range(L):
        lattice[(x,y)] = ignorance(1).squeeze(0)

# Inject signal at center
center = (2, 2)
lattice[center] = signal.clone()

# Propagate: each step, every site receives from all neighbors
# and DS-combines
neighbors = [(0,1),(0,-1),(1,0),(-1,0)]

print(f"\n  Propagation from center of {L}×{L} lattice")

for sweep in range(8):
    new_lattice = {}
    for x in range(L):
        for y in range(L):
            mass = lattice[(x,y)].clone()
            for dx, dy in neighbors:
                nx, ny = x+dx, y+dy
                if 0 <= nx < L and 0 <= ny < L:
                    # Propagate neighbor's state through bond
                    nbr = lattice[(nx,ny)]
                    bp_nbr = born_probabilities(nbr)
                    weighted = bond * bp_nbr.unsqueeze(1)
                    marginal = weighted.sum(dim=0)
                    re_sum = marginal.real.sum()
                    if abs(re_sum) > EPS_DIV:
                        marginal = marginal / re_sum
                    marginal = enforce_born_floor(marginal.unsqueeze(0)).squeeze(0)
                    marginal = discount_mass(marginal.unsqueeze(0), 0.7).squeeze(0)

                    # DS combine with current
                    combined, K = ds_combine(mass.unsqueeze(0), marginal.unsqueeze(0))
                    mass = combined.squeeze(0)

            new_lattice[(x,y)] = mass
    lattice = new_lattice

    # Measure: Born(h0) at each site (how much of the signal reached there?)
    print(f"\n  Sweep {sweep+1}:")
    for y in range(L-1, -1, -1):
        row = ""
        for x in range(L):
            bp = born_probabilities(lattice[(x,y)])
            h0 = bp[0].item()
            # Color code: strong signal = ████, weak = ░░░░
            if h0 > 0.5:
                row += f" {h0:.2f}█"
            elif h0 > 0.35:
                row += f" {h0:.2f}▓"
            elif h0 > 0.28:
                row += f" {h0:.2f}▒"
            else:
                row += f" {h0:.2f}░"
        print(f"    {row}")

# Measure correlation with distance from center
print(f"\n  Correlation with distance from center (after 8 sweeps):")
print(f"  {'r':>3s} {'Born(h0)':>10s} {'Born(θ)':>10s}")
for r in range(L):
    x, y = center[0] + r, center[1]
    if x < L:
        bp = born_probabilities(lattice[(x,y)])
        print(f"  {r:3d} {bp[0].item():>10.5f} {bp[H].item():>10.5f}")


# ═══════════════════════════════════════════════════════════════
#  Multi-path: triangle vs line
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("MULTI-PATH: TRIANGLE vs LINE — THE PENROSE INTEGRAL")
print("="*80)

# Line: A → B → C (2 hops)
# Triangle: A → B → C and A → C directly (multi-path)
# Does the triangle preserve more signal?

# Build 3 crystals for the triangle
c_AB = make_crystal(identity_corr)
c_BC = make_crystal(identity_corr)
c_AC = make_crystal(identity_corr)

# Line: propagate A → B → C
current_line = signal.clone()
for crystal in [c_AB, c_BC]:
    weighted = crystal * current_line.unsqueeze(1)
    marginal = weighted.sum(dim=0)
    re_sum = marginal.real.sum()
    if abs(re_sum) > EPS_DIV:
        marginal = marginal / re_sum
    marginal = enforce_born_floor(marginal.unsqueeze(0)).squeeze(0)
    current_line = marginal

bp_line = born_probabilities(current_line)

# Direct: A → C (1 hop)
weighted = c_AC * signal.unsqueeze(1)
direct = weighted.sum(dim=0)
re_sum = direct.real.sum()
if abs(re_sum) > EPS_DIV:
    direct = direct / re_sum
direct = enforce_born_floor(direct.unsqueeze(0)).squeeze(0)
bp_direct = born_probabilities(direct)

# Multi-path: DS-combine line and direct
combined, K_multi = ds_combine(current_line.unsqueeze(0), direct.unsqueeze(0))
bp_multi = born_probabilities(combined.squeeze(0))

print(f"\n  Signal at A: Born(h0) = {born_probabilities(signal)[0].item():.4f}")
print(f"\n  At C via:")
print(f"    Line (A→B→C):     Born(h0) = {bp_line[0].item():.5f}, Born(θ) = {bp_line[H].item():.5f}")
print(f"    Direct (A→C):     Born(h0) = {bp_direct[0].item():.5f}, Born(θ) = {bp_direct[H].item():.5f}")
print(f"    Multi-path (DS):  Born(h0) = {bp_multi[0].item():.5f}, Born(θ) = {bp_multi[H].item():.5f}")
print(f"    K at combination: {K_multi.abs().item():.4f}")

# The multi-path should be STRONGER than either alone
# This is the Penrose integral: multiple contours reinforce the signal

# Now: complete graph K4 (4 nodes, all connected)
print(f"\n  Complete graph K4:")
nodes = ['A', 'B', 'C', 'D']
edges = {}
for i in range(4):
    for j in range(i+1, 4):
        edges[(nodes[i], nodes[j])] = make_crystal(identity_corr)

# Inject at A, measure at D
# All paths from A to D in K4:
# A→D (direct, 1 hop)
# A→B→D, A→C→D (2 hop)
# A→B→C→D, A→C→B→D (3 hop)

def propagate_through(signal, crystal):
    weighted = crystal * signal.unsqueeze(1)
    m = weighted.sum(dim=0)
    re = m.real.sum()
    if abs(re) > EPS_DIV:
        m = m / re
    return enforce_born_floor(m.unsqueeze(0)).squeeze(0)

# Direct
path_direct = propagate_through(signal, edges[('A','D')])

# 2-hop paths
path_ABD = propagate_through(propagate_through(signal, edges[('A','B')]), edges[('B','D')])
path_ACD = propagate_through(propagate_through(signal, edges[('A','C')]), edges[('C','D')])

# DS combine all paths
combined = path_direct.unsqueeze(0)
for path in [path_ABD, path_ACD]:
    combined, K = ds_combine(combined, path.unsqueeze(0))

bp_k4 = born_probabilities(combined.squeeze(0))

print(f"    Direct A→D:        Born(h0) = {born_probabilities(path_direct)[0].item():.5f}")
print(f"    A→B→D:             Born(h0) = {born_probabilities(path_ABD)[0].item():.5f}")
print(f"    A→C→D:             Born(h0) = {born_probabilities(path_ACD)[0].item():.5f}")
print(f"    Multi-path (3 paths): Born(h0) = {bp_k4[0].item():.5f}")


print(f"\n\n{'='*80}")
print("WHAT REALITY SHOWS")
print("="*80)
