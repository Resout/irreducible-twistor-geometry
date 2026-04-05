"""
The achievable crystal eigenvalue space.

Since composition = matrix multiplication by B (per row), and the
spectral gap of B determines all dynamics, the question becomes:
what eigenvalues can a crystal have?

A crystal B is a [4,4] cfloat matrix with:
1. Real L1 norm = 1 (re-anchoring)
2. Born(θ) ≥ 1/27 (Born floor on marginals)
3. Built by 50 DS iterations from a [3,3] correlation

The eigenvalues {λ₀, λ₁, λ₂, λ₃} encode:
- |λ₀| = overall mass (close to 1 by L1 normalization)
- |λ₁|/|λ₀| = participation ratio ↔ Schmidt
- arg(λᵢ) = complex phases ↔ relationship type
- rank = effective number of nonzero eigenvalues

Questions:
1. What is the boundary of achievable (|λ₀|, |λ₁/λ₀|, |λ₂/λ₀|)?
2. Does the Born floor constrain eigenvalues?
3. What eigenvalue configurations are forbidden?
4. Is there a spectral characterization of "resonant" crystals?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, schmidt_number
from solver.crystals import Entangler, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)


def eigendata(joint):
    ev, _ = torch.linalg.eig(joint)
    idx = ev.abs().argsort(descending=True)
    return ev[idx]


# ═══════════════════════════════════════════════════════════════
#  Survey: eigenvalues of all crystal types × many seeds
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("CRYSTAL EIGENVALUE SPACE")
print("=" * 80)

# Collect eigenvalue data across types and seeds
all_data = []

# All 6 canonical types × 50 seeds
for name, corr in RELATIONSHIP_SIGNATURES.items():
    for seed in range(50):
        ent = Entangler(corr, seed=seed).build()
        ev = eigendata(ent.joint)
        sc = schmidt_number(ent.joint)

        all_data.append({
            "name": name,
            "seed": seed,
            "schmidt": sc,
            "lambda0": ev[0].abs().item(),
            "lambda1": ev[1].abs().item(),
            "lambda2": ev[2].abs().item(),
            "lambda3": ev[3].abs().item(),
            "ratio1": (ev[1].abs() / ev[0].abs()).item(),
            "ratio2": (ev[2].abs() / ev[0].abs()).item(),
            "ratio3": (ev[3].abs() / ev[0].abs()).item(),
            "phase0": math.atan2(ev[0].imag.item(), ev[0].real.item()) / math.pi,
            "phase1": math.atan2(ev[1].imag.item(), ev[1].real.item()) / math.pi,
            "phase2": math.atan2(ev[2].imag.item(), ev[2].real.item()) / math.pi,
            "phase3": math.atan2(ev[3].imag.item(), ev[3].real.item()) / math.pi,
        })

# Summary by type
print(f"\n  {'Type':>15s} {'|λ₀|':>8s} {'|λ₁/λ₀|':>10s} {'|λ₂/λ₀|':>10s} {'|λ₃/λ₀|':>10s} {'Schmidt':>8s}")
print("  " + "-" * 70)

for name in RELATIONSHIP_SIGNATURES:
    subset = [d for d in all_data if d["name"] == name]
    l0 = sum(d["lambda0"] for d in subset) / len(subset)
    r1 = sum(d["ratio1"] for d in subset) / len(subset)
    r2 = sum(d["ratio2"] for d in subset) / len(subset)
    r3 = sum(d["ratio3"] for d in subset) / len(subset)
    sc = sum(d["schmidt"] for d in subset) / len(subset)
    print(f"  {name:>15s} {l0:>8.4f} {r1:>10.4f} {r2:>10.4f} {r3:>10.4f} {sc:>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  Eigenvalue ratios vs Schmidt number
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("EIGENVALUE RATIOS vs SCHMIDT")
print("="*80)

# Correlation between |λ₁/λ₀| and Schmidt
r1_vals = [d["ratio1"] for d in all_data]
sc_vals = [d["schmidt"] for d in all_data]

def corr_coef(xs, ys):
    n = len(xs)
    mx, my = sum(xs)/n, sum(ys)/n
    c = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx = sum((x-mx)**2 for x in xs)**0.5
    sy = sum((y-my)**2 for y in ys)**0.5
    return c/(sx*sy) if sx>1e-10 and sy>1e-10 else 0

print(f"  corr(|λ₁/λ₀|, Schmidt) = {corr_coef(r1_vals, sc_vals):.4f}")
print(f"  corr(|λ₂/λ₀|, Schmidt) = {corr_coef([d['ratio2'] for d in all_data], sc_vals):.4f}")
print(f"  corr(|λ₃/λ₀|, Schmidt) = {corr_coef([d['ratio3'] for d in all_data], sc_vals):.4f}")

# Is Schmidt exactly a function of eigenvalue ratios?
# Schmidt = 1/Σ(Born_i²) = 1/Σ(|λ_i|⁴/Σ|λ_j|⁴)
# For a 4×4 matrix with eigenvalues λ₀...λ₃:
# The JOINT mass has Born probabilities determined by the SVD, not eigenvalues.
# But let's check if there's a simple relationship.

# Also: is there a FORBIDDEN region?
print(f"\n  Eigenvalue ratio bounds:")
print(f"    |λ₁/λ₀|: [{min(r1_vals):.4f}, {max(r1_vals):.4f}]")
print(f"    |λ₂/λ₀|: [{min(d['ratio2'] for d in all_data):.4f}, {max(d['ratio2'] for d in all_data):.4f}]")
print(f"    |λ₃/λ₀|: [{min(d['ratio3'] for d in all_data):.4f}, {max(d['ratio3'] for d in all_data):.4f}]")


# ═══════════════════════════════════════════════════════════════
#  Phase structure: how many distinct phase families?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("PHASE STRUCTURE BY TYPE")
print("="*80)

for name in RELATIONSHIP_SIGNATURES:
    subset = [d for d in all_data if d["name"] == name]
    p0s = [d["phase0"] for d in subset]
    p1s = [d["phase1"] for d in subset]
    p2s = [d["phase2"] for d in subset]
    p3s = [d["phase3"] for d in subset]

    p0_avg = sum(p0s) / len(p0s)
    p1_avg = sum(p1s) / len(p1s)
    p2_avg = sum(p2s) / len(p2s)
    p3_avg = sum(p3s) / len(p3s)

    print(f"  {name:>15s}: φ₀/π={p0_avg:+.3f}, φ₁/π={p1_avg:+.3f}, φ₂/π={p2_avg:+.3f}, φ₃/π={p3_avg:+.3f}")


# ═══════════════════════════════════════════════════════════════
#  The trace constraint
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("TRACE = λ₀ + λ₁ + λ₂ + λ₃")
print("="*80)

# The trace of the joint mass should relate to the marginal Born(θ)
# or to some conservation law

for name in RELATIONSHIP_SIGNATURES:
    subset = [d for d in all_data if d["name"] == name][:5]  # just 5 seeds
    for d in subset:
        ent = Entangler(RELATIONSHIP_SIGNATURES[name], seed=d["seed"]).build()
        tr = ent.joint.diagonal().sum()
        tr_born = ent.joint.abs().pow(2).diagonal().sum().item()
        total_born = ent.joint.abs().pow(2).sum().item()
        print(f"  {name} s{d['seed']}: tr = {tr.real.item():.4f}+{tr.imag.item():.4f}i, "
              f"tr(|m|²)/total(|m|²) = {tr_born/total_born:.4f}")
    break  # just one type for now


# ═══════════════════════════════════════════════════════════════
#  The Born floor constraint on eigenvalues
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("BORN FLOOR AND EIGENVALUE STRUCTURE")
print("="*80)

# A crystal at the Born floor has Born(θ) = 1/27 exactly.
# What does this constrain about the eigenvalues?
# Born(θ) = |m[θ]|² / Σ|m[j]|² for each row.

# Build a crystal at various Born(θ) levels
ent_id = Entangler(torch.eye(H, dtype=torch.float32), seed=42).build()
joint = ent_id.joint

# What are the marginal Born(θ) and eigenvalue ratios?
from solver.algebra import born_probabilities
marginal = joint.abs().pow(2).sum(dim=1)
marginal_bp = marginal / marginal.sum()

ev = eigendata(joint)
print(f"  Identity crystal (seed 42):")
print(f"    Marginal Born(θ) = {marginal_bp[H].item():.6f}")
print(f"    Born floor = {BORN_FLOOR:.6f}")
print(f"    |λ₀| = {ev[0].abs().item():.6f}")
print(f"    |λ₁/λ₀| = {(ev[1].abs()/ev[0].abs()).item():.6f}")
print(f"    |λ₂/λ₀| = {(ev[2].abs()/ev[0].abs()).item():.6f}")
print(f"    |λ₃/λ₀| = {(ev[3].abs()/ev[0].abs()).item():.6f}")

# How do eigenvalue ratios relate to the trace of Born probabilities?
# For the joint mass M, the Born probability matrix is P = |M|²/sum(|M|²)
# Schmidt = 1 / sum(diag(P)²) if we work with the SVD interpretation...
# Actually let me be more precise

# SVD of joint mass
U, S, Vh = torch.linalg.svd(joint)
print(f"\n  Singular values: {S.tolist()}")
print(f"  SV ratios: σ₁/σ₀ = {S[1]/S[0]:.4f}, σ₂/σ₀ = {S[2]/S[0]:.4f}, σ₃/σ₀ = {S[3]/S[0]:.4f}")

# Compare eigenvalue magnitudes to singular values
print(f"  Eigenvalue magnitudes: {[ev[i].abs().item() for i in range(4)]}")
print(f"  Are they the same? max diff = {max(abs(ev[i].abs().item() - S[i].item()) for i in range(4)):.6f}")

# Check: for symmetric matrices eigenvalues = singular values
# But our matrix is complex and not Hermitian, so they differ
sym_check = (joint - joint.T).abs().max().item()
herm_check = (joint - joint.conj().T).abs().max().item()
print(f"\n  Is joint symmetric? max|M-M^T| = {sym_check:.6f}")
print(f"  Is joint Hermitian? max|M-M†| = {herm_check:.6f}")


# ═══════════════════════════════════════════════════════════════
#  The spectral gap as a function of Schmidt
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SPECTRAL GAP vs SCHMIDT: IS THERE A FORMULA?")
print("="*80)

# For each crystal, compute the composition spectral gap
# and plot against Schmidt
gap_data = []

for d in all_data[:100]:  # First 100 (50 seeds × 2 types)
    name = d["name"]
    seed = d["seed"]
    ent = Entangler(RELATIONSHIP_SIGNATURES[name], seed=seed).build()
    ev = eigendata(ent.joint)

    # Composition spectral gap = -ln(|λ₁/λ₀|)
    if ev[1].abs().item() > 1e-12 and ev[0].abs().item() > 1e-12:
        gap = -math.log(ev[1].abs().item() / ev[0].abs().item())
        gap_data.append({
            "name": name,
            "schmidt": d["schmidt"],
            "gap": gap,
            "gap_over_delta": gap / DELTA,
        })

print(f"  corr(gap/Δ, Schmidt) = {corr_coef([d['gap_over_delta'] for d in gap_data], [d['schmidt'] for d in gap_data]):.4f}")

# Group by Schmidt ranges
for sc_lo, sc_hi in [(1.5, 2.0), (2.0, 2.5), (2.5, 3.0), (3.0, 3.5)]:
    subset = [d for d in gap_data if sc_lo <= d["schmidt"] < sc_hi]
    if subset:
        avg_gap = sum(d["gap_over_delta"] for d in subset) / len(subset)
        print(f"  Schmidt [{sc_lo:.1f}, {sc_hi:.1f}): avg gap/Δ = {avg_gap:.3f} (n={len(subset)})")


# Check: is gap = some function of Schmidt?
# For a rank-r matrix with equal singular values: Schmidt = r
# gap = -ln(σ₁/σ₀) and if σ₁ = σ₀ then gap = 0 and Schmidt = rank
# So resonant (high Schmidt) → small gap → slow decay
# This is what we see: bijections have gap ≈ 2Δ, modular has gap ≈ 3.6Δ

# But is gap = f(Schmidt) a FUNCTION or just a correlation?
print(f"\n  Schmidt vs gap/Δ (sample):")
for d in sorted(gap_data, key=lambda x: x["schmidt"])[::-1][:15]:
    print(f"    Schmidt={d['schmidt']:.3f}, gap/Δ={d['gap_over_delta']:.3f} [{d['name']}]")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
