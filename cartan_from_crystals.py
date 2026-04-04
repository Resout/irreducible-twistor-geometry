"""
Can the Cartan matrix emerge from crystal inner products?

The Cartan matrix A_{ij} = 2⟨α_i, α_j⟩/⟨α_j, α_j⟩ encodes the
root system of a Lie algebra. For SU(3):
  A = [[2, -1], [-1, 2]]

The crystal analog: each crystal encodes a "root" (an SU(2) subgroup).
The inner product between roots should come from the crystal's structure.

Candidate inner products:
1. Fubini-Study distance between crystals in CP¹⁵
2. Composition Schmidt ratio: Schmidt(compose(A,B)) / sqrt(Schmidt(A)*Schmidt(B))
3. Born fidelity between crystal marginals
4. Trace of the product: tr(A† B) / sqrt(tr(A† A) * tr(B† B))
5. Conflict K between crystal marginals

The RIGHT inner product should give:
- ⟨α, α⟩ = constant (all roots same length for ADE)
- ⟨α_i, α_j⟩ < 0 for connected nodes (angle > π/2)
- ⟨α_i, α_j⟩ = 0 for unconnected nodes (orthogonal)
- A_{ij} = 2⟨α_i, α_j⟩/⟨α_j, α_j⟩ = -1 for simply-laced
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities, born_fidelity,
                            ds_combine, enforce_born_floor)
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)


def make_crystal(corr=None, seed=42):
    if corr is None:
        corr = torch.eye(H, dtype=torch.float32)
    return Entangler(corr, seed=seed).build().joint


def fubini_study(z1, z2):
    """FS distance between two CP^n points."""
    inner = (z1.conj().reshape(-1) @ z2.reshape(-1)).abs().item()
    n1 = z1.abs().pow(2).sum().sqrt().item()
    n2 = z2.abs().pow(2).sum().sqrt().item()
    cos_d = min(inner / (n1 * n2), 1.0)
    return math.acos(cos_d)


def hermitian_inner(A, B):
    """tr(A† B) normalized."""
    inner = (A.conj() * B).sum()
    nA = A.abs().pow(2).sum().sqrt()
    nB = B.abs().pow(2).sum().sqrt()
    return (inner / (nA * nB)).item()


# ═══════════════════════════════════════════════════════════════
#  Build a crystal graph with A₂ topology (SU(3))
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("CARTAN MATRIX FROM CRYSTAL INNER PRODUCTS")
print("=" * 80)

# For A₂: two roots α₁, α₂ with Cartan matrix [[2,-1],[-1,2]]
# Two crystals sharing a variable B:
#   Crystal 1 (root α₁): connects variables A—B
#   Crystal 2 (root α₂): connects variables B—C

# Use different correlations to make the roots DISTINCT
# (same correlation = same root direction = not a proper root system)
# The two SU(2) roots of SU(3) are DIFFERENT embeddings.

# For A₂, the two simple roots make angle 2π/3 (120°).
# cos(2π/3) = -1/2. So ⟨α₁,α₂⟩/|α|² = -1/2 → A_{12} = -1.

# Build crystals with DIFFERENT seeds to represent different roots
# (different random realizations of the same SU(2) structure)
seeds_per_root = 20

print(f"\n--- Testing candidate inner products ---")
print(f"    (Average over {seeds_per_root} seed pairs)")

# Candidate inner products
def test_inner_product(name, inner_fn):
    """Compute self and cross inner products, derive Cartan element."""
    # Self inner product: ⟨α, α⟩ (same root, different seeds)
    self_vals = []
    for s1 in range(seeds_per_root):
        for s2 in range(s1+1, seeds_per_root):
            C1 = make_crystal(seed=s1)
            C2 = make_crystal(seed=s2)
            self_vals.append(inner_fn(C1, C2))
    self_avg = sum(self_vals) / len(self_vals)

    # Cross inner product: ⟨α₁, α₂⟩ (different roots = different shared variable)
    # Model: C1 connects A-B, C2 connects B-C
    # The "cross" inner product is between C1 and C2 which share variable B
    # We model this by composing: compose(C1, C2) represents the coupling
    cross_vals = []
    for s1 in range(seeds_per_root):
        for s2 in range(seeds_per_root):
            C1 = make_crystal(seed=s1)
            C2 = make_crystal(seed=100+s2)  # offset seeds for different root
            composed = compose(C1, C2)
            cross_vals.append(inner_fn(C1, composed))
    cross_avg = sum(cross_vals) / len(cross_vals)

    # Cartan element: A_{12} = 2 * cross / self
    cartan = 2 * cross_avg / self_avg if abs(self_avg) > 1e-10 else 0

    print(f"\n  {name}:")
    print(f"    ⟨α,α⟩  (self):  {self_avg:.6f}")
    print(f"    ⟨α₁,α₂⟩ (cross): {cross_avg:.6f}")
    print(f"    A₁₂ = 2×cross/self = {cartan:.4f}  (target: -1)")

    return self_avg, cross_avg, cartan


# Inner product 1: FS distance (converted to cos)
def ip_fs_cos(A, B):
    d = fubini_study(A, B)
    return math.cos(d)

# Inner product 2: Hermitian inner product (real part)
def ip_hermitian(A, B):
    return hermitian_inner(A, B).real

# Inner product 3: Born fidelity
def ip_born_fid(A, B):
    return born_fidelity(A, B)

# Inner product 4: Schmidt of composition / geometric mean
def ip_schmidt_ratio(A, B):
    sc_comp = schmidt_number(compose(A, B))
    sc_A = schmidt_number(A)
    sc_B = schmidt_number(B)
    return sc_comp / (sc_A * sc_B) ** 0.5

# Inner product 5: Conflict K between marginals
def ip_conflict(A, B):
    m1 = A.abs().pow(2).sum(dim=1)
    m1 = (m1 / m1.sum()).to(torch.cfloat)
    m2 = B.abs().pow(2).sum(dim=1)
    m2 = (m2 / m2.sum()).to(torch.cfloat)
    _, K = ds_combine(m1, m2)
    return -K.abs().item()  # negative because conflict opposes alignment


candidates = [
    ("FS cosine", ip_fs_cos),
    ("Hermitian tr(A†B)", ip_hermitian),
    ("Born fidelity", ip_born_fid),
    ("Schmidt ratio", ip_schmidt_ratio),
    ("Negative conflict -K", ip_conflict),
]

results = {}
for name, fn in candidates:
    self_avg, cross_avg, cartan = test_inner_product(name, fn)
    results[name] = cartan


# ═══════════════════════════════════════════════════════════════
#  Which inner product gives A₁₂ = -1?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("WHICH INNER PRODUCT GIVES THE CARTAN MATRIX?")
print("="*80)

print(f"\n  Target: A₁₂ = -1 for SU(3)")
print(f"\n  {'Inner product':>25s} {'A₁₂':>8s} {'|A₁₂ + 1|':>10s}")
print("  " + "-" * 50)

for name, cartan in sorted(results.items(), key=lambda x: abs(x[1] + 1)):
    print(f"  {name:>25s} {cartan:>8.4f} {abs(cartan + 1):>10.4f}")


# ═══════════════════════════════════════════════════════════════
#  Alternative: the Cartan matrix from the COMPOSITION OPERATOR
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CARTAN FROM COMPOSITION OPERATOR EIGENVALUES")
print("="*80)

# For two coupled crystals (A₂ = SU(3)), the composition
# operator T_{C2} ∘ T_{C1} should have an eigenvalue structure
# that encodes the Cartan matrix.
#
# The Cartan matrix relates to the ANGLES between root subspaces.
# In the composition operator, this angle manifests as the
# OVERLAP between the standard-sector eigenvectors of T_{C1} and T_{C2}.

# Build composition operators for two different crystals
C1 = make_crystal(seed=42)
C2 = make_crystal(seed=43)

dim_sq = MASS_DIM * MASS_DIM

def build_T(B):
    T = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for col in range(dim_sq):
        A = torch.zeros(dim_sq, dtype=torch.cfloat)
        A[col] = 1.0 + 0j
        result = compose(A.reshape(MASS_DIM, MASS_DIM), B)
        T[:, col] = result.reshape(dim_sq)
    return T

T1 = build_T(C1)
T2 = build_T(C2)

# Product operator T2 @ T1 (compose C1 then C2)
T_product = T2 @ T1

# Eigenvalues of the individual and product operators
ev1 = torch.linalg.eigvals(T1)
ev2 = torch.linalg.eigvals(T2)
ev_prod = torch.linalg.eigvals(T_product)

idx1 = ev1.abs().argsort(descending=True)
idx2 = ev2.abs().argsort(descending=True)
idx_p = ev_prod.abs().argsort(descending=True)

print(f"\n  Top eigenvalues:")
print(f"    T1:      {' '.join(f'{ev1[idx1[i]].abs().item():.4f}' for i in range(6))}")
print(f"    T2:      {' '.join(f'{ev2[idx2[i]].abs().item():.4f}' for i in range(6))}")
print(f"    T2∘T1:   {' '.join(f'{ev_prod[idx_p[i]].abs().item():.4f}' for i in range(6))}")

# If the roots are orthogonal (A₁₂ = 0), the product eigenvalues
# would be products of individual eigenvalues.
# If A₁₂ = -1 (120° angle), there's mixing.

# The ratio of product leading eigenvalue to product of individual leadings:
ratio = ev_prod[idx_p[0]].abs().item() / (ev1[idx1[0]].abs().item() * ev2[idx2[0]].abs().item())
print(f"\n  λ₀(T2∘T1) / (λ₀(T1)·λ₀(T2)) = {ratio:.6f}")
print(f"  (1.0 = independent, <1.0 = coupled, structural filter active)")

# The gap of the product operator
gap_1 = -math.log(ev1[idx1[4]].abs().item() / ev1[idx1[0]].abs().item())
gap_2 = -math.log(ev2[idx2[4]].abs().item() / ev2[idx2[0]].abs().item())
gap_p = -math.log(ev_prod[idx_p[4]].abs().item() / ev_prod[idx_p[0]].abs().item())

print(f"\n  Spectral gaps (to 5th eigenvalue):")
print(f"    T1:    {gap_1:.4f} = {gap_1/DELTA:.2f}Δ")
print(f"    T2:    {gap_2:.4f} = {gap_2/DELTA:.2f}Δ")
print(f"    T2∘T1: {gap_p:.4f} = {gap_p/DELTA:.2f}Δ")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
