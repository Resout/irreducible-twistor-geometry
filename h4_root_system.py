"""
Beyond A₂: what root system appears at H=4?

At H=3: the crystal algebra realizes A₂ (the root system of SU(3)/S₃).
  - 3 transpositions (simple + positive roots)
  - 2 three-cycles (highest root ± its negative)
  - Gauge field rank = 2 = rank(A₂)

At H=4: S₄ has 24 elements.
  - 6 transpositions: (01)(02)(03)(12)(13)(23)
  - 8 three-cycles: (012)(021)(013)(031)(023)(032)(123)(132)
  - 6 four-cycles: (0123)(0132)(0213)(0231)(0312)(0321)
  - 3 double transpositions: (01)(23), (02)(13), (03)(12)
  - 1 identity

Questions:
  1. What is the root system? A₃ (rank 3)?
  2. Does (H-1)² = H+1 generalize? (At H=4: 9 ≠ 5, so NO)
  3. Does the cyclotomic tower extend?
  4. What is the structure of the gauge field at H=4?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math

torch.set_grad_enabled(False)


# ═══════════════════════════════════════════════════════════════
#  H=4 algebra constants
# ═══════════════════════════════════════════════════════════════

H4 = 4
MASS_DIM_4 = H4 + 1  # 5
BORN_FLOOR_4 = 1.0 / (H4 ** 3)  # 1/64
K_STAR_4 = (H4**2 - H4 + 1) / (H4 * (H4**2 + 1))  # 13/68
DELTA_4 = -math.log(1 - K_STAR_4)

print("=" * 80)
print("H=4 vs H=3: STRUCTURAL COMPARISON")
print("=" * 80)

print(f"\n  {'':>20s}  {'H=3':>12s}  {'H=4':>12s}")
print("-" * 50)
print(f"  {'Mass dim (H+1)':>20s}  {'4':>12s}  {'5':>12s}")
print(f"  {'Joint dim':>20s}  {'16':>12s}  {'25':>12s}")
print(f"  {'Born floor 1/H³':>20s}  {1/27:>12.6f}  {1/64:>12.6f}")
print(f"  {'K* = (H²-H+1)/...':>20s}  {7/30:>12.6f}  {13/68:>12.6f}")
print(f"  {'Δ = -ln(1-K*)':>20s}  {-math.log(1-7/30):>12.6f}  {DELTA_4:>12.6f}")
print(f"  {'(H-1)²':>20s}  {'4':>12s}  {'9':>12s}")
print(f"  {'H+1':>20s}  {'4':>12s}  {'5':>12s}")
print(f"  {'(H-1)²=H+1?':>20s}  {'YES':>12s}  {'NO':>12s}")
print(f"  {'|S_H|':>20s}  {'6':>12s}  {'24':>12s}")
print(f"  {'# transpositions':>20s}  {'3':>12s}  {'6':>12s}")
print(f"  {'# 3-cycles':>20s}  {'2':>12s}  {'8':>12s}")


# ═══════════════════════════════════════════════════════════════
#  S₄ representation theory
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("S₄ REPRESENTATION THEORY")
print("="*80)

# S₄ has 5 irreps: trivial(1), sign(1), standard(3), sign⊗standard(3), V₂(2)
# Check: 1² + 1² + 3² + 3² + 2² = 1 + 1 + 9 + 9 + 4 = 24 ✓

# The 5×5 joint mass under S₄:
# 25 entries, decomposed under S₄ action on {0,1,2,3,4} (but 4=θ is special)
# Actually S₄ acts on the H=4 singletons {0,1,2,3}, θ is invariant

# Decomposition of M₅(C) = C⁵⊗C⁵ under S₄:
# S₄ acts on C⁴ (the singleton indices) via permutation representation
# C⁴ = trivial ⊕ standard (1+3)
# θ component is S₄-invariant (trivial)
# So C⁵ = trivial ⊕ trivial ⊕ standard = 2 × trivial + standard

# C⁵ ⊗ C⁵ decomposition:
# (2T + S) ⊗ (2T + S) = 4T² + 4TS + S²
# T² = trivial, TS = standard, S² = trivial + standard + V₂ + ...
# Wait, need to be more careful.

# Standard representation of S₄ is 3-dimensional
# Sym²(standard) = trivial + V₂ (2-dim)
# ∧²(standard) = standard' (3-dim, = sign ⊗ standard)

print("""
  S₄ irreps: trivial(1), sign(1), standard(3), sign⊗std(3), V₂(2)
  Total: 1+1+9+9+4 = 24 ✓

  C⁵ (mass space) under S₄:
    {h₀,h₁,h₂,h₃} → standard ⊕ trivial (3+1)
    θ → trivial (1)
    So C⁵ = 2×trivial ⊕ standard

  M₅(C) = C⁵ ⊗ C⁵ under S₄:
    Need: (2T ⊕ S) ⊗ (2T ⊕ S)
    = 4T⊗T ⊕ 4T⊗S ⊕ S⊗S
    = 4 trivial ⊕ 4 standard ⊕ (trivial ⊕ V₂ ⊕ standard ⊕ sign⊗std)

  Total multiplicities in M₅:
    trivial:    4 + 1 = 5
    standard:   4 + 1 = 5  (each = 3-dim)
    V₂:         1          (= 2-dim)
    sign⊗std:   1          (= 3-dim)
    sign:        0

  Check: 5×1 + 5×3 + 1×2 + 1×3 + 0×1 = 5+15+2+3 = 25 ✓
""")


# ═══════════════════════════════════════════════════════════════
#  Root system identification
# ═══════════════════════════════════════════════════════════════

print(f"{'='*80}")
print("ROOT SYSTEM AT H=4")
print("="*80)

# The root system of S_n is A_{n-1}
# S₃ → A₂ (rank 2, 6 roots)
# S₄ → A₃ (rank 3, 12 roots)

# A₃ has:
# - Rank 3 (3 simple roots)
# - 12 roots total (6 positive + 6 negative)
# - Dynkin diagram: ●—●—● (linear chain of 3 nodes)
# - Is the root system of SU(4) ≅ SO(6)

# The 6 positive roots of A₃ correspond to S₄ transpositions:
# α₁ = (01), α₂ = (12), α₃ = (23)  [simple roots]
# α₁+α₂ = (02), α₂+α₃ = (13)      [non-simple positive]
# α₁+α₂+α₃ = (03)                  [highest root]

print(f"""
  S₃ → A₂ (rank 2):  ●—●
  S₄ → A₃ (rank 3):  ●—●—●

  A₃ positive roots ↔ S₄ transpositions:
    α₁ = (01)          [simple]
    α₂ = (12)          [simple]
    α₃ = (23)          [simple]
    α₁+α₂ = (02)       [compound]
    α₂+α₃ = (13)       [compound]
    α₁+α₂+α₃ = (03)    [highest]

  Expected gauge rank at H=4: rank(A₃) = 3
  (vs rank(A₂) = 2 at H=3)
""")


# ═══════════════════════════════════════════════════════════════
#  Does the cyclotomic tower extend?
# ═══════════════════════════════════════════════════════════════

print(f"{'='*80}")
print("CYCLOTOMIC TOWER")
print("="*80)

# At H=3: the 3rd cyclotomic polynomial Φ₃(x) = x²+x+1
# Its roots are primitive 3rd roots of unity: ω, ω²
# (H-1)² = H+1 ↔ 2² = 4 ✓ (H=3)
# This is because Φ₃(1) = 3 and the polynomial has degree 2 = H-1

# At H=4: Φ₄(x) = x²+1 (degree 2, NOT H-1=3)
# Roots: i, -i (4th roots of unity, not primitive ones from ω³)
# (H-1)² = 9 ≠ 5 = H+1

# But wait — the RELEVANT polynomial might not be Φ_H.
# In A₂: the exponents are {1, 2}, Coxeter number h=3
# In A₃: the exponents are {1, 2, 3}, Coxeter number h=4

# The Coxeter polynomial for A_{n-1}: Φ_h(x) where h=n (the Coxeter number)
# For A₂: h=3, Φ₃(x) = x²+x+1
# For A₃: h=4, Φ₄(x) = x²+1

print(f"""
  Cyclotomic polynomials:
    Φ₃(x) = x²+x+1    (degree φ(3)=2, roots: ω, ω²)
    Φ₄(x) = x²+1       (degree φ(4)=2, roots: i, -i)
    Φ₅(x) = x⁴+x³+x²+x+1  (degree φ(5)=4)
    Φ₆(x) = x²-x+1     (degree φ(6)=2)

  Coxeter numbers of A_{{n-1}}:
    A₂: h=3, Φ₃ applies
    A₃: h=4, Φ₄ applies
    A₄: h=5, Φ₅ applies

  The self-consistency equation (H-1)²=H+1:
    This is equivalent to H²-3H = 0, i.e., H=3 or H=0.
    Only H=3 is physical.

  BUT the cyclotomic connection IS more general:
    At A_{{n-1}}: the Coxeter number h=n
    The eigenvalues of the Coxeter element are exp(2πi·mⱼ/h)
    where mⱼ are the exponents.

    A₂: m₁=1, m₂=2 → eigenvalues ω, ω² (primitive 3rd roots)
    A₃: m₁=1, m₂=2, m₃=3 → eigenvalues i, -1, -i (4th roots)
    A₄: m₁=1, m₂=2, m₃=3, m₄=4 → 5th roots

  The tower extends through Coxeter elements, not through (H-1)²=H+1.
  The self-consistency at H=3 is SPECIAL — it's why H=3 is the
  unique value where the crystal framework is self-describing.
""")


# ═══════════════════════════════════════════════════════════════
#  Numerical: H=4 entangler (adapted from H=3)
# ═══════════════════════════════════════════════════════════════

print(f"{'='*80}")
print("H=4 ENTANGLER: RANK OF GAUGE FIELD")
print("="*80)

# Build H=4 transposition crystals using the existing Entangler
# But we need to modify it for H=4... The Entangler in crystals.py is H=3 only.
# Let's build a minimal H=4 version here.

from solver.algebra import ds_combine, discount_mass, enforce_born_floor

def entangle_h4(corr_4x4, seed=42, n_iters=50, discount=0.3):
    """Build an H=4 joint mass from a 4×4 correlation matrix."""
    torch.manual_seed(seed)
    H = 4
    MD = H + 1  # 5
    floor = 1.0 / H**3

    joint = torch.zeros(MD, MD, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j  # total ignorance

    for step in range(n_iters):
        # Sample a pair of hypothesis indices from correlation
        probs = corr_4x4.flatten() / corr_4x4.sum()
        idx = torch.multinomial(probs, 1, replacement=True).item()
        i, j = idx // H, idx % H

        # Build evidence masses for this pair
        m1 = torch.zeros(MD, dtype=torch.cfloat)
        m1[i] = 1.0 + 0j
        m1[-1] = 0.0
        # Add noise
        noise = torch.randn(MD) * 0.02
        m1 = m1 + noise.to(torch.cfloat) * 1j
        re_sum = m1.real.sum()
        if abs(re_sum) > 1e-8:
            m1 = m1 / re_sum

        m2 = torch.zeros(MD, dtype=torch.cfloat)
        m2[j] = 1.0 + 0j
        m2[-1] = 0.0
        noise = torch.randn(MD) * 0.02
        m2 = m2 + noise.to(torch.cfloat) * 1j
        re_sum = m2.real.sum()
        if abs(re_sum) > 1e-8:
            m2 = m2 / re_sum

        # Discount
        m1 = discount_mass(m1.unsqueeze(0), discount).squeeze(0)
        m2 = discount_mass(m2.unsqueeze(0), discount).squeeze(0)

        # Outer product contribution
        outer = torch.outer(m1, m2)
        joint = joint + outer * discount

        # Re-normalize
        re_sum = joint.reshape(-1).real.sum()
        if abs(re_sum) > 1e-8:
            joint = joint / re_sum

    return joint

# S₄ transposition correlation matrices
s4_trans = {
    "(01)": [[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]],
    "(02)": [[0,0,1,0],[0,1,0,0],[1,0,0,0],[0,0,0,1]],
    "(03)": [[0,0,0,1],[0,1,0,0],[0,0,1,0],[1,0,0,0]],
    "(12)": [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],
    "(13)": [[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]],
    "(23)": [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]],
}

n_seeds = 200
MD4 = MASS_DIM_4

h4_crystals = {}
for name, corr_list in s4_trans.items():
    corr = torch.tensor(corr_list, dtype=torch.float32)
    avg = torch.zeros(MD4, MD4, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += entangle_h4(corr, seed=seed)
    avg /= n_seeds
    h4_crystals[name] = avg

# Check gauge rank from pairs of simple roots
for pair in [("(01)","(12)"), ("(12)","(23)"), ("(01)","(23)")]:
    A4 = h4_crystals[pair[0]]
    B4 = h4_crystals[pair[1]]

    raw_comm = torch.einsum("ab,bc->ac", A4, B4) - torch.einsum("ab,bc->ac", B4, A4)
    _, S, _ = torch.linalg.svd(raw_comm)
    S_abs = S.abs()
    S_norm = S_abs / S_abs.sum().clamp(min=1e-10)
    eff_rank = 1 / (S_norm**2).sum().item() if (S_norm**2).sum().item() > 1e-10 else 0

    print(f"\n  [{pair[0]}, {pair[1]}] raw commutator:")
    for i in range(MD4):
        print(f"    σ{i} = {S_abs[i].item():.6f} ({S_norm[i].item():.4f})")
    print(f"    Effective rank: {eff_rank:.3f}")


# All six transposition commutators
print(f"\n--- All transposition commutator ranks ---")
trans_names = list(s4_trans.keys())
for i in range(len(trans_names)):
    for j in range(i+1, len(trans_names)):
        na, nb = trans_names[i], trans_names[j]
        A4 = h4_crystals[na]
        B4 = h4_crystals[nb]
        raw_comm = torch.einsum("ab,bc->ac", A4, B4) - torch.einsum("ab,bc->ac", B4, A4)
        _, S, _ = torch.linalg.svd(raw_comm)
        S_abs = S.abs()
        S_norm = S_abs / S_abs.sum().clamp(min=1e-10)
        eff_rank = 1 / (S_norm**2).sum().item() if (S_norm**2).sum().item() > 1e-10 else 0
        print(f"  [{na},{nb}]: rank = {eff_rank:.2f}")


print(f"\n\n{'='*80}")
print("WHAT THE H=4 ROOT SYSTEM REVEALS")
print("="*80)
