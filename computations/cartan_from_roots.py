"""
Cartan matrix from crystal ROOT DIRECTIONS.

The roots of SU(3) are the transpositions of S₃:
  α₁ = (01): swap h₀ ↔ h₁
  α₂ = (12): swap h₁ ↔ h₂
  α₃ = (02): swap h₀ ↔ h₂  (= α₁ + α₂, the non-simple root)

The Cartan matrix of SU(3): A = [[2,-1],[-1,2]]
This means ⟨α₁,α₂⟩/⟨α₁,α₁⟩ = -1/2 (angle 120°).

The ROOT DIRECTION of a crystal lives in the standard sector of S₃.
Two transposition crystals have standard-sector content that points
in different directions. Their inner product IN THE STANDARD SECTOR
should give the Cartan matrix element.

The key insight: the standard representation of S₃ is 2-dimensional.
The 3 transpositions project onto 3 directions in this 2D space,
separated by 120° — EXACTLY the root system of SU(3) (= A₂).
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, schmidt_number
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

dim_sq = MASS_DIM * MASS_DIM

# S₃ elements
S3_PERMS = {
    "e":     [0, 1, 2],
    "(01)":  [1, 0, 2],
    "(02)":  [2, 1, 0],
    "(12)":  [0, 2, 1],
    "(012)": [1, 2, 0],
    "(021)": [2, 0, 1],
}

S3_CORRS = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

# Build S₃ projectors
def perm_mat_16(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return torch.kron(P, P)

reps = {n: perm_mat_16(p) for n, p in S3_PERMS.items()}
chi_std = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}

P_std = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
for g in S3_PERMS:
    P_std += chi_std[g] * reps[g]
P_std *= 2 / 6.0  # dim_chi * |G|^{-1}

print(f"Standard projector rank: {torch.linalg.matrix_rank(P_std.real.float()).item()}")


# ═══════════════════════════════════════════════════════════════
#  The standard-sector projections of transposition crystals
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("ROOT DIRECTIONS IN THE STANDARD SECTOR")
print("=" * 80)

# The 3 transpositions = the 3 roots of A₂ (SU(3))
roots = ["(01)", "(02)", "(12)"]

# Build crystals and project onto standard sector
root_vectors = {}

for name in roots:
    # Average over seeds for stability
    v_std_sum = torch.zeros(dim_sq, dtype=torch.cfloat)
    for seed in range(50):
        ent = Entangler(S3_CORRS[name], seed=seed).build()
        v = ent.joint.reshape(dim_sq)
        v_std = P_std @ v
        # Normalize phase: align to real axis
        if v_std.abs().max() > 1e-10:
            phase = torch.angle(v_std[v_std.abs().argmax()])
            v_std = v_std * torch.exp(-1j * phase)
        v_std_sum += v_std

    v_std_avg = v_std_sum / 50.0
    norm = v_std_avg.abs().pow(2).sum().sqrt().item()
    if norm > 1e-10:
        v_std_avg = v_std_avg / norm

    root_vectors[name] = v_std_avg
    print(f"\n  {name}: ||v_std|| = {norm:.6f} (before normalization)")


# ═══════════════════════════════════════════════════════════════
#  Inner product matrix between root vectors
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("INNER PRODUCT MATRIX (ROOT VECTORS IN STANDARD SECTOR)")
print("="*80)

# Compute ⟨α_i, α_j⟩ = Re(v_i† · v_j)
print(f"\n  {'':>6s}", end="")
for j in roots:
    print(f"  {j:>8s}", end="")
print()

inner_products = {}
for i in roots:
    print(f"  {i:>6s}", end="")
    for j in roots:
        ip = (root_vectors[i].conj() @ root_vectors[j]).real.item()
        inner_products[(i,j)] = ip
        print(f"  {ip:>8.4f}", end="")
    print()


# ═══════════════════════════════════════════════════════════════
#  Derive the Cartan matrix
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CARTAN MATRIX FROM ROOT INNER PRODUCTS")
print("="*80)

# Simple roots of A₂: α₁ = (01), α₂ = (12)
# Cartan: A_{ij} = 2⟨α_i, α_j⟩/⟨α_j, α_j⟩
alpha1 = "(01)"
alpha2 = "(12)"

ip_11 = inner_products[(alpha1, alpha1)]
ip_22 = inner_products[(alpha2, alpha2)]
ip_12 = inner_products[(alpha1, alpha2)]
ip_21 = inner_products[(alpha2, alpha1)]

A_11 = 2 * ip_11 / ip_11  # = 2 by definition
A_12 = 2 * ip_12 / ip_22
A_21 = 2 * ip_21 / ip_11
A_22 = 2 * ip_22 / ip_22  # = 2 by definition

print(f"\n  Simple roots: α₁ = {alpha1}, α₂ = {alpha2}")
print(f"  ⟨α₁,α₁⟩ = {ip_11:.6f}")
print(f"  ⟨α₂,α₂⟩ = {ip_22:.6f}")
print(f"  ⟨α₁,α₂⟩ = {ip_12:.6f}")
print(f"  ⟨α₂,α₁⟩ = {ip_21:.6f}")

print(f"\n  Crystal Cartan matrix:")
print(f"    A = [[{A_11:.4f}, {A_12:.4f}],")
print(f"         [{A_21:.4f}, {A_22:.4f}]]")

print(f"\n  Target (SU(3)):")
print(f"    A = [[ 2, -1],")
print(f"         [-1,  2]]")

# Check the ANGLE between roots
cos_angle = ip_12 / (abs(ip_11) * abs(ip_22)) ** 0.5
angle = math.acos(max(min(cos_angle, 1), -1))
print(f"\n  cos(angle) = {cos_angle:.6f}")
print(f"  angle = {angle:.4f} rad = {angle*180/math.pi:.1f}°")
print(f"  Target: 120° (cos = -0.5)")


# ═══════════════════════════════════════════════════════════════
#  Check: do the 3 transposition roots form a regular triangle?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("DO THE 3 TRANSPOSITIONS FORM A ROOT SYSTEM?")
print("="*80)

# In A₂, the 3 positive roots form angles of 60° with each other.
# But the simple roots (01) and (12) should be at 120°.
# The third root (02) = (01)+(12) in the root lattice.

for i in range(len(roots)):
    for j in range(i+1, len(roots)):
        ip = inner_products[(roots[i], roots[j])]
        self_i = inner_products[(roots[i], roots[i])]
        self_j = inner_products[(roots[j], roots[j])]
        cos_a = ip / (abs(self_i) * abs(self_j)) ** 0.5
        angle = math.acos(max(min(cos_a, 1), -1)) * 180 / math.pi
        print(f"  angle({roots[i]}, {roots[j]}) = {angle:.1f}°  (cos = {cos_a:.4f})")


# ═══════════════════════════════════════════════════════════════
#  The standard representation of S₃ IS the root system of A₂
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE STANDARD REPRESENTATION IS THE ROOT SYSTEM")
print("="*80)

print(f"""
  The standard representation of S₃ is 2-dimensional.
  The 3 transpositions act on this 2D space as REFLECTIONS.
  The reflection axes are separated by 60° (= π/3).
  The ROOTS (perpendicular to the reflection axes) are separated by 120° (= 2π/3).

  This IS the root system of A₂ = SU(3):
    α₁ ⊥ reflection axis of (01)
    α₂ ⊥ reflection axis of (12)
    α₃ = α₁ + α₂ ⊥ reflection axis of (02)

  The Cartan matrix A_{ij} = 2cos(angle between αᵢ and αⱼ):
    A₁₂ = 2cos(120°) = 2×(-1/2) = -1

  So the crystal's S₃ standard representation AUTOMATICALLY
  contains the A₂ root system. The Cartan matrix doesn't need
  to be computed — it's built into the representation theory.

  The crystal measures:
    angle((01), (12)) = {inner_products[('(01)','(12)')]/((inner_products[('(01)','(01)')]*inner_products[('(12)','(12)')])**0.5):.4f} as cos
    vs exact cos(120°) = -0.5000
""")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
