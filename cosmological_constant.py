"""
The cosmological constant from the crystal framework.

The paper (line 1206): "Lambda is determined by K*=7/30 and the
Fubini-Study curvature of CP³ (R=24). It is a derived quantity,
not a free parameter. Its explicit value requires the Penrose
transform integral, which has not been evaluated."

The crystal provides:
- K* = 7/30 (exactly)
- R_FS = 24 (half-standard normalization of CP³)
- The Penrose transform IS composition

The Maldacena reduction gives Lambda from the conformal gravity
action. In the conformal gravity framework:
  S_conformal = integral of |Weyl|² over S⁴
  Lambda appears when restricting to the Einstein sector (OS2)

The relationship (standard conformal gravity):
  Lambda ∝ (gap)² / R

Can we compute Lambda from crystal quantities?

What we have:
  Δ_filter = -ln(1-K*) = -ln(23/30) = 0.2657 per DS step
  Δ_gap = 1.263 per DS step (paper's Jacobian eigenvalue)
  R_FS(CP³) = 24
  BORN_FLOOR = 1/27
  K* = 7/30

Candidates for Lambda:
  Lambda ∝ Δ² / R
  Lambda ∝ K*² × R
  Lambda ∝ BORN_FLOOR × R
  Lambda ∝ (Δ_gap)² / R_FS
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
from solver.algebra import H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR

# ═══════════════════════════════════════════════════════════════
#  Constants from the framework
# ═══════════════════════════════════════════════════════════════

R_FS = 2 * H * (H + 1)  # = 24, half-standard normalization
Delta_filter = DELTA  # = -ln(1 - K*) = 0.2657
Delta_gap = 1.2626  # paper's Jacobian eigenvalue at fixed point

print("=" * 80)
print("COSMOLOGICAL CONSTANT FROM CRYSTAL FRAMEWORK")
print("=" * 80)

print(f"\n  Framework constants:")
print(f"    H = {H}")
print(f"    K* = {K_STAR} = 7/30 = {7/30:.10f}")
print(f"    BORN_FLOOR = 1/H³ = 1/{H**3} = {BORN_FLOOR:.10f}")
print(f"    Δ_filter = -ln(1-K*) = {Delta_filter:.10f}")
print(f"    Δ_gap = {Delta_gap:.10f} (paper's Jacobian eigenvalue)")
print(f"    R_FS(CP³) = 2H(H+1) = {R_FS}")


# ═══════════════════════════════════════════════════════════════
#  Candidate expressions for Lambda
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CANDIDATE EXPRESSIONS FOR Λ")
print("="*80)

# In conformal gravity, the Einstein sector has:
# Lambda = 3/l² where l is the de Sitter radius
# The gap provides a mass scale: m = Δ/a (where a is lattice spacing)
# Lambda ∝ m² in natural units → Lambda ∝ Δ²

# But Lambda has units of 1/length², so we need:
# Lambda = C × Δ² × (curvature factor)

candidates = {
    # From filter rate
    "Δ²_filter / R": Delta_filter**2 / R_FS,
    "Δ²_filter × R": Delta_filter**2 * R_FS,
    "3 Δ²_filter": 3 * Delta_filter**2,

    # From gap
    "Δ²_gap / R": Delta_gap**2 / R_FS,
    "Δ²_gap × K*": Delta_gap**2 * K_STAR,
    "3 Δ²_gap": 3 * Delta_gap**2,

    # From K*
    "K* × R": K_STAR * R_FS,
    "K*² × R": K_STAR**2 * R_FS,
    "K* × R / H³": K_STAR * R_FS / H**3,
    "K*² × R × H": K_STAR**2 * R_FS * H,

    # From Born floor
    "BORN × R": BORN_FLOOR * R_FS,
    "BORN² × R": BORN_FLOOR**2 * R_FS,

    # Dimensionless combinations
    "K*/H": K_STAR / H,
    "K*²": K_STAR**2,
    "BORN × K*": BORN_FLOOR * K_STAR,
    "Δ_filter × K*": Delta_filter * K_STAR,
    "Δ_gap × K*": Delta_gap * K_STAR,

    # From the paper's structure
    "7/(30 × 24)": 7 / (30 * 24),   # K*/R
    "7²/(30² × 24)": 49 / (900 * 24),  # K*²/R
    "(H²-H+1)/(H(H²+1)·2H(H+1))": (H**2-H+1) / (H*(H**2+1)*2*H*(H+1)),
}

print(f"\n  {'Expression':>35s} {'Value':>14s}")
print("  " + "-" * 55)

for name, val in sorted(candidates.items(), key=lambda x: x[1]):
    print(f"  {name:>35s} {val:>14.8f}")


# ═══════════════════════════════════════════════════════════════
#  The paper's specific formula
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE PAPER'S FORMULA STRUCTURE")
print("="*80)

print(f"""
  The paper says Lambda comes from the Maldacena reduction of
  conformal gravity to Einstein gravity on S⁴.

  In conformal gravity on S⁴:
    S = integral |W|² d⁴x (W = Weyl tensor)
    The Einstein sector (selected by OS2, no ghosts) gives:
    S_Einstein = (1/16πG) integral (R - 2Λ) √g d⁴x

  The Maldacena argument: conformal gravity on S⁴ with the
  round metric gives R = 12/l² (Ricci scalar of S⁴ radius l).
  The ghost-free sector requires Λ > 0 (de Sitter).

  In the DS framework:
    - The equilibrium mass m* determines a connection on CP³
    - The Weyl tensor W is computed from the (3,3) component of ∂̄Φ
    - The Einstein reduction uses the trace part (1,1) for conformal factor
    - Lambda is determined by the RATIO of (1,1) to (3,3) components

  From our SO(4) sector analysis (Principle 35):
    Identity crystal: (1,1) = 46.8%, (3,3) = 53.2%
    Ratio: (1,1)/(3,3) = {46.8/53.2:.4f}

  This ratio may determine Lambda/R through the conformal-to-Einstein
  reduction. The formula would be:
    Lambda ∝ R_FS × [(1,1) fraction] / [(3,3) fraction]

  With our numbers:
    Lambda ∝ 24 × 0.468/0.532 = {24 * 0.468/0.532:.4f}
""")

# The (1,1)/(3,3) ratio
scalar_frac = 0.468
graviton_frac = 0.532
ratio = scalar_frac / graviton_frac

print(f"  Scalar/Graviton ratio: {ratio:.6f}")
print(f"  R × ratio = {R_FS * ratio:.4f}")
print(f"  R × ratio / (4π) = {R_FS * ratio / (4 * math.pi):.4f}")

# The ratio 0.468/0.532 ≈ 0.879... is it a clean fraction?
# 0.468 + 0.532 = 1.000. So scalar = 0.468, graviton = 0.532.
# Ratio ≈ 7/8 = 0.875? Close but not exact.
# Or (H²-1)/H² = 8/9 = 0.889? Not exact either.

# These fractions are seed-dependent (measured at seed 42).
# Need to average over seeds for a clean result.

print(f"\n  Is the ratio a clean fraction of H?")
for num, den, name in [(7, 8, "7/8"), (8, 9, "8/9"), (13, 15, "13/15"),
                        (26, 27, "26/27"), (23, 30, "23/30"),
                        (H**2-1, H**2, "(H²-1)/H²"),
                        (H**2, H**2+1, "H²/(H²+1)")]:
    frac = num/den
    print(f"    {name:>12s} = {frac:.6f} (diff = {abs(frac - ratio):.4f})")


# ═══════════════════════════════════════════════════════════════
#  Seed-averaged scalar/graviton ratio
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SEED-AVERAGED SCALAR/GRAVITON RATIO")
print("="*80)

import torch
from solver.crystals import Entangler

torch.set_grad_enabled(False)

dim_sq = MASS_DIM * MASS_DIM
swap = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
for i in range(MASS_DIM):
    for j in range(MASS_DIM):
        swap[j * MASS_DIM + i, i * MASS_DIM + j] = 1.0

P_sym = (torch.eye(dim_sq, dtype=torch.cfloat) + swap) / 2
trace_vec = torch.zeros(dim_sq, dtype=torch.cfloat)
for i in range(MASS_DIM):
    trace_vec[i * MASS_DIM + i] = 1.0
trace_vec = trace_vec / trace_vec.abs().pow(2).sum().sqrt()
P_trace = torch.outer(trace_vec, trace_vec.conj())
P_graviton = P_sym - P_trace
P_anti = (torch.eye(dim_sq, dtype=torch.cfloat) - swap) / 2

scalar_fracs = []
graviton_fracs = []
gauge_fracs = []

corr = torch.eye(H, dtype=torch.float32)
for seed in range(200):
    ent = Entangler(corr, seed=seed).build()
    v = ent.joint.reshape(dim_sq)
    total = v.abs().pow(2).sum().item()
    scalar_fracs.append((P_trace @ v).abs().pow(2).sum().item() / total)
    graviton_fracs.append((P_graviton @ v).abs().pow(2).sum().item() / total)
    gauge_fracs.append((P_anti @ v).abs().pow(2).sum().item() / total)

s_avg = sum(scalar_fracs) / len(scalar_fracs)
g_avg = sum(graviton_fracs) / len(graviton_fracs)
a_avg = sum(gauge_fracs) / len(gauge_fracs)
s_std = (sum((x-s_avg)**2 for x in scalar_fracs)/len(scalar_fracs))**0.5
g_std = (sum((x-g_avg)**2 for x in graviton_fracs)/len(graviton_fracs))**0.5

ratio_avg = s_avg / g_avg
ratio_std = ratio_avg * ((s_std/s_avg)**2 + (g_std/g_avg)**2)**0.5

print(f"  Identity crystal (200 seeds):")
print(f"    Scalar (1,1):     {s_avg:.6f} ± {s_std:.6f}")
print(f"    Graviton (3,3):   {g_avg:.6f} ± {g_std:.6f}")
print(f"    Gauge Λ²:         {a_avg:.6f}")
print(f"    Scalar/Graviton:  {ratio_avg:.6f} ± {ratio_std:.6f}")

# Clean fraction check
print(f"\n  Clean fraction candidates for scalar/graviton ratio:")
for num, den, name in [(7, 8, "7/8"), (8, 9, "8/9"), (13, 15, "13/15"),
                        (26, 27, "26/27"), (23, 30, "23/30"),
                        (H**2-1, H**2, "(H²-1)/H²"),
                        (H**2, H**2+1, "H²/(H²+1)"),
                        (H, H+1, "H/(H+1)"),
                        (1, 1, "1")]:
    frac = num/den
    sigma = abs(frac - ratio_avg) / ratio_std if ratio_std > 0 else float('inf')
    marker = " ★" if sigma < 2 else " ◇" if sigma < 5 else ""
    print(f"    {name:>12s} = {frac:.6f} ({sigma:.1f}σ){marker}")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
