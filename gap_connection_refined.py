"""
Refined analysis: the connection between crystal spectral gap and DS mass gap.

Crystal: |λ₁/λ₀| = 1/√H, gap = ln(H)/2 = 0.54931
DS gap: Δ_gap = -ln(0.28291035) = 1.26263

Near-miss: (23/20) × ln(3) = 1.26340 (0.062% off)

What IS the exact ratio?
"""

import math

H = 3
K_STAR = 7/30
DELTA_FILTER = -math.log(1 - K_STAR)  # 0.26570
BORN_FLOOR = 1/H**3  # 1/27

# Exact values
lambda_0_DS = 0.282910346606210  # from spectral_gap_crystal_graph.py
delta_gap = -math.log(lambda_0_DS)  # 1.262625228

crystal_gap = math.log(H) / 2  # ln(3)/2 = 0.54931

print("=" * 80)
print("THE GAP CONNECTION")
print("=" * 80)

# The ratio
ratio = delta_gap / crystal_gap
print(f"\n  Δ_gap = {delta_gap:.12f}")
print(f"  crystal_gap = ln(H)/2 = {crystal_gap:.12f}")
print(f"  ratio = Δ_gap / crystal_gap = {ratio:.12f}")

# Test candidates for the ratio
print(f"\n  Candidates for ratio = {ratio:.10f}:")
candidates = {
    "H(1-K*) = 23/10": H * (1-K_STAR),
    "H-K*": H - K_STAR,
    "(H³-1)/(H(H-1))": (H**3-1)/(H*(H-1)),
    "(H³-1)/(H²+1)×H/(H-1)": (H**3-1)/(H**2+1) * H/(H-1),
    "13/5 × 23/26": (13/5) * (23/26),
    "(H²+H+1)/(H+1)": (H**2+H+1)/(H+1),
    "13/4": 13/4,
    "(H²-1)/H × H/(H-1)": (H**2-1)/H * H/(H-1),
}

for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - ratio)):
    diff = ratio - val
    pct = abs(diff)/ratio * 100
    print(f"    {name:>35s} = {val:.10f} (diff = {diff:+.6e}, {pct:.4f}%)")

# The power: λ₀ = (1-K*)^r
r = delta_gap / DELTA_FILTER
print(f"\n  Power: λ₀ = (1-K*)^r where r = {r:.12f}")

r_candidates = {
    "19/4": 19/4,
    "H² + H - 1/(H+1)": H**2 + H - 1/(H+1),
    "(2H²+1)/(H-1)": (2*H**2+1)/(H-1),
    "(H³+H-1)/(H²-H+1)": (H**3+H-1)/(H**2-H+1),
    "(H²+1)/(H-1)×(H-1)/1": (H**2+1)/(H-1),
    "H²+H-1/4": H**2+H-1/4,
}

for name, val in sorted(r_candidates.items(), key=lambda x: abs(x[1] - r)):
    diff = r - val
    pct = abs(diff)/r * 100
    print(f"    {name:>35s} = {val:.10f} (diff = {diff:+.6e}, {pct:.4f}%)")


# The deep structure: what does the chain rule actually give?
print(f"\n\n{'='*80}")
print("THE CHAIN RULE DECOMPOSITION")
print("="*80)

# λ₀ = dR/ds × ds_new/dR
# where:
#   R = s/w (ratio on state surface)
#   dR/ds = derivative through DS combination on the state surface
#   ds_new/dR = derivative of Born floor projection

# From the computation:
dR_ds = 163.7781877176
ds_new_dR = 0.0017273994

print(f"\n  λ₀ = dR/ds × ds_new/dR")
print(f"     = {dR_ds:.10f} × {ds_new_dR:.10f}")
print(f"     = {dR_ds * ds_new_dR:.10f}")

# What is dR/ds in terms of known quantities?
# R = s/w = 26.873 at the fixed point
R_fp = 26.87287276
print(f"\n  R (fixed point) = s/w = {R_fp:.8f}")
print(f"  R + 2 = {R_fp + 2:.8f}")
print(f"  (R+2) ≈ 28.87 ≈ 30-1.13 ≈ 30 × 26/27")
print(f"  30 × 26/27 = {30*26/27:.8f}")

# The Born floor projection involves c = (R²+2)/(26(R+2)²)
c_fp = (R_fp**2 + 2) / (26 * (R_fp + 2)**2)
print(f"\n  c = (R²+2)/(26(R+2)²) = {c_fp:.10f}")
print(f"  √c = {math.sqrt(c_fp):.10f}")
print(f"  1/(H³ - 1) = {1/(H**3-1):.10f} = 1/26")
print(f"  c × (R+2)² = (R²+2)/26 = {(R_fp**2+2)/26:.6f}")

# Is c related to BORN_FLOOR?
print(f"  BORN_FLOOR = {BORN_FLOOR:.10f}")
print(f"  c / BORN_FLOOR = {c_fp / BORN_FLOOR:.10f}")

# The key ratio: is c ≈ BORN_FLOOR × (something)?
print(f"  c ≈ 1/(H+1)² × correction? 1/16 = {1/16:.6f} vs c = {c_fp:.6f}")


# ═══════════════════════════════════════════════════════════════
# The exact λ₀ in terms of the fixed point
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("ALGEBRAIC STRUCTURE OF λ₀")
print("="*80)

# λ₀ is determined by 6 constraint equations at H=3
# It cannot be expressed as a simple function of K* and H alone
# because it depends on the full nonlinear dynamics (Born floor)

# But we can ask: what polynomial does λ₀ satisfy?
# From the paper: the minimal polynomial has degree 32

# The relationship to crystal quantities:
# Crystal |λ₀| ≈ 0.244, DS λ₀ ≈ 0.283
# Ratio ≈ 1.16 ≈ N (structure constant)

# Let me check this precisely:
lambda_crystal = 0.24422  # from 2000-seed average
N_struct = 1.16198  # from 2000-seed average

print(f"\n  DS λ₀ = {lambda_0_DS:.10f}")
print(f"  Crystal λ₀ = {lambda_crystal:.10f}")
print(f"  Ratio DS/crystal = {lambda_0_DS / lambda_crystal:.10f}")
print(f"  Structure constant N = {N_struct:.10f}")
print(f"  Difference: {lambda_0_DS/lambda_crystal - N_struct:.6f}")
print(f"  Relative: {abs(lambda_0_DS/lambda_crystal - N_struct)/N_struct:.4f}")

# What if: DS_λ₀ = crystal_λ₀ × N ?
predicted = lambda_crystal * N_struct
print(f"\n  Crystal λ₀ × N = {predicted:.10f}")
print(f"  DS λ₀           = {lambda_0_DS:.10f}")
print(f"  Relative error   = {abs(predicted - lambda_0_DS)/lambda_0_DS:.4f}")

# Or: DS_λ₀ = crystal_λ₀ / L1  (where L1 ≈ 0.242 from raw product)
L1_raw = 0.24211  # from structure_constant_exact.py
predicted2 = lambda_crystal / L1_raw
print(f"\n  Crystal λ₀ / L1 = {predicted2:.10f}")
print(f"  DS λ₀             = {lambda_0_DS:.10f}")
print(f"  Relative error     = {abs(predicted2 - lambda_0_DS)/lambda_0_DS:.4f}")


# ═══════════════════════════════════════════════════════════════
# The approximate formula as physical insight
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE APPROXIMATE FORMULA: PHYSICAL INSIGHT")
print("="*80)

print(f"""
  Δ_gap ≈ (H/2) × ln(H) × (1-K*)    [0.06% error]
        = crystal_gap × H(1-K*)
        = (ln(H)/2) × (effective channels)

  The 0.06% correction comes from the nonlinear Born floor dynamics.
  The leading relationship IS:

    mass gap ≈ crystal spectral gap × effective dimension

  where:
    crystal spectral gap = ln(H)/2 (geometric, from the H-simplex)
    effective dimension = H(1-K*) (information-theoretic, non-conflicting channels)

  Numerically:
    crystal gap = {crystal_gap:.8f}
    H(1-K*) = {H*(1-K_STAR):.8f}
    product = {crystal_gap * H*(1-K_STAR):.8f}
    actual Δ = {delta_gap:.8f}
    ratio = {delta_gap / (crystal_gap * H*(1-K_STAR)):.8f}

  The ratio {delta_gap / (crystal_gap * H*(1-K_STAR)):.6f} ≈ 1 - 0.0006 is
  the Born floor correction: the nonlinear dynamics makes the gap
  SLIGHTLY smaller than the product formula predicts.
""")
