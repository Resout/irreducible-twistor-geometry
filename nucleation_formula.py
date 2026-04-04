"""
K_peak = 13/60 = (H²+H+1) / (2H(H²+1))

Verify this formula and explore its meaning.

K* = 7/30 = (H²-H+1) / (H(H²+1))

So:
  K* = (H²-H+1) / (H(H²+1))
  K_peak = (H²+H+1) / (2H(H²+1))

The numerators:
  K*:     H²-H+1 = 7   (counts disagreement channels)
  K_peak: H²+H+1 = 13  (counts... what?)

The denominators:
  K*:     H(H²+1) = 30
  K_peak: 2H(H²+1) = 60

So K_peak = K* × (H²+H+1) / (2(H²-H+1))
         = K* × 13/14
         = K* × (H²+H+1)/(2H²-2H+2)

Note: H²+H+1 = (H³-1)/(H-1) — the number of points in the projective
plane PG(2,H) = the number of elements in the field extension F_{H²}/F_H.
At H=3: 13 = |PG(2,3)|.

And H²-H+1 = (H³+1)/(H+1) — appears in cyclotomic structure.

K_peak/K* = (H²+H+1)/(2(H²-H+1)) = 13/14 at H=3.

The ratio 13/14: entropy peaks when K is 13/14 of the way to K*.
The crystal begins crystallizing one "step" below equilibrium conflict.

Physical meaning: at K_peak, the disorder-generating force of fresh
evidence (driving S up) exactly balances the order-generating force
of accumulated correlation (driving S down). Below K_peak, fresh
evidence dominates → entropy rises. Above K_peak, accumulated
structure dominates → entropy falls (crystallization).

Does K_peak = (H²+H+1)/(2H(H²+1)) hold at other H values?
"""

import math

print("=" * 80)
print("THE NUCLEATION FORMULA")
print("=" * 80)

for H in range(2, 8):
    K_star = (H**2 - H + 1) / (H * (H**2 + 1))
    K_peak = (H**2 + H + 1) / (2 * H * (H**2 + 1))
    DELTA = -math.log(1 - K_star)
    ratio = K_peak / K_star

    print(f"\n  H = {H}:")
    print(f"    K* = {K_star:.6f} = {H**2-H+1}/{H*(H**2+1)}")
    print(f"    K_peak = {K_peak:.6f} = {H**2+H+1}/{2*H*(H**2+1)}")
    print(f"    Ratio K_peak/K* = {ratio:.6f} = {H**2+H+1}/{2*(H**2-H+1)}")
    print(f"    Δ = {DELTA:.6f}")

print(f"\n\n{'='*80}")
print("THE TWO NUMERATORS")
print("="*80)

print(f"\n  H²-H+1: the DISAGREEMENT channels")
print(f"  H²+H+1: the TOTAL channels (agreement + disagreement)")
print(f"  Difference: 2H (the AGREEMENT channels)")
print()
for H in range(2, 8):
    disagreement = H**2 - H + 1
    total = H**2 + H + 1
    agreement = total - disagreement
    print(f"  H={H}: disagree={disagreement}, agree={agreement}=2H, total={total}")

print(f"\n  K* = disagree / (H × (H²+1))")
print(f"  K_peak = total / (2H × (H²+1))")
print(f"  K_peak = (K* × disagree + 2H) / (2 × H × (H²+1))")
print(f"         = K*/2 + 1/(H²+1)")
print()

for H in range(2, 8):
    K_star = (H**2 - H + 1) / (H * (H**2 + 1))
    K_peak = (H**2 + H + 1) / (2 * H * (H**2 + 1))
    alt = K_star / 2 + 1 / (H**2 + 1)
    print(f"  H={H}: K_peak={K_peak:.6f}, K*/2 + 1/(H²+1) = {alt:.6f}, match={abs(K_peak-alt)<1e-10}")


print(f"\n\n{'='*80}")
print("THE NUCLEATION BARRIER AS A BALANCE CONDITION")
print("="*80)

print(f"""
  K_peak = K*/2 + 1/(H²+1)

  The entropy peaks when the conflict rate K equals the arithmetic mean
  of two quantities:
    - K*/2: half the equilibrium conflict (the "background" conflict rate)
    - 1/(H²+1): the natural scale of a single mass coordinate

  At H=3: K*/2 = 7/60, 1/(H²+1) = 1/10 = 6/60
  Sum: 7/60 + 6/60 = 13/60 ✓

  The nucleation happens when K reaches the point where background
  conflict (K*/2) plus one mass coordinate's worth of information
  (1/(H²+1)) exactly balances. Below this, the system is still
  filling up mass coordinates. Above, the correlations start dominating.
""")

# Verify the decomposition identity
print("Verification of K_peak = K*/2 + 1/(H²+1):")
for H in range(2, 8):
    K_star = (H**2 - H + 1) / (H * (H**2 + 1))
    lhs = (H**2 + H + 1) / (2 * H * (H**2 + 1))
    rhs = K_star / 2 + 1 / (H**2 + 1)
    print(f"  H={H}: LHS={lhs:.10f}, RHS={rhs:.10f}, exact={abs(lhs-rhs)<1e-14}")
