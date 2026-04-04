"""
What is K_peak = 0.217?

The entropy peaks at K ≈ 0.217, before K reaches K* = 0.233.
This is the nucleation barrier — the K-value where entropy begins falling,
meaning the crystal starts crystallizing.

Is 0.217 a simple function of H?
"""

import math

H = 3
K_STAR = (H**2 - H + 1) / (H * (H**2 + 1))  # 7/30
DELTA = -math.log(1 - K_STAR)
BORN_FLOOR = 1 / H**3

measured = 0.2167

candidates = {
    "K*": K_STAR,
    "K* - BORN_FLOOR": K_STAR - BORN_FLOOR,
    "K* × (1 - BORN_FLOOR)": K_STAR * (1 - BORN_FLOOR),
    "K* × (H/(H+1))": K_STAR * H / (H+1),
    "K* × (1 - 1/H²)": K_STAR * (1 - 1/H**2),
    "(H-1)/H²": (H-1) / H**2,
    "7/H(H²+1) - 1/H³": K_STAR - BORN_FLOOR,
    "Δ/(Δ+1)": DELTA / (DELTA + 1),
    "1 - exp(-Δ/2)": 1 - math.exp(-DELTA/2),
    "K* × exp(-Δ/2)": K_STAR * math.exp(-DELTA/2),
    "2K*/e": 2*K_STAR/math.e,
    "K*²/K*": K_STAR,  # same
    "7/32": 7/32,
    "7/33": 7/33,
    "13/60": 13/60,
    "(H²-H+1)/(H³+H)": (H**2-H+1)/(H**3+H),
    "K* - K*/H²": K_STAR - K_STAR/H**2,
    "(H²-H)/H(H²+1)": (H**2-H)/(H*(H**2+1)),
    "6/H(H²+1)": 6/(H*(H**2+1)),
    "1/H + 1/H³": 1/H + 1/H**3,
    "(2H-1)/(H(H+1)²)": (2*H-1)/(H*(H+1)**2),
    "exp(-Δ) × K*/(1-K*)": math.exp(-DELTA) * K_STAR / (1-K_STAR),
    "K* × (1-K*)": K_STAR * (1-K_STAR),
    "K* - Δ/H²": K_STAR - DELTA/H**2,
    "(H²-1)/(H(H²+1))": (H**2-1)/(H*(H**2+1)),
    "8/H(H²+1)": 8/(H*(H**2+1)),
    "Δ - BORN_FLOOR": DELTA - BORN_FLOOR,
    "K* × exp(-BORN_FLOOR)": K_STAR * math.exp(-BORN_FLOOR),
    "K* - 1/H(H+1)": K_STAR - 1/(H*(H+1)),
    "H/(H²+1) - 1/H(H²+1)": H/(H**2+1) - 1/(H*(H**2+1)),
    "(H-1)²/H(H²+1)": (H-1)**2 / (H*(H**2+1)),
    "4/H(H²+1)": 4/(H*(H**2+1)),
    "1/H - 1/(H²+1)": 1/H - 1/(H**2+1),
    "(H²-2)/(H(H²+1))": (H**2-2)/(H*(H**2+1)),
    # Deeper: what if K_peak is where dS/dK = 0?
    # S_peak ≈ 2.15, near ln(H² - 1) = ln(8)?
    "ln(H²-1) = ln(8)": math.log(H**2 - 1),
}

print("K_peak = 0.2167 ± 0.008\n")
print(f"{'expression':>30s}  {'value':>8s}  {'diff':>8s}  {'σ from meas':>12s}")
print("-" * 65)
sigma = 0.008
for expr, val in sorted(candidates.items(), key=lambda x: abs(x[1] - measured)):
    diff = val - measured
    nsig = abs(diff) / sigma
    marker = "  ✓" if nsig < 1 else "  ~" if nsig < 2 else ""
    print(f"  {expr:>28s}  {val:>8.5f}  {diff:>+8.5f}  {nsig:>8.1f}σ{marker}")


# What about the ENTROPY at the peak?
print(f"\n\nEntropy at peak ≈ 2.15")
s_peak = 2.15
s_candidates = {
    "ln(H²-1) = ln(8)": math.log(H**2 - 1),
    "ln(H² + H - 1) = ln(11)": math.log(H**2 + H - 1),
    "H × ln(H-1) + Δ": H * math.log(H-1) + DELTA,
    "(H+1)/2 × ln(H)": (H+1)/2 * math.log(H),
    "S_max × K*/(1-K*)": math.log(16) * K_STAR / (1-K_STAR),
    "π²/H²": math.pi**2 / H**2,
    "H - Δ": H - DELTA,
    "ln(2) + ln(H)": math.log(2) + math.log(H),
    "2.0 + Δ/2": 2.0 + DELTA/2,
    "2 + 7/H(H²+1)": 2 + K_STAR,
}
print(f"\n{'expression':>30s}  {'value':>8s}  {'diff':>8s}")
for expr, val in sorted(s_candidates.items(), key=lambda x: abs(x[1] - s_peak)):
    print(f"  {expr:>28s}  {val:>8.5f}  {val-s_peak:>+8.5f}")
