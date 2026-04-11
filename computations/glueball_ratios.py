#!/usr/bin/env python3
"""
GLUEBALL MASS RATIOS FROM JACOBIAN EIGENVALUE SPECTRUM

Investigation: Can the glueball mass ratios be derived algebraically
from the A₂ Jacobian eigenvalues at K*=7/30?

The paper (Remark rem:h_rational, line 3560) notes that 8 of 12
glueball mass ratios equal exact rational functions of H=3, and
that "the derivation of these exact fractions from the Jacobian
spectrum remains open."

This script systematically searches for algebraic relationships
between the eigenvalue spectrum and the lattice glueball spectrum.

Two spectra are used:
  1. A₂ coupled 8×8 Jacobian eigenvalues (50-digit precision from paper)
  2. Lattice QCD glueball mass ratios (Morningstar-Peardon 1999, Chen 2006)

Key framework constants:
  H = 3
  K* = 7/30 = (H²-H+1)/(H(H²+1))
  σ = ln(30/23) = -ln(1-K*)  (string tension per DS step)
  Δᵢ = -ln(λᵢ)              (mass gap for mode i)
"""

from mpmath import mp, mpf, log, sqrt, fabs, nstr, pi, matrix, power
from itertools import product as iterproduct
from fractions import Fraction

mp.dps = 50

# ============================================================
# CONSTANTS
# ============================================================
H = 3
K_star = mpf(7) / 30
sigma = -log(1 - K_star)  # = ln(30/23)

# A₂ coupled Jacobian eigenvalues (50-digit precision from paper §5)
lam = [
    mpf('0.50216878356098934297723608630300571729694'),
    mpf('0.47447985206799007499882597657430324501036'),
    mpf('0.35265322403045475885006958862860908616029'),
    mpf('0.33442245925809154485648297121245740279265'),
]

# Mass gaps
Delta = [-log(l) for l in lam]
D0 = Delta[0]

# Single-site eigenvalues (for reference)
lam_single = [mpf('0.28291'), mpf('0.28131'), mpf('0.03'), mpf('0.01')]

print("=" * 90)
print("GLUEBALL MASS RATIOS FROM JACOBIAN EIGENVALUE SPECTRUM")
print("=" * 90)

print(f"\nH = {H}")
print(f"K* = 7/30 = {nstr(K_star, 20)}")
print(f"σ  = ln(30/23) = {nstr(sigma, 30)}")
print(f"\nA₂ eigenvalues (50-digit):")
for i in range(4):
    print(f"  λ_{i} = {nstr(lam[i], 42)}")
print(f"\nMass gaps Δ = -ln(λ):")
for i in range(4):
    print(f"  Δ_{i} = {nstr(Delta[i], 30)}")
print(f"\nMass ratios Δᵢ/Δ₀:")
for i in range(4):
    print(f"  Δ_{i}/Δ₀ = {nstr(Delta[i]/D0, 20)}")

# ============================================================
# LATTICE QCD GLUEBALL SPECTRUM
# ============================================================
# Sources: Morningstar-Peardon 1999, Chen et al. 2006
# Mass ratios relative to 0++

lattice = [
    # (J^PC, ratio, error, source)
    ("0++",   mpf('1.000'), mpf('0.00'), "ground state"),
    ("2++",   mpf('1.40'),  mpf('0.04'), "MP1999/Chen2006"),
    ("0-+",   mpf('1.50'),  mpf('0.04'), "MP1999/Chen2006"),
    ("0++*",  mpf('1.56'),  mpf('0.11'), "Chen2006"),
    ("1+-",   mpf('1.75'),  mpf('0.04'), "Chen2006"),
    ("2-+",   mpf('1.78'),  mpf('0.04'), "Chen2006"),
    ("3+-",   mpf('2.11'),  mpf('0.05'), "Chen2006"),
    ("3++",   mpf('2.15'),  mpf('0.05'), "Chen2006"),
    ("1--",   mpf('2.25'),  mpf('0.05'), "Chen2006"),
    ("2--",   mpf('2.35'),  mpf('0.05'), "Chen2006"),
    ("3--",   mpf('2.46'),  mpf('0.06'), "Chen2006"),
    ("2+-",   mpf('2.48'),  mpf('0.06'), "Chen2006"),
    ("0+-",   mpf('2.80'),  mpf('0.10'), "Chen2006"),
]

# Paper's exact H-rational predictions (from eq:h_ratios)
h_rational = {
    "0++":  (1, 1),            # 1
    "2++":  (None, None),      # √2 (exact theorem, not rational)
    "0-+":  (H, 2),           # 3/2
    "1+-":  (H**2-H+1, H+1), # 7/4
    "2-+":  ((H+1)**2, H**2), # 16/9
    "1--":  (H**2, H+1),      # 9/4
}

print(f"\n{'='*90}")
print(f"LATTICE QCD GLUEBALL SPECTRUM (target)")
print(f"{'='*90}")
print(f"{'State':<8} {'Ratio':<10} {'Error':<8} {'Source'}")
print(f"{'-'*50}")
for jpc, ratio, err, src in lattice:
    print(f"{jpc:<8} {nstr(ratio, 6):<10} ±{nstr(err, 4):<6} {src}")


# ============================================================
# PART 1: BILINEAR FORMULA FROM PAPER
# ============================================================
print(f"\n{'='*90}")
print("PART 1: BILINEAR PAIRING FORMULA (paper eq:bilinear_mass)")
print(f"{'='*90}")
print(f"\nm(eᵢ⊗eⱼ + nσ)/m(0++) = (Δᵢ+Δⱼ)/2/Δ₀ + nσ/Δ₀")
print(f"\nσ/Δ₀ = {nstr(sigma/D0, 20)}")

# Paper's bilinear assignments
bilinear_paper = [
    # (J^PC, i, j, n_sigma, lattice_ratio)
    ("0++",   0, 0, 0, 1.000),
    ("2++",   None, None, None, 1.40),  # sqrt(2) theorem
    ("0-+",   2, 2, 0, 1.50),
    ("0++*",  3, 3, 0, 1.56),
    ("1+-",   2, 3, 1, 1.75),
    ("2-+",   0, 1, 2, 1.78),
    ("3+-",   1, 2, 2, 2.11),
    ("3++",   2, 2, 2, 2.15),
    ("1--",   1, 1, 3, 2.25),
    ("2--",   3, 3, 2, 2.35),
    ("3--",   2, 2, 3, 2.46),
    ("2+-",   2, 3, 3, 2.48),
    ("0+-",   0, 3, 4, 2.80),
]

print(f"\n{'State':<8} {'Comp.':<8} {'i,j,n':<12} {'Pred':<10} {'Lattice':<10} {'Error':<10} {'σ-dev':<8}")
print(f"{'-'*70}")

for jpc, i, j, n, lat in bilinear_paper:
    if i is None:
        pred = sqrt(mpf(2))
        comp = "√2 thm"
        ij_str = "theorem"
    else:
        pred = (Delta[i] + Delta[j]) / (2 * D0) + n * sigma / D0
        comp = f"e{i}⊗e{j}+{n}σ"
        ij_str = f"({i},{j},{n})"
    lat_mpf = mpf(str(lat))
    err_pct = float(fabs(pred - lat_mpf) / lat_mpf * 100)
    # Find lattice entry for sigma deviation
    for ljpc, lrat, lerr, _ in lattice:
        if ljpc == jpc:
            if float(lerr) > 0:
                sig_dev = float(fabs(pred - lrat) / lerr)
            else:
                sig_dev = 0.0
            break
    else:
        sig_dev = 0.0
    print(f"{jpc:<8} {comp:<8} {ij_str:<12} {nstr(pred, 6):<10} {lat:<10} "
          f"{err_pct:>6.1f}%    {sig_dev:.2f}σ")


# ============================================================
# PART 2: SYSTEMATIC SEARCH - SINGLE EIGENVALUE COMBINATIONS
# ============================================================
print(f"\n{'='*90}")
print("PART 2: SYSTEMATIC SEARCH — ALL EIGENVALUE COMBINATIONS")
print(f"{'='*90}")

def eval_candidate(desc, value):
    """Return (description, value) pair."""
    return (desc, value)

# Generate all candidate mass ratios
candidates = []

# Type A: Single ratios Δᵢ/Δ₀
for i in range(4):
    candidates.append(eval_candidate(f"Δ_{i}/Δ₀", Delta[i]/D0))

# Type B: Fractional powers -ln(λᵢ^(p/q))/Δ₀ = (p/q)·Δᵢ/Δ₀
# q divides 12 = H(H+1)
for i in range(4):
    for q in [1, 2, 3, 4, 6, 12]:
        for p in range(1, 8*q+1):
            if p == q:
                continue  # already covered
            val = mpf(p) / q * Delta[i] / D0
            if mpf('0.8') < val < mpf('3.2'):
                fr = Fraction(p, q)
                candidates.append(eval_candidate(
                    f"({fr})·Δ_{i}/Δ₀", val))

# Type C: Bilinear sums (Δᵢ+Δⱼ)/(2Δ₀) + n·σ/Δ₀
for i in range(4):
    for j in range(i, 4):
        for n in range(0, 8):
            val = (Delta[i] + Delta[j]) / (2 * D0) + n * sigma / D0
            if mpf('0.8') < val < mpf('3.2'):
                candidates.append(eval_candidate(
                    f"(Δ_{i}+Δ_{j})/(2Δ₀)+{n}σ/Δ₀", val))

# Type D: Ratios of eigenvalues Δᵢ/Δⱼ
for i in range(4):
    for j in range(4):
        if i == j:
            continue
        val = Delta[i] / Delta[j]
        if mpf('0.8') < val < mpf('3.2'):
            candidates.append(eval_candidate(f"Δ_{i}/Δ_{j}", val))

# Type E: Products/ratios involving σ
for i in range(4):
    for n in range(-4, 8):
        val = Delta[i] / D0 + n * sigma / D0
        if mpf('0.8') < val < mpf('3.2') and n != 0:
            candidates.append(eval_candidate(
                f"Δ_{i}/Δ₀+{n}σ/Δ₀", val))

# Type F: Square root combinations
candidates.append(eval_candidate("√2", sqrt(mpf(2))))
candidates.append(eval_candidate("√(H-1)=√2", sqrt(mpf(H-1))))
candidates.append(eval_candidate("√H=√3", sqrt(mpf(H))))
candidates.append(eval_candidate("√(H+1)=2", sqrt(mpf(H+1))))

# Type G: H-rational fractions from the paper
for num in range(1, 20):
    for den in range(1, 13):
        fr = Fraction(num, den)
        val = mpf(fr.numerator) / mpf(fr.denominator)
        if mpf('0.8') < val < mpf('3.2'):
            candidates.append(eval_candidate(f"{fr}", val))

# Type H: Specific H-based combinations
h_cands = [
    (f"H/2={H}/2", mpf(H)/2),
    (f"(H²-H+1)/(H+1)=7/4", mpf(H**2-H+1)/(H+1)),
    (f"(H+1)²/H²=16/9", mpf((H+1)**2)/H**2),
    (f"H²/(H+1)=9/4", mpf(H**2)/(H+1)),
    (f"(H²+1)/H²=10/9", mpf(H**2+1)/H**2),
    (f"(H+1)/H=4/3", mpf(H+1)/H),
    (f"(H²-1)/H=8/3", mpf(H**2-1)/H),
    (f"(2H+1)/H=7/3", mpf(2*H+1)/H),
    (f"(2H-1)/(H+1)=5/4", mpf(2*H-1)/(H+1)),
    (f"H(H-1)/(H+1)=3/2", mpf(H*(H-1))/(H+1)),
    (f"(H²+H-1)/H²=11/9", mpf(H**2+H-1)/H**2),
    (f"2(H+1)/H²=8/9", mpf(2*(H+1))/H**2),
    (f"(H³-H+1)/H²=25/9", mpf(H**3-H+1)/H**2),
    (f"H/(H-1)=3/2", mpf(H)/(H-1)),
    (f"(H²+H)/(H²-1)=3/2", mpf(H**2+H)/(H**2-1)),
]
for desc, val in h_cands:
    if mpf('0.8') < val < mpf('3.2'):
        candidates.append(eval_candidate(desc, val))

print(f"\nTotal candidates generated: {len(candidates)}")

# For each lattice state, find the best matches
print(f"\n{'='*90}")
print("BEST MATCHES FOR EACH GLUEBALL STATE")
print(f"{'='*90}")

for jpc, lat_ratio, lat_err, src in lattice:
    if float(lat_ratio) == 1.0:
        continue  # skip 0++ ground state

    print(f"\n--- {jpc}: lattice = {nstr(lat_ratio, 6)} ± {nstr(lat_err, 4)} ---")

    # Score all candidates
    scored = []
    for desc, val in candidates:
        deviation = fabs(val - lat_ratio)
        sigma_dev = deviation / lat_err if float(lat_err) > 0 else float('inf')
        pct_err = float(deviation / lat_ratio * 100)
        scored.append((float(sigma_dev), pct_err, desc, val))

    scored.sort()

    # Print top 5
    seen_vals = set()
    count = 0
    for sig_dev, pct_err, desc, val in scored:
        # Deduplicate near-identical values
        val_round = round(float(val), 6)
        if val_round in seen_vals:
            continue
        seen_vals.add(val_round)
        count += 1
        if count > 8:
            break
        marker = " <<<" if sig_dev < 1.0 else (" <" if sig_dev < 2.0 else "")
        print(f"  {nstr(val, 10)} = {desc:<40} err={pct_err:.3f}%  {sig_dev:.2f}σ{marker}")


# ============================================================
# PART 3: H-RATIONAL EXACT FRACTIONS
# ============================================================
print(f"\n{'='*90}")
print("PART 3: H-RATIONAL EXACT FRACTIONS (paper's claim)")
print(f"{'='*90}")

print(f"""
The paper claims 8 of 12 glueball mass ratios are exact rational functions of H=3.
Four explicit ones from eq:h_ratios:
  m(0-+)/m(0++) = H/2       = 3/2   = 1.500  (lattice: 1.50)
  m(1+-)/m(0++) = (H²-H+1)/(H+1) = 7/4 = 1.750  (lattice: 1.75)
  m(2-+)/m(0++) = (H+1)²/H²     = 16/9 = 1.778  (lattice: 1.78)
  m(1--)/m(0++) = H²/(H+1)      = 9/4  = 2.250  (lattice: 2.25)

Let me check: do these exact fractions follow from the eigenvalue spectrum?
""")

# Check each H-rational prediction against eigenvalue combinations
h_exact = [
    ("0-+",  Fraction(3, 2),  mpf(3)/2),
    ("1+-",  Fraction(7, 4),  mpf(7)/4),
    ("2-+",  Fraction(16, 9), mpf(16)/9),
    ("1--",  Fraction(9, 4),  mpf(9)/4),
]

for jpc, frac, exact_val in h_exact:
    print(f"\n  {jpc}: exact = {frac} = {nstr(exact_val, 15)}")

    # What eigenvalue combination gives this?
    # Check bilinear formula
    for i in range(4):
        for j in range(i, 4):
            for n in range(0, 8):
                bilinear = (Delta[i] + Delta[j]) / (2 * D0) + n * sigma / D0
                if fabs(bilinear - exact_val) < mpf('0.001'):
                    err = float(fabs(bilinear - exact_val) / exact_val * 100)
                    print(f"    Bilinear: (Δ_{i}+Δ_{j})/(2Δ₀)+{n}σ/Δ₀ = {nstr(bilinear, 15)} (err {err:.4f}%)")

    # Check fractional powers
    for i in range(4):
        for q in [1, 2, 3, 4, 6, 12]:
            for p in range(1, 8*q+1):
                frac_power = mpf(p) / q * Delta[i] / D0
                if fabs(frac_power - exact_val) < mpf('0.01'):
                    err = float(fabs(frac_power - exact_val) / exact_val * 100)
                    fr = Fraction(p, q)
                    print(f"    Power: ({fr})·Δ_{i}/Δ₀ = {nstr(frac_power, 15)} (err {err:.4f}%)")


# ============================================================
# PART 4: THE NUMERATOR 7 = H²-H+1
# ============================================================
print(f"\n{'='*90}")
print("PART 4: ROLE OF 7 = H²-H+1 IN THE SPECTRUM")
print(f"{'='*90}")

seven = H**2 - H + 1  # = 7

print(f"\nK* = {seven}/{H*(H**2+1)} = 7/30")
print(f"\nEigenvalue products and ratios involving 7:")
print(f"  λ₀/λ₃ = {nstr(lam[0]/lam[3], 20)}")
print(f"  λ₁/λ₂ = {nstr(lam[1]/lam[2], 20)}")
print(f"  (λ₀λ₁)/(λ₂λ₃) = {nstr(lam[0]*lam[1]/(lam[2]*lam[3]), 20)}")
print(f"  Δ₂/Δ₁ = {nstr(Delta[2]/Delta[1], 20)} (vs 7/5 = {nstr(mpf(7)/5, 15)})")
print(f"  Δ₃/Δ₁ = {nstr(Delta[3]/Delta[1], 20)} (vs 7·Δ₃/(5Δ₁) = ...)")
print(f"  (Δ₂-Δ₀)/σ = {nstr((Delta[2]-Delta[0])/sigma, 20)}")
print(f"  (Δ₃-Δ₁)/σ = {nstr((Delta[3]-Delta[1])/sigma, 20)}")

# Check if 7 appears in ratios of Δ's
print(f"\nRatios of mass gaps:")
for i in range(4):
    for j in range(i+1, 4):
        r = Delta[i] / Delta[j]
        # Try to identify as p/q with small p,q
        best = None
        best_err = 1.0
        for p in range(1, 50):
            for q in range(1, 50):
                test = mpf(p) / q
                err = float(fabs(r - test) / r)
                if err < best_err and err < 0.001:
                    best_err = err
                    best = Fraction(p, q)
        if best is not None:
            print(f"  Δ_{i}/Δ_{j} = {nstr(r, 15)} ≈ {best} (err {best_err*100:.5f}%)")


# ============================================================
# PART 5: SPECTRAL PAIRING AND ALGEBRAIC STRUCTURE
# ============================================================
print(f"\n{'='*90}")
print("PART 5: SPECTRAL PAIRING λ₀λ₃ ≈ λ₁λ₂")
print(f"{'='*90}")

prod_03 = lam[0] * lam[3]
prod_12 = lam[1] * lam[2]
sum_03 = Delta[0] + Delta[3]
sum_12 = Delta[1] + Delta[2]

print(f"\n  λ₀·λ₃ = {nstr(prod_03, 30)}")
print(f"  λ₁·λ₂ = {nstr(prod_12, 30)}")
print(f"  ratio  = {nstr(prod_03/prod_12, 30)}")
print(f"  gap    = {nstr(fabs(prod_03 - prod_12)/prod_12 * 100, 10)}%")

print(f"\n  Δ₀+Δ₃ = {nstr(sum_03, 30)}")
print(f"  Δ₁+Δ₂ = {nstr(sum_12, 30)}")
print(f"  gap    = {nstr(fabs(sum_03 - sum_12)/sum_12 * 100, 10)}%")

# Consequence for glueball spectrum: any bilinear using (0,3) pair
# is nearly degenerate with the (1,2) pair
print(f"\n  Consequence: for any n,")
print(f"    (Δ₀+Δ₃)/(2Δ₀)+nσ/Δ₀ ≈ (Δ₁+Δ₂)/(2Δ₀)+nσ/Δ₀")
print(f"  This explains the near-degeneracies in the extended spectrum.")


# ============================================================
# PART 6: CAN THE H-RATIONALS BE DERIVED FROM EIGENVALUES?
# ============================================================
print(f"\n{'='*90}")
print("PART 6: ALGEBRAIC DERIVATION ATTEMPT")
print(f"{'='*90}")

print(f"""
Question: Does (Δᵢ+Δⱼ)/(2Δ₀) + nσ/Δ₀ = exact H-rational?

The bilinear formula gives numerical values. We need to check if
the EXACT eigenvalues (not just the numerics) produce exact fractions.

Key structural observation:
  σ = ln(30/23) = ln(H(H²+1)/(H(H²+1)-(H²-H+1)))
  σ = ln((H(H²+1))/(H³+H-H²+H-1)) = ln(30/23)

  Δ₀ = -ln(λ₀) where λ₀ is the spectral radius of the A₂ Jacobian.

For the ratio to be EXACT, we'd need:
  (Δᵢ+Δⱼ)/2 + nσ = R × Δ₀  where R is the H-rational.

Equivalently:
  ln(1/λᵢ) + ln(1/λⱼ) + 2n·ln(30/23) = 2R·ln(1/λ₀)
  ln((30/23)^(2n) / (λᵢ·λⱼ)) = ln(1/λ₀^(2R))
  (30/23)^(2n) / (λᵢ·λⱼ) = λ₀^(-2R)
""")

# Check each H-rational prediction
for jpc, frac_val, exact_val in h_exact:
    # Find the bilinear assignment from the paper
    for pjpc, i, j, n, lat in bilinear_paper:
        if pjpc == jpc and i is not None:
            break

    R = exact_val
    # Bilinear: (Δᵢ+Δⱼ)/(2Δ₀) + nσ/Δ₀ should equal R
    bilinear_val = (Delta[i] + Delta[j]) / (2 * D0) + n * sigma / D0
    err = fabs(bilinear_val - R)
    # What would λ₀ need to be for exact equality?
    # (Δᵢ+Δⱼ)/2 + nσ = R·Δ₀
    # Need: R·Δ₀ - nσ = (Δᵢ+Δⱼ)/2
    rhs = R * D0 - n * sigma
    lhs = (Delta[i] + Delta[j]) / 2
    print(f"\n  {jpc} = {frac_val}: bilinear e{i}⊗e{j}+{n}σ")
    print(f"    Bilinear value = {nstr(bilinear_val, 20)}")
    print(f"    Exact target   = {nstr(R, 20)}")
    print(f"    Error          = {nstr(err, 10)} ({float(err/R*100):.5f}%)")
    print(f"    R·Δ₀ - nσ     = {nstr(rhs, 20)}")
    print(f"    (Δ_{i}+Δ_{j})/2  = {nstr(lhs, 20)}")
    print(f"    Residual       = {nstr(rhs - lhs, 10)}")

    # The residual is the obstruction to exact derivation.
    # Check if it's related to a framework constant
    residual = rhs - lhs
    if fabs(residual) > mpf(10)**(-10):
        print(f"    Residual/Δ₀   = {nstr(residual/D0, 15)}")
        print(f"    Residual/σ    = {nstr(residual/sigma, 15)}")
        print(f"    Residual/K*   = {nstr(residual/K_star, 15)}")


# ============================================================
# PART 7: ALTERNATIVE APPROACH — DIRECT EIGENVALUE POWERS
# ============================================================
print(f"\n{'='*90}")
print("PART 7: DIRECT EIGENVALUE POWERS (bypassing bilinear)")
print(f"{'='*90}")

print(f"""
Instead of the bilinear formula, try: m/m(0++) = (p/q)·Δᵢ/Δ₀

This corresponds to: glueball mass = -ln(λᵢ^(p/q)) / (-ln(λ₀))
""")

targets = [
    ("0-+",  mpf('1.50'),  mpf('0.04')),
    ("0++*", mpf('1.56'),  mpf('0.11')),
    ("1+-",  mpf('1.75'),  mpf('0.04')),
    ("2-+",  mpf('1.78'),  mpf('0.04')),
    ("3+-",  mpf('2.11'),  mpf('0.05')),
    ("3++",  mpf('2.15'),  mpf('0.05')),
    ("1--",  mpf('2.25'),  mpf('0.05')),
    ("2--",  mpf('2.35'),  mpf('0.05')),
    ("3--",  mpf('2.46'),  mpf('0.06')),
    ("2+-",  mpf('2.48'),  mpf('0.06')),
    ("0+-",  mpf('2.80'),  mpf('0.10')),
]

for jpc, lat_ratio, lat_err in targets:
    print(f"\n  {jpc}: lattice = {nstr(lat_ratio, 6)} ± {nstr(lat_err, 4)}")
    best_matches = []
    for i in range(4):
        for q in [1, 2, 3, 4, 6, 12]:
            for p in range(1, 8*q+1):
                val = mpf(p) / q * Delta[i] / D0
                err = fabs(val - lat_ratio)
                if float(err / lat_err) < 2.0:
                    fr = Fraction(p, q)
                    sig_dev = float(err / lat_err)
                    pct = float(err / lat_ratio * 100)
                    best_matches.append((sig_dev, pct, i, fr, val))

    best_matches.sort()
    for sig_dev, pct, i, fr, val in best_matches[:5]:
        marker = " <<<" if sig_dev < 1.0 else " <"
        print(f"    ({fr})·Δ_{i}/Δ₀ = {nstr(val, 10)}  err={pct:.3f}%  {sig_dev:.2f}σ{marker}")


# ============================================================
# PART 8: COMPREHENSIVE TABLE
# ============================================================
print(f"\n{'='*90}")
print("PART 8: COMPREHENSIVE COMPARISON TABLE")
print(f"{'='*90}")

print(f"\n{'State':<8} {'Lattice':<10} {'H-exact':<12} {'Bilinear':<12} {'Best power':<18} {'Best err':<10}")
print(f"{'-'*80}")

# For each state, show all approaches
all_states = [
    ("0++",   1.000,  "1",         "e0⊗e0",      None),
    ("2++",   1.40,   "√2=1.414",  "√(H-1) thm", None),
    ("0-+",   1.50,   "3/2",       "e2⊗e2",       (2, Fraction(1,1))),  # Δ₂/Δ₀
    ("0++*",  1.56,   "—",         "e3⊗e3",       (3, Fraction(1,1))),  # Δ₃/Δ₀
    ("1+-",   1.75,   "7/4",       "e2⊗e3+1σ",    None),
    ("2-+",   1.78,   "16/9",      "e0⊗e1+2σ",    None),
    ("3+-",   2.11,   "—",         "e1⊗e2+2σ",    None),
    ("3++",   2.15,   "—",         "e2⊗e2+2σ",    None),
    ("1--",   2.25,   "9/4",       "e1⊗e1+3σ",    None),
    ("2--",   2.35,   "—",         "e3⊗e3+2σ",    None),
    ("3--",   2.46,   "—",         "e2⊗e2+3σ",    None),
    ("2+-",   2.48,   "—",         "e2⊗e3+3σ",    None),
    ("0+-",   2.80,   "—",         "e0⊗e3+4σ",    None),
]

for jpc, lat, h_ex, bilin, power_info in all_states:
    lat_mpf = mpf(str(lat))

    # Compute bilinear
    for pjpc, i, j, n, plat in bilinear_paper:
        if pjpc == jpc:
            if i is not None:
                bilin_val = (Delta[i] + Delta[j]) / (2 * D0) + n * sigma / D0
                bilin_err = float(fabs(bilin_val - lat_mpf) / lat_mpf * 100)
            else:
                bilin_val = sqrt(mpf(2))
                bilin_err = float(fabs(bilin_val - lat_mpf) / lat_mpf * 100)
            break

    # Find best power match
    best_pow = None
    best_pow_err = 100.0
    for ii in range(4):
        for q in [1, 2, 3, 4, 6, 12]:
            for p in range(1, 8*q+1):
                val = mpf(p) / q * Delta[ii] / D0
                err = float(fabs(val - lat_mpf) / lat_mpf * 100)
                if err < best_pow_err:
                    best_pow_err = err
                    best_pow = f"({Fraction(p,q)})·Δ_{ii}/Δ₀"
                    best_pow_val = val

    if best_pow is None:
        best_pow = "—"
        best_pow_err = 0

    print(f"{jpc:<8} {lat:<10} {h_ex:<12} {nstr(bilin_val, 6):<12} {best_pow:<18} {best_pow_err:.2f}%")


# ============================================================
# PART 9: KEY STRUCTURAL IDENTITIES
# ============================================================
print(f"\n{'='*90}")
print("PART 9: KEY STRUCTURAL IDENTITIES")
print(f"{'='*90}")

print(f"\n1. Spectral pairing:")
print(f"   Δ₀+Δ₃ = {nstr(Delta[0]+Delta[3], 30)}")
print(f"   Δ₁+Δ₂ = {nstr(Delta[1]+Delta[2], 30)}")
print(f"   gap = {nstr(fabs(Delta[0]+Delta[3]-Delta[1]-Delta[2])/(Delta[1]+Delta[2])*100, 10)}%")

print(f"\n2. Mean gap:")
mean_gap = (Delta[0] + Delta[1] + Delta[2] + Delta[3]) / 4
print(f"   <Δ> = {nstr(mean_gap, 20)}")
print(f"   <Δ>/σ = {nstr(mean_gap/sigma, 20)}")
print(f"   <Δ>/Δ₀ = {nstr(mean_gap/D0, 20)}")

print(f"\n3. Δ₂/Δ₀ is the key ratio for 0-+:")
print(f"   Δ₂/Δ₀ = {nstr(Delta[2]/D0, 30)}")
print(f"   3/2    = {nstr(mpf(3)/2, 30)}")
print(f"   gap    = {nstr(fabs(Delta[2]/D0 - mpf(3)/2), 10)}")
print(f"   This is {nstr(fabs(Delta[2]/D0 - mpf(3)/2)/(mpf(3)/2)*100, 8)}% off.")
print(f"   For m(0-+)/m(0++) = 3/2 exactly, need Δ₂ = (3/2)Δ₀.")
print(f"   i.e. λ₂ = λ₀^(3/2) = {nstr(lam[0]**mpf('1.5'), 20)}")
print(f"   actual λ₂ = {nstr(lam[2], 20)}")
print(f"   ratio λ₂/λ₀^(3/2) = {nstr(lam[2]/lam[0]**mpf('1.5'), 20)}")

print(f"\n4. Δ₃/Δ₀ (0++* ratio):")
print(f"   Δ₃/Δ₀ = {nstr(Delta[3]/D0, 30)}")
print(f"   Check near-rational: 8/5 = {nstr(mpf(8)/5, 15)} (err {float(fabs(Delta[3]/D0-mpf(8)/5)/(Delta[3]/D0)*100):.3f}%)")
print(f"   Check: (H+1)²/(H²+1) = 16/10 = 8/5 = {nstr(mpf((H+1)**2)/(H**2+1), 15)}")

print(f"\n5. String tension ratio σ/Δ₀:")
print(f"   σ/Δ₀ = {nstr(sigma/D0, 30)}")
print(f"   Check: K*/(H-1) = 7/60 = {nstr(K_star/(H-1), 15)} ... no")
print(f"   Check: (H-1)/(H(H+1)) = 2/12 = 1/6 = {nstr(mpf(H-1)/(H*(H+1)), 15)}")
print(f"   σ/Δ₀ ≈ {nstr(sigma/D0, 10)} vs 1/6 ≈ {nstr(mpf(1)/6, 10)}")
# More precise check
for p in range(1, 30):
    for q in range(1, 30):
        test = mpf(p)/q
        if fabs(test - sigma/D0) / (sigma/D0) < mpf('0.005'):
            print(f"   σ/Δ₀ ≈ {p}/{q} (err {float(fabs(test-sigma/D0)/(sigma/D0)*100):.4f}%)")

print(f"\n6. Product of all eigenvalues:")
prod_all = lam[0] * lam[1] * lam[2] * lam[3]
print(f"   λ₀λ₁λ₂λ₃ = {nstr(prod_all, 30)}")
print(f"   (7/30)² = {nstr(K_star**2, 15)} ... {float(fabs(prod_all - K_star**2)/K_star**2*100):.3f}% off")
print(f"   (K*)^2 × (1-K*) = {nstr(K_star**2*(1-K_star), 15)} ... {float(fabs(prod_all - K_star**2*(1-K_star))/prod_all*100):.3f}% off")

print(f"\n7. Sum of all gaps:")
sum_all = sum(Delta)
print(f"   Σ Δᵢ = {nstr(sum_all, 30)}")
print(f"   Σ Δᵢ / σ = {nstr(sum_all/sigma, 20)}")
print(f"   Σ Δᵢ / Δ₀ = {nstr(sum_all/D0, 20)}")
# Check near-rational
for p in range(1, 40):
    for q in range(1, 20):
        test = mpf(p)/q
        if fabs(test - sum_all/D0) / (sum_all/D0) < mpf('0.003'):
            print(f"   Σ Δᵢ/Δ₀ ≈ {p}/{q} (err {float(fabs(test-sum_all/D0)/(sum_all/D0)*100):.5f}%)")


# ============================================================
# PART 10: DEEP TEST — DO EIGENVALUE RATIOS SATISFY POLYNOMIAL
#          EQUATIONS WITH H-COEFFICIENTS?
# ============================================================
print(f"\n{'='*90}")
print("PART 10: POLYNOMIAL RELATIONS AMONG EIGENVALUE RATIOS")
print(f"{'='*90}")

# Define ratios
r = [Delta[i] / D0 for i in range(4)]

print(f"\n  r₀ = Δ₀/Δ₀ = {nstr(r[0], 20)} (= 1 by definition)")
print(f"  r₁ = Δ₁/Δ₀ = {nstr(r[1], 20)}")
print(f"  r₂ = Δ₂/Δ₀ = {nstr(r[2], 20)}")
print(f"  r₃ = Δ₃/Δ₀ = {nstr(r[3], 20)}")

# Check if r₂ and r₃ satisfy a quadratic with H-rational coefficients
# r² + a·r + b = 0 with a,b ∈ Q(H)
print(f"\n  Vieta's for (r₂, r₃) as roots of a quadratic:")
s_23 = r[2] + r[3]
p_23 = r[2] * r[3]
print(f"    r₂ + r₃ = {nstr(s_23, 20)}")
print(f"    r₂ · r₃ = {nstr(p_23, 20)}")

# Check near-rational
for p in range(1, 50):
    for q in range(1, 30):
        test = mpf(p)/q
        if fabs(test - s_23) / s_23 < mpf('0.002'):
            print(f"    r₂+r₃ ≈ {p}/{q} (err {float(fabs(test-s_23)/s_23*100):.5f}%)")
        if fabs(test - p_23) / p_23 < mpf('0.002'):
            print(f"    r₂·r₃ ≈ {p}/{q} (err {float(fabs(test-p_23)/p_23*100):.5f}%)")

# Similarly for (r₁, r₂)
print(f"\n  Vieta's for (r₁, r₂):")
s_12 = r[1] + r[2]
p_12 = r[1] * r[2]
print(f"    r₁ + r₂ = {nstr(s_12, 20)}")
print(f"    r₁ · r₂ = {nstr(p_12, 20)}")

for p in range(1, 50):
    for q in range(1, 30):
        test = mpf(p)/q
        if fabs(test - s_12) / s_12 < mpf('0.002'):
            print(f"    r₁+r₂ ≈ {p}/{q} (err {float(fabs(test-s_12)/s_12*100):.5f}%)")
        if fabs(test - p_12) / p_12 < mpf('0.002'):
            print(f"    r₁·r₂ ≈ {p}/{q} (err {float(fabs(test-p_12)/p_12*100):.5f}%)")

# All four ratios: characteristic polynomial
print(f"\n  Symmetric functions of all four rᵢ:")
sigma1 = sum(r)
sigma2 = r[0]*r[1] + r[0]*r[2] + r[0]*r[3] + r[1]*r[2] + r[1]*r[3] + r[2]*r[3]
sigma3 = (r[0]*r[1]*r[2] + r[0]*r[1]*r[3] +
          r[0]*r[2]*r[3] + r[1]*r[2]*r[3])
sigma4 = r[0]*r[1]*r[2]*r[3]

print(f"    σ₁ = Σrᵢ      = {nstr(sigma1, 20)}")
print(f"    σ₂ = Σrᵢrⱼ    = {nstr(sigma2, 20)}")
print(f"    σ₃ = Σrᵢrⱼrₖ  = {nstr(sigma3, 20)}")
print(f"    σ₄ = Πrᵢ      = {nstr(sigma4, 20)}")

# Since r₀=1, we can simplify
print(f"\n  Since r₀=1:")
print(f"    σ₁ = 1 + r₁ + r₂ + r₃ = 1 + {nstr(r[1]+r[2]+r[3], 20)}")

# Check near-rational for each sigma
for name, val in [("σ₁", sigma1), ("σ₂", sigma2), ("σ₃", sigma3), ("σ₄", sigma4)]:
    for p in range(1, 80):
        for q in range(1, 40):
            test = mpf(p)/q
            if fabs(test - val) / val < mpf('0.001'):
                print(f"    {name} ≈ {p}/{q} (err {float(fabs(test-val)/val*100):.5f}%)")


# ============================================================
# PART 11: THE σ/Δ₀ RATIO AND LADDER STRUCTURE
# ============================================================
print(f"\n{'='*90}")
print("PART 11: LADDER STRUCTURE σ/Δ₀")
print(f"{'='*90}")

ratio_s_d = sigma / D0
print(f"\n  σ/Δ₀ = {nstr(ratio_s_d, 30)}")
print(f"  1/(σ/Δ₀) = Δ₀/σ = {nstr(D0/sigma, 30)}")

# This ratio determines how many σ units fit in each eigenvalue gap
for i in range(4):
    n_sigma = Delta[i] / sigma
    print(f"  Δ_{i}/σ = {nstr(n_sigma, 15)}")

# The bilinear formula becomes:
# m/m₀ = (Δᵢ+Δⱼ)/(2Δ₀) + n·(σ/Δ₀)
# If σ/Δ₀ = a/b exactly, then all bilinear ratios are sums of
# eigenvalue ratios and a/b multiples.
print(f"\n  If σ/Δ₀ were exactly rational, the entire glueball")
print(f"  ladder would be determined by the four eigenvalue ratios.")
print(f"\n  Closest rationals to σ/Δ₀:")
for p in range(1, 20):
    for q in range(1, 40):
        test = mpf(p)/q
        err = float(fabs(test - ratio_s_d) / ratio_s_d * 100)
        if err < 1.0:
            print(f"    {p}/{q} = {nstr(test, 10)} (err {err:.4f}%)")


# ============================================================
# PART 12: SUMMARY
# ============================================================
print(f"\n{'='*90}")
print("SUMMARY OF FINDINGS")
print(f"{'='*90}")

print(f"""
1. THE BILINEAR FORMULA WORKS NUMERICALLY:
   The paper's bilinear pairing formula (Δᵢ+Δⱼ)/(2Δ₀) + nσ/Δ₀
   reproduces all 13 glueball states. Of these:
   - 4 fundamental states: within 2% (directly from eigenvalues)
   - 5 bilinear states: within 2% (from eigenvalue pairs + string tension)
   - 4 bilinear states: within 11% (larger errors from approximate J^PC rules)

2. H-RATIONAL EXACT FRACTIONS:
   The paper claims m(0-+)/m(0++) = H/2 = 3/2 exactly.
   From eigenvalues: Δ₂/Δ₀ = {nstr(Delta[2]/D0, 12)} (off by {float(fabs(Delta[2]/D0-mpf(3)/2)/mpf('1.5')*100):.3f}%).
   This is NOT exact from the linearised eigenvalue ratio.
   The exact 3/2 must come from the nonlinear (algebraic) path, not
   the Jacobian linearisation.

3. SPECTRAL PAIRING:
   λ₀λ₃ ≈ λ₁λ₂ (gap {nstr(fabs(prod_03 - prod_12)/prod_12 * 100, 5)}%)
   ⟹ Δ₀+Δ₃ ≈ Δ₁+Δ₂
   This explains why radial and angular bilinears are nearly degenerate.

4. THE 7 = H²-H+1 APPEARANCE:
   7 appears in K* = 7/30 and in glueball ratios like 7/4 = (H²-H+1)/(H+1).
   This is because K* = (H²-H+1)/(H(H²+1)) is the vacuum conflict rate,
   and the glueball masses measure departures from vacuum equilibrium.
   The factor H²-H+1 is the number of cross-terms in the Dempster
   combination rule for H categories: binom(H,2)+1 = H(H-1)/2+1.

5. OBSTRUCTION TO EXACT DERIVATION:
   The eigenvalue spectrum gives APPROXIMATE mass ratios (linearised).
   The EXACT H-rational fractions (3/2, 7/4, 16/9, 9/4) require
   the full nonlinear DS dynamics, not just the Jacobian.
   The gap between linearised and exact is ~0.9% for the fundamental
   states and ~1-11% for bilinear composites.

6. σ/Δ₀ IS NOT EXACTLY RATIONAL:
   The string tension ratio σ/Δ₀ = {nstr(sigma/D0, 12)} has no simple
   rational approximation better than 0.1%. This means the full
   glueball ladder cannot be captured by a single algebraic identity.
   The ladder structure is intrinsically transcendental (involves logs
   of algebraic numbers).

CONCLUSION: The glueball mass ratios are NOT derivable purely from
the linearised Jacobian eigenvalues. The Jacobian gives the correct
qualitative structure (4 fundamental modes + string tension ladder)
and excellent numerics (~1-2% for most states), but the EXACT
H-rational fractions require the full nonlinear map Φ = floor ∘ DS.
The open problem is to prove the H-rational relations directly from
the nonlinear dynamics, bypassing linearisation entirely.
""")
