#!/usr/bin/env python3
"""
VERIFICATION: Hadron spectrum predictions from the paper.

Reproduces ALL 21 entries from the main hadron table (Remark rem:hadrons, lines 4419-4447)
and ALL 5 entries from the colour Regge trajectory (Remark rem:colour_regge, lines 4460-4471).

Formula:  m = -ln(lambda_i^{p/q}) * m(0++) / Delta_0
For entries with +n*sigma:  m = [-ln(lambda_i^{p/q}) + n*sigma] * m(0++) / Delta_0
Special case (pion): m = (Delta_1 - Delta_0) * m(0++) / Delta_0

Constants:
  lambda = [0.5022, 0.4745, 0.3527, 0.3344]
  Delta_0 = -ln(0.5022) = 0.6888
  m(0++) = 1710 MeV
  sigma = ln(30/23) = 0.2657
  scale = m(0++) / Delta_0 = 2482.6 MeV
"""

import numpy as np

# ============================================================
# CONSTANTS
# ============================================================
lam = np.array([0.5022, 0.4745, 0.3527, 0.3344])
Delta = -np.log(lam)
D0 = Delta[0]
m0pp = 1710.0  # MeV
sigma = np.log(30.0 / 23.0)
scale = m0pp / D0

print("=" * 90)
print("HADRON SPECTRUM VERIFICATION")
print("=" * 90)
print(f"\nEigenvalues:  lambda = {lam}")
print(f"Deltas:       Delta  = {Delta}")
print(f"Delta_0 = {D0:.6f}")
print(f"m(0++) = {m0pp} MeV")
print(f"sigma = ln(30/23) = {sigma:.6f}")
print(f"scale = m(0++)/Delta_0 = {scale:.2f} MeV")
print()

# ============================================================
# MAIN HADRON TABLE (21 entries from paper lines 4425-4445)
# ============================================================
# Format: (name, PDG_MeV, formula_description, compute_function, paper_predicted, paper_error)

def mass_from_eigenpower(i, p, q, n_sigma=0):
    """m = [-ln(lambda_i^{p/q}) + n*sigma] * scale"""
    val = -np.log(lam[i] ** (p / q)) + n_sigma * sigma
    return val * scale

def mass_pion():
    """m = (Delta_1 - Delta_0) * scale"""
    return (Delta[1] - Delta[0]) * scale

entries = [
    # (name, PDG, formula_str, lambda_idx, p, q, n_sigma, paper_pred, paper_err)
    # Special entries use lambda_idx = -1 for custom formulas
    ("pi+/-",      139.6,  "(Delta1-Delta0)*scale",    -1, 0, 0, 0, 140.9, "0.9%"),
    ("pi0",        135.0,  "(Delta1-Delta0)*scale",    -1, 0, 0, 0, 140.9, "4.3%"),
    ("eta",        547.9,  "lambda_0^{1/3}",            0, 1, 3, 0, 570.0, "4.0%"),
    ("N(938)",     938.3,  "lambda_1^{1/2}",            1, 1, 2, 0, 925.4, "1.4%"),
    ("eta'",       957.8,  "lambda_1^{1/2}",            1, 1, 2, 0, 925.4, "3.4%"),
    ("Lambda",    1115.7,  "lambda_0^{2/3}",            0, 2, 3, 0, 1139.9, "2.2%"),
    ("Delta(1232)",1232.0, "lambda_1^{2/3}",            1, 2, 3, 0, 1233.8, "0.1%"),
    ("Sigma",     1189.4,  "lambda_1^{2/3}",            1, 2, 3, 0, 1233.8, "3.7%"),
    ("Xi",        1314.9,  "lambda_2^{1/2}",            2, 1, 2, 0, 1293.6, "1.6%"),
    ("N(1440)",   1440.0,  "lambda_0^{5/6}",            0, 5, 6, 0, 1424.9, "1.0%"),
    ("Omega-",    1672.5,  "lambda_0^1",                0, 1, 1, 0, 1709.9, "2.2%"),
    ("D0",        1864.8,  "lambda_1^1",                1, 1, 1, 0, 1850.7, "0.8%"),
    ("D+/-",      1869.7,  "lambda_1^1",                1, 1, 1, 0, 1850.7, "1.0%"),
    ("Ds",        1968.3,  "lambda_2^{1/2}+sigma",      2, 1, 2, 1, 1953.2, "0.8%"),
    ("J/psi",     3096.9,  "lambda_1^{5/3}",            1, 5, 3, 0, 3084.6, "0.4%"),
    ("psi(2S)",   3686.1,  "lambda_0^1+3*sigma",        0, 1, 1, 3, 3688.8, "0.1%"),
    ("B+/-",      5279.3,  "lambda_2^1+4*sigma",        2, 1, 1, 4, 5225.7, "1.0%"),
    ("Bs",        5366.9,  "lambda_3^1+4*sigma",        3, 1, 1, 4, 5358.0, "0.2%"),
    ("Ups(1S)",   9460.3,  "lambda_3^{7/2}",            3, 7, 2, 0, 9518.1, "0.6%"),
    ("Ups(2S)",  10023.3,  "lambda_0^6",                0, 6, 1, 0, 10259.4, "2.4%"),
    ("Ups(3S)",  10355.2,  "lambda_2^4",                2, 4, 1, 0, 10348.8, "0.1%"),
]

print("=" * 90)
print(f"{'MAIN HADRON TABLE (21 entries)':^90}")
print("=" * 90)
print(f"{'Hadron':<14} {'PDG':>8} {'Formula':<26} {'Computed':>9} {'Paper':>9} "
      f"{'Comp-PDG':>8} {'Paper-PDG':>9} {'Match?':>7}")
print("-" * 90)

all_ok = True
for (name, pdg, formula, li, p, q, ns, paper_pred, paper_err) in entries:
    if li == -1:
        # Special: pion formula
        computed = mass_pion()
    else:
        computed = mass_from_eigenpower(li, p, q, ns)

    err_computed = abs(computed - pdg) / pdg * 100
    err_paper = abs(paper_pred - pdg) / pdg * 100

    # Check if our computation matches the paper's predicted value
    match = abs(computed - paper_pred) < 1.0  # within 1 MeV
    flag = "OK" if match else "** MISMATCH **"
    if not match:
        all_ok = False

    print(f"{name:<14} {pdg:>8.1f} {formula:<26} {computed:>9.1f} {paper_pred:>9.1f} "
          f"{err_computed:>7.1f}% {err_paper:>8.1f}%  {flag}")

print()

# ============================================================
# COLOUR REGGE TRAJECTORY (5 entries from paper lines 4464-4469)
# ============================================================

regge_entries = [
    # (name, power_str, p, q, paper_pred_MeV, PDG_MeV, paper_error_str)
    ("proton",    "1/2", 1, 2, 925,  938,  "1.4%"),
    ("Delta(1232)", "2/3", 2, 3, 1234, 1232, "0.2%"),
    ("D0",        "1",   1, 1, 1851, 1865, "0.8%"),
    ("Xi_c+",     "4/3", 4, 3, 2468, 2468, "0.004%"),
    ("J/psi",     "5/3", 5, 3, 3085, 3097, "0.4%"),
]

print("=" * 90)
print(f"{'COLOUR REGGE TRAJECTORY (lambda_1 powers)':^90}")
print("=" * 90)
print(f"{'Particle':<14} {'Power':<6} {'Computed':>9} {'Paper':>9} {'PDG':>8} "
      f"{'Comp-PDG':>8} {'Paper-PDG':>9} {'Match?':>7}")
print("-" * 90)

for (name, power_str, p, q, paper_pred, pdg, paper_err) in regge_entries:
    computed = mass_from_eigenpower(1, p, q, 0)  # all lambda_1
    err_computed = abs(computed - pdg) / pdg * 100
    err_paper = abs(paper_pred - pdg) / pdg * 100

    match = abs(computed - paper_pred) < 2.0  # within 2 MeV (paper rounds to integers)
    flag = "OK" if match else "** MISMATCH **"
    if not match:
        all_ok = False

    print(f"{name:<14} {power_str:<6} {computed:>9.1f} {paper_pred:>9.0f} {pdg:>8.0f} "
          f"{err_computed:>7.2f}% {err_paper:>8.3f}%  {flag}")

print()

# ============================================================
# DETAILED BREAKDOWN OF EACH COMPUTATION
# ============================================================
print("=" * 90)
print("DETAILED COMPUTATION LOG")
print("=" * 90)

print(f"\nPion: (Delta1-Delta0)*scale = ({Delta[1]:.6f} - {D0:.6f}) * {scale:.2f}")
print(f"    = {Delta[1]-D0:.6f} * {scale:.2f} = {(Delta[1]-D0)*scale:.2f} MeV")

print(f"\neta: -ln(0.5022^(1/3)) * scale = -ln({lam[0]**(1/3):.6f}) * {scale:.2f}")
print(f"    = {-np.log(lam[0]**(1/3)):.6f} * {scale:.2f} = {-np.log(lam[0]**(1/3))*scale:.2f} MeV")

print(f"\nN(938): -ln(0.4745^(1/2)) * scale = -ln({lam[1]**(1/2):.6f}) * {scale:.2f}")
print(f"    = {-np.log(lam[1]**(1/2)):.6f} * {scale:.2f} = {-np.log(lam[1]**(1/2))*scale:.2f} MeV")

print(f"\nXi_c+(2468): -ln(0.4745^(4/3)) * scale = -ln({lam[1]**(4/3):.6f}) * {scale:.2f}")
print(f"    = {-np.log(lam[1]**(4/3)):.6f} * {scale:.2f} = {-np.log(lam[1]**(4/3))*scale:.2f} MeV")

print(f"\nDs: [-ln(0.3527^(1/2)) + sigma] * scale = [{-np.log(lam[2]**(1/2)):.6f} + {sigma:.6f}] * {scale:.2f}")
print(f"    = {-np.log(lam[2]**(1/2)) + sigma:.6f} * {scale:.2f} = {(-np.log(lam[2]**(1/2)) + sigma)*scale:.2f} MeV")

print(f"\npsi(2S): [-ln(0.5022^1) + 3*sigma] * scale = [{-np.log(lam[0]):.6f} + {3*sigma:.6f}] * {scale:.2f}")
print(f"    = {-np.log(lam[0]) + 3*sigma:.6f} * {scale:.2f} = {(-np.log(lam[0]) + 3*sigma)*scale:.2f} MeV")

print(f"\nB+/-: [-ln(0.3527^1) + 4*sigma] * scale = [{-np.log(lam[2]):.6f} + {4*sigma:.6f}] * {scale:.2f}")
print(f"    = {-np.log(lam[2]) + 4*sigma:.6f} * {scale:.2f} = {(-np.log(lam[2]) + 4*sigma)*scale:.2f} MeV")

print(f"\nBs: [-ln(0.3344^1) + 4*sigma] * scale = [{-np.log(lam[3]):.6f} + {4*sigma:.6f}] * {scale:.2f}")
print(f"    = {-np.log(lam[3]) + 4*sigma:.6f} * {scale:.2f} = {(-np.log(lam[3]) + 4*sigma)*scale:.2f} MeV")

print(f"\nUps(1S): -ln(0.3344^(7/2)) * scale = -ln({lam[3]**(7/2):.10f}) * {scale:.2f}")
print(f"    = {-np.log(lam[3]**(7/2)):.6f} * {scale:.2f} = {-np.log(lam[3]**(7/2))*scale:.2f} MeV")

print(f"\nUps(2S): -ln(0.5022^6) * scale = -ln({lam[0]**6:.10f}) * {scale:.2f}")
print(f"    = {-np.log(lam[0]**6):.6f} * {scale:.2f} = {-np.log(lam[0]**6)*scale:.2f} MeV")

print(f"\nUps(3S): -ln(0.3527^4) * scale = -ln({lam[2]**4:.10f}) * {scale:.2f}")
print(f"    = {-np.log(lam[2]**4):.6f} * {scale:.2f} = {-np.log(lam[2]**4)*scale:.2f} MeV")

# ============================================================
# SUMMARY
# ============================================================
print()
print("=" * 90)
print("SUMMARY")
print("=" * 90)
if all_ok:
    print("ALL computed values match the paper's predicted values (within rounding).")
else:
    print("SOME values do NOT match the paper. See ** MISMATCH ** entries above.")

# Count entries within various error thresholds
within_2pct = 0
within_5pct = 0
for (name, pdg, formula, li, p, q, ns, paper_pred, paper_err) in entries:
    if li == -1:
        computed = mass_pion()
    else:
        computed = mass_from_eigenpower(li, p, q, ns)
    err = abs(computed - pdg) / pdg * 100
    if err < 2.0:
        within_2pct += 1
    if err < 5.0:
        within_5pct += 1

print(f"\nOf 21 hadrons:")
print(f"  Within 2% of PDG: {within_2pct}")
print(f"  Within 5% of PDG: {within_5pct}")
print(f"  Outside 5%:       {21 - within_5pct}")
