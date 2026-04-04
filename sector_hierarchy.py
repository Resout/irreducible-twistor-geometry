"""
The sector hierarchy: all physical scales from S/dim.

From Principle 91, the instanton action S = 810/7 divided by sector
dimension gives:
  dim=1  (sign):     exp(-S)    = 10⁻⁵⁰  → CC
  dim=5  (trivial):  exp(-S/5)  = 10⁻¹⁰  → θ_QCD
  dim=10 (standard): exp(-S/10) = 10⁻⁵   → m_e/m_W
  dim=16 (total):    exp(-S/16) = 10⁻³   → QCD/EW ratio

But there are MORE natural dimensions: H, H+1, H², (H-1)², etc.
Each gives a different physical scale. How many match?

Also: η_baryon = 7 × exp(-S/5) to 2%. What other physical constants
involve sector dimensions multiplied by crystal invariants?
"""

import math
from fractions import Fraction

H = 3
K_STAR = Fraction(7, 30)
BORN_FLOOR = Fraction(1, 27)
MASS_DIM = H + 1
S_FRAC = Fraction(810, 7)
S = float(S_FRAC)
ln10 = math.log(10)

# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("ALL NATURAL DIVISORS OF S AND THEIR PHYSICAL SCALES")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Natural divisors: sector dimensions, H-powers, combinations
divisors = {}

def add(name, val, origin=""):
    divisors[name] = (val, origin)

# Sector dimensions
add("dim(sign)", 1, "sign irrep")
add("dim(triv)", 5, "trivial irrep")
add("dim(std)", 10, "standard irrep")
add("dim(total)", 16, "MASS_DIM²")

# H-powers
add("H", 3, "hypotheses")
add("H²", 9, "H squared")
add("H³", 27, "H cubed = 1/BF")

# MASS_DIM powers
add("MASS_DIM", 4, "H+1")
add("MASS_DIM²", 16, "= dim(total)")

# Combined
add("2H", 6, "two hypotheses")
add("H(H+1)", 12, "H × MASS_DIM")
add("2^H", 8, "= 2×MASS_DIM")
add("H+2", 5, "= dim(triv)")

# From K*
add("1/K*", float(1/K_STAR), "inverse conflict")
add("H/K*", float(H/K_STAR), "= 90/7")
add("H²/K*", float(H**2/K_STAR), "= S/H = S_gauge")

# Known from hierarchy
add("S/H=S_gauge", S/H, "single gauge instanton")

# Other combinations
add("5!/S", 120/S, "= 1/K*... no, (S+1/K*)/S")
add("H(1-K*)", float(H * (1-K_STAR)), "effective channels from Pr.86")
add("2H(H+1)", 24, "= R_FS")

# Physical constants for comparison
constants = {
    "CC (Λ)":                      (-122, "Λ/m_Pl⁴, cosmological constant"),
    "CC (dimensionless)":          (-50.3, "Λ in natural units"),
    "EW hierarchy":                (-16.75, "v/m_Pl"),
    "θ_QCD":                       (-10.0, "strong CP angle upper bound"),
    "η_baryon":                    (-9.21, "baryon-to-photon ratio"),
    "α_EM":                        (-2.14, "1/137 fine structure"),
    "α_s":                         (-0.93, "strong coupling"),
    "α_W":                         (-1.47, "weak coupling"),
    "sin²θ_W":                     (-0.63, "Weinberg angle"),
    "m_e/m_W":                     (-5.20, "electron/W mass ratio"),
    "m_e/m_Z":                     (-5.25, "electron/Z mass ratio"),
    "m_μ/m_W":                     (-2.88, "muon/W mass ratio"),
    "m_τ/m_W":                     (-1.65, "tau/W mass ratio"),
    "m_t/v":                       (-0.15, "top/Higgs vev"),
    "QCD/EW ratio":                (-3.14, "Λ_QCD/v"),
    "m_e/m_p":                     (-3.26, "electron/proton"),
    "Λ_QCD/m_Pl":                  (-19.9, "QCD scale"),
    "V_us (Cabibbo)":              (-0.65, "CKM element"),
    "V_cb":                        (-1.38, "CKM element"),
    "V_ub":                        (-2.43, "CKM element"),
    "m_ν/m_τ (seesaw)":           (-10, "seesaw neutrino mass"),
    "Jarlskog J":                  (-4.64, "CP violation invariant"),
    "m_u/m_t":                     (-4.90, "up/top quark ratio"),
    "m_d/m_b":                     (-2.95, "down/bottom quark ratio"),
    "m_e/m_τ":                     (-3.54, "electron/tau ratio"),
}

print(f"\n  {'Divisor':>20s} {'value':>8s} {'S/div':>8s} {'exp(-S/div)':>12s} {'log₁₀':>8s}  Best match (≤0.3 orders)")
print("  " + "-" * 95)

matches = []

for name, (val, origin) in sorted(divisors.items(), key=lambda x: x[1][0]):
    if val <= 0:
        continue
    s_d = S / val
    exp_val = math.exp(-s_d)
    log_val = -s_d / ln10

    # Find closest physical constant
    best_name = ""
    best_err = float('inf')
    for cname, (clog, _) in constants.items():
        err = abs(log_val - clog)
        if err < best_err:
            best_err = err
            best_name = cname

    match_str = f"{best_name} (Δ={best_err:.2f})" if best_err < 0.3 else ""
    if best_err < 0.3:
        matches.append((name, val, log_val, best_name, best_err))

    print(f"  {name:>20s} {val:>8.3f} {s_d:>8.3f} {exp_val:>12.4e} {log_val:>8.4f}  {match_str}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("MATCHES WITHIN 0.3 ORDERS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

matches.sort(key=lambda x: x[4])
for name, val, log_val, cname, err in matches:
    clog = constants[cname][0]
    print(f"  exp(-S/{val:.3f}) = 10^({log_val:.4f})  ↔  {cname} = 10^({clog:.2f})  Δ={err:.4f} orders")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("COMBINED: SECTOR × CRYSTAL INVARIANT")
print("Constants = (crystal factor) × exp(-S/dim)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# η_baryon = 7 × exp(-S/5). What others?
# Try: constant = (integer or simple fraction) × exp(-S/dim)

crystal_factors = {
    "1": 1,
    "H": H,
    "H²-H+1": H**2-H+1,  # = 7
    "H²": H**2,
    "1/H": 1/H,
    "K*": float(K_STAR),
    "1-K*": float(1-K_STAR),
    "1/K*": float(1/K_STAR),
    "BF": float(BORN_FLOOR),
    "1/BF": float(1/BORN_FLOOR),
    "H+1": H+1,
    "MASS_DIM²": MASS_DIM**2,
}

sector_dims = [1, 5, 10, 16, H, H**2, MASS_DIM]

print(f"\n  Looking for (factor) × exp(-S/dim) = physical constant:")
print(f"  {'factor':>12s} × {'exp(-S/dim)':>12s} {'= value':>12s}  {'log₁₀':>8s}  Match (≤0.15 orders)")
print("  " + "-" * 80)

found = []
for fname, fval in crystal_factors.items():
    for dim in sector_dims:
        if dim <= 0:
            continue
        product = fval * math.exp(-S/dim)
        if product <= 0 or product >= 1:
            continue
        log_p = math.log10(product)

        for cname, (clog, cdesc) in constants.items():
            err = abs(log_p - clog)
            if err < 0.15:
                found.append((fname, dim, fval, product, log_p, cname, clog, err))

found.sort(key=lambda x: x[7])
seen = set()
for fname, dim, fval, product, log_p, cname, clog, err in found:
    key = (fname, dim, cname)
    if key in seen:
        continue
    seen.add(key)
    print(f"  {fname:>12s} × exp(-S/{dim}) = {product:.4e}  {log_p:>8.4f}  {cname} ({clog:.2f}) Δ={err:.4f}")

# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("THE COMPLETE SECTOR HIERARCHY")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print(f"""
  Instanton action S = 810/7 ≈ 116.

  SECTOR TOWER (S/dim of irrep sector):
  ┌─────────────────────────────────────────────────────────┐
  │ dim=1  (sign/chirality): exp(-S)    = 10⁻⁵⁰  → CC     │
  │ dim=5  (trivial/gravity): exp(-S/5) = 10⁻¹⁰  → θ_QCD  │
  │ dim=10 (standard/gauge): exp(-S/10) = 10⁻⁵   → m_e/m_W│
  │ dim=16 (total):          exp(-S/16) = 10⁻³   → QCD/EW  │
  └─────────────────────────────────────────────────────────┘

  INSTANTON TOWER (S × n/H):
  ┌─────────────────────────────────────────────────────────┐
  │ n=1: exp(-S/H)   = 10⁻¹⁷  → EW hierarchy (v/m_Pl)     │
  │ n=H: exp(-S)     = 10⁻⁵⁰  → CC (same as dim=1)        │
  └─────────────────────────────────────────────────────────┘

  COMBINED (factor × exp):
  ┌─────────────────────────────────────────────────────────┐
  │ (H²-H+1) × exp(-S/5) = 6.23×10⁻¹⁰ → η_baryon (2.1%)  │
  └─────────────────────────────────────────────────────────┘

  The CC appears at the intersection: dim=1 AND n=H.
  The chiral sector IS the multi-instanton level.
""")
