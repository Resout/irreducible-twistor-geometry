import numpy as np

# ============================================================
# FRAMEWORK CONSTANTS
# ============================================================

# A₂ eigenvalues
lam = [0.5022, 0.4745, 0.3527, 0.3344]
Delta = [-np.log(l) for l in lam]
Delta_0 = Delta[0]  # 0.6888

# Single-site eigenvalues
lam_s = 0.28291  # radial (scalar)
lam_f = 0.28131  # angular (spinor)

# Other framework quantities
K_star = 7/30
sigma = -np.log(1 - K_star)  # string tension = 0.2657
sqrt2 = np.sqrt(2)
H = 3

# Scale: 0++ glueball = 1710 MeV
m_0pp = 1710  # MeV
scale = m_0pp / Delta_0  # MeV per DS unit

print("="*70)
print("FRAMEWORK CONSTANTS")
print("="*70)
print(f"A₂ eigenvalues: {[f'{l:.4f}' for l in lam]}")
print(f"A₂ gaps: {[f'{d:.4f}' for d in Delta]}")
print(f"Single-site: λ_s={lam_s:.5f}, λ_f={lam_f:.5f}")
print(f"K* = {K_star:.6f}, σ = {sigma:.4f}")
print(f"Scale: {scale:.1f} MeV per DS unit")
print(f"Δ₀ = {Delta_0:.4f}")

# ============================================================
# KNOWN MASSES: REVERSE SOLVE
# ============================================================

print("\n" + "="*70)
print("REVERSE SOLVE: KNOWN MASSES → FRAMEWORK EIGENVALUES")
print("="*70)

# All known masses in MeV
particles = {
    # Glueball spectrum (Chen 2006)
    "0++ glueball":     1710,
    "2++ glueball":     2390,
    "0-+ glueball":     2560,
    "1+- glueball":     2980,
    "2-+ glueball":     3040,
    "3+- glueball":     3600,
    "3++ glueball":     3670,
    "1-- glueball":     3830,
    "2-- glueball":     4010,
    "3-- glueball":     4200,
    "2+- glueball":     4230,
    "0+- glueball":     4780,
    
    # Constituent quarks
    "u constituent":    336,
    "d constituent":    340,
    "s constituent":    486,
    
    # Current quarks (MS bar, 2 GeV)
    "u current":        2.16,
    "d current":        4.70,
    "s current":        93.5,
    "c current":        1270,
    "b current":        4180,
    "t current":        173000,
    
    # Leptons
    "electron":         0.511,
    "muon":             105.7,
    "tau":              1777,
    
    # Gauge bosons
    "W":                80400,
    "Z":                91200,
    "Higgs":            125000,
    
    # Hadrons
    "pion":             140,
    "kaon":             494,
    "rho":              775,
    "proton":           938,
    "neutron":          940,
    "Delta baryon":     1232,
    "Omega baryon":     1672,
    
    # QCD scale
    "Lambda_QCD":       332,
    "sqrt(sigma)":      430,
}

# For each particle, compute:
# 1. Δ = m / scale (gap in DS units)
# 2. λ = exp(-Δ) (eigenvalue)
# 3. Δ/Δ₀ (ratio to 0++ glueball)
# 4. Try to express as combinations of framework quantities

print(f"\n{'Particle':<20} {'MeV':>8} {'Δ':>8} {'Δ/Δ₀':>7} {'λ':>10} {'Possible framework expression':>40}")
print("-"*100)

# Framework building blocks for pattern matching
building_blocks = {
    "Δ₀": Delta_0,
    "Δ₁": Delta[1],
    "Δ₂": Delta[2],
    "Δ₃": Delta[3],
    "σ": sigma,
    "K*": K_star,
    "1/H": 1/3,
    "1/H²": 1/9,
    "1/H³": 1/27,
    "ln2": np.log(2),
    "√2·Δ₀": sqrt2*Delta_0,
}

# Generate combinations
combos = {}

# Single quantities
for name, val in building_blocks.items():
    combos[name] = val

# Integer multiples of eigenvalues (Koopman products)
for i in range(4):
    for n in range(1, 6):
        combos[f"{n}·Δ_{i}"] = n * Delta[i]

# Sums of eigenvalues
for i in range(4):
    for j in range(i, 4):
        combos[f"Δ_{i}+Δ_{j}"] = Delta[i] + Delta[j]
        for k in range(j, 4):
            combos[f"Δ_{i}+Δ_{j}+Δ_{k}"] = Delta[i] + Delta[j] + Delta[k]

# Eigenvalue + string tension
for i in range(4):
    for n in range(1, 4):
        combos[f"Δ_{i}+{n}σ"] = Delta[i] + n*sigma
        combos[f"Δ_{i}-{n}σ"] = Delta[i] - n*sigma

# Products with √2
for i in range(4):
    combos[f"√2·Δ_{i}"] = sqrt2 * Delta[i]

# Ratios involving H
for i in range(4):
    combos[f"Δ_{i}/H"] = Delta[i] / H
    combos[f"Δ_{i}/(H+1)"] = Delta[i] / (H+1)
    combos[f"Δ_{i}·H"] = Delta[i] * H

# String tension combos
for n in range(1, 8):
    combos[f"{n}σ"] = n * sigma

# Single-site eigenvalue combos
for n in range(1, 6):
    combos[f"{n}·Δ_s"] = -n * np.log(lam_s)
    combos[f"{n}·Δ_f"] = -n * np.log(lam_f)

# Mixed
combos["Δ_s+Δ_f"] = -np.log(lam_s) - np.log(lam_f)
combos["Δ₀-Δ₁"] = Delta[0] - Delta[1]
combos["Δ₂-Δ₃"] = Delta[2] - Delta[3]
combos["(Δ₁-Δ₀)"] = Delta[1] - Delta[0]
combos["Δ₀·(Δ₁-Δ₀)/Δ₁"] = Delta[0] * (Delta[1]-Delta[0])/Delta[1]

# H-based ratios
combos["Δ₀/27"] = Delta_0/27
combos["Δ₀·K*"] = Delta_0 * K_star
combos["σ/H"] = sigma/H
combos["σ·K*"] = sigma * K_star
combos["Δ₀·σ"] = Delta_0 * sigma

for pname, mass in sorted(particles.items(), key=lambda x: x[1]):
    delta_p = mass / scale
    lambda_p = np.exp(-delta_p)
    ratio = delta_p / Delta_0
    
    # Find best matching combination
    best_match = ""
    best_err = float('inf')
    for cname, cval in combos.items():
        if cval <= 0: continue
        err = abs(delta_p - cval) / delta_p
        if err < best_err:
            best_err = err
            best_match = f"{cname} = {cval:.4f}"
    
    match_str = f"{best_match} ({best_err*100:.1f}%)" if best_err < 0.15 else "no match <15%"
    
    print(f"{pname:<20} {mass:>8.1f} {delta_p:8.4f} {ratio:7.3f} {lambda_p:10.6f} {match_str:>45}")

# ============================================================
# FOCUSED: GLUEBALL SPECTRUM RATIOS
# ============================================================

print("\n" + "="*70)
print("GLUEBALL SPECTRUM: RATIOS AND FRAMEWORK MATCHES")
print("="*70)

glueballs = [
    ("0++",  1.000), ("2++",  1.40), ("0-+",  1.50),
    ("1+-",  1.75),  ("2-+",  1.78), ("3+-",  2.11),
    ("3++",  2.15),  ("1--",  2.25), ("2--",  2.35),
    ("3--",  2.46),  ("2+-",  2.48), ("0+-",  2.80),
]

print(f"\n{'J^PC':<8} {'m/m₀':>6} {'Δ':>8} {'Best match':>35} {'Error':>8}")
print("-"*75)

for jpc, ratio in glueballs:
    delta_g = ratio * Delta_0
    
    best_match = ""
    best_err = float('inf')
    for cname, cval in combos.items():
        if cval <= 0: continue
        err = abs(delta_g - cval) / delta_g
        if err < best_err:
            best_err = err
            best_match = cname
    
    print(f"{jpc:<8} {ratio:6.2f} {delta_g:8.4f} {best_match:>35} {best_err*100:7.1f}%")

# ============================================================
# FOCUSED: MASS RATIOS BETWEEN PARTICLE CLASSES
# ============================================================

print("\n" + "="*70)
print("INTER-PARTICLE RATIOS")
print("="*70)

# Key ratios
ratios_to_check = [
    ("proton/0++", 938/1710),
    ("pion/0++", 140/1710),
    ("electron/0++", 0.511/1710),
    ("muon/0++", 105.7/1710),
    ("tau/0++", 1777/1710),
    ("W/0++", 80400/1710),
    ("Z/0++", 91200/1710),
    ("Higgs/0++", 125000/1710),
    ("u_const/0++", 336/1710),
    ("Lambda_QCD/0++", 332/1710),
    ("sqrt(sigma)/0++", 430/1710),
    ("proton/sqrt(sigma)", 938/430),
    ("0++/sqrt(sigma)", 1710/430),
    ("W/proton", 80400/938),
    ("tau/muon", 1777/105.7),
    ("muon/electron", 105.7/0.511),
    ("t/b", 173000/4180),
    ("b/c", 4180/1270),
    ("c/s", 1270/93.5),
    ("s/d", 93.5/4.70),
    ("d/u", 4.70/2.16),
]

print(f"\n{'Ratio':<25} {'Value':>10} {'ln(ratio)':>10} {'Framework?':>30}")
print("-"*80)

for name, val in ratios_to_check:
    ln_val = np.log(val) if val > 0 else 0
    
    # Check if ln(ratio) matches a framework quantity
    best = ""
    best_err = float('inf')
    
    # Check against ratios of framework quantities
    framework_ratios = {
        "1": 1.0,
        "H": 3.0, "H+1": 4.0, "H²": 9.0, "H²-1": 8.0,
        "H³": 27.0, "1/H": 1/3, "1/(H+1)": 0.25, "1/H³": 1/27,
        "K*": 7/30, "1-K*": 23/30, "σ/Δ₀": sigma/Delta_0,
        "Δ₁/Δ₀": Delta[1]/Delta[0], "Δ₂/Δ₀": Delta[2]/Delta[0],
        "Δ₃/Δ₀": Delta[3]/Delta[0],
        "√2": sqrt2, "√3": np.sqrt(3), "√(H-1)": np.sqrt(2),
        "2π": 2*np.pi, "π": np.pi,
        "λ₀": lam[0], "λ₁": lam[1], "λ₂": lam[2], "λ₃": lam[3],
        "λ₀²": lam[0]**2, "λ₁²": lam[1]**2,
        "λ₀³": lam[0]**3, "λ₁³": lam[1]**3,
    }
    
    for fname, fval in framework_ratios.items():
        err = abs(val - fval) / val if val > 0 else float('inf')
        if err < best_err:
            best_err = err
            best = f"{fname}={fval:.4f}"
        # Also check ln
        if ln_val != 0:
            err_ln = abs(ln_val - fval) / abs(ln_val)
            if err_ln < best_err:
                best_err = err_ln
                best = f"ln=>{fname}={fval:.4f}"
    
    match_str = f"{best} ({best_err*100:.1f}%)" if best_err < 0.20 else ""
    print(f"{name:<25} {val:10.4f} {ln_val:10.4f} {match_str:>30}")

# ============================================================
# REVERSE SOLVE: WHAT λ GIVES EACH KNOWN MASS?
# ============================================================

print("\n" + "="*70)
print("REVERSE EIGENVALUES: WHAT λ PRODUCES EACH MASS?")
print("="*70)

print(f"\n{'Particle':<20} {'MeV':>8} {'Required λ':>12} {'Nearest framework λ':>25}")
print("-"*70)

framework_lambdas = {
    "λ₀(A₂)": lam[0],
    "λ₁(A₂)": lam[1],
    "λ₂(A₂)": lam[2],
    "λ₃(A₂)": lam[3],
    "λ_s(single)": lam_s,
    "λ_f(single)": lam_f,
    "λ₀²": lam[0]**2,
    "λ₁²": lam[1]**2,
    "λ₀·λ₁": lam[0]*lam[1],
    "λ₀³": lam[0]**3,
    "λ₁³": lam[1]**3,
    "λ₀²·λ₁": lam[0]**2*lam[1],
    "λ₀·λ₁²": lam[0]*lam[1]**2,
    "λ₂²": lam[2]**2,
    "λ₂·λ₃": lam[2]*lam[3],
    "λ₃²": lam[3]**2,
    "λ₀⁴": lam[0]**4,
    "λ₁⁴": lam[1]**4,
    "√λ₀": np.sqrt(lam[0]),
    "√λ₁": np.sqrt(lam[1]),
    "√λ₂": np.sqrt(lam[2]),
    "√λ₃": np.sqrt(lam[3]),
    "λ₀^(1/3)": lam[0]**(1/3),
    "λ₁^(1/3)": lam[1]**(1/3),
    "λ₂^(1/3)": lam[2]**(1/3),
    "e^(-σ)": np.exp(-sigma),
    "e^(-2σ)": np.exp(-2*sigma),
    "λ₀·e^(-σ)": lam[0]*np.exp(-sigma),
    "λ₁·e^(-σ)": lam[1]*np.exp(-sigma),
    "λ₀·e^(-2σ)": lam[0]*np.exp(-2*sigma),
}

for pname, mass in sorted(particles.items(), key=lambda x: x[1]):
    delta_p = mass / scale
    lambda_p = np.exp(-delta_p)
    
    best = ""
    best_err = float('inf')
    for fname, fval in framework_lambdas.items():
        if fval <= 0 or fval >= 1: continue
        err = abs(lambda_p - fval) / lambda_p if lambda_p > 1e-15 else float('inf')
        if err < best_err:
            best_err = err
            best = f"{fname}={fval:.6f}"
    
    if lambda_p > 1e-15:
        match_str = f"{best} ({best_err*100:.1f}%)" if best_err < 0.30 else "no match"
    else:
        match_str = "λ ≈ 0 (too heavy)"
    
    print(f"{pname:<20} {mass:>8.1f} {lambda_p:>12.8f} {match_str}")

# ============================================================
# PATTERN SEARCH: LOGARITHMIC MASS RATIOS
# ============================================================

print("\n" + "="*70)
print("LOGARITHMIC MASS SPECTRUM")
print("="*70)

print(f"\nAll masses as multiples of Δ₀ (= m(0++)):")
print(f"Looking for integer, half-integer, or H-related patterns\n")

for pname, mass in sorted(particles.items(), key=lambda x: x[1]):
    ratio = mass / m_0pp
    if ratio > 0:
        log_ratio = np.log(ratio)
        delta_ratio = -np.log(np.exp(-Delta_0 * ratio)) / Delta_0  # = ratio
        
        # Check if ratio is close to n/H, n/(H+1), n/H², etc.
        for denom in [1, 2, 3, 4, 6, 8, 9, 27]:
            for numer in range(1, 300):
                test = numer / denom
                if abs(ratio - test) / ratio < 0.02:
                    print(f"  {pname:<20} m/m₀ = {ratio:.4f} ≈ {numer}/{denom} = {test:.4f} ({abs(ratio-test)/ratio*100:.1f}%)")
                    break
            else:
                continue
            break

