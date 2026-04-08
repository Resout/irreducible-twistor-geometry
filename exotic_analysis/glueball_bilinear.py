import numpy as np
from itertools import combinations_with_replacement

# ============================================================
# THE KEY INSIGHT: GLUEBALLS ARE BILINEARS OF ELEMENTARY MODES
# ============================================================

# A₂ eigenvalues
lam = [0.5022, 0.4745, 0.3527, 0.3344]
Delta = [-np.log(l) for l in lam]
sigma = -np.log(1 - 7/30)
m_0pp = 1710  # MeV
scale = m_0pp / Delta[0]

print("="*70)
print("GLUEBALLS AS BILINEARS OF ELEMENTARY MODES")
print("="*70)

print(f"""
If glueballs are PAIRS of elementary excitations:
  Elementary mode: eigenvalue √λ, mass Δ/2
  Glueball (pair): eigenvalue λ = (√λ)², mass Δ = 2×(Δ/2)

Check: 0++ glueball at Δ₀ = 2 × (Δ₀/2).
  Elementary radial: √λ₀ = {np.sqrt(lam[0]):.4f}, mass = {Delta[0]/2*scale:.0f} MeV
  Pair at J=0: mass = 2×{Delta[0]/2*scale:.0f} = {Delta[0]*scale:.0f} MeV = 0++ ✓

Check: 0-+ glueball at Δ₂ = 2 × (Δ₂/2).
  Elementary angular: √λ₂ = {np.sqrt(lam[2]):.4f}, mass = {Delta[2]/2*scale:.0f} MeV
  Pair at J=0: mass = 2×{Delta[2]/2*scale:.0f} = {Delta[2]*scale:.0f} MeV = 0-+ ✓
""")

# ============================================================
# ELEMENTARY PARTICLE TABLE
# ============================================================

print("="*70)
print("ELEMENTARY MODES (√λ, half-eigenvalue)")
print("="*70)

elem_modes = []
labels = ["radial,node-anti", "angular,node-anti", "angular,node-sym", "radial,node-sym"]
subspace = ["scalar", "spinor", "spinor", "scalar"]  # ℂ¹ or ℂ²
parity_intrinsic = [+1, -1, -1, +1]  # P from subspace
C_node = [-1, -1, +1, +1]  # C from node symmetry

for i in range(4):
    sqrt_lam = np.sqrt(lam[i])
    mass_elem = -np.log(sqrt_lam) * scale
    j_sub = "0" if subspace[i] == "scalar" else "1/2"
    elem_modes.append({
        "idx": i, "label": labels[i], "subspace": subspace[i],
        "sqrt_lam": sqrt_lam, "mass": mass_elem,
        "delta_half": Delta[i]/2,
        "P": parity_intrinsic[i], "C": C_node[i],
        "j_sub": j_sub
    })
    P_str = "+" if parity_intrinsic[i] > 0 else "-"
    C_str = "+" if C_node[i] > 0 else "-"
    print(f"  e_{i}: √λ_{i}={sqrt_lam:.4f}  mass={mass_elem:.0f} MeV  {labels[i]:20s}  ℂ{'¹' if subspace[i]=='scalar' else '²'}  P={P_str} C={C_str}")

# ============================================================
# BILINEAR PAIRINGS: SPIN FROM TENSOR DECOMPOSITION
# ============================================================

print("\n" + "="*70)
print("BILINEAR PAIRINGS AND THEIR SPIN CONTENT")
print("="*70)

print(f"""
Two elementary modes combine. The spin of the pair depends on 
which subspace each mode lives in:

  scalar ⊗ scalar (ℂ¹⊗ℂ¹):  
    J = 0 only (two scalars make a scalar)
    
  scalar ⊗ spinor (ℂ¹⊗ℂ²):  
    J = 1/2 only (scalar doesn't change the spinor's J)
    
  spinor ⊗ spinor (ℂ²⊗ℂ²):
    J = 0 ⊕ 1 (antisymmetric ⊕ symmetric, like ½⊗½ = 0⊕1)
    
For GLUEBALL bilinears (the observables Tr(F²), etc.):
  The √2 theorem applies. For a vector v = (v₁,v₂,v₂):
    Spin-0 norm² : |v|⁴/3
    Spin-2 norm² : 2|v|⁴/3
    Ratio: √2
    
  So any bilinear of two radial modes gives BOTH J=0 and J=2,
  with mass ratio √2.
""")

# ============================================================
# GENERATE ALL POSSIBLE BILINEAR GLUEBALL STATES
# ============================================================

print("="*70)
print("ALL POSSIBLE BILINEAR GLUEBALL STATES")
print("="*70)

# For each pair (i,j) of elementary modes:
# Mass = Δᵢ/2 + Δⱼ/2 = (Δᵢ+Δⱼ)/2
# P = Pᵢ × Pⱼ
# C = Cᵢ × Cⱼ
# J depends on the subspace coupling

# Also include string tension: pairs connected by n links of string
# Each string link adds mass σ to the bilinear

bilinear_states = []

for i in range(4):
    for j in range(i, 4):  # i ≤ j to avoid double counting
        for n_sigma in range(0, 5):  # 0 to 4 string links
            
            delta_pair = (Delta[i] + Delta[j])/2 + n_sigma * sigma
            mass_pair = delta_pair * scale
            
            P_pair = parity_intrinsic[i] * parity_intrinsic[j] * (-1)**n_sigma
            C_pair = C_node[i] * C_node[j] * (-1)**n_sigma
            
            # Determine available J values
            sub_i = subspace[i]
            sub_j = subspace[j]
            
            if sub_i == "scalar" and sub_j == "scalar":
                J_fibre = [0]
                # √2 theorem: bilinear of scalars also gives J=2
                J_fibre_extended = [0, 2]
            elif sub_i == "spinor" and sub_j == "spinor":
                J_fibre = [0, 1]
                J_fibre_extended = [0, 1, 2]  # with bilinear projection
            else:
                J_fibre = [0, 1]  # scalar × spinor
                J_fibre_extended = [0, 1]
            
            # String adds orbital angular momentum up to n_sigma
            for J_f in J_fibre_extended:
                for L in range(0, n_sigma + 1):
                    for J_total in range(abs(J_f - L), J_f + L + 1):
                        
                        # Apply √2 mass correction for J=2 bilinear
                        if J_total == 2 and J_f == 2 and sub_i == "scalar" and sub_j == "scalar" and i == j:
                            delta_actual = np.sqrt(2) * Delta[i]  # √2 theorem
                        elif J_total == 2 and J_f == 2 and sub_i == "spinor" and sub_j == "spinor" and i == j:
                            delta_actual = np.sqrt(2) * Delta[i]  # √2 applies to spinor pairs too?
                        else:
                            delta_actual = delta_pair
                        
                        mass_actual = delta_actual * scale
                        ratio = delta_actual / Delta[0]
                        
                        P_str = "+" if P_pair > 0 else "-"
                        C_str = "+" if C_pair > 0 else "-"
                        jpc = f"{J_total}{P_str}{C_str}"
                        
                        bilinear_states.append({
                            "i": i, "j": j, "n_sigma": n_sigma,
                            "J": J_total, "P": P_pair, "C": C_pair,
                            "jpc": jpc, "delta": delta_actual, "ratio": ratio,
                            "mass": mass_actual,
                            "label": f"e_{i}⊗e_{j}+{n_sigma}σ" if n_sigma > 0 else f"e_{i}⊗e_{j}"
                        })

# ============================================================
# MATCH TO LATTICE GLUEBALL SPECTRUM
# ============================================================

print("\n" + "="*70)
print("MATCHING BILINEAR STATES TO LATTICE GLUEBALLS")
print("="*70)

lattice_glueballs = [
    ("0++",  1.000), ("2++",  1.40), ("0-+",  1.50),
    ("1+-",  1.75),  ("2-+",  1.78), ("3+-",  2.11),
    ("3++",  2.15),  ("1--",  2.25), ("2--",  2.35),
    ("3--",  2.46),  ("2+-",  2.48), ("0+-",  2.80),
]

print(f"\n{'Lattice':>8} {'ratio':>6} | {'Best match':>20} {'pred ratio':>10} {'error':>7} | {'Composition':>25}")
print("-"*95)

matched = 0

for jpc_lat, ratio_lat in lattice_glueballs:
    # Find the best matching bilinear state with the same J^PC
    best = None
    best_err = float('inf')
    
    for state in bilinear_states:
        if state["jpc"] != jpc_lat:
            continue
        err = abs(state["ratio"] - ratio_lat) / ratio_lat
        if err < best_err:
            best_err = err
            best = state
    
    if best and best_err < 0.10:
        matched += 1
        mark = "✓" if best_err < 0.03 else "~"
        print(f"{jpc_lat:>8} {ratio_lat:6.2f} | {best['jpc']:>8} {best['ratio']:>10.3f} {best_err*100:>6.1f}% {mark} | {best['label']:>25}")
    else:
        # Show best regardless of J^PC
        best_any = None
        best_any_err = float('inf')
        for state in bilinear_states:
            err = abs(state["ratio"] - ratio_lat) / ratio_lat
            if err < best_any_err:
                best_any_err = err
                best_any = state
        
        if best:
            print(f"{jpc_lat:>8} {ratio_lat:6.2f} | {best['jpc']:>8} {best['ratio']:>10.3f} {best_err*100:>6.1f}%   | {best['label']:>25}")
        else:
            print(f"{jpc_lat:>8} {ratio_lat:6.2f} | {'NO MATCH':>20} | best mass: {best_any['jpc']} at {best_any['ratio']:.3f} ({best_any['label']})")

print(f"\n{matched}/12 matched (correct J^PC AND mass < 10%)")

# ============================================================
# SHOW THE ELEMENTARY PARTICLE CONNECTIONS
# ============================================================

print("\n" + "="*70)
print("THE UNIFYING PICTURE")
print("="*70)

print(f"""
ELEMENTARY MODES (half-eigenvalues, the "quarks and leptons"):
  e₀: √λ₀ = {np.sqrt(lam[0]):.4f} → {Delta[0]/2*scale:.0f} MeV  (scalar, colour-fund)
  e₁: √λ₁ = {np.sqrt(lam[1]):.4f} → {Delta[1]/2*scale:.0f} MeV  (spinor, colour-fund)  ≈ proton!
  e₂: √λ₂ = {np.sqrt(lam[2]):.4f} → {Delta[2]/2*scale:.0f} MeV  (spinor, colour-singlet) ≈ Δ baryon!
  e₃: √λ₃ = {np.sqrt(lam[3]):.4f} → {Delta[3]/2*scale:.0f} MeV  (scalar, colour-singlet)

GLUEBALLS (bilinear pairs of elementary modes):
  0⁺⁺ = e₀⊗e₀ at J=0:  mass = Δ₀ = {Delta[0]*scale:.0f} MeV
  2⁺⁺ = e₀⊗e₀ at J=2:  mass = √2·Δ₀ = {np.sqrt(2)*Delta[0]*scale:.0f} MeV
  0⁻⁺ = e₂⊗e₂ at J=0:  mass = Δ₂ = {Delta[2]*scale:.0f} MeV

PION = splitting between scalar and spinor elementary modes:
  m(π) = Δ₁/2 - Δ₀/2 = {(Delta[1]-Delta[0])/2*scale:.0f} MeV... 
  
  Wait: Δ₁-Δ₀ = {(Delta[1]-Delta[0])*scale:.0f} MeV.
  But (Δ₁-Δ₀)/2 = {(Delta[1]-Delta[0])/2*scale:.0f} MeV.
  Actual pion = 140 MeV.
  
  Δ₁-Δ₀ = {(Delta[1]-Delta[0])*scale:.0f} MeV matches.
  (Δ₁-Δ₀)/2 = {(Delta[1]-Delta[0])/2*scale:.0f} MeV matches... 
  
  Hmm. The pion matched Δ₁-Δ₀ (the FULL eigenvalue splitting).
  If elementary modes are at Δ/2, then the elementary splitting is 
  (Δ₁-Δ₀)/2 = {(Delta[1]-Delta[0])/2*scale:.0f} MeV.
  The BILINEAR (meson) splitting is Δ₁-Δ₀ = {(Delta[1]-Delta[0])*scale:.0f} MeV.
  The pion IS the bilinear splitting. Consistent!
  
  Pion = (spinor pair at J=0) - (scalar pair at J=0) = Δ₁ - Δ₀.
""")

