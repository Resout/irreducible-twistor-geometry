"""
Thread 98: Standard sector 5-doublet structure → quarks?

The 10D standard sector has 5 degenerate doublets with a 1+3+1 band.
Exponents are symmetric around a₀ = 1/3. Pairs sum to 2/3.

Questions:
1. Does the 1+3+1 structure persist across ALL conjugacy classes?
2. What determines the band width (0.08 in exponent space)?
3. Do the standard sector PHASES carry CKM information?
4. Does the leading standard eigenvalue differ between classes
   in a way that matches quark mass ratios?
5. The 1+3+1 band: is it rank(A₂) + |Φ⁺(A₂)| = 2 + 3?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import numpy as np
from solver.algebra import H, MASS_DIM, K_STAR, BORN_FLOOR, DELTA
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

dim_sq = MASS_DIM ** 2

S3_PERMS = {
    "e":     [0, 1, 2],
    "(01)":  [1, 0, 2],
    "(02)":  [2, 1, 0],
    "(12)":  [0, 2, 1],
    "(012)": [1, 2, 0],
    "(021)": [2, 0, 1],
}

S3_CORR = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

CONJUGACY_CLASSES = {
    "identity": ["e"],
    "3-cycles": ["(012)", "(021)"],
    "transpositions": ["(01)", "(02)", "(12)"],
}

def perm_mat_4x4(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return P

reps = {n: torch.kron(perm_mat_4x4(p), perm_mat_4x4(p)) for n, p in S3_PERMS.items()}

chi_d = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}
chi_s = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}
chi_t = {"e": 1, "(01)": 1, "(02)": 1, "(12)": 1, "(012)": 1, "(021)": 1}

def proj(chi, dim_chi):
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for n in S3_PERMS:
        P += chi[n] * reps[n]
    return P * dim_chi / 6.0

P_std = proj(chi_d, 2)
P_sign = proj(chi_s, 1)
P_triv = proj(chi_t, 1)

def build_T_B(B_joint):
    T = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for col in range(dim_sq):
        A = torch.zeros(dim_sq, dtype=torch.cfloat)
        A[col] = 1.0
        result = compose(A.reshape(MASS_DIM, MASS_DIM), B_joint)
        T[:, col] = result.reshape(dim_sq)
    return T

N_SEEDS = 200


# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("PART 1: FULL STANDARD SECTOR SPECTRUM FOR ALL CLASSES")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

class_std_spectra = {}

for class_name, members in CONJUGACY_CLASSES.items():
    all_evals = []  # list of sorted magnitude arrays

    for name in members:
        for seed in range(N_SEEDS):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T = build_T_B(ent.joint)
            restricted = P_std @ T @ P_std
            evals = torch.linalg.eigvals(restricted)
            mags = evals.abs().sort(descending=True).values
            # Take top 10 (the standard sector should have 10 significant eigenvalues)
            above = mags[mags > 1e-6]
            all_evals.append(above[:10].tolist())

    # Average
    max_len = max(len(s) for s in all_evals)
    avg_spec = []
    std_spec = []
    for i in range(min(max_len, 10)):
        vals = [s[i] for s in all_evals if len(s) > i]
        avg_spec.append(np.mean(vals))
        std_spec.append(np.std(vals) / np.sqrt(len(vals)))

    class_std_spectra[class_name] = avg_spec

    print(f"\n  {class_name} ({len(all_evals)} measurements):")
    print(f"  {'idx':>4s}  {'|λ|':>10s}  {'±':>8s}  {'exponent':>10s}  {'pair sum':>10s}")
    print(f"  {'-'*50}")

    for i, (val, err) in enumerate(zip(avg_spec, std_spec)):
        a = -math.log(val) / math.log(6) if val > 1e-20 else float('inf')
        # Pair sum with symmetric partner
        j = len(avg_spec) - 1 - i
        if j > i and j < len(avg_spec):
            pair_sum = a + (-math.log(avg_spec[j]) / math.log(6) if avg_spec[j] > 1e-20 else float('inf'))
            print(f"  {i:>4d}  {val:>10.6f}  {err:>8.6f}  {a:>10.6f}  {pair_sum:>10.6f}")
        else:
            print(f"  {i:>4d}  {val:>10.6f}  {err:>8.6f}  {a:>10.6f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 2: DOUBLET STRUCTURE AND BAND COMPARISON")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print("\n  Doublet exponents by class (pairing consecutive eigenvalues):")
print(f"  {'Doublet':>8s}  {'identity':>12s}  {'3-cycles':>12s}  {'transpos.':>12s}")
print(f"  {'-'*50}")

for d in range(5):
    vals = []
    for cn in ["identity", "3-cycles", "transpositions"]:
        spec = class_std_spectra[cn]
        if 2*d + 1 < len(spec):
            avg = (spec[2*d] + spec[2*d+1]) / 2
            a = -math.log(avg) / math.log(6)
            vals.append(a)
        else:
            vals.append(float('nan'))
    print(f"  {d:>8d}  {vals[0]:>12.6f}  {vals[1]:>12.6f}  {vals[2]:>12.6f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 3: CLASS DIFFERENCES IN STANDARD EIGENVALUES")
print("Do the standard eigenvalues shift between classes like quark masses?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# For each of the 5 doublets, compare the eigenvalue across classes
print("\n  Doublet eigenvalue shifts (ratio to identity):")
print(f"  {'Doublet':>8s}  {'3c/id':>10s}  {'tr/id':>10s}  {'3c−id(exp)':>12s}  {'tr−id(exp)':>12s}")
print(f"  {'-'*60}")

id_spec = class_std_spectra["identity"]
cyc_spec = class_std_spectra["3-cycles"]
tr_spec = class_std_spectra["transpositions"]

for d in range(5):
    if 2*d + 1 < min(len(id_spec), len(cyc_spec), len(tr_spec)):
        id_avg = (id_spec[2*d] + id_spec[2*d+1]) / 2
        cyc_avg = (cyc_spec[2*d] + cyc_spec[2*d+1]) / 2
        tr_avg = (tr_spec[2*d] + tr_spec[2*d+1]) / 2

        r_cyc = cyc_avg / id_avg
        r_tr = tr_avg / id_avg

        a_id = -math.log(id_avg) / math.log(6)
        a_cyc = -math.log(cyc_avg) / math.log(6)
        a_tr = -math.log(tr_avg) / math.log(6)

        print(f"  {d:>8d}  {r_cyc:>10.6f}  {r_tr:>10.6f}  {a_cyc-a_id:>12.6f}  {a_tr-a_id:>12.6f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 4: STANDARD SECTOR COMPLEX EIGENVALUES")
print("Extract phases, not just magnitudes")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# For the identity crystal, get the FULL complex eigenvalues of P_std @ T @ P_std
print("\n  Identity crystal: standard sector full complex eigenvalues (100 seeds):")

all_complex_evals = []
for seed in range(100):
    ent = Entangler(S3_CORR["e"], seed=seed).build()
    T = build_T_B(ent.joint)
    restricted = P_std @ T @ P_std
    evals = torch.linalg.eigvals(restricted)

    # Sort by magnitude
    mags = evals.abs()
    idx = mags.argsort(descending=True)
    sorted_evals = evals[idx][:10]  # top 10
    all_complex_evals.append(sorted_evals)

# Average the complex eigenvalues (accounting for phase alignment)
# Actually, average magnitudes and phases separately
avg_mags = []
avg_phases = []
for i in range(10):
    vals = [s[i] for s in all_complex_evals if len(s) > i]
    mags = [v.abs().item() for v in vals]
    phases = [math.atan2(v.imag.item(), v.real.item()) for v in vals]
    avg_mags.append(np.mean(mags))
    avg_phases.append(np.mean(phases))

print(f"\n  {'idx':>4s}  {'|λ|':>10s}  {'φ':>10s}  {'φ/π':>8s}  {'φ/(2π/3)':>10s}")
print(f"  {'-'*50}")
for i, (m, p) in enumerate(zip(avg_mags, avg_phases)):
    print(f"  {i:>4d}  {m:>10.6f}  {p:>10.4f}  {p/math.pi:>8.4f}  {p/(2*math.pi/3):>10.4f}")

# Do the phases show 2π/3 structure?
print(f"\n  Phase differences (consecutive):")
for i in range(min(len(avg_phases)-1, 9)):
    diff = avg_phases[i+1] - avg_phases[i]
    print(f"    φ_{i+1} − φ_{i} = {diff:>8.4f} rad = {diff/math.pi:>8.4f}π")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 5: SIGN × STANDARD CROSS-COUPLING (YUKAWA)")
print("The off-diagonal sector coupling sign→standard")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# The Yukawa coupling connects the sign sector to the standard sector.
# Compute: P_std @ T @ P_sign  (sign → standard transition)

print("\n  Sign → Standard coupling (per conjugacy class, 100 seeds):")

for class_name, members in CONJUGACY_CLASSES.items():
    all_couplings = []
    for name in members:
        for seed in range(100):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T = build_T_B(ent.joint)
            cross = P_std @ T @ P_sign
            # Frobenius norm of the cross-coupling
            fnorm = cross.abs().pow(2).sum().sqrt().item()
            all_couplings.append(fnorm)

    avg_c = np.mean(all_couplings)
    print(f"    {class_name:>15s}: ||P_std T P_sign|| = {avg_c:.6f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 6: TRIVIAL SECTOR STRUCTURE")
print("5 trivial eigenvalues — do they pair with standard?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print("\n  Identity crystal: trivial sector eigenvalues (200 seeds):")

all_triv = []
for seed in range(200):
    ent = Entangler(S3_CORR["e"], seed=seed).build()
    T = build_T_B(ent.joint)
    restricted = P_triv @ T @ P_triv
    evals = torch.linalg.eigvals(restricted)
    mags = evals.abs().sort(descending=True).values
    above = mags[mags > 1e-6]
    all_triv.append(above[:5].tolist())

max_len = max(len(s) for s in all_triv)
avg_triv = []
for i in range(min(max_len, 5)):
    vals = [s[i] for s in all_triv if len(s) > i]
    avg_triv.append(np.mean(vals))

print(f"  {'idx':>4s}  {'|λ_triv|':>12s}  {'exponent':>10s}")
print(f"  {'-'*30}")
for i, val in enumerate(avg_triv):
    a = -math.log(val) / math.log(6)
    print(f"  {i:>4d}  {val:>12.6f}  {a:>10.6f}")

# Compare trivial spectrum to standard doublet centroids
print(f"\n  Comparison: trivial vs standard doublet centroids:")
std_centroids = []
spec = class_std_spectra["identity"]
for d in range(5):
    if 2*d+1 < len(spec):
        std_centroids.append((spec[2*d] + spec[2*d+1]) / 2)

for i in range(min(len(avg_triv), len(std_centroids))):
    a_t = -math.log(avg_triv[i]) / math.log(6)
    a_s = -math.log(std_centroids[i]) / math.log(6)
    print(f"    Doublet {i}: triv={avg_triv[i]:.6f}(a={a_t:.4f}), std={std_centroids[i]:.6f}(a={a_s:.4f}), ratio={avg_triv[i]/std_centroids[i]:.4f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 7: THE BAND WIDTH — WHAT SETS THE 0.08 SPREAD?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

spec = class_std_spectra["identity"]
if len(spec) >= 10:
    # Doublet exponents
    doublet_a = []
    for d in range(5):
        avg = (spec[2*d] + spec[2*d+1]) / 2
        doublet_a.append(-math.log(avg) / math.log(6))

    band = doublet_a[-1] - doublet_a[0]
    mid_spread = doublet_a[3] - doublet_a[1]
    gap_low = doublet_a[1] - doublet_a[0]
    gap_high = doublet_a[4] - doublet_a[3]
    mid_center = doublet_a[2]

    print(f"\n  Band parameters (identity crystal):")
    print(f"    Total band: {band:.6f}  (2/3 = {2/3:.6f}, dev = {abs(band-2/3)/(2/3)*100:.2f}%)")
    print(f"    Middle spread: {mid_spread:.6f}")
    print(f"    Gap low→mid: {gap_low:.6f}")
    print(f"    Gap mid→high: {gap_high:.6f}")
    print(f"    Middle center: {mid_center:.6f}  (a₀ = {1/3:.6f})")
    print(f"    Symmetry: low gap/high gap = {gap_low/gap_high:.4f}  (1 = perfect)")

    # Middle spread candidates
    print(f"\n  Middle spread {mid_spread:.6f} candidates:")
    for name, val in [("Δ/H", DELTA/H), ("K*/H", K_STAR/H),
                       ("1/h(E₈)", 1/30), ("BORN_FLOOR", BORN_FLOOR),
                       ("Δ/3", DELTA/3), ("(a₂-a₀)/h(E₈)", (59/30-1/3)/30),
                       ("1/(H²+1)", 1/(H**2+1)), ("2Δ/H²", 2*DELTA/H**2),
                       ("K*/3", K_STAR/3), ("1/12", 1/12),
                       ("Δ²", DELTA**2), ("K*²", K_STAR**2)]:
        dev = abs(val - mid_spread) / mid_spread * 100
        if dev < 30:
            print(f"      {name:>15s} = {val:.6f} ({dev:.1f}% off)")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 8: INDIVIDUAL TRANSPOSITION STANDARD SPECTRA")
print("Do the THREE transpositions give different standard spectra?")
print("(If so, this could be the quark mass splitting)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

for name in ["(01)", "(02)", "(12)"]:
    all_evals = []
    for seed in range(N_SEEDS):
        ent = Entangler(S3_CORR[name], seed=seed).build()
        T = build_T_B(ent.joint)
        restricted = P_std @ T @ P_std
        evals = torch.linalg.eigvals(restricted)
        mags = evals.abs().sort(descending=True).values
        above = mags[mags > 1e-6]
        all_evals.append(above[:10].tolist())

    avg = []
    for i in range(10):
        vals = [s[i] for s in all_evals if len(s) > i]
        avg.append(np.mean(vals))

    # Show doublet centroids
    doublets = []
    for d in range(5):
        if 2*d+1 < len(avg):
            doublets.append((avg[2*d] + avg[2*d+1]) / 2)

    print(f"\n  {name}: doublets = [{', '.join(f'{d:.6f}' for d in doublets)}]")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("SUMMARY")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════
print(f"""
  Results above. Looking for:
  1. Universal 1+3+1 band across all classes
  2. Class-dependent shifts → quark masses
  3. Phase structure → CKM
  4. Sign×standard coupling → Yukawa
""")
