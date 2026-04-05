"""
Thread 98b: The Yukawa coupling from crystal sector structure.

The sign→standard coupling is 8.4× stronger for transpositions
than for identity/3-cycles. What is this ratio exactly?

Also: the standard→sign coupling (reverse direction), and the
trivial sector couplings. Full inter-sector coupling matrix.
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

chi_s = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}
chi_t = {"e": 1, "(01)": 1, "(02)": 1, "(12)": 1, "(012)": 1, "(021)": 1}
chi_d = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}

def proj(chi, dim_chi):
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for n in S3_PERMS:
        P += chi[n] * reps[n]
    return P * dim_chi / 6.0

P_triv = proj(chi_t, 1)
P_sign = proj(chi_s, 1)
P_std = proj(chi_d, 2)

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
print("FULL INTER-SECTOR COUPLING MATRIX")
print("3×3 matrix: ||P_α T_B P_β||_F for α,β ∈ {triv, sign, std}")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

sectors = {"triv": P_triv, "sign": P_sign, "std": P_std}

for class_name, members in CONJUGACY_CLASSES.items():
    print(f"\n  {class_name} ({len(members)} elements, {N_SEEDS} seeds each):")

    coupling_matrix = {}
    for alpha in ["triv", "sign", "std"]:
        for beta in ["triv", "sign", "std"]:
            all_norms = []
            for name in members:
                for seed in range(N_SEEDS):
                    ent = Entangler(S3_CORR[name], seed=seed).build()
                    T = build_T_B(ent.joint)
                    cross = sectors[alpha] @ T @ sectors[beta]
                    fnorm = cross.abs().pow(2).sum().sqrt().item()
                    all_norms.append(fnorm)
            coupling_matrix[(alpha, beta)] = np.mean(all_norms)

    # Print as matrix
    print(f"  {'→':>8s}", end="")
    for beta in ["triv", "sign", "std"]:
        print(f"  {beta:>10s}", end="")
    print()
    for alpha in ["triv", "sign", "std"]:
        print(f"  {alpha:>8s}", end="")
        for beta in ["triv", "sign", "std"]:
            print(f"  {coupling_matrix[(alpha,beta)]:>10.6f}", end="")
        print()


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("YUKAWA RATIO ANALYSIS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Re-measure with higher precision for just sign→std
print("\n  High-precision sign→standard coupling:")

class_yukawa = {}
for class_name, members in CONJUGACY_CLASSES.items():
    all_norms = []
    for name in members:
        for seed in range(N_SEEDS):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T = build_T_B(ent.joint)
            cross = P_std @ T @ P_sign
            fnorm = cross.abs().pow(2).sum().sqrt().item()
            all_norms.append(fnorm)

    mean = np.mean(all_norms)
    sem = np.std(all_norms) / np.sqrt(len(all_norms))
    class_yukawa[class_name] = (mean, sem)
    print(f"    {class_name:>15s}: {mean:.8f} ± {sem:.8f}")

y_id = class_yukawa["identity"][0]
y_3c = class_yukawa["3-cycles"][0]
y_tr = class_yukawa["transpositions"][0]

ratio = y_tr / y_id
print(f"\n  Ratio transposition/identity: {ratio:.6f}")
print(f"  Ratio transposition/3-cycles: {y_tr / y_3c:.6f}")

# Crystal expressions for the ratio
print(f"\n  Crystal expression search for ratio {ratio:.6f}:")
S = H**3 / K_STAR  # 810/7

for name, val in [("H²", H**2),
                   ("2^H", 2**H),
                   ("H!", math.factorial(H)),
                   ("H(H+1)/2", H*(H+1)/2),
                   ("1/K*", 1/K_STAR),
                   ("H/K*", H/K_STAR),
                   ("dim_std/dim_sign", 10/1),
                   ("dim_std/dim_triv", 10/5),
                   ("√dim_std", math.sqrt(10)),
                   ("(1-K*)/K*", (1-K_STAR)/K_STAR),
                   ("H/(1-K*)", H/(1-K_STAR)),
                   ("1/BORN_FLOOR", 1/BORN_FLOOR),
                   ("H³", H**3),
                   ("√(H³/K*)", math.sqrt(S)),
                   ("H²/(1-K*)", H**2/(1-K_STAR)),
                   ("(H²+1)/(H-1)", (H**2+1)/(H-1)),
                   ("dim_std-dim_sign", 10-1),
                   ("H² + K*_num", H**2 + 7),
                   ("h(E₈)/H", 30/H),
                   ("(H+1)²/H", (H+1)**2/H)]:
    dev = abs(val - ratio) / ratio * 100
    if dev < 20:
        print(f"    {name:>20s} = {val:>10.4f}  ({dev:.2f}% off)")

# The sign-character weighted version
# For transpositions χ_sign = -1, so the sign→standard coupling
# through the sign sector is NEGATED. The Frobenius norm doesn't
# see the sign, but the coupling strength might be |χ_sign| × something.
# For identity/3-cycles: χ = +1, coupling ≈ 0.066
# For transpositions: χ = -1, coupling ≈ 0.554
# The transposition coupling is |χ| × base, but the base is different!

# Maybe: transposition coupling = identity_coupling / |λ_sign(trans)|?
lambda_sign_trans = 0.02944  # from earlier measurement
lambda_sign_id = 0.55008
print(f"\n  Identity Yukawa / |λ_sign(id)| = {y_id / lambda_sign_id:.6f}")
print(f"  Transposition Yukawa / |λ_sign(tr)| = {y_tr / lambda_sign_trans:.6f}")
print(f"  Ratio: {(y_tr/lambda_sign_trans) / (y_id/lambda_sign_id):.6f}")

# Or: the Yukawa coupling IS the sign coupling expressed in the standard basis
# For the identity: weak coupling (sign and standard are "aligned")
# For transpositions: strong coupling (sign and standard "mix")
print(f"\n  Yukawa ≈ |χ_sign| × (dim_std/H!) × |λ_sign(class)|^{-1}?")
for cn, yval in class_yukawa.items():
    lam = {"identity": lambda_sign_id, "3-cycles": 0.27665, "transpositions": lambda_sign_trans}[cn]
    chi_val = {"identity": 1, "3-cycles": 1, "transpositions": -1}[cn]
    pred = abs(chi_val) * (10 / math.factorial(H)) * (1 / lam)
    print(f"    {cn:>15s}: pred = {pred:.6f}, measured = {yval[0]:.6f}, ratio = {yval[0]/pred:.4f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("DIAGONAL SECTOR COUPLINGS — THE EIGENVALUE STRUCTURE")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Diagonal sector couplings are the eigenvalues within each sector.
# Compare the LEADING eigenvalue of each sector across classes.

print("\n  Leading eigenvalues per sector per class:")
print(f"  {'class':>15s}  {'|λ_triv|':>10s}  {'|λ_sign|':>10s}  {'|λ_std|':>10s}  {'t/s ratio':>10s}")
print(f"  {'-'*60}")

for class_name, members in CONJUGACY_CLASSES.items():
    triv_vals = []
    sign_vals = []
    std_vals = []
    for name in members:
        for seed in range(100):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T = build_T_B(ent.joint)

            # Trivial
            r = P_triv @ T @ P_triv
            evals = torch.linalg.eigvals(r)
            triv_vals.append(evals.abs().max().item())

            # Sign
            r = P_sign @ T @ P_sign
            evals = torch.linalg.eigvals(r)
            sign_vals.append(evals.abs().max().item())

            # Standard
            r = P_std @ T @ P_std
            evals = torch.linalg.eigvals(r)
            std_vals.append(evals.abs().max().item())

    t = np.mean(triv_vals)
    s = np.mean(sign_vals)
    st = np.mean(std_vals)
    print(f"  {class_name:>15s}  {t:>10.6f}  {s:>10.6f}  {st:>10.6f}  {t/s:>10.4f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("THE 8:1 RATIO — IS IT H³/K*_num = 27/7 ≈ 3.86?")
print("Or dim_std/dim_sign = 10? Or something else?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# The measured ratio is y_tr/y_id
print(f"\n  Measured Yukawa ratio: {ratio:.6f}")
print(f"\n  Key crystal quantities near {ratio:.2f}:")

candidates = sorted([
    ("H²", H**2),
    ("2^H", 2**H),
    ("H!", math.factorial(H)),
    ("H(H+1)/2", H*(H+1)//2),
    ("1/K*", 1/K_STAR),
    ("dim_std", 10),
    ("dim_std/dim_sign", 10),
    ("dim_std/dim_triv", 2),
    ("(1-K*)/K*", (1-K_STAR)/K_STAR),
    ("H/K*_num", H/7),
    ("h(E₈)/H²", 30/9),
    ("23/H", 23/3),
    ("H²+1", H**2+1),
    ("H²-1", H**2-1),
    ("S/h(E₈)", (810/7)/30),
    ("2H+1", 2*H+1),
    ("2H+2", 2*H+2),
    ("3H-1", 3*H-1),
    ("3H", 3*H),
    ("H×(H-1)+2", H*(H-1)+2),
], key=lambda x: abs(x[1] - ratio))

for name, val in candidates[:10]:
    dev = abs(val - ratio) / ratio * 100
    print(f"    {name:>20s} = {val:>10.4f}  ({dev:.2f}% off)")


print(f"\n\n{'='*80}")
print("SUMMARY")
print("=" * 80)
print(f"""
  The Yukawa coupling (sign→standard) is:
    identity/3-cycles: ~0.066
    transpositions:    ~{y_tr:.3f}

  Ratio: {ratio:.4f}

  This is the crystal's Yukawa mechanism: the gauge sector (transpositions)
  couples sign to standard {ratio:.1f}× more strongly than the gravity/chirality
  sectors. The sign representation character χ(σ) = -1 activates the mixing.
""")
