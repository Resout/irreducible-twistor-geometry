"""
Gauge-invariant Jacobi: does the trivial sector form a closed subalgebra?

The S₃-symmetrized (gauge-invariant) crystals have 4 parameters.
Under normalized composition, they form a 4-parameter family.

Questions:
1. Does the trivial sector close under normalized composition?
2. Does it close under the bracket [A,B] = A∘B - B∘A?
3. Does the Jacobi identity hold within the trivial sector?
4. If so, what Lie algebra is it?

The trivial sector projection: P_triv(M) = (1/6) Σ_{σ∈S₃} σ·M·σ⁻¹
This should commute with composition IF normalization respects symmetry.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, EPS_LOG,
                            schmidt_number, born_fidelity)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


def s3_project(joint):
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    signs = [1, -1, -1, -1, 1, 1]
    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T
    triv = sum(permute(joint, p) for p in perms) / 6
    sign = sum(s * permute(joint, p) for p, s in zip(perms, signs)) / 6
    std = joint - triv - sign
    return triv, sign, std


n_seeds = 500

# Build all S₃ element crystals and their trivial projections
corrs = {
    "e": torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

crystals = {}
for name, corr in corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    crystals[name] = avg


# ═══════════════════════════════════════════════════════════════
#  Trivial projections of all crystals
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("TRIVIAL SECTOR PROJECTIONS")
print("=" * 80)

triv_crystals = {}
for name in crystals:
    triv, _, _ = s3_project(crystals[name])
    triv_crystals[name] = triv

    # How many independent parameters in the trivial projection?
    # For a 4×4 matrix invariant under S₃ (acting on first 3 indices),
    # the structure should be: a on diagonal[0:3], b on off-diagonal[0:3],
    # c on row/col 3, d on (3,3)

    print(f"\n  Trivial projection of crystal({name}):")
    for i in range(MASS_DIM):
        row = ""
        for j in range(MASS_DIM):
            row += f"  {triv[i,j].real.item():>9.6f}"
        print(f"    [{row} ]")

    # Verify it's S₃-invariant
    _, s, st = s3_project(triv)
    print(f"    sign residual: {s.abs().pow(2).sum().sqrt().item():.2e}")
    print(f"    std residual:  {st.abs().pow(2).sum().sqrt().item():.2e}")


# ═══════════════════════════════════════════════════════════════
#  Parameter extraction: the 4 gauge-invariant parameters
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("GAUGE-INVARIANT PARAMETERS (a, b, c, d)")
print("="*80)

print("""
  A 4×4 S₃-invariant matrix (S₃ acts on indices 0,1,2; θ=3 fixed):
    a = diagonal entry M[i,i] for i ∈ {0,1,2}
    b = off-diagonal entry M[i,j] for i≠j ∈ {0,1,2}
    c = θ-coupling M[i,3] = M[3,i] for i ∈ {0,1,2}
    d = θ-self M[3,3]
  4 real parameters (if we restrict to real part).
""")

for name in crystals:
    triv = triv_crystals[name]
    a = triv[0, 0].real.item()
    b = triv[0, 1].real.item()
    c = triv[0, 3].real.item()
    d = triv[3, 3].real.item()
    c2 = triv[3, 0].real.item()  # should equal c for symmetric

    print(f"  crystal({name:>5s}): a={a:>9.6f}, b={b:>9.6f}, c={c:>9.6f}, d={d:>9.6f}  [c_row={c2:.6f}]")


# ═══════════════════════════════════════════════════════════════
#  Closure under normalized composition
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CLOSURE: does compose(P_triv(A), P_triv(B)) = P_triv(compose(A,B))?")
print("="*80)

names5 = ["e", "(01)", "(02)", "(12)", "(012)", "(021)"]

# Test 1: compose trivial projections
print(f"\n--- Method 1: compose(P(A), P(B)) vs P(compose(A,B)) ---")

for na in ["(01)", "(012)"]:
    for nb in ["(02)", "(021)"]:
        # Path 1: compose trivial projections
        comp1 = compose(triv_crystals[na], triv_crystals[nb])
        triv_comp1, _, _ = s3_project(comp1)

        # Path 2: project the composition of full crystals
        comp2 = compose(crystals[na], crystals[nb])
        triv_comp2, _, _ = s3_project(comp2)

        # Path 3: compose full, then project
        diff = (triv_comp1 - triv_comp2).abs().pow(2).sum().sqrt().item()
        ref_norm = triv_comp2.abs().pow(2).sum().sqrt().item()

        # Is comp1 itself S₃-invariant?
        _, s_res, st_res = s3_project(comp1)
        s_frac = s_res.abs().pow(2).sum().item() / comp1.abs().pow(2).sum().item()
        st_frac = st_res.abs().pow(2).sum().item() / comp1.abs().pow(2).sum().item()

        print(f"\n  compose(P({na}), P({nb})):")
        print(f"    ||P(comp(P,P)) - P(comp(full,full))||/||ref|| = {diff/ref_norm:.6f}")
        print(f"    non-trivial leakage: sign={100*s_frac:.2f}%, std={100*st_frac:.2f}%")


# ═══════════════════════════════════════════════════════════════
#  Bracket closure in trivial sector
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("BRACKET IN TRIVIAL SECTOR")
print("="*80)

# The bracket [P(A), P(B)] = compose(P(A), P(B)) - compose(P(B), P(A))
# If both P(A) and P(B) are S₃-invariant, their raw product AB is S₃-invariant
# But the normalized composition may break this

for na, nb in [("(01)","(02)"), ("(01)","(012)"), ("(012)","(021)")]:
    PA = triv_crystals[na]
    PB = triv_crystals[nb]

    AB = compose(PA, PB)
    BA = compose(PB, PA)
    bracket = AB - BA

    br_norm = bracket.abs().pow(2).sum().sqrt().item()

    triv_br, sign_br, std_br = s3_project(bracket)
    nt = triv_br.abs().pow(2).sum().item()
    ns = sign_br.abs().pow(2).sum().item()
    nst = std_br.abs().pow(2).sum().item()
    total = nt + ns + nst

    print(f"\n  [P({na}), P({nb})]:")
    print(f"    ||bracket|| = {br_norm:.8f}")
    if total > 1e-15:
        print(f"    triv={100*nt/total:.2f}%, sign={100*ns/total:.2f}%, std={100*nst/total:.2f}%")
    else:
        print(f"    bracket ≈ 0 (norm < 1e-15)")

    # Raw bracket for comparison
    raw_br = torch.einsum("ab,bc->ac", PA, PB) - torch.einsum("ab,bc->ac", PB, PA)
    raw_norm = raw_br.abs().pow(2).sum().sqrt().item()
    triv_raw, sign_raw, std_raw = s3_project(raw_br)
    nt_r = triv_raw.abs().pow(2).sum().item()
    ns_r = sign_raw.abs().pow(2).sum().item()
    nst_r = std_raw.abs().pow(2).sum().item()
    total_r = nt_r + ns_r + nst_r
    if total_r > 1e-15:
        print(f"    Raw bracket: triv={100*nt_r/total_r:.2f}%, sign={100*ns_r/total_r:.2f}%, std={100*nst_r/total_r:.2f}%")


# ═══════════════════════════════════════════════════════════════
#  Jacobi in trivial sector
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("JACOBI IDENTITY IN TRIVIAL SECTOR")
print("="*80)

PA = triv_crystals["(01)"]
PB = triv_crystals["(02)"]
PC = triv_crystals["(12)"]

def norm(X):
    return X.abs().pow(2).sum().sqrt().item()

# Raw brackets (should satisfy Jacobi exactly)
raw_AB = torch.einsum("ab,bc->ac", PA, PB) - torch.einsum("ab,bc->ac", PB, PA)
raw_BC = torch.einsum("ab,bc->ac", PB, PC) - torch.einsum("ab,bc->ac", PC, PB)
raw_CA = torch.einsum("ab,bc->ac", PC, PA) - torch.einsum("ab,bc->ac", PA, PC)

raw_AB_C = torch.einsum("ab,bc->ac", raw_AB, PC) - torch.einsum("ab,bc->ac", PC, raw_AB)
raw_BC_A = torch.einsum("ab,bc->ac", raw_BC, PA) - torch.einsum("ab,bc->ac", PA, raw_BC)
raw_CA_B = torch.einsum("ab,bc->ac", raw_CA, PB) - torch.einsum("ab,bc->ac", PB, raw_CA)

jacobi_raw = raw_AB_C + raw_BC_A + raw_CA_B
print(f"\n  Raw Jacobi in trivial sector:")
print(f"    ||Jacobi|| = {norm(jacobi_raw):.2e}")
print(f"    ||double brackets|| = {max(norm(raw_AB_C), norm(raw_BC_A), norm(raw_CA_B)):.2e}")
if max(norm(raw_AB_C), norm(raw_BC_A), norm(raw_CA_B)) > 1e-15:
    print(f"    relative = {norm(jacobi_raw)/max(norm(raw_AB_C), norm(raw_BC_A), norm(raw_CA_B)):.2e}")

# Normalized brackets
def bracket_norm(X, Y):
    return compose(X, Y) - compose(Y, X)

norm_AB = bracket_norm(PA, PB)
norm_BC = bracket_norm(PB, PC)
norm_CA = bracket_norm(PC, PA)

norm_AB_C = bracket_norm(norm_AB, PC)
norm_BC_A = bracket_norm(norm_BC, PA)
norm_CA_B = bracket_norm(norm_CA, PB)

jacobi_norm = norm_AB_C + norm_BC_A + norm_CA_B
max_db = max(norm(norm_AB_C), norm(norm_BC_A), norm(norm_CA_B))
print(f"\n  Normalized Jacobi in trivial sector:")
print(f"    ||Jacobi|| = {norm(jacobi_norm):.6f}")
print(f"    ||double brackets|| = {max_db:.6f}")
if max_db > 1e-15:
    print(f"    relative = {norm(jacobi_norm)/max_db:.6f}")

# Compare with full crystal Jacobi
A_full = crystals["(01)"]
B_full = crystals["(02)"]
C_full = crystals["(12)"]

full_AB = bracket_norm(A_full, B_full)
full_BC = bracket_norm(B_full, C_full)
full_CA = bracket_norm(C_full, A_full)

full_AB_C = bracket_norm(full_AB, C_full)
full_BC_A = bracket_norm(full_BC, A_full)
full_CA_B = bracket_norm(full_CA, B_full)

jacobi_full = full_AB_C + full_BC_A + full_CA_B
max_db_full = max(norm(full_AB_C), norm(full_BC_A), norm(full_CA_B))
print(f"\n  Full crystal Jacobi (for comparison):")
print(f"    ||Jacobi|| = {norm(jacobi_full):.6f}")
print(f"    relative = {norm(jacobi_full)/max_db_full:.6f}")

print(f"\n  Trivial/Full Jacobi ratio: {norm(jacobi_norm)/norm(jacobi_full) if norm(jacobi_full) > 1e-10 else 0:.6f}")


# ═══════════════════════════════════════════════════════════════
#  The 4-parameter algebra: explicit structure
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE 4-PARAMETER TRIVIAL ALGEBRA")
print("="*80)

# A basis for the trivial sector:
# B1 = diag entries [0:3], i.e., proportional to diag(1,1,1,0)
# B2 = off-diag [0:3]×[0:3], proportional to [[0,1,1],[1,0,1],[1,1,0],0]
# B3 = coupling entries [i,3] and [3,i], proportional to [[0,0,0,1]...]
# B4 = (3,3) entry, proportional to [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]]

B1 = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for i in range(3):
    B1[i, i] = 1.0
B1 = B1 / 3

B2 = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for i in range(3):
    for j in range(3):
        if i != j:
            B2[i, j] = 1.0
B2 = B2 / 6

B3 = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for i in range(3):
    B3[i, 3] = 1.0
    B3[3, i] = 1.0
B3 = B3 / 6

B4 = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
B4[3, 3] = 1.0

basis = [B1, B2, B3, B4]
basis_names = ["diag", "off-diag", "coupling", "θ-self"]

# Raw multiplication table
print(f"\n  Multiplication table (raw: B_i · B_j expanded in basis):")
for i, ni in enumerate(basis_names):
    for j, nj in enumerate(basis_names):
        prod = torch.einsum("ab,bc->ac", basis[i], basis[j])
        # Project onto basis
        coeffs = []
        for k in range(4):
            # Use Frobenius inner product
            ip = torch.sum(prod.conj() * basis[k]).real.item()
            norm_k = torch.sum(basis[k].conj() * basis[k]).real.item()
            coeffs.append(ip / norm_k if norm_k > 1e-10 else 0)
        # Residual
        reconstructed = sum(c * basis[k] for k, c in enumerate(coeffs))
        res = (prod - reconstructed).abs().pow(2).sum().sqrt().item()

        nonzero = [(basis_names[k], coeffs[k]) for k in range(4) if abs(coeffs[k]) > 1e-6]
        res_str = f" + res={res:.2e}" if res > 1e-6 else ""
        terms = " + ".join(f"{c:.4f}·{n}" for n, c in nonzero) if nonzero else "0"
        print(f"    {ni:>9s} · {nj:<9s} = {terms}{res_str}")


# Raw bracket table
print(f"\n  Bracket table [B_i, B_j]:")
for i in range(4):
    for j in range(i+1, 4):
        br = torch.einsum("ab,bc->ac", basis[i], basis[j]) - \
             torch.einsum("ab,bc->ac", basis[j], basis[i])
        br_norm_val = br.abs().pow(2).sum().sqrt().item()
        if br_norm_val < 1e-10:
            print(f"    [{basis_names[i]}, {basis_names[j]}] = 0")
        else:
            coeffs = []
            for k in range(4):
                ip = torch.sum(br.conj() * basis[k]).real.item()
                norm_k = torch.sum(basis[k].conj() * basis[k]).real.item()
                coeffs.append(ip / norm_k if norm_k > 1e-10 else 0)
            nonzero = [(basis_names[k], coeffs[k]) for k in range(4) if abs(coeffs[k]) > 1e-6]
            reconstructed = sum(c * basis[k] for k, c in enumerate(coeffs))
            res = (br - reconstructed).abs().pow(2).sum().sqrt().item()
            terms = " + ".join(f"{c:.6f}·{n}" for n, c in nonzero) if nonzero else "0"
            res_str = f" + res={res:.2e}" if res > 1e-6 else ""
            print(f"    [{basis_names[i]}, {basis_names[j]}] = {terms}{res_str}")


print(f"\n\n{'='*80}")
print("WHAT THE GAUGE-INVARIANT JACOBI REVEALS")
print("="*80)
