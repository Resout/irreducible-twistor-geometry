"""
eigenvalue_spectrum.py — Eigenvalue structure of [4,4] joint mass crystals.

Investigates:
1. Eigenvalues of each canonical relationship signature crystal
2. Eigenvalue evolution under iterated self-composition
3. Spectral radius vs Schmidt number
4. Eigenvalue ratios and condition number under composition
5. Normality: A*A† vs A†*A
6. Eigenvalues of the composition operator T_B (16×16 linear map)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch
import math
from solver.algebra import (
    H, MASS_DIM, schmidt_number, DELTA, K_STAR,
    EPS_LOG, EPS_DIV,
)
from solver.crystals import (
    Entangler, compose, RELATIONSHIP_SIGNATURES,
)


def eigenvalues(M):
    """Eigenvalues of a [4,4] complex matrix, sorted by descending |λ|."""
    vals, _ = torch.linalg.eig(M)
    order = vals.abs().argsort(descending=True)
    return vals[order]


def singular_values(M):
    """Singular values of a [4,4] complex matrix, descending."""
    _, S, _ = torch.linalg.svd(M)
    return S


def spectral_radius(M):
    """Largest |eigenvalue|."""
    vals = eigenvalues(M)
    return vals[0].abs().item()


def build_all_signatures():
    """Build crystal joints for all 6 relationship signatures."""
    crystals = {}
    for name, corr in RELATIONSHIP_SIGNATURES.items():
        ent = Entangler(corr, seed=42).build()
        crystals[name] = ent.joint.clone()
    return crystals


# ═══════════════════════════════════════════════════════════════
#  1. Eigenvalues of each signature crystal
# ═══════════════════════════════════════════════════════════════

def investigate_signature_eigenvalues(crystals):
    print("=" * 72)
    print("  1. EIGENVALUES OF CANONICAL RELATIONSHIP SIGNATURE CRYSTALS")
    print("=" * 72)

    for name, joint in crystals.items():
        eigs = eigenvalues(joint)
        svs = singular_values(joint)
        sn = schmidt_number(joint)
        det_val = torch.linalg.det(joint)

        print(f"\n  {name.upper()}")
        print(f"    Schmidt number: {sn:.4f}")
        print(f"    det = {det_val.item():.6e}")
        print(f"    |det| = {abs(det_val.item()):.6e}")
        print(f"    det phase = {math.atan2(det_val.imag, det_val.real)/math.pi:.4f}π")

        print(f"    Eigenvalues (sorted by |λ|):")
        for i, e in enumerate(eigs):
            print(f"      λ_{i} = {e.real:.6f} + {e.imag:.6f}i"
                  f"   |λ_{i}| = {e.abs():.6f}"
                  f"   phase = {math.atan2(e.imag, e.real)/math.pi:.4f}π")

        print(f"    Singular values:")
        for i, s in enumerate(svs):
            print(f"      σ_{i} = {s:.6f}")

        # Check if eigenvalue magnitudes match singular values
        eig_mags = sorted([e.abs().item() for e in eigs], reverse=True)
        sv_list = svs.tolist()
        max_diff = max(abs(eig_mags[i] - sv_list[i]) for i in range(4))
        print(f"    Max |  |λ_i| - σ_i  | = {max_diff:.6f}"
              f"  {'(NORMAL)' if max_diff < 0.01 else '(NOT NORMAL)'}")

        # Eigenvalue product = det
        eig_product = torch.prod(eigs)
        print(f"    Product of eigenvalues: {eig_product.item():.6e}")


# ═══════════════════════════════════════════════════════════════
#  2. Eigenvalue evolution under iterated self-composition
# ═══════════════════════════════════════════════════════════════

def investigate_composition_eigenvalues(crystals):
    print("\n" + "=" * 72)
    print("  2. EIGENVALUE EVOLUTION UNDER ITERATED SELF-COMPOSITION")
    print("=" * 72)

    # Use proportional (identity) crystal
    joint = crystals["proportional"]
    N_ITER = 12

    print(f"\n  Starting with PROPORTIONAL crystal, composing {N_ITER} times:")
    print(f"  {'Step':>4}  {'|λ₀|':>10}  {'|λ₁|':>10}  {'|λ₂|':>10}  {'|λ₃|':>10}"
          f"  {'|det|':>12}  {'Schmidt':>8}  {'ρ(A)':>8}")

    current = joint.clone()
    spectral_radii = []
    det_magnitudes = []
    eig0_magnitudes = []

    for step in range(N_ITER + 1):
        eigs = eigenvalues(current)
        det_val = torch.linalg.det(current)
        sn = schmidt_number(current)
        rho = eigs[0].abs().item()

        mags = [e.abs().item() for e in eigs]
        spectral_radii.append(rho)
        det_magnitudes.append(abs(det_val.item()))
        eig0_magnitudes.append(mags[0])

        print(f"  {step:>4}  {mags[0]:>10.6f}  {mags[1]:>10.6f}  "
              f"{mags[2]:>10.6f}  {mags[3]:>10.6f}  "
              f"{abs(det_val.item()):>12.4e}  {sn:>8.4f}  {rho:>8.6f}")

        if step < N_ITER:
            current = compose(current, joint)

    # Decay rate analysis
    print(f"\n  Decay rate analysis:")
    for i in range(1, min(N_ITER + 1, len(eig0_magnitudes))):
        if eig0_magnitudes[i] > 1e-15 and eig0_magnitudes[i-1] > 1e-15:
            ratio = eig0_magnitudes[i] / eig0_magnitudes[i-1]
            log_ratio = -math.log(ratio) if ratio > 0 else float('inf')
            print(f"    Step {i-1}→{i}: |λ₀| ratio = {ratio:.6f},"
                  f"  -ln(ratio) = {log_ratio:.6f}")

    if len(det_magnitudes) > 2 and det_magnitudes[1] > 1e-30:
        det_ratio = det_magnitudes[2] / det_magnitudes[1] if det_magnitudes[1] > 1e-30 else 0
        print(f"\n  |det| contraction ratio (step 1→2): {det_ratio:.6f}")

    print(f"\n  Reference: Δ = {DELTA:.6f}, K* = {K_STAR:.6f}")

    # Also try with inverse crystal
    print(f"\n  --- Now with INVERSE crystal ---")
    joint_inv = crystals["inverse"]
    current = joint_inv.clone()
    print(f"  {'Step':>4}  {'|λ₀|':>10}  {'|λ₁|':>10}  {'|λ₂|':>10}  {'|λ₃|':>10}")

    for step in range(N_ITER + 1):
        eigs = eigenvalues(current)
        mags = [e.abs().item() for e in eigs]
        print(f"  {step:>4}  {mags[0]:>10.6f}  {mags[1]:>10.6f}  "
              f"{mags[2]:>10.6f}  {mags[3]:>10.6f}")
        if step < N_ITER:
            current = compose(current, joint_inv)


# ═══════════════════════════════════════════════════════════════
#  3. Spectral radius vs Schmidt number
# ═══════════════════════════════════════════════════════════════

def investigate_spectral_vs_schmidt(crystals):
    print("\n" + "=" * 72)
    print("  3. SPECTRAL RADIUS vs SCHMIDT NUMBER")
    print("=" * 72)

    # Build many crystals with different correlations
    test_correlations = {
        **{name: corr for name, corr in RELATIONSHIP_SIGNATURES.items()},
        "strong_diag": torch.tensor([[2, 0, 0], [0, 2, 0], [0, 0, 2]], dtype=torch.float32),
        "off_diag": torch.tensor([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=torch.float32),
        "uniform": torch.tensor([[1, 1, 1], [1, 1, 1], [1, 1, 1]], dtype=torch.float32),
        "corner": torch.tensor([[3, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=torch.float32),
        "L_shape": torch.tensor([[1, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=torch.float32),
        "checkerboard": torch.tensor([[1, 0, 1], [0, 1, 0], [1, 0, 1]], dtype=torch.float32),
        "band": torch.tensor([[1, .5, 0], [.5, 1, .5], [0, .5, 1]], dtype=torch.float32),
        "skew": torch.tensor([[0, 0, 1], [0, 0, 0], [0, 0, 0]], dtype=torch.float32),
        "triangular": torch.tensor([[1, 1, 1], [0, 1, 1], [0, 0, 1]], dtype=torch.float32),
    }

    # Also compose pairs to get more variety
    print(f"\n  {'Name':>20}  {'ρ(A)':>10}  {'Schmidt':>10}  {'|det|':>12}  {'κ(A)':>10}")
    print(f"  {'-'*20}  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*10}")

    results = []
    for name, corr in test_correlations.items():
        ent = Entangler(corr, seed=42).build()
        joint = ent.joint
        eigs = eigenvalues(joint)
        sn = schmidt_number(joint)
        rho = eigs[0].abs().item()
        det_val = abs(torch.linalg.det(joint).item())
        cond = eigs[0].abs().item() / max(eigs[-1].abs().item(), 1e-12)

        results.append((name, rho, sn, det_val, cond))
        print(f"  {name:>20}  {rho:>10.6f}  {sn:>10.4f}  {det_val:>12.4e}  {cond:>10.2f}")

    # Also some composed crystals
    print(f"\n  Composed crystals:")
    composed_pairs = [
        ("prop∘inv", "proportional", "inverse"),
        ("prop∘mod", "proportional", "modular"),
        ("inv∘exp", "inverse", "exponential"),
        ("exp∘log", "exponential", "logarithmic"),
        ("quad∘quad", "quadratic", "quadratic"),
    ]
    for label, n1, n2 in composed_pairs:
        j1 = Entangler(RELATIONSHIP_SIGNATURES[n1], seed=42).build().joint
        j2 = Entangler(RELATIONSHIP_SIGNATURES[n2], seed=42).build().joint
        comp = compose(j1, j2)
        eigs = eigenvalues(comp)
        sn = schmidt_number(comp)
        rho = eigs[0].abs().item()
        det_val = abs(torch.linalg.det(comp).item())
        cond = eigs[0].abs().item() / max(eigs[-1].abs().item(), 1e-12)
        results.append((label, rho, sn, det_val, cond))
        print(f"  {label:>20}  {rho:>10.6f}  {sn:>10.4f}  {det_val:>12.4e}  {cond:>10.2f}")

    # Correlation analysis
    rhos = [r[1] for r in results]
    schmidts = [r[2] for r in results]
    n = len(rhos)
    mean_r = sum(rhos) / n
    mean_s = sum(schmidts) / n
    cov = sum((rhos[i] - mean_r) * (schmidts[i] - mean_s) for i in range(n)) / n
    std_r = (sum((r - mean_r) ** 2 for r in rhos) / n) ** 0.5
    std_s = (sum((s - mean_s) ** 2 for s in schmidts) / n) ** 0.5
    corr = cov / (std_r * std_s) if std_r * std_s > 1e-12 else 0
    print(f"\n  Pearson correlation(ρ, Schmidt) = {corr:.4f}")


# ═══════════════════════════════════════════════════════════════
#  4. Eigenvalue ratios under composition
# ═══════════════════════════════════════════════════════════════

def investigate_eigenvalue_ratios(crystals):
    print("\n" + "=" * 72)
    print("  4. EIGENVALUE RATIOS UNDER COMPOSITION")
    print("=" * 72)

    N_ITER = 12
    for name in ["proportional", "inverse", "modular", "exponential"]:
        joint = crystals[name]
        current = joint.clone()

        print(f"\n  {name.upper()} crystal:")
        print(f"  {'Step':>4}  {'λ₀/λ₁':>10}  {'λ₀/λ₂':>10}  {'λ₀/λ₃':>10}  {'λ₁/λ₂':>10}")

        for step in range(N_ITER + 1):
            eigs = eigenvalues(current)
            mags = [e.abs().item() for e in eigs]

            r01 = mags[0] / max(mags[1], 1e-15)
            r02 = mags[0] / max(mags[2], 1e-15)
            r03 = mags[0] / max(mags[3], 1e-15)
            r12 = mags[1] / max(mags[2], 1e-15)

            print(f"  {step:>4}  {r01:>10.4f}  {r02:>10.4f}  {r03:>10.4f}  {r12:>10.4f}")

            if step < N_ITER:
                current = compose(current, joint)


# ═══════════════════════════════════════════════════════════════
#  5. Normality check: A A† vs A† A
# ═══════════════════════════════════════════════════════════════

def investigate_normality(crystals):
    print("\n" + "=" * 72)
    print("  5. NORMALITY CHECK: A*A† vs A†*A")
    print("=" * 72)

    for name, joint in crystals.items():
        AAt = joint @ joint.conj().T
        AtA = joint.conj().T @ joint
        diff = (AAt - AtA).norm().item()
        norm_A = joint.norm().item()

        # Relative normality defect
        rel_defect = diff / max(norm_A ** 2, 1e-12)

        # Compare eigenvalue magnitudes to singular values
        eigs = eigenvalues(joint)
        svs = singular_values(joint)
        eig_mags = sorted([e.abs().item() for e in eigs], reverse=True)
        sv_list = svs.tolist()
        max_eig_sv_diff = max(abs(eig_mags[i] - sv_list[i]) for i in range(4))

        print(f"\n  {name.upper()}")
        print(f"    ||A A† - A† A||_F = {diff:.8f}")
        print(f"    ||A||_F = {norm_A:.6f}")
        print(f"    Relative normality defect = {rel_defect:.8f}")
        print(f"    Max ||λ_i| - σ_i| = {max_eig_sv_diff:.8f}")
        print(f"    {'NORMAL' if rel_defect < 0.01 else 'NOT NORMAL'}"
              f" (threshold: 0.01)")

        # For non-normal: show the departure structure
        if rel_defect >= 0.01:
            print(f"    Structure of A*A†:")
            for i in range(4):
                row = [f"{AAt[i,j].real:.4f}+{AAt[i,j].imag:.4f}i" for j in range(4)]
                print(f"      [{', '.join(row)}]")


# ═══════════════════════════════════════════════════════════════
#  6. Eigenvalues of the composition operator T_B
# ═══════════════════════════════════════════════════════════════

def investigate_composition_operator(crystals):
    print("\n" + "=" * 72)
    print("  6. EIGENVALUES OF THE COMPOSITION OPERATOR T_B")
    print("=" * 72)
    print("  T_B(A) = compose(A, B) = (A @ B) / L1(A @ B)")
    print("  This is a 16→16 map. We linearize around the identity crystal.")

    for name in ["proportional", "inverse", "modular"]:
        joint_B = crystals[name]

        # Build the 16x16 Jacobian of T_B at the fixed point (joint_B itself)
        # Since L1 normalization makes this nonlinear, we use numerical differentiation
        # T_B(A) = A @ B / L1(A @ B)
        # Linearize: for small perturbation δA around A₀:
        #   T_B(A₀ + εδA) ≈ T_B(A₀) + ε * J * vec(δA)

        A0 = joint_B.clone()
        eps = 1e-5
        J = torch.zeros(16, 16, dtype=torch.cfloat)

        T0 = compose(A0, joint_B).reshape(16)

        for k in range(16):
            dA = torch.zeros(16, dtype=torch.cfloat)
            dA[k] = eps + 0j
            A_plus = (A0.reshape(16) + dA).reshape(4, 4)
            T_plus = compose(A_plus, joint_B).reshape(16)

            dA[k] = eps * 1j
            A_plus_i = (A0.reshape(16) + dA).reshape(4, 4)
            T_plus_i = compose(A_plus_i, joint_B).reshape(16)

            # Complex derivative: df/dz = (df/dx - i df/dy) / 2
            # But since compose is not holomorphic (due to L1 norm),
            # just use real Jacobian approach: df/d(Re z) and df/d(Im z)
            J[:, k] = (T_plus - T0) / eps

        # Eigenvalues of the Jacobian
        eigs_J, _ = torch.linalg.eig(J)
        order = eigs_J.abs().argsort(descending=True)
        eigs_J = eigs_J[order]

        print(f"\n  {name.upper()} — Jacobian of T_B at fixed point:")
        print(f"    Top eigenvalues of J (16×16):")
        for i, e in enumerate(eigs_J[:8]):
            print(f"      λ_{i:2d} = {e.real:>10.6f} + {e.imag:>10.6f}i"
                  f"   |λ| = {e.abs():>10.6f}")
        print(f"      ... (remaining 8 eigenvalues)")
        for i, e in enumerate(eigs_J[8:], start=8):
            print(f"      λ_{i:2d} = {e.real:>10.6f} + {e.imag:>10.6f}i"
                  f"   |λ| = {e.abs():>10.6f}")

        # The spectral gap is 1 - |λ₁| (since |λ₀| should be 1 for fixed point)
        rho_0 = eigs_J[0].abs().item()
        rho_1 = eigs_J[1].abs().item()
        gap = rho_0 - rho_1
        print(f"    Spectral radius ρ₀ = {rho_0:.6f}")
        print(f"    Second eigenvalue |λ₁| = {rho_1:.6f}")
        print(f"    Spectral gap (ρ₀ - |λ₁|) = {gap:.6f}")
        print(f"    Gap ratio |λ₁|/|λ₀| = {rho_1/max(rho_0, 1e-12):.6f}")
        print(f"    -ln(|λ₁|/|λ₀|) = {-math.log(rho_1/max(rho_0, 1e-12)) if rho_1 > 1e-12 else float('inf'):.6f}")
        print(f"    Reference Δ = {DELTA:.6f}")


# ═══════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("╔" + "═" * 70 + "╗")
    print("║  EIGENVALUE SPECTRUM OF [4,4] JOINT MASS CRYSTALS" + " " * 19 + "║")
    print("║  H=3, MASS_DIM=4, dtype=cfloat" + " " * 39 + "║")
    print("╚" + "═" * 70 + "╝")

    crystals = build_all_signatures()

    investigate_signature_eigenvalues(crystals)
    investigate_composition_eigenvalues(crystals)
    investigate_spectral_vs_schmidt(crystals)
    investigate_eigenvalue_ratios(crystals)
    investigate_normality(crystals)
    investigate_composition_operator(crystals)

    # ═══════════════════════════════════════════════════════════
    #  7. Sub-dominant eigenvalue decay rate (where Δ lives)
    # ═══════════════════════════════════════════════════════════
    print("\n" + "=" * 72)
    print("  7. SUB-DOMINANT EIGENVALUE DECAY RATE vs Δ")
    print("=" * 72)
    print("  Under self-composition, |λ₁| decays while |λ₀| stays ~constant.")
    print("  The decay rate of |λ₁|/|λ₀| per step should relate to Δ.")

    for name in ["proportional", "inverse", "modular", "exponential", "logarithmic", "quadratic"]:
        joint = crystals[name]
        current = joint.clone()
        ratios = []

        for step in range(13):
            eigs = eigenvalues(current)
            mags = [e.abs().item() for e in eigs]
            ratio = mags[1] / max(mags[0], 1e-15)
            ratios.append(ratio)
            if step < 12:
                current = compose(current, joint)

        print(f"\n  {name.upper()}:")
        print(f"    {'Step':>4}  {'|λ₁/λ₀|':>12}  {'ln(|λ₁/λ₀|)':>14}  {'step decay':>12}")
        for i, r in enumerate(ratios):
            ln_r = math.log(r) if r > 1e-15 else float('-inf')
            if i > 0 and ratios[i-1] > 1e-15 and r > 1e-15:
                step_decay = math.log(ratios[i-1] / r)
            else:
                step_decay = 0
            print(f"    {i:>4}  {r:>12.8f}  {ln_r:>14.6f}  {step_decay:>12.6f}")

        # Average decay rate from steps 2-10 (skip transient)
        if len(ratios) > 4:
            decay_rates = []
            for i in range(2, min(10, len(ratios))):
                if ratios[i] > 1e-15 and ratios[i-1] > 1e-15:
                    decay_rates.append(math.log(ratios[i-1] / ratios[i]))
            if decay_rates:
                avg_decay = sum(decay_rates) / len(decay_rates)
                print(f"    Average per-step decay rate (steps 2-9): {avg_decay:.6f}")
                print(f"    exp(-Δ) = {math.exp(-DELTA):.6f}")
                print(f"    Ratio to Δ: {avg_decay / DELTA:.6f}")

    # ═══════════════════════════════════════════════════════════
    #  8. Eigenvalue phases: do they lock?
    # ═══════════════════════════════════════════════════════════
    print("\n" + "=" * 72)
    print("  8. EIGENVALUE PHASE EVOLUTION UNDER SELF-COMPOSITION")
    print("=" * 72)

    for name in ["proportional", "inverse"]:
        joint = crystals[name]
        current = joint.clone()

        print(f"\n  {name.upper()}:")
        print(f"    {'Step':>4}  {'φ(λ₀)/π':>10}  {'φ(λ₁)/π':>10}  {'φ(λ₂)/π':>10}  {'φ(λ₃)/π':>10}")

        for step in range(13):
            eigs = eigenvalues(current)
            phases = [math.atan2(e.imag, e.real) / math.pi for e in eigs]
            print(f"    {step:>4}  {phases[0]:>10.6f}  {phases[1]:>10.6f}  "
                  f"{phases[2]:>10.6f}  {phases[3]:>10.6f}")
            if step < 12:
                current = compose(current, joint)

    print("\n" + "=" * 72)
    print("  SUMMARY")
    print("=" * 72)
    print(f"  Δ (structural filter rate) = {DELTA:.6f}")
    print(f"  K* = {K_STAR:.6f}")
    print(f"  exp(-Δ) = {math.exp(-DELTA):.6f}")
    print(f"  H = {H}, MASS_DIM = {MASS_DIM}")
    print("  Done.")
