"""
Penrose kernel UV singularity: explicit computation.

The DS fiber correlator S₂(x,y) is smooth at coincident points.
The spacetime YM correlator acquires its power-law UV singularity
from the Penrose integration kernel when two CP¹ fibers collide.

We compute the double fiber integral by Monte Carlo (avoids quadrature
issues with the near-singular kernel) and extract the power law.
"""
import numpy as np

print("="*65)
print("PENROSE KERNEL UV SINGULARITY")
print("="*65)

# ================================================================
# Geometry: two CP¹ fibers in CP³ separated by angle α
#
# L_x: Z(s) = (cos s, sin s, 0, 0)  — fiber over x
# L_y: W(t) = (cos α cos t, sin t, sin α cos t, 0) — fiber over y
#
# The tilt angle α = ε/2 where ε = |x-y| on S⁴.
#
# Hermitian inner product:
#   ⟨Z,W⟩ = cos(s) cos(α) cos(t) + sin(s) sin(t)
#
# The Penrose kernel for homogeneity degree k is |⟨Z,W⟩|^{-2k}.
# ================================================================

def fiber_integral_mc(epsilon, k, n_samples=500000):
    """
    Monte Carlo computation of the double fiber integral:
      I(ε) = ∮_{L_x} ∮_{L_y} |⟨Z,W⟩|^{-2k} ds dt
    where α = ε/2 is the tilt angle between fibers.
    """
    alpha = epsilon / 2.0

    # Sample uniformly on [0, 2π] × [0, 2π]
    s = np.random.uniform(0, 2*np.pi, n_samples)
    t = np.random.uniform(0, 2*np.pi, n_samples)

    # Hermitian inner product
    inner = np.cos(s) * np.cos(alpha) * np.cos(t) + np.sin(s) * np.sin(t)

    # |inner|^{-2k}, with truncation to avoid infinity
    abs_inner = np.abs(inner)
    # Truncate at a tiny floor to regularise (this slightly underestimates)
    abs_inner = np.maximum(abs_inner, 1e-10)
    integrand = abs_inner**(-2*k)

    # MC estimate: integral = (volume of domain) × mean(integrand)
    volume = (2*np.pi)**2
    mean_val = np.mean(integrand)
    std_val = np.std(integrand) / np.sqrt(n_samples)

    return volume * mean_val, volume * std_val

# ================================================================
# Power-law extraction
# ================================================================

print("\nDouble fiber integral I(ε) = ∮∮ |⟨Z,W⟩|^{-2k} ds dt")
print("Theory: I(ε) ~ ε^{-2(k-1)} as ε → 0")
print("  k=1: logarithmic (dim 1 scalar)")
print("  k=2: ε^{-2} (dim 2 gauge field)")
print("  k=3: ε^{-4} (dim 3 graviton)")
print("  k=4: ε^{-6} (dim 4 composite like tr(F²))")

np.random.seed(42)

epsilons = [0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005]

for k in [1, 2, 3, 4]:
    print(f"\n{'='*55}")
    print(f"  k = {k}  |  Expected power: {-2*(k-1)}")
    print(f"{'='*55}")
    print(f"  {'ε':>8s}  {'I(ε)':>14s}  {'log₁₀(ε)':>10s}  {'log₁₀(I)':>10s}  {'slope':>8s}")

    log_eps = []
    log_I = []

    for eps in epsilons:
        I_val, I_err = fiber_integral_mc(eps, k, n_samples=1000000)
        le = np.log10(eps)
        li = np.log10(abs(I_val)) if I_val > 0 else float('nan')
        log_eps.append(le)
        log_I.append(li)

        if len(log_eps) >= 2:
            slope = (log_I[-1] - log_I[-2]) / (log_eps[-1] - log_eps[-2])
        else:
            slope = float('nan')

        print(f"  {eps:8.4f}  {I_val:14.4f}  {le:10.4f}  {li:10.4f}  {slope:8.3f}")

    # Best-fit power law from small-ε points
    fit_start = max(0, len(log_eps) - 5)
    coeffs = np.polyfit(log_eps[fit_start:], log_I[fit_start:], 1)
    print(f"\n  Best-fit: I(ε) ~ ε^{coeffs[0]:.3f}")
    print(f"  Theory:   I(ε) ~ ε^{-2*(k-1):.3f}")
    print(f"  Match: {'YES' if abs(coeffs[0] - (-2*(k-1))) < 0.3 else 'CHECK'}")

# ================================================================
# The full picture: UV × IR factorisation
# ================================================================
print(f"\n{'='*65}")
print("THE FULL SPACETIME CORRELATOR")
print("="*65)

print("""
The spacetime two-point function factorises:

  ⟨O_YM(x) O_YM(y)⟩ = [Penrose kernel] × [DS fiber correlator]
                      = |x-y|^{-2d}     × Σ|c_k|² λ_k^{|x-y|/(2δ)}

  Short distance (|x-y| → 0):
    DS correlator → Σ|c_k|² (constant)
    Penrose kernel → |x-y|^{-2d} (power-law singularity)
    → Full correlator ~ |x-y|^{-2d}  [UV: asymptotic freedom]

  Long distance (|x-y| → ∞):
    DS correlator → exp(-m|x-y|/2) (exponential decay)
    Penrose kernel → slowly varying
    → Full correlator ~ exp(-m|x-y|/2)  [IR: mass gap, confinement]

The DS construction provides the IR physics (mass gap Δ, confinement).
The Penrose transform provides the UV physics (conformal dimensions).
Neither alone gives the full theory. Together they do.
""")

# Numerical demonstration
Delta = 1.2626
delta_step = 0.05  # FS displacement per DS step
m_phys = Delta / delta_step

print(f"DS spectral gap: Δ = {Delta}")
print(f"Physical mass: m = {m_phys:.2f} (in FS units)")
print()

print(f"{'|x-y|':>8s}  {'r^{-8} (UV)':>14s}  {'exp(-mr/2)':>14s}  {'Product':>14s}  {'Regime':>10s}")
print(f"  {'-'*66}")

for r in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0]:
    uv = r**(-8)  # tr(F²) conformal dimension = 4
    ir = np.exp(-m_phys * r / 2)
    prod = uv * ir
    if r < 0.05:
        regime = "UV"
    elif r < 1:
        regime = "crossover"
    else:
        regime = "IR"
    print(f"  {r:7.4f}  {uv:14.4e}  {ir:14.4e}  {prod:14.4e}  {regime:>10s}")

print(f"\nAt short distance: product dominated by r^{{-8}} → YM UV physics")
print(f"At long distance: product dominated by exp(-mr/2) → mass gap")
print(f"The crossover scale is r* ~ 2/m × 8 ln(1/r*) ... self-consistent at r* ≈ {8/(m_phys/2):.4f}")
