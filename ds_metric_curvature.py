"""
Compute the curvature of the DS emergent metric.

g = diag(θ², -(s₁+θ)², -(s₂+θ)², -(s₃+θ)²)

This is a diagonal metric on the 4D mass function space.
We compute Christoffel symbols, Riemann tensor, Ricci tensor,
Ricci scalar, and Einstein tensor.

We evaluate at:
  - Complete ignorance: θ=1, sᵢ=0
  - Equilibrium: θ*=0.5959, sᵢ=0.1347
  - Generic points along the DS trajectory
"""
import numpy as np
from itertools import product

def ds_metric(theta, s):
    """Return the diagonal metric components g_μμ."""
    s1, s2, s3 = s
    return np.array([
        theta**2,
        -(s1 + theta)**2,
        -(s2 + theta)**2,
        -(s3 + theta)**2
    ])

def christoffel(theta, s, eps=1e-7):
    """Compute Christoffel symbols Γ^λ_μν numerically."""
    n = 4
    coords = np.array([theta, s[0], s[1], s[2]])

    # Compute metric and its derivatives
    def metric_at(c):
        return np.diag(ds_metric(c[0], c[1:]))

    g = metric_at(coords)
    g_inv = np.diag(1.0 / np.diag(g))

    # Partial derivatives of metric: dg[μν, α] = ∂g_μν/∂x^α
    dg = np.zeros((n, n, n))
    for a in range(n):
        coords_plus = coords.copy()
        coords_minus = coords.copy()
        coords_plus[a] += eps
        coords_minus[a] -= eps
        dg[:, :, a] = (metric_at(coords_plus) - metric_at(coords_minus)) / (2 * eps)

    # Γ^λ_μν = (1/2) g^λσ (∂_μ g_σν + ∂_ν g_σμ - ∂_σ g_μν)
    Gamma = np.zeros((n, n, n))
    for lam in range(n):
        for mu in range(n):
            for nu in range(n):
                total = 0.0
                for sig in range(n):
                    total += g_inv[lam, sig] * (
                        dg[sig, nu, mu] + dg[sig, mu, nu] - dg[mu, nu, sig]
                    )
                Gamma[lam, mu, nu] = 0.5 * total

    return Gamma

def riemann(theta, s, eps=1e-5):
    """Compute Riemann tensor R^ρ_σμν numerically."""
    n = 4
    coords = np.array([theta, s[0], s[1], s[2]])

    # Get Christoffel symbols at nearby points for derivatives
    def gamma_at(c):
        return christoffel(c[0], c[1:], eps=eps*0.1)

    G = gamma_at(coords)

    # ∂_μ Γ^ρ_νσ
    dG = np.zeros((n, n, n, n))  # dG[ρ, ν, σ, μ] = ∂_μ Γ^ρ_νσ
    for mu in range(n):
        cp = coords.copy(); cp[mu] += eps
        cm = coords.copy(); cm[mu] -= eps
        dG[:, :, :, mu] = (gamma_at(cp) - gamma_at(cm)) / (2 * eps)

    # R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ
    R = np.zeros((n, n, n, n))
    for rho in range(n):
        for sig in range(n):
            for mu in range(n):
                for nu in range(n):
                    R[rho, sig, mu, nu] = (
                        dG[rho, nu, sig, mu] - dG[rho, mu, sig, nu]
                    )
                    for lam in range(n):
                        R[rho, sig, mu, nu] += (
                            G[rho, mu, lam] * G[lam, nu, sig]
                            - G[rho, nu, lam] * G[lam, mu, sig]
                        )
    return R

def ricci_tensor(R):
    """Contract Riemann to get Ricci tensor R_μν = R^λ_μλν."""
    n = 4
    Ric = np.zeros((n, n))
    for mu in range(n):
        for nu in range(n):
            for lam in range(n):
                Ric[mu, nu] += R[lam, mu, lam, nu]
    return Ric

def ricci_scalar(Ric, theta, s):
    """Contract Ricci tensor with inverse metric."""
    g = ds_metric(theta, s)
    g_inv = 1.0 / g
    R_scalar = 0.0
    for mu in range(4):
        R_scalar += g_inv[mu] * Ric[mu, mu]  # diagonal metric
    return R_scalar

def einstein_tensor(Ric, R_scalar, theta, s):
    """G_μν = R_μν - (1/2) g_μν R."""
    g = np.diag(ds_metric(theta, s))
    G = Ric - 0.5 * R_scalar * g
    return G

# ============================================================
print("=" * 60)
print("COMPLETE IGNORANCE: θ=1, sᵢ=0")
print("=" * 60)

theta0, s0 = 1.0, [0.0, 0.0, 0.0]
g0 = ds_metric(theta0, s0)
print(f"Metric: {g0}")
print(f"Signature: ({'+' if g0[0]>0 else '-'}, {'+' if g0[1]>0 else '-'}, {'+' if g0[2]>0 else '-'}, {'+' if g0[3]>0 else '-'})")

# Problem: sᵢ=0 means ∂g/∂sᵢ involves (sᵢ+θ) = θ
# The metric is well-defined but let's check curvature
G0 = christoffel(theta0, s0)
print(f"\nNon-zero Christoffel symbols:")
for i, j, k in product(range(4), repeat=3):
    if abs(G0[i, j, k]) > 1e-8:
        labels = ['θ', 's₁', 's₂', 's₃']
        print(f"  Γ^{labels[i]}_{labels[j]}{labels[k]} = {G0[i,j,k]:.6f}")

R0 = riemann(theta0, s0)
Ric0 = ricci_tensor(R0)
RS0 = ricci_scalar(Ric0, theta0, s0)
print(f"\nRicci tensor (diagonal):")
for i in range(4):
    labels = ['θ', 's₁', 's₂', 's₃']
    print(f"  R_{labels[i]}{labels[i]} = {Ric0[i,i]:.8f}")
print(f"\nRicci scalar R = {RS0:.8f}")

E0 = einstein_tensor(Ric0, RS0, theta0, s0)
print(f"\nEinstein tensor (diagonal):")
for i in range(4):
    labels = ['θ', 's₁', 's₂', 's₃']
    print(f"  G_{labels[i]}{labels[i]} = {E0[i,i]:.8f}")

# ============================================================
print("\n" + "=" * 60)
print("EQUILIBRIUM: θ*=0.5959, sᵢ=0.1347")
print("=" * 60)

theta_eq, s_eq = 0.5959, [0.1347, 0.1347, 0.1347]
g_eq = ds_metric(theta_eq, s_eq)
print(f"Metric: {g_eq}")
ratio = g_eq[0] / g_eq[1]
print(f"g₀₀/g₁₁ = {ratio:.6f}")
print(f"(H-1)/H = {2/3:.6f}")

G_eq = christoffel(theta_eq, s_eq)
print(f"\nNon-zero Christoffel symbols:")
for i, j, k in product(range(4), repeat=3):
    if abs(G_eq[i, j, k]) > 1e-8:
        labels = ['θ', 's₁', 's₂', 's₃']
        print(f"  Γ^{labels[i]}_{labels[j]}{labels[k]} = {G_eq[i,j,k]:.6f}")

R_eq = riemann(theta_eq, s_eq)
Ric_eq = ricci_tensor(R_eq)
RS_eq = ricci_scalar(Ric_eq, theta_eq, s_eq)
print(f"\nRicci tensor (diagonal):")
for i in range(4):
    labels = ['θ', 's₁', 's₂', 's₃']
    print(f"  R_{labels[i]}{labels[i]} = {Ric_eq[i,i]:.8f}")
print(f"\nRicci scalar R = {RS_eq:.8f}")

E_eq = einstein_tensor(Ric_eq, RS_eq, theta_eq, s_eq)
print(f"\nEinstein tensor (diagonal):")
for i in range(4):
    labels = ['θ', 's₁', 's₂', 's₃']
    print(f"  G_{labels[i]}{labels[i]} = {E_eq[i,i]:.8f}")

# ============================================================
print("\n" + "=" * 60)
print("ALONG THE DS TRAJECTORY")
print("=" * 60)

# Simple DS step
def ds_step(theta, s, e_theta=0.5, e_s=None):
    """One DS combination step with evidence."""
    if e_s is None:
        e_s = [0.2, 0.15, 0.15]
    # Unnormalised
    m = [theta, s[0], s[1], s[2]]
    e = [e_theta, e_s[0], e_s[1], e_s[2]]

    # DS combination (simplified): multiply and renormalise
    raw = [m[i] * e[i] for i in range(4)]
    K = 1 - sum(raw)  # conflict
    if abs(K) < 1e-15:
        return theta, s, 0
    combined = [r / (1 - K) for r in raw]  # Dempster normalisation

    # This is a simplification - real DS is more complex
    # but captures the contraction
    total = sum(combined)
    normalised = [c / total for c in combined]

    return normalised[0], normalised[1:], K

# Track curvature along trajectory
print(f"{'Step':>4s} {'θ':>8s} {'s₁':>8s} {'R':>12s} {'g00/g11':>10s}")
print("-" * 50)

theta_t, s_t = 0.9, [0.03, 0.03, 0.04]
for step in range(15):
    g_t = ds_metric(theta_t, s_t)
    if abs(g_t[1]) > 1e-10:  # avoid division by zero
        try:
            R_t = riemann(theta_t, s_t)
            Ric_t = ricci_tensor(R_t)
            RS_t = ricci_scalar(Ric_t, theta_t, s_t)
            ratio_t = g_t[0] / g_t[1]
            print(f"{step:4d} {theta_t:8.4f} {s_t[0]:8.4f} {RS_t:12.6f} {ratio_t:10.6f}")
        except:
            print(f"{step:4d} {theta_t:8.4f} {s_t[0]:8.4f} {'error':>12s}")

    theta_t, s_t, K_t = ds_step(theta_t, s_t)

# ============================================================
print("\n" + "=" * 60)
print("KRETSCHNER SCALAR (curvature invariant)")
print("=" * 60)

def kretschner(R, theta, s):
    """K = R_αβμν R^αβμν - the curvature invariant."""
    g = ds_metric(theta, s)
    g_inv = 1.0 / g  # diagonal

    K = 0.0
    for a, b, m, n in product(range(4), repeat=4):
        # Lower all indices of R^ρ_σμν to get R_ρσμν
        R_lower = g[a] * R[a, b, m, n]  # R_aσμν = g_aρ R^ρ_σμν (diagonal)
        # Raise all indices
        R_upper = g_inv[a] * g_inv[b] * g_inv[m] * g_inv[n] * R_lower
        # Actually need R_abmn R^abmn
        # R^abmn = g^aα g^bβ g^mμ g^nν R_αβμν
        # For diagonal metric: R^abmn = R_abmn / (g_a g_b g_m g_n) ... but signs matter
        # Let me just compute R_abmn R_abmn with appropriate index raising
        pass

    # Simpler: K = R_ρσμν R^ρσμν
    # R^ρσμν = g^σα R^ρ_αμν (one index already up)
    # Then raise μ,ν with g^μβ, g^νγ
    K = 0.0
    for rho, sig, mu, nu in product(range(4), repeat=4):
        R_up = R[rho, sig, mu, nu]  # R^ρ_σμν
        # Raise σ: g^σσ (diagonal)
        R_up_2 = g_inv[sig] * R_up  # R^ρσ_μν  (with sign from g_inv)
        # Raise μ:
        R_up_3 = g_inv[mu] * R_up_2
        # Raise ν:
        R_up_4 = g_inv[nu] * R_up_3

        # Lower version: R_ρσμν = g_ρρ R^ρ_σμν (diagonal)
        R_low = g[rho] * R[rho, sig, mu, nu]

        # But we need R_ρσμν with all lower = g_ρα R^α_σμν
        # = g_ρρ R^ρ_σμν (diagonal, no sum)
        # Hmm this only works for the first index
        # Actually R_ρσμν = g_ρα R^α_σμν, and for diagonal g, = g_ρρ R^ρ_σμν (no sum on ρ)

        # K = R^ρσμν R_ρσμν = Σ (g^σσ g^μμ g^νν) R^ρ_σμν × g_ρρ R^ρ_σμν
        pass

    # Let me just do this cleanly
    K = 0.0
    for rho, sig, mu, nu in product(range(4), repeat=4):
        # R_ρσμν = g_ρρ R^ρ_σμν
        R_lower_all = g[rho] * R[rho, sig, mu, nu]
        # R^ρσμν = g^σσ g^μμ g^νν R^ρ_σμν
        R_upper_all = g_inv[sig] * g_inv[mu] * g_inv[nu] * R[rho, sig, mu, nu]
        K += R_lower_all * R_upper_all  # Hmm, double counting ρ

    # Actually the clean way:
    # K = g^αρ g^βσ g^γμ g^δν R_αβγδ R_ρσμν
    # First compute R_αβγδ = g_αα R^α_βγδ
    K = 0.0
    for a, b, c, d in product(range(4), repeat=4):
        R1 = g[a] * R[a, b, c, d]  # R_abcd
        R2 = g[a] * R[a, b, c, d]  # R_abcd again
        # Raise all of R2: g^aa g^bb g^cc g^dd R_abcd
        R2_up = g_inv[a] * g_inv[b] * g_inv[c] * g_inv[d] * R2
        K += R1 * R2_up

    return K

for name, th, ss in [
    ("Ignorance", 1.0, [0.0, 0.0, 0.0]),
    ("Equilibrium", 0.5959, [0.1347, 0.1347, 0.1347]),
    ("Mid-trajectory", 0.7, [0.1, 0.1, 0.1]),
]:
    R_k = riemann(th, ss)
    K_k = kretschner(R_k, th, ss)
    print(f"{name:20s}: Kretschner = {K_k:.8f}")
