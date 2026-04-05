"""
What does the Einstein equation look like in the DS framework?

We have:
  - Emergent metric g from mass function
  - Born floor as the source of non-holomorphicity
  - Ricci/Weyl ratio fixed at ~2/3

Question: what effective stress-energy tensor does the floor produce?
Does G_μν = 8πG T_μν take a recognisable form?
"""
import numpy as np
from scipy.optimize import brentq

H = 3
FLOOR = 1.0 / H**3

# ============================================================
# DS machinery (same as before)
# ============================================================
def ds_combine(m, e, apply_floor=True):
    s, theta = m[:3], m[3]
    se, theta_e = e[:3], e[3]
    s_new = s * se + s * theta_e + theta * se
    theta_new = theta * theta_e
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    if abs(denom) < 1e-15:
        return m.copy(), K
    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    total = np.sum(np.abs(m_out))
    if total > 0:
        m_out = m_out / total
    if apply_floor:
        b = np.abs(m_out[3])**2 / np.sum(np.abs(m_out)**2)
        if b < FLOOR:
            m_out = enforce_floor(m_out)
    return m_out, K

def enforce_floor(m):
    s = m[:3].copy()
    theta = m[3]
    ph_t = theta / np.abs(theta) if np.abs(theta) > 1e-15 else 1.0
    ph_s = np.array([si / np.abs(si) if np.abs(si) > 1e-15 else 1.0 for si in s])
    abs_s = np.abs(s)
    lo, hi = np.abs(theta), 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        ss = np.sum(abs_s)
        sc = (1.0 - mid) / ss if ss > 0 else 0
        st = abs_s * sc
        b = mid**2 / (np.sum(st**2) + mid**2)
        if b < FLOOR:
            lo = mid
        else:
            hi = mid
    tn = (lo + hi) / 2
    ss = np.sum(abs_s)
    sc = (1.0 - tn) / ss if ss > 0 else 0
    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = ph_s * abs_s * sc
    m_out[3] = ph_t * tn
    return m_out

def born(m):
    return np.abs(m[3])**2 / np.sum(np.abs(m)**2)

# ============================================================
# Find exact K*=7/30 equilibrium
# ============================================================
def K_at_eq(p_dom_val):
    p_w = (1 - p_dom_val) / 2
    sc = 1 - FLOOR
    raw = np.array([np.sqrt(p_dom_val*sc), np.sqrt(p_w*sc), np.sqrt(p_w*sc), np.sqrt(FLOOR)])
    e = raw / np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(1000):
        m_new, _ = ds_combine(m, e.astype(complex))
        if np.max(np.abs(m_new - m)) < 1e-15:
            break
        m = np.real(m_new)
    _, K = ds_combine(m_new, e.astype(complex), apply_floor=False)
    return np.real(K)

p_exact = brentq(lambda p: K_at_eq(p) - 7/30, 0.92, 0.94, xtol=1e-14)
p_w = (1 - p_exact) / 2
sc = 1 - FLOOR
raw = np.array([np.sqrt(p_exact*sc), np.sqrt(p_w*sc), np.sqrt(p_w*sc), np.sqrt(FLOOR)])
e_star = (raw / np.sum(raw)).astype(complex)

m = np.array([0.4, 0.2, 0.2, 0.2], dtype=complex)
for _ in range(1000):
    m_new, _ = ds_combine(m, e_star)
    if np.max(np.abs(m_new - m)) < 1e-15:
        break
    m = m_new
m_star = m_new

print("=" * 60)
print("K*=7/30 EQUILIBRIUM")
print("=" * 60)
print(f"m* = {np.real(m_star)}")
print(f"K* = {np.real(sum(m_star[i]*e_star[j] for i in range(3) for j in range(3) if i!=j)):.10f}")
print()

# ============================================================
# The emergent metric and its curvature
# ============================================================
# g_μν = diag(θ², -(s₁+θ)², -(s₂+θ)², -(s₃+θ)²)

def metric_from_mass(m):
    """Emergent metric from mass function."""
    th = np.real(m[3])
    s = np.real(m[:3])
    return np.array([th**2, -(s[0]+th)**2, -(s[1]+th)**2, -(s[2]+th)**2])

def christoffel_symbols(m, eps=1e-6):
    """Γ^λ_μν from the emergent metric."""
    n = 4
    coords = np.array([np.real(m[3]), np.real(m[0]), np.real(m[1]), np.real(m[2])])
    # coords = (θ, s₁, s₂, s₃)

    def g_at(c):
        # c = (θ, s₁, s₂, s₃)
        th = c[0]
        return np.diag([th**2, -(c[1]+th)**2, -(c[2]+th)**2, -(c[3]+th)**2])

    g = g_at(coords)
    g_inv = np.diag(1.0 / np.diag(g))

    dg = np.zeros((n, n, n))
    for a in range(n):
        cp = coords.copy(); cp[a] += eps
        cm = coords.copy(); cm[a] -= eps
        dg[:, :, a] = (g_at(cp) - g_at(cm)) / (2 * eps)

    Gamma = np.zeros((n, n, n))
    for lam in range(n):
        for mu in range(n):
            for nu in range(n):
                total = 0.0
                for sig in range(n):
                    total += g_inv[lam, sig] * (dg[sig, nu, mu] + dg[sig, mu, nu] - dg[mu, nu, sig])
                Gamma[lam, mu, nu] = 0.5 * total
    return Gamma

def riemann_tensor(m, eps=1e-4):
    """R^ρ_σμν"""
    n = 4
    coords = np.array([np.real(m[3]), np.real(m[0]), np.real(m[1]), np.real(m[2])])

    def gamma_at(c):
        m_temp = np.array([c[1], c[2], c[3], c[0]], dtype=complex)
        return christoffel_symbols(m_temp, eps=eps*0.1)

    G = gamma_at(coords)
    dG = np.zeros((n, n, n, n))
    for mu in range(n):
        cp = coords.copy(); cp[mu] += eps
        cm = coords.copy(); cm[mu] -= eps
        dG[:, :, :, mu] = (gamma_at(cp) - gamma_at(cm)) / (2 * eps)

    R = np.zeros((n, n, n, n))
    for rho in range(n):
        for sig in range(n):
            for mu in range(n):
                for nu in range(n):
                    R[rho, sig, mu, nu] = dG[rho, nu, sig, mu] - dG[rho, mu, sig, nu]
                    for lam in range(n):
                        R[rho, sig, mu, nu] += G[rho, mu, lam] * G[lam, nu, sig] - G[rho, nu, lam] * G[lam, mu, sig]
    return R

def ricci_from_riemann(R):
    n = 4
    Ric = np.zeros((n, n))
    for mu in range(n):
        for nu in range(n):
            for lam in range(n):
                Ric[mu, nu] += R[lam, mu, lam, nu]
    return Ric

def scalar_curvature(Ric, m):
    g = metric_from_mass(m)
    g_inv = 1.0 / g
    return sum(g_inv[i] * Ric[i, i] for i in range(4))

def einstein_tensor(Ric, R_scalar, m):
    g = np.diag(metric_from_mass(m))
    return Ric - 0.5 * R_scalar * g

def weyl_tensor(R, Ric, R_scalar, m):
    """C^ρ_σμν = R^ρ_σμν - (Ricci terms) + (scalar terms)"""
    n = 4
    g = np.diag(metric_from_mass(m))
    g_inv = np.diag(1.0 / np.diag(g))

    # C_ρσμν = R_ρσμν - (2/(n-2))(g_ρ[μ R_ν]σ - g_σ[μ R_ν]ρ) + (2/((n-1)(n-2))) R g_ρ[μ g_ν]σ
    # For n=4:
    # C_ρσμν = R_ρσμν - (g_ρμ R_νσ - g_ρν R_μσ - g_σμ R_νρ + g_σν R_μρ)/2
    #          + R (g_ρμ g_νσ - g_ρν g_μσ)/6

    C = np.zeros((n, n, n, n))
    for rho in range(n):
        for sig in range(n):
            for mu in range(n):
                for nu in range(n):
                    # R_ρσμν = g_ρα R^α_σμν
                    R_lower = g[rho, rho] * R[rho, sig, mu, nu]

                    # Ricci terms
                    ricci_part = (
                        g[rho, rho]*(1 if rho==mu else 0) * Ric[nu, sig]
                        - g[rho, rho]*(1 if rho==nu else 0) * Ric[mu, sig]
                        - g[sig, sig]*(1 if sig==mu else 0) * Ric[nu, rho]
                        + g[sig, sig]*(1 if sig==nu else 0) * Ric[mu, rho]
                    ) / 2

                    # Scalar terms
                    scalar_part = R_scalar * (
                        g[rho, rho]*(1 if rho==mu else 0) * g[nu, nu]*(1 if nu==sig else 0)
                        - g[rho, rho]*(1 if rho==nu else 0) * g[mu, mu]*(1 if mu==sig else 0)
                    ) / 6

                    C[rho, sig, mu, nu] = R_lower - ricci_part + scalar_part
    return C

# ============================================================
print("=" * 60)
print("CURVATURE AT K*=7/30 EQUILIBRIUM")
print("=" * 60)

g_eq = metric_from_mass(m_star)
print(f"Metric: {g_eq}")
print()

R = riemann_tensor(m_star)
Ric = ricci_from_riemann(R)
R_s = scalar_curvature(Ric, m_star)
G = einstein_tensor(Ric, R_s, m_star)

labels = ['θ', 's₁', 's₂', 's₃']
print("Ricci tensor (diagonal):")
for i in range(4):
    print(f"  R_{labels[i]}{labels[i]} = {Ric[i,i]:+.6f}")
print(f"\nScalar curvature R = {R_s:.6f}")
print()

print("Einstein tensor (diagonal):")
for i in range(4):
    print(f"  G_{labels[i]}{labels[i]} = {G[i,i]:+.6f}")
print()

# ============================================================
# Effective stress-energy: T_μν = G_μν / (8πG)
# We don't know G, but we can look at the STRUCTURE of T
# ============================================================
print("=" * 60)
print("EFFECTIVE STRESS-ENERGY TENSOR")
print("=" * 60)

# T_μν = G_μν / (8πG). Setting 8πG = 1 for now.
T = G.copy()

# For a perfect fluid: T_μν = diag(ρ, -p, -p, -p) in the rest frame
# (with our signature (+,-,-,-))
# ρ = energy density, p = pressure

rho = T[0, 0] / g_eq[0]  # T^0_0 = ρ for diagonal metric
p1 = -T[1, 1] / (-g_eq[1])  # T^1_1 = -p
p2 = -T[2, 2] / (-g_eq[2])
p3 = -T[3, 3] / (-g_eq[3])

# Actually for diagonal metric: T^μ_ν = g^μμ T_μν
# T^0_0 = G^0_0 = g^00 G_00 = G_00/g_00
# T^i_i = G^i_i = g^ii G_ii = G_ii/g_ii

T_mixed = np.array([G[i,i] / g_eq[i] for i in range(4)])
print(f"T^μ_ν (mixed, diagonal):")
print(f"  T^0_0 = ρ  = {T_mixed[0]:+.6f}")
print(f"  T^1_1 = -p₁ = {T_mixed[1]:+.6f}")
print(f"  T^2_2 = -p₂ = {T_mixed[2]:+.6f}")
print(f"  T^3_3 = -p₃ = {T_mixed[3]:+.6f}")
print()

rho_eff = T_mixed[0]
p_dom = -T_mixed[1]
p_weak = -T_mixed[2]

print(f"Energy density: ρ = {rho_eff:+.6f}")
print(f"Pressure (dominant): p₁ = {p_dom:+.6f}")
print(f"Pressure (weak): p₂ = p₃ = {p_weak:+.6f}")
print()

# Equation of state
if abs(rho_eff) > 1e-10:
    w_dom = p_dom / rho_eff
    w_weak = p_weak / rho_eff
    w_avg = (p_dom + 2*p_weak) / (3 * rho_eff)
    print(f"Equation of state:")
    print(f"  w_dom = p₁/ρ = {w_dom:.6f}")
    print(f"  w_weak = p₂/ρ = {w_weak:.6f}")
    print(f"  w_avg = <p>/ρ = {w_avg:.6f}")
    print()
    print(f"Reference values:")
    print(f"  w = 0: dust (non-relativistic matter)")
    print(f"  w = 1/3: radiation")
    print(f"  w = -1: cosmological constant")
    print(f"  w = 1: stiff matter")
    print()

# Trace of stress-energy
T_trace = sum(T_mixed)
print(f"Trace T = ρ - p₁ - p₂ - p₃ = {T_trace:.6f}")
print(f"(Traceless = conformal matter, e.g. radiation with T=0)")
print()

# ============================================================
print("=" * 60)
print("KEY RATIOS")
print("=" * 60)

print(f"ρ/R = {rho_eff/R_s:.6f}")
print(f"K* = {7/30:.6f}")
print(f"1/H = {1/3:.6f}")
print(f"(H-1)/H = {2/3:.6f}")
print(f"1/H³ = {1/27:.6f}")
print()

# Check: is ρ related to K*?
print("Is ρ related to framework constants?")
for name, val in [
    ("K*", 7/30),
    ("1-K*", 23/30),
    ("K*²", (7/30)**2),
    ("1/H", 1/3),
    ("H/(H²+1)", 3/10),
    ("(H²-1)/H²", 8/9),
    ("Born floor", 1/27),
    ("2/3", 2/3),
    ("Δ", 1.263),
    ("λ₀", 0.2829),
]:
    if abs(val) > 1e-10:
        ratio = rho_eff / val
        print(f"  ρ/{name} = {ratio:.6f}")

print()

# ============================================================
print("=" * 60)
print("WEYL TENSOR NORM")
print("=" * 60)

C = weyl_tensor(R, Ric, R_s, m_star)

# |C|² = C_ρσμν C^ρσμν
C_sq = 0.0
g_inv = 1.0 / g_eq
for rho in range(4):
    for sig in range(4):
        for mu in range(4):
            for nu in range(4):
                # C^ρσμν = g^σσ g^μμ g^νν C_ρσμν / g_ρρ...
                # Actually C has first index up already from Riemann
                # Let me compute C_lower then raise
                C_lower = C[rho, sig, mu, nu]  # this is already C_ρσμν
                C_upper = C_lower * g_inv[rho] * g_inv[sig] * g_inv[mu] * g_inv[nu]
                C_sq += C_lower * C_upper

print(f"|C|² (Weyl squared) = {C_sq:.6f}")
print(f"|R_μν|² (Ricci squared) = {sum(Ric[i,i]**2 * g_inv[i]**2 for i in range(4)):.6f}")
print(f"R² (scalar squared) = {R_s**2:.6f}")
print()

Ric_sq = sum(Ric[i,i]**2 * abs(g_inv[i]) for i in range(4))  # |Ric|²
print(f"Ricci/Weyl ratio: {np.sqrt(Ric_sq)/np.sqrt(abs(C_sq)) if abs(C_sq)>1e-10 else 'inf'}")

# ============================================================
print("\n" + "=" * 60)
print("COSMOLOGICAL CONSTANT?")
print("=" * 60)

# Einstein equations with Λ: G_μν + Λg_μν = 8πG T_μν
# If T=0 (vacuum): G_μν = -Λg_μν
# Check: is G_μν proportional to g_μν?

ratios = [G[i,i] / g_eq[i] for i in range(4)]
print(f"G_μμ/g_μμ for each direction:")
for i in range(4):
    print(f"  {labels[i]}: {ratios[i]:+.6f}")

print(f"\nMean = {np.mean(ratios):.6f}")
print(f"Spread = {np.std(ratios):.6f}")
print()

if np.std(ratios) / abs(np.mean(ratios)) < 0.1:
    Lambda_eff = -np.mean(ratios)
    print(f"G_μν ≈ -Λ g_μν with Λ = {Lambda_eff:.6f}")
    print(f"This is a VACUUM + COSMOLOGICAL CONSTANT solution!")
    print(f"Λ/K* = {Lambda_eff/(7/30):.6f}")
    print(f"Λ × H² = {Lambda_eff * 9:.6f}")
else:
    print("G_μν is NOT proportional to g_μν.")
    print("There is non-trivial matter content (anisotropic stress).")
    print()
    # Decompose: G_μν = -Λg_μν + 8πG T^matter_μν
    Lambda_eff = -np.mean(ratios)
    T_matter = np.array([G[i,i] + Lambda_eff * g_eq[i] for i in range(4)])
    print(f"Decomposition: G = -Λg + T_matter")
    print(f"  Λ_eff = {Lambda_eff:.6f}")
    print(f"  T_matter (diagonal): {T_matter}")
    print(f"  Anisotropy: T_matter[1] vs T_matter[2] = {T_matter[1]/T_matter[2] if abs(T_matter[2])>1e-10 else 'inf':.4f}")
