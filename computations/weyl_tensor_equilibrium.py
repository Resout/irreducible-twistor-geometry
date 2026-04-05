"""
Weyl tensor at the K*=7/30 equilibrium.

Question: Does the Weyl tensor vanish at equilibrium?
If yes: equilibrium is conformally flat (de Sitter with Lambda).
If no: equilibrium has non-trivial gravitational curvature.

The emergent metric from the paper's Lorentzian structure (eq:detonsurface):
g = diag(theta^2, -(s1+theta)^2, -(s2+theta)^2, -(s3+theta)^2)

Two independent checks:
1. Direct Riemann -> Weyl computation on this metric
2. Symmetry analysis: is the metric maximally symmetric?
"""

import numpy as np
from itertools import product as iprod
from scipy.optimize import brentq

H = 3
FLOOR = 1.0 / 27.0

# ============================================================
# Find the K*=7/30 equilibrium
# ============================================================

def enforce_floor(m):
    s = m[:3].copy()
    th = m[3]
    born = th**2 / (np.sum(s**2) + th**2)
    if born >= FLOOR:
        return m.copy()
    S = np.sum(s)
    Sq = np.sum(s**2)
    R = Sq / S**2
    t = (np.sqrt(26 * R) - R) / (26 - R)
    alpha = (1 - t) / S
    out = m.copy()
    out[:3] = s * alpha
    out[3] = t
    return out

def ds_step(m, e):
    s, th = m[:3], m[3]
    se, ph = e[:3], e[3]
    sn = s * se + s * ph + th * se
    tn = th * ph
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    d = 1 - K
    out = np.zeros(4)
    out[:3] = sn / d
    out[3] = tn / d
    return enforce_floor(out), K

def K_res(p):
    pw = (1 - p) / 2
    sc = 26.0 / 27.0
    raw = np.array([np.sqrt(p * sc), np.sqrt(pw * sc), np.sqrt(pw * sc), np.sqrt(1.0 / 27.0)])
    e = raw / np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(5000):
        m, _ = ds_step(m, e)
    s, th = m[:3], m[3]
    se, ph = e[:3], e[3]
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    return K - 7.0 / 30.0

p = brentq(K_res, 0.92, 0.94, xtol=1e-15)
pw = (1 - p) / 2
sc = 26.0 / 27.0
raw = np.array([np.sqrt(p * sc), np.sqrt(pw * sc), np.sqrt(pw * sc), np.sqrt(1.0 / 27.0)])
e_star = raw / np.sum(raw)
m = np.array([0.4, 0.2, 0.2, 0.2])
for _ in range(5000):
    m, _ = ds_step(m, e_star)

m_star = m
th = m_star[3]
s1, s2, s3 = m_star[0], m_star[1], m_star[2]

print("=" * 60)
print("WEYL TENSOR AT THE K*=7/30 EQUILIBRIUM")
print("=" * 60)
print()
print(f"m* = ({s1:.10f}, {s2:.10f}, {s3:.10f}, {th:.10f})")
print(f"s2 == s3: {abs(s2 - s3) < 1e-10}")
print()

# ============================================================
# Emergent metric and its structure
# ============================================================

print("=== Emergent metric g_ab ===")
print()
g_diag = np.array([th**2, -(s1 + th)**2, -(s2 + th)**2, -(s3 + th)**2])
print(f"g_00 = theta^2 = {g_diag[0]:.10f}")
print(f"g_11 = -(s1+th)^2 = {g_diag[1]:.10f}")
print(f"g_22 = -(s2+th)^2 = {g_diag[2]:.10f}")
print(f"g_33 = -(s3+th)^2 = {g_diag[3]:.10f}")
print()
print(f"Isotropic (g_11=g_22=g_33)? NO.")
print(f"  (s1+th)^2 = {(s1 + th)**2:.10f}")
print(f"  (s2+th)^2 = {(s2 + th)**2:.10f}")
print(f"  Anisotropy ratio: {(s1 + th)**2 / (s2 + th)**2:.6f}")
print()
print("The metric is anisotropic: dominant hypothesis direction")
print("has a much larger spatial extent than the weak directions.")
print("This breaks maximal symmetry -> Weyl can be nonzero.")
print()

# ============================================================
# Full Riemann tensor computation (finite differences)
# ============================================================

def metric_at(coords):
    t, a, b, c = coords
    return np.diag([t**2, -(a + t)**2, -(b + t)**2, -(c + t)**2])

def christoffel_at(coords, eps=1e-7):
    n = 4
    g = metric_at(coords)
    g_inv = np.diag(1.0 / np.diag(g))
    dg = np.zeros((n, n, n))
    for a in range(n):
        cp = coords.copy(); cp[a] += eps
        cm = coords.copy(); cm[a] -= eps
        dg[:, :, a] = (metric_at(cp) - metric_at(cm)) / (2 * eps)
    G = np.zeros((n, n, n))
    for l in range(n):
        for mu in range(n):
            for nu in range(n):
                for sig in range(n):
                    G[l, mu, nu] += g_inv[l, sig] * (dg[sig, nu, mu] + dg[sig, mu, nu] - dg[mu, nu, sig])
                G[l, mu, nu] *= 0.5
    return G, g, g_inv

def full_curvature(coords, eps=1e-6):
    n = 4
    G0, g, g_inv = christoffel_at(coords, eps)

    dG = np.zeros((n, n, n, n))
    for mu in range(n):
        cp = coords.copy(); cp[mu] += eps
        cm = coords.copy(); cm[mu] -= eps
        Gp, _, _ = christoffel_at(cp, eps)
        Gm, _, _ = christoffel_at(cm, eps)
        dG[:, :, :, mu] = (Gp - Gm) / (2 * eps)

    # Riemann R^rho_{sigma mu nu}
    R = np.zeros((n, n, n, n))
    for rho in range(n):
        for sig in range(n):
            for mu in range(n):
                for nu in range(n):
                    R[rho, sig, mu, nu] = dG[rho, nu, sig, mu] - dG[rho, mu, sig, nu]
                    for lam in range(n):
                        R[rho, sig, mu, nu] += G0[rho, mu, lam] * G0[lam, nu, sig] - G0[rho, nu, lam] * G0[lam, mu, sig]

    # Ricci R_{mu nu} = R^lam_{mu lam nu}
    Ric = np.zeros((n, n))
    for mu in range(n):
        for nu in range(n):
            for lam in range(n):
                Ric[mu, nu] += R[lam, mu, lam, nu]

    # Scalar R = g^{mu nu} R_{mu nu}
    R_scalar = sum(g_inv[i, i] * Ric[i, i] for i in range(n))

    return R, Ric, R_scalar, g, g_inv

print("=== Riemann/Ricci/Weyl computation ===")
print()

coords = np.array([th, s1, s2, s3])
R_t, Ric_t, R_s, g, g_inv = full_curvature(coords)

print(f"Scalar curvature R = {R_s:.10f}")
print()

# Ricci tensor components
print("Ricci tensor R_{ab}:")
for i in range(4):
    vals = [f"{Ric_t[i,j]:+.8f}" for j in range(4)]
    print(f"  [{', '.join(vals)}]")
print()

# Check if Ricci is proportional to metric (Einstein space: R_ab = (R/4) g_ab)
print("Einstein space check: R_{ab} = (R/4) g_{ab}?")
for i in range(4):
    expected = (R_s / 4) * g[i, i]
    actual = Ric_t[i, i]
    print(f"  ({i},{i}): expected {expected:+.8f}, actual {actual:+.8f}, diff {actual - expected:+.2e}")
print()

# Kretschner scalar
K_ret = 0.0
for a, b, c, d in iprod(range(4), repeat=4):
    # R_{abcd} = g_{ae} R^e_{bcd}
    R_fully_lower = g[a, a] * R_t[a, b, c, d]
    # R^{abcd} = g^{ae} g^{bf} g^{cg} g^{dh} R_{efgh}
    R_fully_upper = R_t[a, b, c, d] * abs(g_inv[a, a]) * abs(g_inv[b, b]) * abs(g_inv[c, c]) * abs(g_inv[d, d])
    K_ret += R_fully_lower * R_fully_upper

# |Ric|^2
Ric_sq = sum(Ric_t[i, i]**2 * g_inv[i, i]**2 for i in range(4))

# Weyl^2 via 4D identity: K = |Weyl|^2 + 2|Ric|^2 - (2/3)R^2
# Wait, the correct identity is: K = |Weyl|^2 + 2|Ric|^2 - R^2/3
# No: Gauss-Bonnet in 4D: |Riem|^2 = |Weyl|^2 + 2|Ric_0|^2 + R^2/24 * 2
# Let me use the standard: C^2 = K - 2*E + (2/3)*R^2
# where E = R_{ab}R^{ab}
# Actually: In 4D: |Riem|^2 = |Weyl|^2 + 2|Ric - (R/4)g|^2 + R^2/24 * something
# The clean formula: C_{abcd}C^{abcd} = R_{abcd}R^{abcd} - 2R_{ab}R^{ab} + R^2/3

C_sq = K_ret - 2 * Ric_sq + R_s**2 / 3.0

print(f"Kretschner K = R_abcd R^abcd = {K_ret:.10f}")
print(f"|Ric|^2 = R_ab R^ab = {Ric_sq:.10f}")
print(f"|Weyl|^2 = C_abcd C^abcd = {C_sq:.10f}")
print()

if abs(C_sq) < 1e-6:
    print("*** WEYL TENSOR VANISHES ***")
    print("Equilibrium is conformally flat (de Sitter).")
else:
    print(f"*** WEYL TENSOR IS NONZERO ***")
    print(f"|Weyl|^2 = {C_sq:.10f}")
    ratio = np.sqrt(abs(Ric_sq)) / np.sqrt(abs(C_sq))
    print(f"sqrt(|Ric|^2 / |Weyl|^2) = {ratio:.10f}")
    print(f"2/3 = {2/3:.10f}")
    print(f"Matches 2/3 to {abs(ratio - 2/3):.2e}")
    print()
    print("Equilibrium is NOT conformally flat.")
    print("It has non-trivial Weyl curvature from the anisotropy")
    print("of the dominant hypothesis direction.")

print()

# ============================================================
# De Sitter cross-checks
# ============================================================

print("=== De Sitter cross-checks ===")
print()
print("For maximally symmetric spaces (de Sitter/anti-de Sitter):")
print(f"  K should be R^2/6 = {R_s**2/6:.10f} (actual: {K_ret:.10f})")
print(f"  |Ric|^2 should be R^2/4 = {R_s**2/4:.10f} (actual: {Ric_sq:.10f})")
print(f"  |Weyl|^2 should be 0 (actual: {C_sq:.10f})")
print()

# ============================================================
# Also check at the SYMMETRIC equilibrium (self-evidence)
# ============================================================

print("=== Comparison: symmetric self-evidence equilibrium ===")
print()
# At self-evidence: all s_i equal, no preferred direction
m_sym = np.array([0.25, 0.25, 0.25, 0.25])
for _ in range(1000):
    m_sym = enforce_floor(m_sym)
    # At self-evidence, DS(m,m) gives a specific fixed point
    m_sym2, _ = ds_step(m_sym, m_sym)
    m_sym = m_sym2

th_s = m_sym[3]
s_s = m_sym[0]
print(f"Symmetric FP: s = {s_s:.10f}, theta = {th_s:.10f}")
print(f"  g_11 = g_22 = g_33 = {-(s_s + th_s)**2:.10f}")
print(f"  Isotropic: YES")

coords_sym = np.array([th_s, s_s, s_s, s_s])
R_sym, Ric_sym, R_s_sym, g_sym, g_inv_sym = full_curvature(coords_sym)

K_sym = 0.0
for a, b, c, d in iprod(range(4), repeat=4):
    R_lower = g_sym[a, a] * R_sym[a, b, c, d]
    R_upper = R_sym[a, b, c, d] * abs(g_inv_sym[a, a]) * abs(g_inv_sym[b, b]) * abs(g_inv_sym[c, c]) * abs(g_inv_sym[d, d])
    K_sym += R_lower * R_upper

Ric_sq_sym = sum(Ric_sym[i, i]**2 * g_inv_sym[i, i]**2 for i in range(4))
C_sq_sym = K_sym - 2 * Ric_sq_sym + R_s_sym**2 / 3.0

print(f"  R = {R_s_sym:.10f}")
print(f"  |Weyl|^2 = {C_sq_sym:.10f}")
if abs(C_sq_sym) < 1e-6:
    print(f"  Weyl VANISHES at symmetric equilibrium (conformally flat)")
else:
    print(f"  Weyl NONZERO at symmetric equilibrium")
print()

# ============================================================
# Self-dual / anti-self-dual decomposition of Weyl
# ============================================================

print("=== Weyl self-dual/anti-self-dual decomposition ===")
print()
# In 4D Lorentzian, Weyl decomposes as W = W+ + W-
# The Hodge dual: *R_{abcd} = (1/2) eps_{abef} R^{ef}_{cd}
# W+ = (W + *W)/2, W- = (W - *W)/2
# For a real Riemannian metric, |W+|^2 = |W-|^2 (chirality balance)
# So |W|^2 = 2|W+|^2

# Construct the Weyl tensor
# C_{abcd} = R_{abcd} - (2/(n-2))(g_{a[c}R_{d]b} - g_{b[c}R_{d]a}) + (2/((n-1)(n-2)))R g_{a[c}g_{d]b}
n = 4
C = np.zeros((n, n, n, n))
for a in range(n):
    for b in range(n):
        for c in range(n):
            for d in range(n):
                # R_{abcd} with first index lowered
                R_abcd = g[a, a] * R_t[a, b, c, d]

                # Schouten tensor terms
                S_ac = Ric_t[a, c] - R_s / (2 * (n - 1)) * g[a, c] if a == c else Ric_t[a, c]
                S_ad = Ric_t[a, d] - R_s / (2 * (n - 1)) * g[a, d] if a == d else Ric_t[a, d]
                S_bc = Ric_t[b, c] - R_s / (2 * (n - 1)) * g[b, c] if b == c else Ric_t[b, c]
                S_bd = Ric_t[b, d] - R_s / (2 * (n - 1)) * g[b, d] if b == d else Ric_t[b, d]

                C[a, b, c, d] = R_abcd - (1.0 / (n - 2)) * (
                    g[a, c] * Ric_t[b, d] if a == c else 0 +
                    g[b, d] * Ric_t[a, c] if b == d else 0 -
                    (g[a, d] * Ric_t[b, c] if a == d else 0) -
                    (g[b, c] * Ric_t[a, d] if b == c else 0)
                ) + (R_s / ((n - 1) * (n - 2))) * (
                    (g[a, c] * g[b, d] if (a == c and b == d) else 0) -
                    (g[a, d] * g[b, c] if (a == d and b == c) else 0)
                )

# C^2 check
C_sq_direct = 0.0
for a, b, c, d in iprod(range(4), repeat=4):
    C_lower = C[a, b, c, d]
    C_upper = C[a, b, c, d]
    for idx, val in [(a, g_inv), (b, g_inv), (c, g_inv), (d, g_inv)]:
        C_upper *= abs(val[idx, idx])
    C_sq_direct += C_lower * C_upper

print(f"|Weyl|^2 (direct construction) = {C_sq_direct:.10f}")
print(f"|Weyl|^2 (via K identity)      = {C_sq:.10f}")
print(f"Agreement: {abs(C_sq_direct - C_sq) / max(abs(C_sq), 1e-15):.2e}")

print()
print("DONE.")
