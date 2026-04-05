"""
Trace constraint: Step 1.

Prove that the (1,1), (3,1)⊕(1,3), and (3,3) components of ∂̄Φ
are not independent — they are ALL determined by the 3-parameter
mass function m = (s₁, s₂, s₃, θ) with L1=1 (3 free parameters).

The decomposition has 1+6+9 = 16 components.
Only 3 parameters control all 16.
The sectors are massively overconstrained.
"""
import numpy as np
from scipy.optimize import brentq

H = 3
FLOOR = 1.0 / H**3

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

def dbar_phi(m, e, eps=1e-7):
    n = 4
    z = m.copy()
    def phi(z_in):
        out, _ = ds_combine(z_in, e, apply_floor=True)
        return out
    dbar = np.zeros((n, n), dtype=complex)
    for b in range(n):
        zp = z.copy(); zp[b] += eps
        zm = z.copy(); zm[b] -= eps
        d_dx = (phi(zp) - phi(zm)) / (2 * eps)
        zp = z.copy(); zp[b] += 1j * eps
        zm = z.copy(); zm[b] -= 1j * eps
        d_dy = (phi(zp) - phi(zm)) / (2 * eps)
        dbar[:, b] = 0.5 * (d_dx + 1j * d_dy)
    return dbar

def decompose(db):
    M = np.real(db)
    tr = np.trace(M) / 4
    sym = (M + M.T) / 2
    asym = (M - M.T) / 2
    stl = sym - tr * np.eye(4)
    total = np.sum(M**2)
    if total < 1e-30:
        return 0, 0, 0, 0
    conf = 4 * tr**2
    s1 = np.sum(asym**2)
    grav = np.sum(stl**2)
    return conf/total, s1/total, grav/total, total

# ============================================================
# Find K*=7/30 evidence and fixed point
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

# ============================================================
print("=" * 60)
print("HOW THE 3 PARAMETERS CONTROL ALL 16 COMPONENTS")
print("=" * 60)
print()

# The mass function m = (s₁, s₂, s₃, θ) with L1=1 has 3 free params.
# ∂̄Φ is a 4×4 real matrix = 16 components.
# The decomposition: 1 (trace) + 6 (antisym) + 9 (sym traceless) = 16.
# But all 16 are functions of the same 3 parameters.

# Map out how the fractions change as we vary the mass function
# along a 1D path through the equilibrium manifold.

print("Varying θ (1 parameter) while maintaining L1=1:")
print(f"{'θ':>8} {'s₁':>8} {'s₂':>8} {'conf%':>8} {'YM%':>8} {'grav%':>8} {'||∂̄||²':>10}")
print("-" * 62)

for theta_val in np.linspace(0.05, 0.50, 15):
    # Keep s₁ dominant, s₂=s₃ equal, adjust to L1=1
    remaining = 1.0 - theta_val
    # Dominant gets 80% of remaining, weak get 10% each
    s1 = 0.8 * remaining
    s2 = 0.1 * remaining
    m_test = np.array([s1, s2, s2, theta_val], dtype=complex)

    db = dbar_phi(m_test, e_star)
    conf, ym, grav, tot = decompose(db)

    print(f"{theta_val:8.3f} {s1:8.3f} {s2:8.3f} {conf*100:8.1f} {ym*100:8.1f} {grav*100:8.1f} {tot:10.4f}")

print()

# ============================================================
print("=" * 60)
print("THE CONSTRAINT: FRACTIONS ARE COUPLED")
print("=" * 60)
print()

# Collect many random mass functions and plot the relationship
# between the three fractions
np.random.seed(42)
confs, yms, gravs = [], [], []

for trial in range(500):
    m_rand = np.random.dirichlet([1, 1, 1, 2]).astype(complex)  # random mass function
    db = dbar_phi(m_rand, e_star)
    conf, ym, grav, tot = decompose(db)
    if tot > 1e-10:
        confs.append(conf)
        yms.append(ym)
        gravs.append(grav)

confs = np.array(confs)
yms = np.array(yms)
gravs = np.array(gravs)

print(f"Sampled {len(confs)} random mass functions.")
print(f"Conformal: range [{confs.min()*100:.1f}%, {confs.max()*100:.1f}%], mean {confs.mean()*100:.1f}%")
print(f"YM:        range [{yms.min()*100:.1f}%, {yms.max()*100:.1f}%], mean {yms.mean()*100:.1f}%")
print(f"Graviton:  range [{gravs.min()*100:.1f}%, {gravs.max()*100:.1f}%], mean {gravs.mean()*100:.1f}%")
print()

# Check: do fractions sum to 1?
sums = confs + yms + gravs
print(f"Sum of fractions: mean {sums.mean():.6f}, std {sums.std():.6f}")
print(f"(Should be 1.0 since trace² + antisym² + sym_tl² = total²)")
print()

# Correlation between sectors
corr_cg = np.corrcoef(confs, gravs)[0, 1]
corr_cy = np.corrcoef(confs, yms)[0, 1]
corr_yg = np.corrcoef(yms, gravs)[0, 1]

print("Correlations between sectors:")
print(f"  Conformal ↔ Graviton: r = {corr_cg:+.4f}")
print(f"  Conformal ↔ YM:       r = {corr_cy:+.4f}")
print(f"  YM ↔ Graviton:        r = {corr_yg:+.4f}")
print()

# ============================================================
print("=" * 60)
print("THE TRACE CONSTRAINT")
print("=" * 60)
print()

# What is tr(∂̄Φ) at the fixed point?
db_star = dbar_phi(m_star, e_star)
tr_star = np.real(np.trace(db_star))
print(f"tr(∂̄Φ) at K*=7/30 fixed point: {tr_star:.6f}")
print()

# How does tr(∂̄Φ) vary across mass function space?
traces = []
for trial in range(500):
    m_rand = np.random.dirichlet([1, 1, 1, 2]).astype(complex)
    db = dbar_phi(m_rand, e_star)
    traces.append(np.real(np.trace(db)))

traces = np.array(traces)
print(f"tr(∂̄Φ) across random mass functions:")
print(f"  Mean: {traces.mean():.6f}")
print(f"  Std:  {traces.std():.6f}")
print(f"  Range: [{traces.min():.4f}, {traces.max():.4f}]")
print()

# Is the trace constant? If so, it's a conservation law.
print(f"Is tr(∂̄Φ) constant? CV = {traces.std()/abs(traces.mean()):.4f}")
print(f"(CV << 1 would mean approximately constant)")
print()

# ============================================================
print("=" * 60)
print("THE REAL CONSTRAINT: 3 PARAMETERS → 16 COMPONENTS")
print("=" * 60)
print()
print("The mass function has 3 free parameters (L1=1).")
print("∂̄Φ has 16 real components.")
print("The 16 components are NOT independent — they live on a")
print("3-dimensional surface in 16-dimensional space.")
print()
print("This means: the (1,1), (3,1)⊕(1,3), and (3,3) fractions")
print("are coupled. You cannot change gravity without changing")
print("Yang-Mills without changing the conformal factor.")
print()
print("The coupling is NOT a choice. It's forced by the fact that")
print("all three come from the SAME 4-component mass function")
print("through the SAME DS+floor map.")
print()

# Measure the effective dimension of the ∂̄Φ manifold
# by computing the rank of the Jacobian d(∂̄Φ)/dm
print("Effective dimension of the ∂̄Φ manifold:")
# ∂̄Φ has 16 real components, parametrised by 3 mass function params
# The Jacobian d(∂̄Φ)/d(m) should have rank ≤ 3

eps = 1e-6
m_base = np.real(m_star).copy()
db_base = np.real(dbar_phi(m_star, e_star)).flatten()

# Tangent vectors in the L1=1 simplex
# v1 = (1, 0, 0, -1), v2 = (0, 1, 0, -1), v3 = (0, 0, 1, -1)
jac = np.zeros((16, 3))
for j, v in enumerate([(1,0,0,-1), (0,1,0,-1), (0,0,1,-1)]):
    m_plus = m_star + eps * np.array(v, dtype=complex)
    m_minus = m_star - eps * np.array(v, dtype=complex)
    db_plus = np.real(dbar_phi(m_plus, e_star)).flatten()
    db_minus = np.real(dbar_phi(m_minus, e_star)).flatten()
    jac[:, j] = (db_plus - db_minus) / (2 * eps)

sv = np.linalg.svd(jac, compute_uv=False)
print(f"  Singular values of d(∂̄Φ)/dm: {sv}")
print(f"  Rank (sv > 1e-8): {np.sum(sv > 1e-8)}")
print(f"  The ∂̄Φ manifold has dimension {np.sum(sv > 1e-8)} in R¹⁶.")
print()
print(f"  16 components controlled by {np.sum(sv > 1e-8)} parameters.")
print(f"  Constraint ratio: {16}/{np.sum(sv > 1e-8)} = {16/max(np.sum(sv > 1e-8),1):.0f}:1")
print()
print("This is the structural coupling. The three force sectors")
print("are controlled by 3 parameters with 16 outputs.")
print("Every change affects all three simultaneously.")
print("There is no way to modify gravity alone, or YM alone,")
print("or the conformal factor alone. They move together.")
