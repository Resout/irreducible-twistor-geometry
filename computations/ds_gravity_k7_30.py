"""
Gravity from DS at EXACT K*=7/30 equilibrium.

Uses the exact evidence and fixed point from spectral_gap_computation.py:
  p_dom = 0.9322756157
  e = (0.6312, 0.1203, 0.1203, 0.1282)
  m* = (0.7869, 0.0293, 0.0293, 0.1545)
  K* = 7/30 exactly
  Born(θ*) = 1/27 exactly
  λ₀ = 0.2829, Δ = 1.263
"""
import numpy as np

H = 3
FLOOR = 1.0 / H**3

# ============================================================
# DS combination (exact)
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

def K_conflict(m, e):
    return sum(m[i] * e[j] for i in range(3) for j in range(3) if i != j)

# ============================================================
# Construct exact K*=7/30 evidence
# ============================================================
from scipy.optimize import brentq

def K_at_eq(p_dom_val):
    p_w = (1 - p_dom_val) / 2
    p_th = FLOOR
    sc = 1 - p_th
    raw = np.array([np.sqrt(p_dom_val*sc), np.sqrt(p_w*sc), np.sqrt(p_w*sc), np.sqrt(p_th)])
    e = raw / np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(1000):
        m_new, _ = ds_combine(m, e.astype(complex))
        if np.max(np.abs(m_new - m)) < 1e-15:
            break
        m = np.real(m_new)
    _, K = ds_combine(m_new, e.astype(complex), apply_floor=False)
    return np.real(K)

p_dom_exact = brentq(lambda p: K_at_eq(p) - 7/30, 0.92, 0.94, xtol=1e-14)

p_w = (1 - p_dom_exact) / 2
sc = 1 - FLOOR
raw = np.array([np.sqrt(p_dom_exact*sc), np.sqrt(p_w*sc), np.sqrt(p_w*sc), np.sqrt(FLOOR)])
e_star = (raw / np.sum(raw)).astype(complex)

# Find exact fixed point
m = np.array([0.4, 0.2, 0.2, 0.2], dtype=complex)
for _ in range(1000):
    m_new, _ = ds_combine(m, e_star)
    if np.max(np.abs(m_new - m)) < 1e-15:
        break
    m = m_new
m_star = m_new

_, K_star_check = ds_combine(m_star, e_star, apply_floor=False)

print("=" * 60)
print("EXACT K*=7/30 EQUILIBRIUM")
print("=" * 60)
print(f"p_dom = {p_dom_exact:.10f}")
print(f"e* = ({np.real(e_star[0]):.10f}, {np.real(e_star[1]):.10f}, {np.real(e_star[2]):.10f}, {np.real(e_star[3]):.10f})")
print(f"m* = ({np.real(m_star[0]):.10f}, {np.real(m_star[1]):.10f}, {np.real(m_star[2]):.10f}, {np.real(m_star[3]):.10f})")
print(f"K* = {np.real(K_star_check):.15f}")
print(f"7/30 = {7/30:.15f}")
print(f"|K*-7/30| = {abs(K_star_check - 7/30):.2e}")
print(f"Born(θ*) = {np.real(born(m_star)):.10f}")
print(f"1/27 = {1/27:.10f}")
print()

# ============================================================
# Anti-holomorphic Jacobian
# ============================================================
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
    tr = np.trace(db) / 4
    sym = (db + db.T) / 2
    asym = (db - db.T) / 2
    stl = sym - tr * np.eye(4)
    frob = np.sqrt(np.sum(np.abs(db)**2))
    return {
        'frob': frob,
        'conf': 4 * np.abs(tr)**2,
        'spin1': np.sum(np.abs(asym)**2),
        'grav': np.sum(np.abs(stl)**2),
        'sv': np.linalg.svd(db, compute_uv=False),
    }

# ============================================================
print("=" * 60)
print("∂̄Φ AT K*=7/30 FIXED POINT")
print("=" * 60)

db = dbar_phi(m_star, e_star)
dec = decompose(db)
print(f"||∂̄Φ|| = {dec['frob']:.8f}")
print(f"Singular values: {dec['sv']}")
tot = dec['conf'] + dec['spin1'] + dec['grav']
if tot > 1e-15:
    print(f"Conformal (spin-0): {dec['conf']/tot*100:.1f}%")
    print(f"Spin-1:             {dec['spin1']/tot*100:.1f}%")
    print(f"Graviton (spin-2):  {dec['grav']/tot*100:.1f}%")
print(f"Born at FP = {np.real(born(m_star)):.8f}, floor = {FLOOR:.8f}")
print(f"Floor marginal: {abs(np.real(born(m_star)) - FLOOR) < 1e-4}")
print()

# ============================================================
print("=" * 60)
print("COMPLEX TRAJECTORY TOWARD K*=7/30")
print("=" * 60)

# Start away from FP with complex components
m_init = m_star + np.array([0.05+0.02j, -0.02+0.01j, -0.01-0.015j, -0.02-0.015j])
m_init = m_init / np.sum(np.abs(m_init))

m = m_init.copy()
print(f"Initial: Born={born(m):.5f}, K={K_conflict(m, e_star):.5f}")
print()
print(f"{'Step':>4} {'Born':>8} {'|K|':>8} {'Im(K)':>10} {'||∂̄||':>8} {'grav':>8} {'conf':>8} {'s1':>6}")
print("-" * 72)

integ = {'conf': 0, 'spin1': 0, 'grav': 0}

for step in range(40):
    b = born(m)
    K_val = K_conflict(m, e_star)
    db = dbar_phi(m, e_star)
    dec = decompose(db)

    integ['conf'] += dec['conf']
    integ['spin1'] += dec['spin1']
    integ['grav'] += dec['grav']

    if step < 20 or step % 5 == 0:
        print(f"{step:4d} {np.real(b):8.5f} {np.abs(K_val):8.5f} {np.imag(K_val):10.6f} "
              f"{dec['frob']:8.4f} {np.sqrt(dec['grav']):8.4f} {np.sqrt(dec['conf']):8.4f} "
              f"{np.abs(dec['sv'][1]/dec['sv'][0]) if dec['sv'][0]>1e-12 else 0:6.3f}")

    m_new, _ = ds_combine(m, e_star)
    m = m_new

tot = integ['conf'] + integ['spin1'] + integ['grav']
print()
print("INTEGRATED TRAJECTORY DECOMPOSITION:")
print(f"  Conformal (spin-0): {integ['conf']/tot*100:.1f}%")
print(f"  Spin-1:             {integ['spin1']/tot*100:.1f}%")
print(f"  Graviton (spin-2):  {integ['grav']/tot*100:.1f}%")
print()

# ============================================================
print("=" * 60)
print("RANK-2 COUPLED DS AT K*=7/30")
print("=" * 60)

# Two sites near equilibrium, complex, coupled
mA = m_star + np.array([0.03+0.01j, -0.01+0.005j, -0.005-0.01j, -0.015-0.005j])
mA = mA / np.sum(np.abs(mA))
mB = m_star + np.array([-0.02+0.015j, 0.008-0.003j, 0.005+0.002j, 0.007-0.014j])
mB = mB / np.sum(np.abs(mB))

print(f"Initial Born_A={born(mA):.5f}, Born_B={born(mB):.5f}")
print()
print(f"{'Step':>4} {'Born_A':>8} {'Born_B':>8} {'|K_AB|':>8} {'||∂̄||':>8} {'grav':>8} {'σ2/σ1':>6} {'Im(K)':>10}")
print("-" * 72)

integ2 = {'conf': 0, 'spin1': 0, 'grav': 0}
n_fl = 0

for step in range(60):
    bA, bB = born(mA), born(mB)
    K_AB = K_conflict(mA, mB)
    db = dbar_phi(mA, mB)
    dec = decompose(db)

    integ2['conf'] += dec['conf']
    integ2['spin1'] += dec['spin1']
    integ2['grav'] += dec['grav']
    if np.real(bA) <= FLOOR + 1e-6:
        n_fl += 1

    sv_ratio = np.abs(dec['sv'][1]/dec['sv'][0]) if dec['sv'][0] > 1e-12 else 0

    if step < 20 or step % 10 == 0:
        print(f"{step:4d} {np.real(bA):8.5f} {np.real(bB):8.5f} {np.abs(K_AB):8.5f} "
              f"{dec['frob']:8.4f} {np.sqrt(dec['grav']):8.4f} {sv_ratio:6.3f} {np.imag(K_AB):10.6f}")

    mA_new, _ = ds_combine(mA, mB)
    mB_new, _ = ds_combine(mB, mA)
    mA, mB = mA_new, mB_new

tot2 = integ2['conf'] + integ2['spin1'] + integ2['grav']
print()
print(f"Floor active/marginal: {n_fl}/60 steps")
print(f"Final Born_A={born(mA):.6f}, Born_B={born(mB):.6f}")
print(f"Final K_AB = {K_conflict(mA, mB):.8f}")
print()
print("INTEGRATED RANK-2 DECOMPOSITION:")
print(f"  Conformal (spin-0): {integ2['conf']/tot2*100:.1f}%")
print(f"  Spin-1:             {integ2['spin1']/tot2*100:.1f}%")
print(f"  Graviton (spin-2):  {integ2['grav']/tot2*100:.1f}%")

# ============================================================
print("\n" + "=" * 60)
print("EMERGENT METRIC AT K*=7/30 EQUILIBRIUM")
print("=" * 60)

th = np.real(m_star[3])
s1 = np.real(m_star[0])
s2 = np.real(m_star[1])
g00 = th**2
g11 = (s1 + th)**2
g22 = (s2 + th)**2

print(f"θ* = {th:.6f}")
print(f"s_dom = {s1:.6f}, s_weak = {s2:.6f}")
print(f"g₀₀ = θ² = {g00:.6f}")
print(f"|g₁₁| = (s_dom+θ)² = {g11:.6f}")
print(f"|g₂₂| = |g₃₃| = (s_weak+θ)² = {g22:.6f}")
print(f"g₀₀/|g₁₁| = {g00/g11:.6f}")
print(f"g₀₀/|g₂₂| = {g00/g22:.6f}")
print()
print(f"Metric: diag({g00:.4f}, -{g11:.4f}, -{g22:.4f}, -{g22:.4f})")
print(f"Anisotropic! Dominant direction has different scale than weak.")
print()

# Key ratios
print("Key ratios:")
print(f"  g₀₀/|g₂₂| = {g00/g22:.6f}")
print(f"  (H-1)/H = {2/3:.6f}")
print(f"  K* = {7/30:.6f}")
print(f"  g₀₀/(g₁₁+2g₂₂) = {g00/(g11+2*g22):.6f}")
print(f"  g₀₀/mean(g_spatial) = {g00/((g11+2*g22)/3):.6f}")
