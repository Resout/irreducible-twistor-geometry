"""
Gravity from DS: self-consistent equilibrium.

K*=7/30 comes from self-evidence: m combines with m.
The fixed point is m* where DS(m*, m*) = m* and K(m*,m*) = 7/30.
"""
import numpy as np

H = 3
FLOOR = 1.0 / H**3

def ds_combine(m, e, apply_floor=True):
    """DS combination: m=(s1,s2,s3,theta), e=(se1,se2,se3,theta_e)."""
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
    phase_theta = theta / np.abs(theta) if np.abs(theta) > 1e-15 else 1.0
    phases_s = np.array([si / np.abs(si) if np.abs(si) > 1e-15 else 1.0 for si in s])
    abs_s = np.abs(s)

    lo, hi = np.abs(theta), 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        s_sum = np.sum(abs_s)
        sc = (1.0 - mid) / s_sum if s_sum > 0 else 0
        st = abs_s * sc
        b = mid**2 / (np.sum(st**2) + mid**2)
        if b < FLOOR:
            lo = mid
        else:
            hi = mid

    tn = (lo + hi) / 2
    s_sum = np.sum(abs_s)
    sc = (1.0 - tn) / s_sum if s_sum > 0 else 0
    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = phases_s * abs_s * sc
    m_out[3] = phase_theta * tn
    return m_out

def born(m):
    return np.abs(m[3])**2 / np.sum(np.abs(m)**2)

def K_conflict(m, e):
    return sum(m[i] * e[j] for i in range(3) for j in range(3) if i != j)

def dbar_phi(m, e, eps=1e-7):
    """∂̄Φ where Φ(m) = DS(m, e) with floor."""
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

def decompose_dbar(db):
    tr = np.trace(db) / 4
    sym = (db + db.T) / 2
    antisym = (db - db.T) / 2
    sym_tl = sym - tr * np.eye(4)
    frob = np.sqrt(np.sum(np.abs(db)**2))
    return {
        'frobenius': frob,
        'conf': 4 * np.abs(tr)**2,
        'spin1': np.sum(np.abs(antisym)**2),
        'grav': np.sum(np.abs(sym_tl)**2),
        'sv': np.linalg.svd(db, compute_uv=False),
    }

# ============================================================
print("=" * 60)
print("SELF-CONSISTENT FIXED POINT: DS(m*,m*) = m*")
print("=" * 60)

# Iterate: m_{n+1} = DS(m_n, m_n)
m = np.array([0.4, 0.2, 0.2, 0.2], dtype=complex)
for i in range(500):
    m_new, K = ds_combine(m, m)
    diff = np.max(np.abs(m_new - m))
    if i < 5 or diff < 1e-14:
        print(f"  Step {i}: m=({np.real(m_new[0]):.8f}, {np.real(m_new[1]):.8f}, "
              f"{np.real(m_new[2]):.8f}, {np.real(m_new[3]):.8f}), K={np.real(K):.8f}")
    if diff < 1e-15:
        print(f"  Converged at step {i}")
        break
    m = m_new

m_star = m_new
K_star = np.real(K_conflict(m_star, m_star))
print(f"\nm* = ({np.real(m_star[0]):.10f}, {np.real(m_star[1]):.10f}, "
      f"{np.real(m_star[2]):.10f}, {np.real(m_star[3]):.10f})")
print(f"K* = {K_star:.10f}")
print(f"7/30 = {7/30:.10f}")
print(f"|K* - 7/30| = {abs(K_star - 7/30):.2e}")
print(f"Born(θ*) = {np.real(born(m_star)):.10f}")
print(f"1/27 = {1/27:.10f}")

# The paper's values
print(f"\nPaper values: θ*=0.5959, s*=0.1347")
print(f"Our values:   θ*={np.real(m_star[3]):.4f}, s_dom={np.real(m_star[0]):.4f}, s_weak={np.real(m_star[1]):.4f}")
print()

# ============================================================
print("=" * 60)
print("SYMMETRIC SELF-EVIDENCE (s1=s2=s3)")
print("=" * 60)

# Try symmetric initial condition
m_sym = np.array([0.2, 0.2, 0.2, 0.4], dtype=complex)
for i in range(500):
    m_new, K = ds_combine(m_sym, m_sym)
    diff = np.max(np.abs(m_new - m_sym))
    if diff < 1e-15:
        print(f"  Converged at step {i}")
        break
    m_sym = m_new

m_star_sym = m_new
K_sym = np.real(K_conflict(m_star_sym, m_star_sym))
print(f"m* = ({np.real(m_star_sym[0]):.10f}, {np.real(m_star_sym[1]):.10f}, "
      f"{np.real(m_star_sym[2]):.10f}, {np.real(m_star_sym[3]):.10f})")
print(f"K* = {K_sym:.10f}")
print(f"7/30 = {7/30:.10f}")
print(f"|K* - 7/30| = {abs(K_sym - 7/30):.2e}")
print(f"Born(θ*) = {np.real(born(m_star_sym)):.10f}")
print()

# Check: m* = (s, s, s, θ) with s+s+s+θ=1 => θ=1-3s
# K(m,m) = s*s + s*s + s*s + s*s + s*s + s*s = 6s²
# At H=3 conservation law: K*(H²+1) - η*H² = 1 where η=(H²-H+1)/(H(H²+1))
# K* = 7/30 => 6s² = 7/30 => s² = 7/180 => s = √(7/180) = 0.19720
# θ = 1 - 3*√(7/180) = 1 - 0.59161 = 0.40839
print("Analytical check:")
s_anal = np.sqrt(7/180)
theta_anal = 1 - 3*s_anal
K_anal = 6 * s_anal**2
print(f"  s = √(7/180) = {s_anal:.10f}")
print(f"  θ = 1 - 3s = {theta_anal:.10f}")
print(f"  K = 6s² = {K_anal:.10f}")
print(f"  7/30 = {7/30:.10f}")
born_anal = theta_anal**2 / (3*s_anal**2 + theta_anal**2)
print(f"  Born(θ) = {born_anal:.10f}")
print(f"  1/27 = {1/27:.10f}")
print()

# So K*=7/30 is satisfied by s=√(7/180), θ=1-3s
# But Born(θ) = 0.3048... >> 1/27 => floor is NOT active here!
# The floor doesn't activate at the symmetric self-evidence fixed point.

# Let's check: does the symmetric iteration converge to this?
m_test = np.array([s_anal, s_anal, s_anal, theta_anal], dtype=complex)
m_out, K_test = ds_combine(m_test, m_test, apply_floor=True)
print(f"DS(m_anal, m_anal):")
print(f"  Input:  ({np.real(m_test[0]):.8f}, {np.real(m_test[1]):.8f}, {np.real(m_test[2]):.8f}, {np.real(m_test[3]):.8f})")
print(f"  Output: ({np.real(m_out[0]):.8f}, {np.real(m_out[1]):.8f}, {np.real(m_out[2]):.8f}, {np.real(m_out[3]):.8f})")
print(f"  K = {np.real(K_test):.10f}")
print(f"  Is fixed point? {np.max(np.abs(m_out - m_test)) < 1e-10}")
print()

# ============================================================
# If it's NOT a fixed point, find the actual symmetric FP
# ============================================================
print("=" * 60)
print("ACTUAL SYMMETRIC SELF-EVIDENCE FIXED POINT")
print("=" * 60)

# Force symmetry: m = (s, s, s, 1-3s)
from scipy.optimize import brentq

def sym_fp_eq(s_val):
    """DS(m,m) = m for m=(s,s,s,1-3s). Return s_out - s_in."""
    if s_val <= 0 or s_val >= 1/3:
        return 1.0
    theta_val = 1 - 3*s_val
    m_in = np.array([s_val, s_val, s_val, theta_val], dtype=complex)
    m_out, _ = ds_combine(m_in, m_in, apply_floor=True)
    return np.real(m_out[0]) - s_val

# Scan for sign change
print("Scanning for fixed point:")
for s_try in np.linspace(0.01, 0.32, 20):
    val = sym_fp_eq(s_try)
    print(f"  s={s_try:.3f}: f(s) = {val:+.6f}")

# Find root
try:
    s_fp = brentq(sym_fp_eq, 0.01, 0.32)
    theta_fp = 1 - 3*s_fp
    m_fp = np.array([s_fp, s_fp, s_fp, theta_fp], dtype=complex)
    K_fp = np.real(K_conflict(m_fp, m_fp))

    print(f"\nFixed point found:")
    print(f"  s* = {s_fp:.10f}")
    print(f"  θ* = {theta_fp:.10f}")
    print(f"  K* = {K_fp:.10f}")
    print(f"  7/30 = {7/30:.10f}")
    print(f"  |K* - 7/30| = {abs(K_fp - 7/30):.2e}")
    print(f"  Born(θ*) = {np.real(born(m_fp)):.10f}")
    print(f"  1/27 = {1/27:.10f}")

    # Verify
    m_check, K_check = ds_combine(m_fp, m_fp)
    print(f"  DS(m*,m*) - m* max diff: {np.max(np.abs(m_check - m_fp)):.2e}")
    print()

    # ============================================================
    print("=" * 60)
    print("GRAVITY AT THE CORRECT FIXED POINT")
    print("=" * 60)

    # ∂̄Φ at fixed point
    db_fp = dbar_phi(m_fp, m_fp)
    dec_fp = decompose_dbar(db_fp)
    print(f"||∂̄Φ|| at m* = {dec_fp['frobenius']:.6f}")
    print(f"Singular values: {dec_fp['sv']}")
    print(f"Graviton: {np.sqrt(dec_fp['grav']):.6f}")
    print(f"Conformal: {np.sqrt(dec_fp['conf']):.6f}")
    print(f"Spin-1: {np.sqrt(dec_fp['spin1']):.6f}")
    floor_active = np.real(born(m_fp)) < FLOOR + 1e-10
    print(f"Born at FP: {np.real(born(m_fp)):.8f}, floor: {FLOOR:.8f}")
    print(f"Floor marginal: {abs(np.real(born(m_fp)) - FLOOR) < 1e-4}")
    print()

    # ============================================================
    print("=" * 60)
    print("COUPLED RANK-2 AT CORRECT EQUILIBRIUM")
    print("=" * 60)

    # Start near FP with complex perturbation
    mA = m_fp + np.array([0.01+0.005j, -0.005+0.003j, -0.003-0.002j, -0.002-0.006j])
    mA = mA / np.sum(np.abs(mA))
    mB = m_fp + np.array([-0.008+0.004j, 0.003-0.002j, 0.004+0.001j, 0.001-0.003j])
    mB = mB / np.sum(np.abs(mB))

    print(f"Initial Born_A={born(mA):.5f}, Born_B={born(mB):.5f}")
    print()

    print(f"{'Step':>4} {'Born_A':>8} {'|K_AB|':>8} {'||∂̄||':>8} {'grav':>8} {'conf':>8} {'Im(K)':>10}")
    print("-" * 62)

    integ = {'conf': 0, 'spin1': 0, 'grav': 0}
    n_fl = 0

    for step in range(60):
        bA = born(mA)
        K_AB = K_conflict(mA, mB)
        db = dbar_phi(mA, mB)
        dec = decompose_dbar(db)

        integ['conf'] += dec['conf']
        integ['spin1'] += dec['spin1']
        integ['grav'] += dec['grav']
        if np.real(bA) < FLOOR + 1e-6:
            n_fl += 1

        if step < 15 or step % 10 == 0:
            print(f"{step:4d} {np.real(bA):8.5f} {np.abs(K_AB):8.5f} "
                  f"{dec['frobenius']:8.4f} {np.sqrt(dec['grav']):8.4f} "
                  f"{np.sqrt(dec['conf']):8.4f} {np.imag(K_AB):10.6f}")

        mA_new, _ = ds_combine(mA, mB)
        mB_new, _ = ds_combine(mB, mA)
        mA, mB = mA_new, mB_new

    tot = integ['conf'] + integ['spin1'] + integ['grav']
    print()
    print(f"Floor active/marginal: {n_fl}/60 steps")
    print(f"Final K_AB = {K_conflict(mA,mB):.8f}")
    print()
    if tot > 1e-15:
        print("INTEGRATED DECOMPOSITION:")
        print(f"  Conformal (spin-0): {integ['conf']/tot*100:.1f}%")
        print(f"  Spin-1:             {integ['spin1']/tot*100:.1f}%")
        print(f"  Graviton (spin-2):  {integ['grav']/tot*100:.1f}%")

    # ============================================================
    print("\n" + "=" * 60)
    print("THE EMERGENT METRIC AT EQUILIBRIUM")
    print("=" * 60)

    # g = diag(θ², -(s₁+θ)², -(s₂+θ)², -(s₃+θ)²)
    th = np.real(theta_fp)
    sv = np.real(s_fp)
    g00 = th**2
    g11 = -(sv + th)**2
    print(f"θ* = {th:.6f}, s* = {sv:.6f}")
    print(f"g₀₀ = θ² = {g00:.6f}")
    print(f"g₁₁ = g₂₂ = g₃₃ = -(s+θ)² = {g11:.6f}")
    print(f"g₀₀/|g₁₁| = {g00/abs(g11):.6f}")
    print(f"(H-1)/H = {2/3:.6f}")
    print(f"Match: {abs(g00/abs(g11) - 2/3) < 0.01}")

except Exception as ex:
    print(f"Error: {ex}")
    import traceback
    traceback.print_exc()
