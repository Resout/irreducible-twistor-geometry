"""
Gravity from DS: using the REAL DS combination.

Follow the framework's own requirements:
- Complex mass functions
- Entangled (rank-2 coupled)
- Born floor active during transients
- On CP³

We track the anti-holomorphic Jacobian along the trajectory
to see where and how gravity is generated.
"""
import numpy as np

H = 3
FLOOR = 1.0 / H**3  # 1/27

# ============================================================
# Real DS combination (from spectral_gap_computation.py)
# Extended to complex mass functions
# ============================================================
def ds_combine(m, e, apply_floor=True):
    """DS combination for complex mass functions.
    m = (s1, s2, s3, theta), e = (se1, se2, se3, theta_e).
    """
    s = m[:3]
    theta = m[3]
    se = e[:3]
    theta_e = e[3]

    # Pre-normalisation
    s_new = s * se + s * theta_e + theta * se
    theta_new = theta * theta_e

    # Conflict
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)

    # Normalise
    denom = 1.0 - K
    if abs(denom) < 1e-15:
        return m.copy(), K
    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom

    # L1 normalise
    total = np.sum(np.abs(m_out))
    if total > 0:
        m_out = m_out / total

    # Born floor enforcement
    if apply_floor:
        born = np.abs(m_out[3])**2 / np.sum(np.abs(m_out)**2)
        if born < FLOOR:
            m_out = enforce_floor_complex(m_out)

    return m_out, K

def enforce_floor_complex(m):
    """Enforce Born(theta) >= 1/27 for complex m, maintaining L1=1."""
    s = m[:3].copy()
    theta = m[3]

    # Keep phases, adjust magnitudes
    phase_theta = theta / np.abs(theta) if np.abs(theta) > 1e-15 else 1.0
    phases_s = np.array([si / np.abs(si) if np.abs(si) > 1e-15 else 1.0 for si in s])
    abs_s = np.abs(s)

    lo, hi = np.abs(theta), 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        s_scale = (1.0 - mid) / np.sum(abs_s) if np.sum(abs_s) > 0 else 0
        s_trial = abs_s * s_scale
        born = mid**2 / (np.sum(s_trial**2) + mid**2)
        if born < FLOOR:
            lo = mid
        else:
            hi = mid

    theta_new_abs = (lo + hi) / 2
    s_scale = (1.0 - theta_new_abs) / np.sum(abs_s) if np.sum(abs_s) > 0 else 0

    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = phases_s * abs_s * s_scale
    m_out[3] = phase_theta * theta_new_abs
    return m_out

def born(m):
    return np.abs(m[3])**2 / np.sum(np.abs(m)**2)

def conflict(m, e):
    s, se = m[:3], e[:3]
    return sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)

# ============================================================
# Anti-holomorphic Jacobian
# ============================================================
def dbar_phi(m, e, eps=1e-7):
    """Compute dbar of one DS step: partial Phi / partial m-bar.
    Phi(m) = DS(m, e) with floor.
    dbar = (1/2)(d/dx + i d/dy) for z = x + iy.
    """
    n = 4
    z = m.copy()

    def phi(z_in):
        out, _ = ds_combine(z_in, e, apply_floor=True)
        return out

    dbar = np.zeros((n, n), dtype=complex)
    for b in range(n):
        # d/dx_b
        zp = z.copy(); zp[b] += eps
        zm = z.copy(); zm[b] -= eps
        d_dx = (phi(zp) - phi(zm)) / (2 * eps)

        # d/dy_b
        zp = z.copy(); zp[b] += 1j * eps
        zm = z.copy(); zm[b] -= 1j * eps
        d_dy = (phi(zp) - phi(zm)) / (2 * eps)

        dbar[:, b] = 0.5 * (d_dx + 1j * d_dy)

    return dbar


# ============================================================
# Build equilibrium evidence
# ============================================================
# At equilibrium: p_dom = 9/10, p_weak = 1/20 each
# Evidence that sustains K* = 7/30
def make_equilibrium_evidence():
    """Construct evidence at K*=7/30 equilibrium."""
    # From the paper: at equilibrium, the dominant hypothesis has
    # Born probability 9/10, weak ones have 1/20 each.
    # The evidence mass function has the same structure.
    p_dom = 0.932  # from spectral_gap_computation.py Brent search
    p_weak = (1 - p_dom) / 2

    # Convert Born probs to mass function
    # m = (s1, s2, s3, theta) where s1 is dominant
    # Simple: evidence strongly supports hypothesis 1
    e = np.array([p_dom, p_weak, p_weak, 0.0], dtype=complex)
    # Normalize L1
    e = e / np.sum(np.abs(e))
    return e


# ============================================================
# MAIN COMPUTATION
# ============================================================
print("=" * 60)
print("DS GRAVITY: Real DS on CP³")
print("=" * 60)

# Evidence at equilibrium
e_eq = make_equilibrium_evidence()
print(f"Evidence: {np.real(e_eq)}")
print(f"Evidence Born: {born(e_eq):.6f}")
print()

# ============================================================
# PART 1: Find the real fixed point
# ============================================================
print("--- Finding fixed point ---")
m = np.array([0.25, 0.25, 0.25, 0.25], dtype=complex)
for i in range(200):
    m_new, K = ds_combine(m, e_eq)
    if np.max(np.abs(m_new - m)) < 1e-14:
        print(f"Converged at step {i}")
        break
    m = m_new

m_star = m
print(f"m* = {np.real(m_star)}")
print(f"K* = {np.real(conflict(m_star, e_eq)):.10f}")
print(f"7/30 = {7/30:.10f}")
print(f"Born(m*) = {born(m_star):.10f}")
print(f"1/27 = {1/27:.10f}")
print()

# ============================================================
# PART 2: dbar at the fixed point (should be ~0)
# ============================================================
print("--- dbar at fixed point ---")
db_star = dbar_phi(m_star, e_eq)
print(f"||dbar|| at m* = {np.sqrt(np.sum(np.abs(db_star)**2)):.2e}")
print(f"(Should be ~0: floor marginally saturated at fixed point)")
print()

# ============================================================
# PART 3: Track dbar along trajectory FROM BELOW
# ============================================================
# Start with low theta (floor will activate)
print("--- Trajectory from floor-active initial condition ---")
m_init = np.array([0.05, 0.35, 0.30, 0.30], dtype=complex)
print(f"Initial: m = {np.real(m_init)}")
print(f"Initial Born = {born(m_init):.6f} (floor = {FLOOR:.6f})")
print(f"Floor active: {born(m_init) < FLOOR}")
print()

m = m_init.copy()
print(f"{'Step':>4} {'Born':>8} {'|K|':>8} {'||dbar||':>10} {'trace':>10} {'graviton':>10} {'floor?':>6}")
print("-" * 62)

for step in range(40):
    b = born(m)
    K_val = conflict(m, e_eq)

    # Compute dbar BEFORE the DS step (what the floor does to current state)
    db = dbar_phi(m, e_eq)
    frob = np.sqrt(np.sum(np.abs(db)**2))
    tr = np.trace(db) / 4
    sym = (db + db.T) / 2
    sym_tl = sym - tr * np.eye(4)
    grav_norm = np.sqrt(np.sum(np.abs(sym_tl)**2))

    floor_active = "YES" if b < FLOOR else "no"

    if step < 20 or step % 5 == 0:
        print(f"{step:4d} {np.abs(b):8.5f} {np.abs(K_val):8.5f} "
              f"{frob:10.6f} {np.abs(tr):10.6f} {grav_norm:10.6f} {floor_active:>6}")

    m_new, _ = ds_combine(m, e_eq)
    m = m_new

print()

# ============================================================
# PART 4: Complex mass functions (both chiralities)
# ============================================================
print("--- Complex trajectory (both chiralities) ---")
m_c = np.array([0.05+0.02j, 0.35+0.01j, 0.30-0.015j, 0.30+0.005j], dtype=complex)
m_c = m_c / np.sum(np.abs(m_c))
print(f"Initial Born = {born(m_c):.6f}")
print()

print(f"{'Step':>4} {'Born':>8} {'|K|':>8} {'||dbar||':>10} {'Im(K)':>10} {'graviton':>10} {'floor?':>6}")
print("-" * 68)

for step in range(40):
    b = born(m_c)
    K_val = conflict(m_c, e_eq)

    db = dbar_phi(m_c, e_eq)
    frob = np.sqrt(np.sum(np.abs(db)**2))
    tr = np.trace(db) / 4
    sym = (db + db.T) / 2
    sym_tl = sym - tr * np.eye(4)
    grav_norm = np.sqrt(np.sum(np.abs(sym_tl)**2))

    floor_active = "YES" if np.real(b) < FLOOR else "no"

    if step < 20 or step % 5 == 0:
        print(f"{step:4d} {np.abs(b):8.5f} {np.abs(K_val):8.5f} "
              f"{frob:10.6f} {np.imag(K_val):10.6f} {grav_norm:10.6f} {floor_active:>6}")

    m_new, _ = ds_combine(m_c, e_eq)
    m_c = m_new

print()

# ============================================================
# PART 5: Rank-2 entangled coupled system
# ============================================================
print("=" * 60)
print("RANK-2 COUPLED DS (entangled)")
print("=" * 60)

# Two sites, each uses the other as evidence
mA = np.array([0.05+0.02j, 0.35+0.01j, 0.30-0.015j, 0.30+0.005j], dtype=complex)
mA = mA / np.sum(np.abs(mA))
mB = np.array([0.08-0.01j, 0.30+0.02j, 0.32-0.01j, 0.30+0.005j], dtype=complex)
mB = mB / np.sum(np.abs(mB))

print(f"Initial Born_A = {born(mA):.6f}, Born_B = {born(mB):.6f}")
print()

print(f"{'Step':>4} {'Born_A':>8} {'Born_B':>8} {'|K_AB|':>8} {'dbar_A':>10} {'grav_A':>10} {'floor_A':>7} {'floor_B':>7}")
print("-" * 75)

for step in range(50):
    bA = born(mA)
    bB = born(mB)

    # Coupled: A uses B as evidence, B uses A
    K_AB = conflict(mA, mB)

    # dbar at site A (using B as evidence)
    db_A = dbar_phi(mA, mB)
    frob_A = np.sqrt(np.sum(np.abs(db_A)**2))
    tr_A = np.trace(db_A) / 4
    sym_A = (db_A + db_A.T) / 2
    sym_tl_A = sym_A - tr_A * np.eye(4)
    grav_A = np.sqrt(np.sum(np.abs(sym_tl_A)**2))

    fA = "YES" if np.real(bA) < FLOOR else "no"
    fB = "YES" if np.real(bB) < FLOOR else "no"

    if step < 25 or step % 5 == 0:
        print(f"{step:4d} {np.abs(bA):8.5f} {np.abs(bB):8.5f} {np.abs(K_AB):8.5f} "
              f"{frob_A:10.6f} {grav_A:10.6f} {fA:>7} {fB:>7}")

    # DS step: coupled
    mA_new, _ = ds_combine(mA, mB)
    mB_new, _ = ds_combine(mB, mA)
    mA, mB = mA_new, mB_new

print()
print(f"Final state A: {mA}")
print(f"Final state B: {mB}")
print(f"Final Born_A = {born(mA):.6f}, Born_B = {born(mB):.6f}")
print(f"Final K_AB = {conflict(mA, mB):.6f}")

# ============================================================
# PART 6: Decompose the gravitational signal
# ============================================================
print()
print("=" * 60)
print("DECOMPOSITION OF GRAVITATIONAL SIGNAL")
print("=" * 60)

# Collect dbar along the coupled trajectory
mA = np.array([0.05+0.02j, 0.35+0.01j, 0.30-0.015j, 0.30+0.005j], dtype=complex)
mA = mA / np.sum(np.abs(mA))
mB = np.array([0.08-0.01j, 0.30+0.02j, 0.32-0.01j, 0.30+0.005j], dtype=complex)
mB = mB / np.sum(np.abs(mB))

max_grav = 0
max_grav_step = 0
total_grav = 0
total_conf = 0
total_spin1 = 0
n_floor_active = 0

for step in range(100):
    bA = born(mA)
    if np.real(bA) < FLOOR:
        n_floor_active += 1

    db = dbar_phi(mA, mB)
    frob = np.sqrt(np.sum(np.abs(db)**2))

    if frob > 1e-10:
        tr = np.trace(db) / 4
        antisym = (db - db.T) / 2
        sym = (db + db.T) / 2
        sym_tl = sym - tr * np.eye(4)

        conf = 4 * np.abs(tr)**2
        spin1 = np.sum(np.abs(antisym)**2)
        grav = np.sum(np.abs(sym_tl)**2)

        total_conf += conf
        total_spin1 += spin1
        total_grav += grav

        if grav > max_grav:
            max_grav = grav
            max_grav_step = step

    mA_new, _ = ds_combine(mA, mB)
    mB_new, _ = ds_combine(mB, mA)
    mA, mB = mA_new, mB_new

total = total_conf + total_spin1 + total_grav
print(f"Floor active in {n_floor_active}/100 steps")
print(f"Peak graviton signal at step {max_grav_step}")
print()
if total > 1e-15:
    print(f"Integrated signal decomposition:")
    print(f"  Conformal (spin-0): {total_conf/total*100:.1f}%")
    print(f"  Spin-1:             {total_spin1/total*100:.1f}%")
    print(f"  Graviton (spin-2):  {total_grav/total*100:.1f}%")
else:
    print("No gravitational signal detected (floor never activated)")
    print("This means the initial conditions don't cross the Born floor.")
    print("Need initial conditions with Born < 1/27 = 0.037")
