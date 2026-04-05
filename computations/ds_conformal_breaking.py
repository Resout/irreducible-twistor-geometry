"""
Does the Born floor break conformal invariance?

Conformal gravity: g → Ω²g is a symmetry. Action = ∫|W|².
Einstein gravity: g → Ω²g is NOT a symmetry. Action = ∫R.

The Born floor 1/27 is a fixed scale. If ∂̄Φ transforms non-covariantly
under rescaling, conformal invariance is broken → Einstein, not Weyl.

Test: scale the mass function by Ω and check how ∂̄Φ transforms.
Under conformal: ∂̄Φ should transform homogeneously.
Under Einstein: it shouldn't.
"""
import numpy as np

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

def born(m):
    return np.abs(m[3])**2 / np.sum(np.abs(m)**2)

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

def dbar_no_floor(m, e, eps=1e-7):
    """Same but without floor enforcement."""
    n = 4
    z = m.copy()
    def phi(z_in):
        out, _ = ds_combine(z_in, e, apply_floor=False)
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

# ============================================================
print("=" * 60)
print("TEST 1: CONFORMAL SCALING OF ∂̄Φ")
print("=" * 60)

# If conformal: ∂̄Φ(Ωm, Ωe) = Ω^k ∂̄Φ(m, e) for some weight k
# On CP³ (projective), scaling m → Ωm is trivial (same point).
# But the Born floor breaks this: Born(Ωm) = Born(m)
# (Born is scale-invariant), so the floor activates at the same
# states regardless of Ω. But the CORRECTION depends on magnitudes.

# Better test: the floor introduces a non-projective element.
# On CP³, everything should depend only on ratios m_i/m_j.
# The Born floor depends on |m_i|² / Σ|m_j|², which IS projective.
# BUT the floor enforcement changes |m_i| not m_i/m_j.
# It changes the POINT on the simplex, not just the ray.

# The real test of conformal breaking:
# DS without floor is polynomial → holomorphic → conformal on CP³
# DS with floor involves |θ|² → not holomorphic → breaks conformal

# Let's just measure it directly.

# Test state where floor is active
m_test = np.array([0.02+0.01j, 0.38+0.02j, 0.32-0.015j, 0.28+0.005j], dtype=complex)
m_test = m_test / np.sum(np.abs(m_test))
e_test = np.array([0.63+0.01j, 0.12+0.005j, 0.12-0.003j, 0.13+0.002j], dtype=complex)
e_test = e_test / np.sum(np.abs(e_test))

print(f"Born(m) = {born(m_test):.6f}, floor = {FLOOR:.6f}")
print(f"Floor active: {born(m_test) < FLOOR}")
print()

# ∂̄Φ with and without floor
db_with = dbar_phi(m_test, e_test)
db_without = dbar_no_floor(m_test, e_test)
db_floor_only = db_with - db_without

frob_with = np.sqrt(np.sum(np.abs(db_with)**2))
frob_without = np.sqrt(np.sum(np.abs(db_without)**2))
frob_floor = np.sqrt(np.sum(np.abs(db_floor_only)**2))

print(f"||∂̄Φ_total|| = {frob_with:.6f}")
print(f"||∂̄Φ_DS||    = {frob_without:.6f}  (without floor, should be ~0 if DS is holomorphic)")
print(f"||∂̄Φ_floor|| = {frob_floor:.6f}  (floor contribution only)")
print()

# Decompose the floor-only part
tr_f = np.trace(db_floor_only) / 4
sym_f = (db_floor_only + db_floor_only.T) / 2
asym_f = (db_floor_only - db_floor_only.T) / 2
stl_f = sym_f - tr_f * np.eye(4)

tot_f = 4*np.abs(tr_f)**2 + np.sum(np.abs(asym_f)**2) + np.sum(np.abs(stl_f)**2)
if tot_f > 1e-15:
    print("Floor-only ∂̄Φ decomposition:")
    print(f"  Conformal: {4*np.abs(tr_f)**2/tot_f*100:.1f}%")
    print(f"  Spin-1:    {np.sum(np.abs(asym_f)**2)/tot_f*100:.1f}%")
    print(f"  Graviton:  {np.sum(np.abs(stl_f)**2)/tot_f*100:.1f}%")

print()

# ============================================================
print("=" * 60)
print("TEST 2: IS DS COMBINATION HOLOMORPHIC?")
print("=" * 60)

# Pure DS (no floor) should be holomorphic (polynomial in m)
# So ∂̄Φ_DS should be 0 at all points
print("Testing ∂̄ of pure DS (no floor) at various states:")
for name, m_t in [
    ("Near ignorance", np.array([0.01+0.001j, 0.01+0.002j, 0.01-0.001j, 0.97+0.003j], dtype=complex)),
    ("Democratic", np.array([0.25+0.01j, 0.25-0.02j, 0.25+0.015j, 0.25+0.005j], dtype=complex)),
    ("Low theta", np.array([0.02+0.01j, 0.38+0.02j, 0.32-0.01j, 0.28+0.005j], dtype=complex)),
]:
    m_t = m_t / np.sum(np.abs(m_t))
    db = dbar_no_floor(m_t, e_test)
    print(f"  {name:20s}: ||∂̄|| = {np.sqrt(np.sum(np.abs(db)**2)):.2e}")

print()

# ============================================================
print("=" * 60)
print("TEST 3: SCALE DEPENDENCE OF FLOOR")
print("=" * 60)

# The floor activates when Born(θ) < 1/27
# Born = |θ|²/Σ|m_i|² is scale-invariant
# But the correction is NOT: it changes the L1-normalised point
# This means the floor introduces a preferred scale on the simplex

# Test: compare floor correction at different "distances" from threshold
print("Floor correction magnitude vs Born deficit:")
print(f"{'Born':>10} {'||correction||':>15} {'graviton frac':>15}")
for born_target in [0.035, 0.030, 0.025, 0.020, 0.015, 0.010, 0.005]:
    # Construct state with given Born
    # Born = θ²/(θ² + Σs²). Set s1=s2=s3=s, θ such that Born = target
    # θ² = target * Σm² / 1, with L1=1: θ + 3s = 1
    # θ²/(θ² + 3s²) = target → θ² = target*(θ² + 3s²)
    # θ²(1-target) = 3*target*s² → θ/s = √(3*target/(1-target))
    # θ + 3s = 1 → s(√(3t/(1-t)) + 3) = 1
    t = born_target
    ratio = np.sqrt(3*t/(1-t))
    s_val = 1.0 / (ratio + 3)
    th_val = ratio * s_val

    m_t = np.array([s_val+0.01j, s_val-0.005j, s_val+0.003j, th_val+0.002j], dtype=complex)
    m_t = m_t / np.sum(np.abs(m_t))

    db_w = dbar_phi(m_t, e_test)
    db_wo = dbar_no_floor(m_t, e_test)
    db_fl = db_w - db_wo
    frob_fl = np.sqrt(np.sum(np.abs(db_fl)**2))

    tr_fl = np.trace(db_fl) / 4
    sym_fl = (db_fl + db_fl.T) / 2
    stl_fl = sym_fl - tr_fl * np.eye(4)
    grav_frac = np.sum(np.abs(stl_fl)**2) / max(np.sum(np.abs(db_fl)**2), 1e-15)

    print(f"{born_target:10.3f} {frob_fl:15.6f} {grav_frac*100:14.1f}%")

print()

# ============================================================
print("=" * 60)
print("TEST 4: THE CONFORMAL BREAKING EQUATION")
print("=" * 60)

# In conformal gravity: W_μν = 0 at equilibrium (conformally flat)
# In Einstein gravity: R_μν = 8πG T_μν (Ricci curvature = matter)
#
# The DS framework gives us both Weyl (traceless) and Ricci (trace)
# The floor contribution to ∂̄Φ has both components.
#
# Key question: what sets the RATIO of Ricci to Weyl?
# If it's determined by the floor structure, we get Einstein.
# If it's freely adjustable, we stay conformal.

# At the K*=7/30 fixed point:
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

# Perturb slightly into complex
m_pert = m_star + np.array([0.001j, -0.0003j, -0.0002j, -0.0005j])
m_pert = m_pert / np.sum(np.abs(m_pert))

db_total = dbar_phi(m_pert, e_star)
db_ds = dbar_no_floor(m_pert, e_star)
db_floor = db_total - db_ds

# The trace of ∂̄Φ is the Ricci part (conformal factor)
# The traceless symmetric part is the Weyl/graviton part
# In Einstein gravity: Ricci is DETERMINED by matter (T_μν)
# In conformal gravity: Ricci is FREELY specifiable

tr_total = np.trace(db_total) / 4
tr_ds = np.trace(db_ds) / 4
tr_floor = np.trace(db_floor) / 4

print(f"At near-equilibrium (K*=7/30):")
print(f"  tr(∂̄Φ_total)/4 = {tr_total:.6f}")
print(f"  tr(∂̄Φ_DS)/4    = {tr_ds:.6f}")
print(f"  tr(∂̄Φ_floor)/4 = {tr_floor:.6f}")
print()

# The floor FIXES the trace — it's not a free parameter.
# The trace is determined by the floor structure (1/27, L1=1).
# This is exactly what Einstein gravity does: Ricci = matter.
# The "matter" here is the floor enforcement.

# Check: is the trace determined by the Born deficit?
born_val = born(m_pert)
born_deficit = FLOOR - np.real(born_val)
print(f"  Born = {np.real(born_val):.8f}")
print(f"  Born deficit = {born_deficit:.8f}")
print(f"  |tr(floor)|/|deficit| = {np.abs(tr_floor)/max(abs(born_deficit), 1e-15):.4f}")
print()

# ============================================================
print("=" * 60)
print("TEST 5: RICCI-TO-WEYL RATIO IS FIXED")
print("=" * 60)

# If the ratio Ricci/Weyl is fixed by the framework (not freely
# adjustable), then we have Einstein gravity, not conformal.

print("Ricci/Weyl ratio at various states near equilibrium:")
print(f"{'perturbation':>15} {'|Ricci|':>10} {'|Weyl|':>10} {'ratio':>10}")

for i, delta in enumerate([
    np.array([0.001j, -0.0003j, -0.0002j, -0.0005j]),
    np.array([0.002j, -0.001j, -0.0005j, -0.0005j]),
    np.array([0.003j, 0.001j, -0.002j, -0.002j]),
    np.array([-0.001j, 0.002j, -0.001j, 0.000j]),
    np.array([0.005j, -0.002j, -0.001j, -0.002j]),
    np.array([0.01j, -0.003j, -0.004j, -0.003j]),
    np.array([0.001+0.002j, -0.0005-0.001j, -0.0005+0.0005j, 0.000-0.0015j]),
    np.array([0.02j, -0.008j, -0.005j, -0.007j]),
]):
    mp = m_star + delta
    mp = mp / np.sum(np.abs(mp))
    db = dbar_phi(mp, e_star)
    tr = np.trace(db) / 4
    sym = (db + db.T) / 2
    stl = sym - tr * np.eye(4)
    ricci = 2 * np.abs(tr)  # proportional to scalar curvature
    weyl = np.sqrt(np.sum(np.abs(stl)**2))
    ratio = ricci / weyl if weyl > 1e-12 else float('inf')
    print(f"  pert {i:2d}:      {ricci:10.6f} {weyl:10.6f} {ratio:10.4f}")

print()
print("If the ratio is approximately constant, Ricci is determined by Weyl")
print("(or vice versa) — conformal symmetry is broken → Einstein gravity.")
print("If it varies freely, conformal symmetry survives → Weyl gravity.")

# ============================================================
print("\n" + "=" * 60)
print("TEST 6: WHAT BREAKS CONFORMAL INVARIANCE")
print("=" * 60)

# Three things could break it:
# 1. The Born floor (introduces absolute scale 1/27)
# 2. The L1=1 constraint (projects to simplex, not full projective space)
# 3. The conflict K (removes probability mass)

# Test: compute ∂̄Φ for DS steps where floor does NOT activate
# If ∂̄ = 0, then ONLY the floor breaks holomorphicity
# If ∂̄ ≠ 0 even without floor, something else is non-holomorphic too

# L1 renormalization involves |m_i|, which is non-holomorphic!
# L1 = Σ|m_i| depends on m̄_i. So even without Born floor,
# the L1 normalization is anti-holomorphic.

print("∂̄ of DS with L1 normalization but no Born floor:")
m_hi_born = np.array([0.1+0.01j, 0.1+0.005j, 0.1-0.003j, 0.7+0.002j], dtype=complex)
m_hi_born = m_hi_born / np.sum(np.abs(m_hi_born))
print(f"  Born = {born(m_hi_born):.4f} (>> 1/27, floor inactive)")

db_nf = dbar_no_floor(m_hi_born, e_star)
print(f"  ||∂̄Φ_no_floor|| = {np.sqrt(np.sum(np.abs(db_nf)**2)):.6f}")
print()

# Now test: DS without ANY normalization (raw polynomial)
def ds_raw(m, e):
    """DS without normalization or floor. Pure polynomial."""
    s, theta = m[:3], m[3]
    se, theta_e = e[:3], e[3]
    s_new = s * se + s * theta_e + theta * se
    theta_new = theta * theta_e
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    if abs(denom) < 1e-15:
        return m.copy()
    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    return m_out

def dbar_raw(m, e, eps=1e-7):
    n = 4
    z = m.copy()
    dbar = np.zeros((n, n), dtype=complex)
    for b in range(n):
        zp = z.copy(); zp[b] += eps
        zm = z.copy(); zm[b] -= eps
        d_dx = (ds_raw(zp, e) - ds_raw(zm, e)) / (2 * eps)
        zp = z.copy(); zp[b] += 1j * eps
        zm = z.copy(); zm[b] -= 1j * eps
        d_dy = (ds_raw(zp, e) - ds_raw(zm, e)) / (2 * eps)
        dbar[:, b] = 0.5 * (d_dx + 1j * d_dy)
    return dbar

print("∂̄ of raw DS (no normalization, no floor):")
db_raw = dbar_raw(m_hi_born, e_star)
print(f"  ||∂̄Φ_raw|| = {np.sqrt(np.sum(np.abs(db_raw)**2)):.2e}")
print(f"  (Should be ~0: raw DS is a rational function of z, hence holomorphic)")
print()

print("∂̄ of DS with L1 norm only:")
print(f"  ||∂̄Φ_L1|| = {np.sqrt(np.sum(np.abs(db_nf)**2)):.6f}")
print()

print("CONCLUSION:")
print("  Raw DS: holomorphic (∂̄ = 0)")
print("  + L1 normalization: non-holomorphic (|m_i| depends on m̄_i)")
print("  + Born floor: additional non-holomorphic contribution")
print()
print("  Both L1 and Born floor break conformal invariance.")
print("  L1 is the projective structure (CP³ constraint).")
print("  Born floor is the scale-setting (gravity).")
