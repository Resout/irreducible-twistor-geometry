"""
The K*=7/30 curve as moduli space of gauge theories.

Investigation:
1. Map the K*=7/30 curve precisely — what is its exact equation?
2. At each point, compute the full spectral data (λ₀, λ₁, Δ)
3. Check: do the 36 gauge groups from the any-G gauntlet land on this curve?
4. What physical quantities vary along the curve and which are constant?
5. Is there a natural coordinate on the curve?
6. What are the ENDPOINTS of the curve? What physics do they represent?
7. The off-curve family — what does K* ≠ 7/30 mean physically?
"""

import numpy as np
from scipy.optimize import brentq, fsolve

H = 3
FLOOR = 1.0 / H**3
ETA = (H-1)**2 / H**3

# ============================================================
# DS machinery (same as before)
# ============================================================
def enforce_floor(m):
    s, theta = m[:3], m[3]
    ssq = sum(si**2 for si in s)
    total_sq = ssq + theta**2
    if theta**2 / total_sq >= FLOOR:
        return m.copy()
    ss = sum(s)
    if ss < 1e-15: return m.copy()
    r = ssq / ss**2
    disc = (2*r)**2 + 4*(26-r)*r
    t = (-2*r + np.sqrt(disc)) / (2*(26-r))
    alpha = (1 - t) / ss
    return np.array([s[0]*alpha, s[1]*alpha, s[2]*alpha, t])

def ds_step(m, e):
    s, theta = m[:3], m[3]
    se, phi = e[:3], e[3]
    s_new = s * se + s * phi + theta * se
    theta_new = theta * phi
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)
    d = 1.0 - K
    if abs(d) < 1e-15: return m.copy(), K
    m_out = np.zeros(4)
    m_out[:3] = s_new / d
    m_out[3] = theta_new / d
    m_out = m_out / np.sum(m_out)
    return enforce_floor(m_out), K

def find_fp(e, n=1000):
    m = np.array([0.7, 0.05, 0.05, 0.2])
    for _ in range(n):
        m2, K = ds_step(m, e)
        if np.max(np.abs(m2 - m)) < 1e-15: break
        m = m2
    return m, sum(m[i]*e[j] for i in range(3) for j in range(3) if i != j)

def jacobian_3x3(m, e, eps=1e-8):
    """Projected Jacobian on L₁=1 tangent space."""
    m0, _ = ds_step(m, e)
    J = np.zeros((3, 3))
    for j in range(3):
        idx = [0, 1, 3][j]
        m_p = m.copy()
        m_p[idx] += eps
        m_p[3 if idx != 3 else 0] -= eps
        mp_out, _ = ds_step(m_p, e)
        for i in range(3):
            J[i, j] = (mp_out[[0,1,3][i]] - m0[[0,1,3][i]]) / eps
    return J

# ============================================================
# PART 1: Precise mapping of the K*=7/30 curve
# ============================================================
print("=" * 80)
print("PART 1: The K*=7/30 curve — precise mapping")
print("=" * 80)

# For each φ, find the p_dom that gives K*=7/30
def K_at_pdom_phi(p_dom, phi):
    p_w = (1 - p_dom) / 2
    e = np.array([p_dom*(1-phi), p_w*(1-phi), p_w*(1-phi), phi])
    if any(ei <= 0 for ei in e): return -1
    e = e / np.sum(e)  # ensure L1
    m, K = find_fp(e)
    return K

curve_data = []
print(f"\n{'phi':>8s} {'p_dom':>8s} {'K*':>12s} {'lam0':>10s} {'lam1':>10s} "
      f"{'Delta':>8s} {'theta':>8s} {'s1':>8s} {'s2':>8s}")
print("-" * 95)

for phi in np.linspace(0.005, 0.65, 130):
    try:
        # Find p_dom where K* = 7/30
        def obj(p):
            return K_at_pdom_phi(p, phi) - 7.0/30

        # Bracket search
        p_lo, p_hi = 0.35, 0.99
        f_lo, f_hi = obj(p_lo), obj(p_hi)
        if f_lo * f_hi > 0:
            continue

        p_dom = brentq(obj, p_lo, p_hi, xtol=1e-12)

        p_w = (1 - p_dom) / 2
        e = np.array([p_dom*(1-phi), p_w*(1-phi), p_w*(1-phi), phi])
        e = e / np.sum(e)
        m, K = find_fp(e)

        if abs(K - 7/30) > 1e-6:
            continue

        # Spectral data
        J = jacobian_3x3(m, e)
        eigs = np.sort(np.abs(np.linalg.eigvals(J)))[::-1]
        lam0, lam1 = eigs[0], eigs[1] if len(eigs) > 1 else 0
        gap = -np.log(lam0) if 0 < lam0 < 1 else 0

        # Agreement/commitment decomposition
        agree = sum(m[i]*e[i] for i in range(3))
        commit = sum(m[i]*e[3] + m[3]*e[i] for i in range(3))
        ignor = m[3]*e[3]

        # Commutator content: |s × e|²
        sx = np.array([m[1]*e[2] - m[2]*e[1],
                       m[2]*e[0] - m[0]*e[2],
                       m[0]*e[1] - m[1]*e[0]])
        comm_sq = np.sum(sx**2)

        curve_data.append({
            'phi': phi, 'p_dom': p_dom, 'K': K,
            'lam0': lam0, 'lam1': lam1, 'gap': gap,
            'm': m, 'e': e,
            'theta': m[3], 's1': m[0], 's2': m[1],
            'agree': agree, 'commit': commit, 'ignor': ignor,
            'comm_sq': comm_sq,
            'e1': e[0], 'e2': e[1], 'phi_e': e[3],
            'splitting': abs(lam0 - lam1) / lam0 if lam0 > 0 else 0,
        })

        if len(curve_data) % 5 == 1:
            d = curve_data[-1]
            print(f"{phi:8.4f} {p_dom:8.4f} {K:12.8f} {lam0:10.6f} {lam1:10.6f} "
                  f"{gap:8.4f} {m[3]:8.5f} {m[0]:8.5f} {m[1]:8.5f}")
    except:
        continue

print(f"\nTotal points on K*=7/30 curve: {len(curve_data)}")

# ============================================================
# PART 2: Invariants along the curve
# ============================================================
print("\n" + "=" * 80)
print("PART 2: What is CONSTANT along the K*=7/30 curve?")
print("=" * 80)

if curve_data:
    quantities = {
        'K*': [d['K'] for d in curve_data],
        'θ*': [d['theta'] for d in curve_data],
        's₁*': [d['s1'] for d in curve_data],
        's₂*': [d['s2'] for d in curve_data],
        'λ₀': [d['lam0'] for d in curve_data],
        'λ₁': [d['lam1'] for d in curve_data],
        'Δ': [d['gap'] for d in curve_data],
        'splitting': [d['splitting'] for d in curve_data],
        'agree': [d['agree'] for d in curve_data],
        'commit': [d['commit'] for d in curve_data],
        'ignore': [d['ignor'] for d in curve_data],
        '|s×e|²': [d['comm_sq'] for d in curve_data],
        'K*+φ': [d['K'] + d['phi'] for d in curve_data],
        'θ*/φ': [d['theta']/d['phi'] for d in curve_data],
        's₁/e₁': [d['s1']/d['e1'] if d['e1'] > 0.001 else 0 for d in curve_data],
        'agree/K': [d['agree']/d['K'] if d['K'] > 0.001 else 0 for d in curve_data],
        'commit/agree': [d['commit']/d['agree'] if d['agree'] > 0.001 else 0
                        for d in curve_data],
        'θ*·(1-K*)': [d['theta']*(1-d['K']) for d in curve_data],
        'agree+ignore': [d['agree']+d['ignor'] for d in curve_data],
        '(agree-K)/commit': [(d['agree']-d['K'])/d['commit']
                            if d['commit'] > 0.001 else 0 for d in curve_data],
        'λ₀·λ₁': [d['lam0']*d['lam1'] for d in curve_data],
        'λ₀+λ₁': [d['lam0']+d['lam1'] for d in curve_data],
        'λ₀²+λ₁²': [d['lam0']**2+d['lam1']**2 for d in curve_data],
    }

    print(f"\n{'Quantity':>20s} {'Mean':>12s} {'Std':>12s} {'CV':>8s} {'Min':>10s} {'Max':>10s}  Status")
    print("-" * 95)

    for name, vals in quantities.items():
        vals = np.array(vals)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0: continue
        mn, sd = np.mean(vals), np.std(vals)
        cv = sd / abs(mn) if abs(mn) > 1e-10 else float('inf')
        status = "*** CONSTANT ***" if cv < 0.01 else ("~ constant" if cv < 0.05 else "varies")
        print(f"{name:>20s} {mn:12.6f} {sd:12.6f} {cv:8.4f} {vals.min():10.6f} {vals.max():10.6f}  {status}")

# ============================================================
# PART 3: Endpoints of the curve
# ============================================================
print("\n" + "=" * 80)
print("PART 3: Endpoints — what are the limits?")
print("=" * 80)

if curve_data:
    cd = sorted(curve_data, key=lambda d: d['phi'])

    print("\nLow-φ endpoint (specific evidence):")
    d = cd[0]
    print(f"  φ = {d['phi']:.4f}, p_dom = {d['p_dom']:.4f}")
    print(f"  λ₀ = {d['lam0']:.6f}, Δ = {d['gap']:.4f}")
    print(f"  m* = [{d['m'][0]:.4f}, {d['m'][1]:.4f}, {d['m'][2]:.4f}, {d['m'][3]:.4f}]")
    print(f"  e* = [{d['e'][0]:.4f}, {d['e'][1]:.4f}, {d['e'][2]:.4f}, {d['e'][3]:.4f}]")
    print(f"  |s×e|² = {d['comm_sq']:.6f}")

    print("\nHigh-φ endpoint (ignorant evidence):")
    d = cd[-1]
    print(f"  φ = {d['phi']:.4f}, p_dom = {d['p_dom']:.4f}")
    print(f"  λ₀ = {d['lam0']:.6f}, Δ = {d['gap']:.4f}")
    print(f"  m* = [{d['m'][0]:.4f}, {d['m'][1]:.4f}, {d['m'][2]:.4f}, {d['m'][3]:.4f}]")
    print(f"  e* = [{d['e'][0]:.4f}, {d['e'][1]:.4f}, {d['e'][2]:.4f}, {d['e'][3]:.4f}]")
    print(f"  |s×e|² = {d['comm_sq']:.6f}")

    # What happens as φ → 0? Evidence becomes completely specific
    # What happens as φ → max? Evidence becomes pure ignorance
    print(f"\n  As φ→0: evidence becomes maximally specific (one hypothesis dominates)")
    print(f"    → λ₀→{cd[0]['lam0']:.4f}, gap→{cd[0]['gap']:.4f} (LARGE gap)")
    print(f"  As φ→max: evidence approaches pure ignorance")
    print(f"    → λ₀→{cd[-1]['lam0']:.4f}, gap→{cd[-1]['gap']:.4f} (SMALL gap)")

# ============================================================
# PART 4: The paper's specific equilibrium — where is it on the curve?
# ============================================================
print("\n" + "=" * 80)
print("PART 4: The paper's equilibrium on the K*=7/30 curve")
print("=" * 80)

# The paper's equilibrium: p_dom ≈ 0.9322, e* ≈ (0.631, 0.120, 0.120, 0.128)
# Find it on our curve
paper_phi = 0.128  # approximate
closest = min(curve_data, key=lambda d: abs(d['phi'] - paper_phi))
print(f"\nPaper's approximate φ = {paper_phi}")
print(f"Closest point on curve: φ = {closest['phi']:.4f}, p_dom = {closest['p_dom']:.4f}")
print(f"  λ₀ = {closest['lam0']:.6f} (paper: 0.28291)")
print(f"  Δ  = {closest['gap']:.4f} (paper: 1.263)")
print(f"  m* = {closest['m']}")
print(f"  e* = {closest['e']}")

# ============================================================
# PART 5: Natural coordinate on the curve
# ============================================================
print("\n" + "=" * 80)
print("PART 5: Natural coordinate — what parametrises the curve?")
print("=" * 80)

if curve_data:
    # The curve is parametrised by φ. But is there a more natural coordinate?
    # Options: δ (diagonal content), commutator content, λ₀, Δ, p_dom

    phis = np.array([d['phi'] for d in curve_data])
    lam0s = np.array([d['lam0'] for d in curve_data])
    pdoms = np.array([d['p_dom'] for d in curve_data])
    comms = np.array([d['comm_sq'] for d in curve_data])
    gaps = np.array([d['gap'] for d in curve_data])
    thetas = np.array([d['theta'] for d in curve_data])

    # Is λ₀ a simple function of φ?
    print("\nλ₀ as a function of φ along the curve:")
    for deg in [1, 2, 3]:
        coeffs = np.polyfit(phis, lam0s, deg)
        fitted = np.polyval(coeffs, phis)
        rmse = np.sqrt(np.mean((fitted - lam0s)**2))
        print(f"  Degree {deg}: RMSE = {rmse:.6f}")

    # Is λ₀ a simple function of p_dom?
    print("\nλ₀ as a function of p_dom along the curve:")
    for deg in [1, 2, 3]:
        coeffs = np.polyfit(pdoms, lam0s, deg)
        fitted = np.polyval(coeffs, pdoms)
        rmse = np.sqrt(np.mean((fitted - lam0s)**2))
        print(f"  Degree {deg}: RMSE = {rmse:.6f}")

    # Is |s×e|² a natural coordinate?
    print("\nλ₀ as a function of |s×e|² along the curve:")
    mask = comms > 1e-10
    if np.sum(mask) > 5:
        for deg in [1, 2]:
            coeffs = np.polyfit(comms[mask], lam0s[mask], deg)
            fitted = np.polyval(coeffs, comms[mask])
            rmse = np.sqrt(np.mean((fitted - lam0s[mask])**2))
            print(f"  Degree {deg}: RMSE = {rmse:.6f}")

    # θ* as coordinate? (it varies only slightly along the curve)
    print(f"\nθ* range along curve: [{thetas.min():.6f}, {thetas.max():.6f}]")
    print(f"θ* CV: {thetas.std()/thetas.mean():.4f}")

# ============================================================
# PART 6: Off-curve interpretation
# ============================================================
print("\n" + "=" * 80)
print("PART 6: What does K* ≠ 7/30 mean?")
print("=" * 80)

print("""
The conservation law K*(H²+1) - η*H² = 1 at the self-consistent equilibrium.

For K* ≠ 7/30:
  CL = K*·10 - 4/3

  CL < 1 (K* < 7/30): drain < fill + budget
    → evidence is "too weak" — doesn't generate enough conflict
    → the system would drift toward higher K* under self-consistent evolution

  CL > 1 (K* > 7/30): drain > fill + budget
    → evidence is "too strong" — generates excess conflict
    → the system would drift toward lower K* under self-consistent evolution

  CL = 1 (K* = 7/30): drain = fill + budget
    → evidence exactly matches the channel geometry
    → self-consistent equilibrium
""")

# Verify: at K* < 7/30, does iterating with self-evidence push K up?
# At K* > 7/30, does it push K down?
print("Verification: does the system self-correct toward K*=7/30?")

for K_target, label in [(0.10, "K*=0.10 (below 7/30)"),
                         (0.40, "K*=0.40 (above 7/30)")]:
    # Find a fixed point with this K*
    for phi_try in np.linspace(0.05, 0.60, 50):
        try:
            def obj(p):
                return K_at_pdom_phi(p, phi_try) - K_target
            p = brentq(obj, 0.35, 0.99, xtol=1e-10)
            p_w = (1-p)/2
            e = np.array([p*(1-phi_try), p_w*(1-phi_try), p_w*(1-phi_try), phi_try])
            e = e / np.sum(e)
            m, K = find_fp(e)
            if abs(K - K_target) < 0.01:
                # Now use m as BOTH state and evidence (self-evidence)
                m_self, K_self = find_fp(m)
                K_self_val = sum(m_self[i]*m[j] for i in range(3) for j in range(3) if i != j)

                # Or: iterate with m generating its own evidence
                # Use the state as evidence for itself
                m_iter = m.copy()
                K_track = []
                for _ in range(50):
                    # Evidence IS the current state
                    m_iter, K_iter = ds_step(m_iter, m_iter)
                    K_actual = sum(m_iter[i]*m_iter[j] for i in range(3) for j in range(3) if i != j)
                    K_track.append(K_actual)

                print(f"\n  {label}:")
                print(f"    Initial K* = {K:.4f}")
                print(f"    After 50 self-evidence steps: K_self = {K_track[-1]:.6f}")
                print(f"    Direction: {'toward' if abs(K_track[-1] - 7/30) < abs(K - 7/30) else 'away from'} 7/30")
                break
        except:
            continue
