"""
Deep investigation of the 2D family of DS fixed points.

The DS map with Born floor has a 2-parameter family of equilibria,
parametrised by the evidence distribution. We map this family completely.

Questions:
  1. What is the geometry of this family? Boundaries?
  2. How do K*, Δ, λ₀, δ vary across it?
  3. What is the K*=7/30 curve? What's its equation?
  4. Are there other distinguished curves?
  5. Is the conservation law value a simple function on this family?
  6. What happens at the boundaries?
  7. Is there a NATURAL parametrisation?
"""

import numpy as np
from scipy.optimize import fsolve

H = 3
FLOOR = 1.0 / H**3
ETA = (H-1)**2 / H**3  # 4/27

# ============================================================
# Core DS machinery
# ============================================================
def enforce_floor(m):
    s, theta = m[:3], m[3]
    ssq = sum(si**2 for si in s)
    total_sq = ssq + theta**2
    born_val = theta**2 / total_sq if total_sq > 0 else 1
    if born_val >= FLOOR:
        return m.copy()
    ss = sum(s)
    if ss < 1e-15:
        return m.copy()
    r = ssq / ss**2
    a_c = 26 - r
    b_c = 2 * r
    c_c = -r
    disc = b_c**2 - 4 * a_c * c_c
    t = (-b_c + np.sqrt(disc)) / (2 * a_c)
    alpha = (1 - t) / ss
    return np.array([s[0]*alpha, s[1]*alpha, s[2]*alpha, t])

def ds_combine_raw(m, e):
    s, theta = m[:3], m[3]
    se, phi = e[:3], e[3]
    s_new = s * se + s * phi + theta * se
    theta_new = theta * phi
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    if abs(denom) < 1e-15:
        return m.copy(), K
    m_out = np.zeros(4)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    total = np.sum(m_out)
    if total > 0:
        m_out = m_out / total
    return m_out, K

def full_step(m, e):
    m_ds, K = ds_combine_raw(m, e)
    return enforce_floor(m_ds), K

def born(m):
    return m[3]**2 / np.sum(m**2) if np.sum(m**2) > 0 else 1

def K_conflict(m, e):
    return sum(m[i] * e[j] for i in range(3) for j in range(3) if i != j)

def find_fixed_point(e, max_iter=800):
    m = np.array([0.7, 0.05, 0.05, 0.2])
    for _ in range(max_iter):
        m_new, K = full_step(m, e)
        if np.max(np.abs(m_new - m)) < 1e-14:
            return m_new, K_conflict(m_new, e), True
        m = m_new
    return m, K_conflict(m, e), False

def compute_eigenvalues(m, e, eps=1e-8):
    """Eigenvalues of 3x3 projected Jacobian."""
    m0, _ = full_step(m, e)
    J = np.zeros((3, 3))
    for j in range(3):
        idx = [0, 1, 3][j]
        m_p = m.copy()
        m_p[idx] += eps
        comp_idx = 3 if idx != 3 else 0
        m_p[comp_idx] -= eps
        m_p_out, _ = full_step(m_p, e)
        for i in range(3):
            idx_i = [0, 1, 3][i]
            J[i, j] = (m_p_out[idx_i] - m0[idx_i]) / eps
    return np.sort(np.abs(np.linalg.eigvals(J)))[::-1]

# ============================================================
# PART 1: Map the full family in natural coordinates
# ============================================================
# Evidence is S₂-symmetric: e = (e₁, e₂, e₂, φ) with e₁+2e₂+φ=1
# Natural parameters: p_dom = e₁/(e₁+2e₂) [singleton concentration]
#                     φ = ignorance fraction
# Then e₁ = p_dom*(1-φ), e₂ = (1-p_dom)/2*(1-φ)

print("=" * 80)
print("PART 1: Full map of the 2D family")
print("=" * 80)

# Dense grid
p_dom_grid = np.linspace(0.35, 0.995, 200)
phi_grid = np.linspace(0.005, 0.70, 200)

family = []
for p_dom in p_dom_grid:
    for phi in phi_grid:
        p_weak = (1 - p_dom) / 2
        e1 = p_dom * (1 - phi)
        e2 = p_weak * (1 - phi)
        if e1 <= 0 or e2 <= 0 or phi <= 0 or phi >= 1:
            continue
        e = np.array([e1, e2, e2, phi])

        m_fp, K_val, converged = find_fixed_point(e)
        if not converged:
            continue
        if K_val <= -0.5 or K_val >= 1:
            continue
        b = born(m_fp)
        if abs(b - FLOOR) > 1e-4:
            continue

        # Compute spectral data
        eigs = compute_eigenvalues(m_fp, e)
        lam0 = eigs[0] if len(eigs) > 0 else 0
        lam1 = eigs[1] if len(eigs) > 1 else 0

        # Conservation law value
        cons = K_val * (H**2 + 1) - ETA * H**2

        # Diagonal content
        diag_prod = sum(m_fp[i]*e[i] for i in range(3))
        total_prod = sum(m_fp[i]*e[j] for i in range(3) for j in range(3))
        delta = diag_prod / total_prod if total_prod > 0 else 0

        # Self-conflict of evidence
        K_self_e = sum(e[i]*e[j] for i in range(3) for j in range(3) if i != j)

        # θ ratio: θ*/φ and (1-K*)
        theta_over_phi_ratio = m_fp[3] / phi if phi > 0 else 0

        family.append({
            'p_dom': p_dom, 'phi': phi,
            'K': K_val, 'cons': cons,
            'lam0': lam0, 'lam1': lam1,
            'gap': -np.log(lam0) if 0 < lam0 < 1 else 0,
            'delta': delta, 'theta': m_fp[3],
            'm': m_fp, 'e': e,
            'K_self_e': K_self_e,
            'theta_phi_ratio': theta_over_phi_ratio,
            'e1': e1, 'e2': e2,
        })

print(f"Total fixed points mapped: {len(family)}")

# ============================================================
# PART 2: What are the boundaries of the family?
# ============================================================
print("\n" + "=" * 80)
print("PART 2: Boundaries and shape")
print("=" * 80)

K_vals = np.array([r['K'] for r in family])
phi_vals = np.array([r['phi'] for r in family])
pdom_vals = np.array([r['p_dom'] for r in family])
cons_vals = np.array([r['cons'] for r in family])
lam0_vals = np.array([r['lam0'] for r in family])
gap_vals = np.array([r['gap'] for r in family])
theta_vals = np.array([r['theta'] for r in family])
delta_vals = np.array([r['delta'] for r in family])

print(f"K*    range: [{K_vals.min():.6f}, {K_vals.max():.6f}]")
print(f"φ     range: [{phi_vals.min():.6f}, {phi_vals.max():.6f}]")
print(f"p_dom range: [{pdom_vals.min():.6f}, {pdom_vals.max():.6f}]")
print(f"Δ     range: [{gap_vals.min():.6f}, {gap_vals.max():.6f}]")
print(f"λ₀    range: [{lam0_vals.min():.6f}, {lam0_vals.max():.6f}]")
print(f"θ*    range: [{theta_vals.min():.6f}, {theta_vals.max():.6f}]")
print(f"δ     range: [{delta_vals.min():.6f}, {delta_vals.max():.6f}]")
print(f"CL    range: [{cons_vals.min():.6f}, {cons_vals.max():.6f}]")

# ============================================================
# PART 3: The K*=7/30 curve
# ============================================================
print("\n" + "=" * 80)
print("PART 3: The K*=7/30 curve within the family")
print("=" * 80)

K_target = 7/30
curve_730 = [r for r in family if abs(r['K'] - K_target) < 0.002]
curve_730.sort(key=lambda r: r['phi'])

print(f"\nPoints on K*=7/30 curve: {len(curve_730)}")
print(f"\n{'phi':>8s} {'p_dom':>8s} {'K*':>12s} {'lambda0':>10s} {'Delta':>8s} "
      f"{'delta':>8s} {'theta':>8s} {'cons':>8s}")
print("-" * 85)

for r in curve_730[::max(1, len(curve_730)//25)]:
    print(f"{r['phi']:8.4f} {r['p_dom']:8.4f} {r['K']:12.8f} {r['lam0']:10.6f} "
          f"{r['gap']:8.4f} {r['delta']:8.4f} {r['theta']:8.5f} {r['cons']:8.5f}")

# What's the equation of this curve? Fit p_dom as function of phi
if len(curve_730) > 5:
    phi_c = np.array([r['phi'] for r in curve_730])
    pdom_c = np.array([r['p_dom'] for r in curve_730])
    lam0_c = np.array([r['lam0'] for r in curve_730])
    gap_c = np.array([r['gap'] for r in curve_730])

    # Try polynomial fit
    for deg in [1, 2, 3]:
        coeffs = np.polyfit(phi_c, pdom_c, deg)
        fitted = np.polyval(coeffs, phi_c)
        rmse = np.sqrt(np.mean((fitted - pdom_c)**2))
        print(f"\n  Degree-{deg} fit of p_dom(φ) on K*=7/30 curve: RMSE = {rmse:.6f}")
        if deg <= 2:
            print(f"    Coefficients: {coeffs}")

    print(f"\n  Along K*=7/30 curve:")
    print(f"    λ₀ range: [{lam0_c.min():.6f}, {lam0_c.max():.6f}]")
    print(f"    Δ  range: [{gap_c.min():.6f}, {gap_c.max():.6f}]")
    print(f"    The gap VARIES along the curve: Δ is NOT fixed by K*=7/30 alone")

# ============================================================
# PART 4: Is K* a simple function of the evidence parameters?
# ============================================================
print("\n" + "=" * 80)
print("PART 4: What determines K*? Functional relationships")
print("=" * 80)

# Test: is K* = f(p_dom, φ)?
# Already know it's a function of evidence. But is there a simple formula?

# Check: K* vs φ at fixed p_dom
print("\nK* vs φ at fixed p_dom:")
for p_target in [0.50, 0.70, 0.90]:
    subset = [r for r in family if abs(r['p_dom'] - p_target) < 0.01]
    if len(subset) > 3:
        phis = np.array([r['phi'] for r in subset])
        Ks = np.array([r['K'] for r in subset])
        order = np.argsort(phis)
        phis, Ks = phis[order], Ks[order]
        # Is K linear in φ?
        if len(phis) > 2:
            coeffs = np.polyfit(phis, Ks, 1)
            fitted = np.polyval(coeffs, phis)
            rmse = np.sqrt(np.mean((fitted - Ks)**2))
            print(f"  p_dom≈{p_target}: K* = {coeffs[0]:.4f}·φ + {coeffs[1]:.4f} "
                  f"(RMSE={rmse:.4f}, n={len(subset)})")

# Check: K* vs p_dom at fixed φ
print("\nK* vs p_dom at fixed φ:")
for phi_target in [0.05, 0.15, 0.30]:
    subset = [r for r in family if abs(r['phi'] - phi_target) < 0.008]
    if len(subset) > 3:
        pdoms = np.array([r['p_dom'] for r in subset])
        Ks = np.array([r['K'] for r in subset])
        order = np.argsort(pdoms)
        pdoms, Ks = pdoms[order], Ks[order]
        if len(pdoms) > 2:
            coeffs = np.polyfit(pdoms, Ks, 2)
            fitted = np.polyval(coeffs, pdoms)
            rmse = np.sqrt(np.mean((fitted - Ks)**2))
            print(f"  φ≈{phi_target}: K* = {coeffs[0]:.4f}·p² + {coeffs[1]:.4f}·p + "
                  f"{coeffs[2]:.4f} (RMSE={rmse:.6f}, n={len(subset)})")

# ============================================================
# PART 5: The conservation law as a function on the family
# ============================================================
print("\n" + "=" * 80)
print("PART 5: Conservation law value across the family")
print("=" * 80)

# CL = K*(H²+1) - η*H² = 10K* - 9η
# At 7/30: CL = 10*(7/30) - 9*(4/27) = 7/3 - 4/3 = 1
# What is CL as a function of (p_dom, φ)?

print(f"\nCL = K*·10 - η·9")
print(f"CL at K*=7/30: {10*7/30 - 9*4/27:.6f}")

# Is CL = simple function of K*?
print(f"\nCL vs K*: is it just 10K* - {9*ETA:.6f}?")
print(f"  Yes, by definition: CL = 10·K* - 1.333...")
print(f"  CL = 1  ⟺  K* = (1 + 9η)/10 = (1 + 4/3)/10 = 7/30")
print(f"  This is tautological — CL is linear in K*")

# The real question: what ELSE is special about K*=7/30?
# Let's look at the relationship between m* and e* at the fixed point

# ============================================================
# PART 6: State-evidence relationship across the family
# ============================================================
print("\n" + "=" * 80)
print("PART 6: State-evidence relationship")
print("=" * 80)

# At the fixed point, what's the relationship between m* and e*?
# From θ equation: θ_out = floor(θ*φ/(1-K*))
# Since m*=Φ(m*,e*), we need floor(θ*φ/(1-K*)) = θ*
# This means θ*φ/(1-K*) < θ* (floor boosts it back up)
# i.e., φ/(1-K*) < 1, i.e., φ < 1-K*

print("\nRelationship: φ vs (1-K*)")
print(f"{'phi':>8s} {'1-K*':>10s} {'phi/(1-K*)':>12s} {'K*':>10s}")
print("-" * 45)
for r in sorted(family, key=lambda r: r['K'])[::len(family)//20]:
    ratio = r['phi'] / (1 - r['K']) if r['K'] < 1 else float('inf')
    print(f"{r['phi']:8.4f} {1-r['K']:10.6f} {ratio:12.6f} {r['K']:10.6f}")

# Check: is there a universal relationship between θ* and the evidence?
print("\n\nθ*/φ ratio across the family:")
ratios = np.array([r['theta']/r['phi'] for r in family if r['phi'] > 0.01])
print(f"  Range: [{ratios.min():.4f}, {ratios.max():.4f}]")
print(f"  Mean:  {ratios.mean():.4f}")
print(f"  Std:   {ratios.std():.4f}")
if ratios.std() / ratios.mean() < 0.1:
    print(f"  >>> θ*/φ is approximately CONSTANT = {ratios.mean():.4f}")
else:
    print(f"  >>> θ*/φ varies significantly (CV = {ratios.std()/ratios.mean():.2f})")

# Check: what about s₁*/e₁* ratio?
print("\ns₁*/e₁* ratio across the family:")
s1_e1_ratios = np.array([r['m'][0]/r['e'][0] for r in family if r['e'][0] > 0.01])
print(f"  Range: [{s1_e1_ratios.min():.4f}, {s1_e1_ratios.max():.4f}]")
print(f"  Mean:  {s1_e1_ratios.mean():.4f}")
print(f"  Std:   {s1_e1_ratios.std():.4f}")

# ============================================================
# PART 7: The deep question — what makes K*=7/30 self-consistent?
# ============================================================
print("\n" + "=" * 80)
print("PART 7: Self-consistency characterisation")
print("=" * 80)

# At the self-consistent equilibrium, the evidence is "produced by" the
# channel geometry. What does this mean concretely?
#
# The pre-normalisation DS step produces:
#   s_new,i = s_i*e_i + s_i*φ + θ*e_i  (agreement + commitment terms)
#   θ_new = θ*φ                          (joint ignorance)
#   K = Σ_{i≠j} s_i*e_j                  (cross-focal, discarded)
#
# The "channel geometry" is: H²=9 singleton-pair channels, H²+1=10 effective.
# The "self-consistent" condition might mean: the RATIO of agreement to
# commitment to ignorance matches the channel dimensions.

print("\nDecomposition of the pre-normalisation output at each fixed point:")
print(f"{'K*':>10s} {'agree':>10s} {'commit':>10s} {'ignor':>10s} "
      f"{'agree/K':>10s} {'commit/K':>10s} {'a/(a+c+i)':>12s}")
print("-" * 80)

for r in sorted(family, key=lambda r: r['K'])[::len(family)//20]:
    m, e = r['m'], r['e']
    agree = sum(m[i]*e[i] for i in range(3))
    commit = sum(m[i]*e[3] + m[3]*e[i] for i in range(3))
    ignor = m[3]*e[3]
    K = r['K']
    total_pre = agree + commit + ignor  # = 1 - K

    agree_frac = agree / total_pre if total_pre > 0 else 0

    print(f"{K:10.6f} {agree:10.6f} {commit:10.6f} {ignor:10.6f} "
          f"{agree/K if K>0.001 else 0:10.4f} {commit/K if K>0.001 else 0:10.4f} "
          f"{agree_frac:12.6f}")

# Check: at K*=7/30, does the ratio agree/commit/ignore match H²:2H:1?
print(f"\nChannel dimension ratios: H²:2H:1 = 9:6:1 (fractions: 0.5625, 0.375, 0.0625)")
print(f"These sum to 16 = (H+1)² = total channels")
print(f"\nAt K*=7/30 points:")
for r in curve_730[::max(1, len(curve_730)//8)]:
    m, e = r['m'], r['e']
    agree = sum(m[i]*e[i] for i in range(3))
    commit = sum(m[i]*e[3] + m[3]*e[i] for i in range(3))
    ignor = m[3]*e[3]
    total = agree + commit + ignor
    print(f"  φ={r['phi']:.3f}: agree={agree/total:.4f} commit={commit/total:.4f} "
          f"ignore={ignor/total:.4f} | "
          f"a/i={agree/ignor:.2f} c/i={commit/ignor:.2f}")

# ============================================================
# PART 8: Can we find a CLOSED FORM for K*(p_dom, φ)?
# ============================================================
print("\n" + "=" * 80)
print("PART 8: Closed form for K*(p_dom, φ)")
print("=" * 80)

# The fixed-point condition WITHOUT floor would give:
# θ(1-K) = θφ → K = 1-φ (exact, no floor)
# s_i(1-K) = s_i*e_i + s_i*φ + θ*e_i → 0 = s_i*e_i + θ*e_i → e_i = 0

# WITH floor, the relationship is more complex.
# Let's check: is K = 1 - φ a good approximation?
print("\nK* vs (1-φ) — the no-floor prediction:")
print(f"{'K*':>10s} {'1-phi':>10s} {'K*-(1-phi)':>12s} {'phi':>8s}")
print("-" * 45)
for r in sorted(family, key=lambda r: r['K'])[::len(family)//15]:
    print(f"{r['K']:10.6f} {1-r['phi']:10.6f} {r['K']-(1-r['phi']):12.6f} {r['phi']:8.4f}")

# The floor modifies this. Let's find the actual relationship.
# At the fixed point WITH floor:
# 1. DS gives θ_raw = θ*φ/(1-K)
# 2. Floor boosts θ_raw → θ*
# 3. So θ* > θ*φ/(1-K), meaning the floor adds Δθ = θ* - θ*φ/(1-K)
# 4. This Δθ comes from redistributing singleton mass

# The floor-corrected θ equation:
# θ* = F(θ*φ/(1-K*)) where F is the floor enforcement
# At Born floor equality: 26θ*² = Σs_i*²
# After DS+normalise: θ_raw = θ*φ/(1-K*), s_raw,i = (s_i*e_i + s_i*φ + θ*e_i)/(1-K*)
# Floor adjusts: θ_raw → θ* with rescaling of singletons

# The key insight: the floor maps (s_raw, θ_raw) → (s_raw*α, t)
# where 26t² = (s_raw*α)² and t + 3*s_raw*α*... = 1
# This gives a specific relationship between θ_raw and θ*

# Let's just measure the floor correction across the family
print("\n\nFloor correction across the family:")
print(f"{'K*':>10s} {'theta_raw':>12s} {'theta*':>10s} {'boost':>10s} {'boost/theta':>12s}")
print("-" * 60)
for r in sorted(family, key=lambda r: r['K'])[::len(family)//15]:
    m, e = r['m'], r['e']
    theta_raw = m[3]*e[3]/(1-r['K']) if r['K'] < 1 else 0
    # After L1 normalise:
    m_ds, _ = ds_combine_raw(m, e)
    theta_after_ds = m_ds[3]
    boost = m[3] - theta_after_ds
    print(f"{r['K']:10.6f} {theta_after_ds:12.8f} {m[3]:10.6f} "
          f"{boost:10.6f} {boost/m[3] if m[3]>0 else 0:12.6f}")

# ============================================================
# PART 9: The 1-K* = φ + floor_correction identity
# ============================================================
print("\n" + "=" * 80)
print("PART 9: Exact relationship at each fixed point")
print("=" * 80)

# At the fixed point: Φ(m*, e*) = m*
# The θ component: after DS, θ_ds = θ*φ/(1-K*), then L1 normalise, then floor
# After L1: θ_ds is part of a vector that sums to 1
# After floor: θ is boosted to satisfy Born = 1/27

# What fraction of the total does the floor redistribute?
print(f"{'K*':>8s} {'φ':>8s} {'θ*':>8s} {'θ_ds':>8s} {'floor_frac':>10s} "
      f"{'(1-K*)/φ':>10s} {'K*+φ':>8s}")
print("-" * 70)

for r in sorted(family, key=lambda r: r['K'])[::len(family)//20]:
    m, e = r['m'], r['e']
    m_ds, K = ds_combine_raw(m, e)
    floor_frac = (m[3] - m_ds[3]) / m[3] if m[3] > 0 else 0
    ratio_1mK_phi = (1-r['K'])/r['phi'] if r['phi'] > 0 else 0
    print(f"{r['K']:8.5f} {r['phi']:8.5f} {m[3]:8.5f} {m_ds[3]:8.5f} "
          f"{floor_frac:10.6f} {ratio_1mK_phi:10.4f} {r['K']+r['phi']:8.5f}")

# Check if K* + φ = constant?
Kphi = np.array([r['K'] + r['phi'] for r in family])
print(f"\nK* + φ: range [{Kphi.min():.4f}, {Kphi.max():.4f}], "
      f"mean {Kphi.mean():.4f}, std {Kphi.std():.4f}")
if Kphi.std() < 0.01:
    print(f"  >>> K* + φ ≈ {Kphi.mean():.4f} is approximately CONSTANT")

# Check K* * φ
Ktphi = np.array([r['K'] * r['phi'] for r in family])
print(f"K* · φ: range [{Ktphi.min():.6f}, {Ktphi.max():.6f}]")

# Check (1-K*)/φ
ratio = np.array([(1-r['K'])/r['phi'] for r in family if r['phi'] > 0.01])
print(f"(1-K*)/φ: range [{ratio.min():.4f}, {ratio.max():.4f}], "
      f"mean {ratio.mean():.4f}, std {ratio.std():.4f}")

# Check floor_frac
ff = np.array([(r['theta'] - ds_combine_raw(r['m'], r['e'])[0][3])/r['theta']
               for r in family if r['theta'] > 0.01])
print(f"floor_frac: range [{ff.min():.4f}, {ff.max():.4f}], "
      f"mean {ff.mean():.4f}, std {ff.std():.4f}")
if ff.std() / ff.mean() < 0.1:
    print(f"  >>> Floor correction fraction ≈ {ff.mean():.4f} is approximately CONSTANT!")
    print(f"  >>> This would mean the floor always corrects by ~{ff.mean()*100:.1f}%")
