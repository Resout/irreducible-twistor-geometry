"""Heavy-sampling computation: does min transverse λ converge to H/(H-1) = 3/2?

Strategy:
  1. Sample 50,000 random det=0 states (symmetric and asymmetric)
  2. For each, compute the transverse amplification ratio |det(F(m+εv))|/|det(m+εv)|
     at multiple ε scales to extract the linearised eigenvalue
  3. Also: follow 1000 dynamical trajectories for 100 steps each,
     measuring λ at each step where det is small
  4. Report distribution of λ values and whether min → 3/2
"""
import numpy as np

# === Core DS functions (H=3) ===

def ds_combine(m1, m2):
    H = 3
    s1, th1 = m1[:H], m1[H]
    s2, th2 = m2[:H], m2[H]
    s_new = s1*s2 + s1*th2 + th1*s2
    th_new = th1*th2
    K = np.sum(s1)*np.sum(s2) - np.sum(s1*s2)
    n = 1.0 - K
    if abs(n) < 1e-15: n = 1e-15
    return np.append(s_new, th_new) / n, K

def det_M(m):
    return 0.5*(m[3]**2 - m[0]**2 - m[1]**2 - m[2]**2)

def construct_R0_evidence(m):
    s = m[:3]
    theta = m[3]
    c = np.zeros(3, dtype=complex)
    for i in range(3):
        others_sq = sum(s[j]**2 for j in range(3) if j != i)
        den = 2*s[i]*(s[i]+theta) + others_sq
        if abs(den) < 1e-15: c[i] = 0
        else: c[i] = s[i]*(s[i]+theta) / den
    denom = 1.0 - 2.0*np.sum(c)
    if abs(denom) < 1e-15: return None
    phi = 1.0 / denom
    e = -2.0 * phi * c
    return np.append(e, phi)

def F_map(m):
    ev = construct_R0_evidence(m)
    if ev is None: return m
    result, K = ds_combine(m, ev)
    return result

def transverse_eigenvalue(m0, eps=1e-8):
    """Compute linearised transverse eigenvalue at a det≈0 point."""
    grad_det = np.array([-m0[0], -m0[1], -m0[2], m0[3]])
    norm = np.linalg.norm(grad_det)
    if norm < 1e-15: return None
    grad_det = grad_det / norm

    m_pert = m0 + eps * grad_det
    m_pert = m_pert / np.sum(m_pert)

    det_in = abs(det_M(m_pert))
    if det_in < 1e-30: return None

    m_out = F_map(m_pert)
    det_out = abs(det_M(m_out))

    return det_out / det_in

# === Part 1: Static sampling at det=0 ===

print("=" * 60)
print("PART 1: Transverse λ at random det=0 states")
print("=" * 60)

np.random.seed(42)
lambdas_sym = []
lambdas_asym = []

# Symmetric (s2=s3)
for trial in range(25000):
    b = (0.05 * np.random.randn()) + (0.05 * np.random.randn() - 0.2) * 1j
    denom = -2 + 4 * b
    if abs(denom) < 0.01: continue
    a = -(1 - 4*b + 2*b**2) / denom
    theta = 1 - a - 2*b
    m = np.array([a, b, b, theta], dtype=complex)
    if abs(np.sum(m) - 1) > 1e-8: continue
    if abs(det_M(m)) > 1e-6: continue

    lam = transverse_eigenvalue(m)
    if lam is not None and np.isfinite(lam) and lam > 0:
        lambdas_sym.append(lam)

# Asymmetric (all s_i different)
for trial in range(25000):
    s1 = 0.3 + 0.15*np.random.randn() + 0.15j*np.random.randn()
    s2 = 0.15 + 0.1*np.random.randn() + 0.1j*np.random.randn()
    s3 = 0.1 + 0.08*np.random.randn() + 0.08j*np.random.randn()
    theta = 1 - s1 - s2 - s3
    m = np.array([s1, s2, s3, theta], dtype=complex)

    # Project onto det=0
    ss = np.sum(m[:3]**2)
    theta_new = np.sqrt(ss)
    m[3] = theta_new
    m = m / np.sum(m)

    if abs(det_M(m)) > 1e-4: continue

    lam = transverse_eigenvalue(m)
    if lam is not None and np.isfinite(lam) and lam > 0:
        lambdas_asym.append(lam)

ls = np.array(lambdas_sym)
la = np.array(lambdas_asym)

print(f"\nSymmetric det=0 states: n={len(ls)}")
print(f"  mean λ = {ls.mean():.6f}")
print(f"  std  λ = {ls.std():.6f}")
print(f"  min  λ = {ls.min():.6f}")
print(f"  max  λ = {ls.max():.6f}")
print(f"  ALL > 1: {(ls > 1).all()}")

print(f"\nAsymmetric det=0 states: n={len(la)}")
print(f"  mean λ = {la.mean():.6f}")
print(f"  std  λ = {la.std():.6f}")
print(f"  min  λ = {la.min():.6f}")
print(f"  max  λ = {la.max():.6f}")
print(f"  ALL > 1: {(la > 1).all()}")

all_static = np.concatenate([ls, la])
print(f"\nAll static: n={len(all_static)}")
print(f"  min λ = {all_static.min():.6f}")

# === Part 2: Along dynamical trajectories ===

print("\n" + "=" * 60)
print("PART 2: λ along dynamical trajectories (near det=0)")
print("=" * 60)

lambdas_traj = []
traj_details = []  # (trial, step, |det|, λ)

np.random.seed(123)
n_traj = 2000
n_steps = 100
det_threshold = 0.05  # measure λ when |det| < this

for trial in range(n_traj):
    # Start from random mass function near ignorance
    s = 0.2 + 0.05*np.random.randn(3) + 0.05j*np.random.randn(3)
    theta = 1 - np.sum(s)
    m = np.array([s[0], s[1], s[2], theta], dtype=complex)

    for step in range(n_steps):
        d = abs(det_M(m))

        if d < det_threshold and d > 1e-10:
            lam = transverse_eigenvalue(m, eps=1e-8)
            if lam is not None and np.isfinite(lam) and lam > 0:
                lambdas_traj.append(lam)
                traj_details.append((trial, step, d, lam))

        # Random evidence
        e_s = 0.1*np.random.rand(3) + 0.05j*np.random.randn(3)
        e_th = 1 - np.sum(e_s)
        e = np.array([e_s[0], e_s[1], e_s[2], e_th], dtype=complex)

        m, K = ds_combine(m, e)

lt = np.array(lambdas_traj) if lambdas_traj else np.array([])

if len(lt) > 0:
    print(f"\nTrajectory measurements (|det|<{det_threshold}): n={len(lt)}")
    print(f"  mean λ = {lt.mean():.6f}")
    print(f"  std  λ = {lt.std():.6f}")
    print(f"  min  λ = {lt.min():.6f}")
    print(f"  max  λ = {lt.max():.6f}")
    print(f"  ALL > 1: {(lt > 1).all()}")

    # Find the minimum λ point
    min_idx = np.argmin(lt)
    trial_m, step_m, det_m, lam_m = traj_details[min_idx]
    print(f"\n  Min λ details: trial={trial_m}, step={step_m}, |det|={det_m:.6e}, λ={lam_m:.6f}")

    # Distribution near the minimum
    bottom_10 = np.sort(lt)[:min(10, len(lt))]
    print(f"  Bottom 10 λ values: {[f'{x:.4f}' for x in bottom_10]}")
else:
    print("No trajectory points with |det| < threshold found")

# === Part 3: Targeted search for minimum λ ===

print("\n" + "=" * 60)
print("PART 3: Targeted search — minimize λ over det=0 surface")
print("=" * 60)

# Parametric sweep: systematically explore the det=0 surface
# det=0 means θ² = s₁² + s₂² + s₃²
# Parametrise: s_i = r_i e^{iφ_i}, θ = √(Σr_i² e^{2iφ_i})

np.random.seed(456)
lambdas_targeted = []

for trial in range(50000):
    # Random magnitudes and phases
    r = np.random.exponential(0.3, 3)
    phi = np.random.uniform(-np.pi, np.pi, 3)
    s = r * np.exp(1j * phi)

    # θ from det=0 condition
    ss = np.sum(s**2)
    theta = np.sqrt(ss)  # one branch
    if np.random.rand() > 0.5:
        theta = -theta  # other branch

    m = np.array([s[0], s[1], s[2], theta], dtype=complex)

    # Normalise L1
    total = np.sum(m)
    if abs(total) < 1e-10: continue
    m = m / total

    if abs(det_M(m)) > 1e-4: continue

    # Check that components aren't too extreme
    if np.max(np.abs(m)) > 10: continue

    lam = transverse_eigenvalue(m, eps=1e-9)
    if lam is not None and np.isfinite(lam) and lam > 0:
        lambdas_targeted.append(lam)

lt2 = np.array(lambdas_targeted)
print(f"\nTargeted det=0 sweep: n={len(lt2)}")
if len(lt2) > 0:
    print(f"  mean λ = {lt2.mean():.6f}")
    print(f"  std  λ = {lt2.std():.6f}")
    print(f"  min  λ = {lt2.min():.6f}")
    print(f"  max  λ = {lt2.max():.6f}")

    bottom_20 = np.sort(lt2)[:min(20, len(lt2))]
    print(f"  Bottom 20: {[f'{x:.4f}' for x in bottom_20]}")

    # Percentiles
    for p in [1, 5, 10, 25, 50]:
        print(f"  {p}th percentile: {np.percentile(lt2, p):.6f}")

# === Part 4: Check against constants ===

print("\n" + "=" * 60)
print("PART 4: Comparison with framework constants")
print("=" * 60)

all_lambdas = []
if len(ls) > 0: all_lambdas.append(ls)
if len(la) > 0: all_lambdas.append(la)
if len(lt) > 0: all_lambdas.append(lt)
if len(lt2) > 0: all_lambdas.append(lt2)
all_lambdas = np.concatenate(all_lambdas)

lmin = all_lambdas.min()
lmean = all_lambdas.mean()

print(f"\nGlobal minimum λ = {lmin:.6f}")
print(f"Global mean λ = {lmean:.6f}")
print(f"Total measurements: {len(all_lambdas)}")

H = 3
candidates = {
    "H/(H-1)": H/(H-1),           # 3/2 = 1.5
    "H/(H-1) squared": (H/(H-1))**2,  # 9/4 = 2.25
    "1 + 1/H": 1 + 1/H,           # 4/3
    "1 + 1/(H-1)": 1 + 1/(H-1),   # 3/2 = 1.5
    "1 + K*": 1 + 7/30,            # 37/30
    "1/(1-K*)": 1/(1-7/30),        # 30/23
    "H²/(H²-1)": H**2/(H**2-1),   # 9/8
    "sqrt(H)": np.sqrt(H),         # 1.732
    "(H+1)/H": (H+1)/H,           # 4/3
    "2H/(H+1)": 2*H/(H+1),        # 3/2 = 1.5
    "golden ratio": (1+np.sqrt(5))/2,  # 1.618
}

print(f"\nCandidate constants (comparing to global min = {lmin:.6f}):")
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - lmin)):
    diff = abs(val - lmin)
    print(f"  {name:25s} = {val:.6f}  (diff = {diff:.6f})")

# Also check: does min λ depend on dimensionality H?
# We're at H=3. Check H=2 and H=4 for comparison.

print("\n" + "=" * 60)
print("PART 5: Dimensionality dependence (H=2 and H=4)")
print("=" * 60)

def ds_combine_H(m1, m2, H):
    s1, th1 = m1[:H], m1[H]
    s2, th2 = m2[:H], m2[H]
    s_new = s1*s2 + s1*th2 + th1*s2
    th_new = th1*th2
    K = np.sum(s1)*np.sum(s2) - np.sum(s1*s2)
    n = 1.0 - K
    if abs(n) < 1e-15: n = 1e-15
    return np.append(s_new, th_new) / n, K

def det_M_H(m, H):
    return 0.5*(m[H]**2 - np.sum(m[:H]**2))

def construct_R0_H(m, H):
    s = m[:H]
    theta = m[H]
    c = np.zeros(H, dtype=complex)
    for i in range(H):
        others_sq = sum(s[j]**2 for j in range(H) if j != i)
        den = 2*s[i]*(s[i]+theta) + others_sq
        if abs(den) < 1e-15: c[i] = 0
        else: c[i] = s[i]*(s[i]+theta) / den
    denom = 1.0 - 2.0*np.sum(c)
    if abs(denom) < 1e-15: return None
    phi = 1.0 / denom
    e = -2.0 * phi * c
    return np.append(e, phi)

def F_map_H(m, H):
    ev = construct_R0_H(m, H)
    if ev is None: return m
    result, K = ds_combine_H(m, ev, H)
    return result

def transverse_eigenvalue_H(m0, H, eps=1e-8):
    grad_det = np.append(-m0[:H], m0[H])
    norm = np.linalg.norm(grad_det)
    if norm < 1e-15: return None
    grad_det = grad_det / norm
    m_pert = m0 + eps * grad_det
    m_pert = m_pert / np.sum(m_pert)
    det_in = abs(det_M_H(m_pert, H))
    if det_in < 1e-30: return None
    m_out = F_map_H(m_pert, H)
    det_out = abs(det_M_H(m_out, H))
    return det_out / det_in

for H_test in [2, 4]:
    np.random.seed(789 + H_test)
    lambdas_H = []

    for trial in range(20000):
        r = np.random.exponential(0.3, H_test)
        phi = np.random.uniform(-np.pi, np.pi, H_test)
        s = r * np.exp(1j * phi)
        ss = np.sum(s**2)
        theta = np.sqrt(ss)
        if np.random.rand() > 0.5: theta = -theta

        m = np.append(s, theta).astype(complex)
        total = np.sum(m)
        if abs(total) < 1e-10: continue
        m = m / total

        if abs(det_M_H(m, H_test)) > 1e-4: continue
        if np.max(np.abs(m)) > 10: continue

        lam = transverse_eigenvalue_H(m, H_test, eps=1e-9)
        if lam is not None and np.isfinite(lam) and lam > 0:
            lambdas_H.append(lam)

    lH = np.array(lambdas_H)
    ratio = H_test / (H_test - 1) if H_test > 1 else float('inf')
    print(f"\nH={H_test}: n={len(lH)}")
    if len(lH) > 0:
        print(f"  min λ  = {lH.min():.6f}")
        print(f"  mean λ = {lH.mean():.6f}")
        print(f"  H/(H-1) = {ratio:.6f}")
        print(f"  |min - H/(H-1)| = {abs(lH.min() - ratio):.6f}")
        bottom_5 = np.sort(lH)[:5]
        print(f"  Bottom 5: {[f'{x:.4f}' for x in bottom_5]}")
