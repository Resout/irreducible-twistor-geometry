"""
Sub-leading eigenvalues of the DS transfer operator.

The transfer operator at K*=7/30 has eigenvalues:
  λ₀ = 0.2829 (leading)
  λ₁ = 0.2813 (sub-leading)
  λ₂ ≈ 0 (negligible)

The near-degeneracy λ₀ ≈ λ₁ (0.6% splitting) predicts
two nearly-degenerate glueball states.

In YM: the 0++ glueball mass is m₀ = -ln(λ₀)/a.
A second state at m₁ = -ln(λ₁)/a would have:
  m₁/m₀ = ln(λ₁)/ln(λ₀) = 1.268/1.263 = 1.004

This means the two lightest excited states are within 0.4% of
each other in mass. In lattice YM, the 0++ and 2++ glueballs
have masses m(0++)/m(2++) ≈ 0.7. This doesn't match the near-degeneracy.

But: our eigenvalues are at the K*=7/30 equilibrium with specific
evidence. Different equilibria might have different splittings.
"""
import numpy as np

H = 3
FLOOR = 1.0 / H**3

def enforce_floor(m):
    s = m[:3].copy()
    lo, hi = m[3], 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        s_scale = (1.0 - mid) / np.sum(s) if np.sum(s) > 0 else 0
        s_trial = s * s_scale
        born = mid**2 / (np.sum(s_trial**2) + mid**2)
        if born < FLOOR: lo = mid
        else: hi = mid
    theta_new = (lo + hi) / 2
    s_scale = (1.0 - theta_new) / np.sum(s) if np.sum(s) > 0 else 0
    m_out = np.zeros(4)
    m_out[:3] = s * s_scale
    m_out[3] = theta_new
    return m_out

def ds_combine(m, e, apply_floor=True):
    s, theta = m[:3], m[3]
    se, theta_e = e[:3], e[3]
    s_new = s * se + s * theta_e + theta * se
    theta_new = theta * theta_e
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    m_out = np.zeros(4)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    if apply_floor:
        born_theta = m_out[3]**2 / np.sum(m_out**2)
        if born_theta < FLOOR:
            m_out = enforce_floor(m_out)
    return m_out, K

# K*=7/30 equilibrium
p_dom_exact = 0.9322756157
p_weak_exact = (1 - p_dom_exact) / 2
raw = np.array([np.sqrt(p_dom_exact*(1-FLOOR)), np.sqrt(p_weak_exact*(1-FLOOR)),
                np.sqrt(p_weak_exact*(1-FLOOR)), np.sqrt(FLOOR)])
e_exact = raw / np.sum(raw)

m = np.array([0.4, 0.2, 0.2, 0.2])
for _ in range(1000):
    m_new, _ = ds_combine(m, e_exact)
    if np.max(np.abs(m_new - m)) < 1e-15: break
    m = m_new
m_star = m_new

# Full eigenvalue analysis
eps = 1e-8
V = np.array([[1,0,0,-1], [0,1,0,-1], [0,0,1,-1]], dtype=float).T
VTV_inv = np.linalg.inv(V.T @ V)

J = np.zeros((4, 4))
f0 = ds_combine(m_star, e_exact)[0]
for j in range(4):
    mp = m_star.copy(); mp[j] += eps
    fp = ds_combine(mp, e_exact)[0]
    J[:, j] = (fp - f0) / eps

J_proj = VTV_inv @ V.T @ J @ V
evals_raw = np.linalg.eigvals(J_proj)
evals = np.sort(np.abs(evals_raw))[::-1]

# Also get eigenvectors
evals_full, evecs = np.linalg.eig(J_proj)
idx = np.argsort(np.abs(evals_full))[::-1]
evals_full = evals_full[idx]
evecs = evecs[:, idx]

print("=" * 60)
print("EXCITED STATE SPECTRUM")
print("=" * 60)
print()
print("Eigenvalues of the projected Jacobian (L₁=1 tangent space):")
for i in range(3):
    lam = evals_full[i]
    abs_lam = abs(lam)
    if abs_lam > 1e-10:
        delta = -np.log(abs_lam)
        print(f"  λ_{i} = {lam:+.10f}  (|λ| = {abs_lam:.10f}, Δ = {delta:.6f})")
    else:
        print(f"  λ_{i} = {lam:+.10f}  (|λ| ≈ 0, Δ → ∞)")

    # Eigenvector
    v = evecs[:, i]
    v_4d = V @ v  # back to 4D
    print(f"         eigenvector (4D): ({v_4d[0]:+.4f}, {v_4d[1]:+.4f}, {v_4d[2]:+.4f}, {v_4d[3]:+.4f})")
    print()

# Mass ratios
print("Mass ratios:")
if evals[0] > 1e-10 and evals[1] > 1e-10:
    m0 = -np.log(evals[0])
    m1 = -np.log(evals[1])
    print(f"  m₀ = -ln(λ₀) = {m0:.6f}")
    print(f"  m₁ = -ln(λ₁) = {m1:.6f}")
    print(f"  m₁/m₀ = {m1/m0:.6f}")
    print(f"  Splitting: {(m1/m0-1)*100:.2f}%")
    print()

# The near-degeneracy is because the two weak hypotheses (s₂=s₃)
# are symmetric. The eigenvectors should reveal this:
print("Symmetry analysis:")
print("  The fixed point has s₂=s₃ (two weak hypotheses symmetric).")
print("  λ₀ eigenvector: perturbation mixing s₁ and (s₂+s₃)/2")
print("  λ₁ eigenvector: perturbation in (s₂-s₃) direction")
print("  λ₂ ≈ 0: perturbation in the θ direction (strongly contracted)")
print()

# What if we break the s₂=s₃ symmetry?
print("=" * 60)
print("SYMMETRY-BROKEN SPECTRUM")
print("=" * 60)
print()

# Use evidence that breaks s₂ ↔ s₃ symmetry
e_broken = e_exact.copy()
e_broken[1] *= 1.3  # boost s₂ evidence
e_broken[2] *= 0.7  # reduce s₃ evidence
e_broken = e_broken / np.sum(e_broken)

m_b = np.array([0.4, 0.2, 0.2, 0.2])
for _ in range(1000):
    m_new, _ = ds_combine(m_b, e_broken)
    if np.max(np.abs(m_new - m_b)) < 1e-15: break
    m_b = m_new

_, K_b = ds_combine(m_b, e_broken, apply_floor=False)
print(f"Broken-symmetry FP: m = ({m_b[0]:.6f}, {m_b[1]:.6f}, {m_b[2]:.6f}, {m_b[3]:.6f})")
print(f"K* = {K_b:.6f}")

J_b = np.zeros((4, 4))
f0_b = ds_combine(m_b, e_broken)[0]
for j in range(4):
    mp = m_b.copy(); mp[j] += eps
    fp = ds_combine(mp, e_broken)[0]
    J_b[:, j] = (fp - f0_b) / eps

J_b_proj = VTV_inv @ V.T @ J_b @ V
evals_b = np.sort(np.abs(np.linalg.eigvals(J_b_proj)))[::-1]

print(f"Eigenvalues: {evals_b}")
if evals_b[0] > 1e-10 and evals_b[1] > 1e-10:
    m0_b = -np.log(evals_b[0])
    m1_b = -np.log(evals_b[1])
    print(f"m₁/m₀ = {m1_b/m0_b:.4f} (splitting: {(m1_b/m0_b-1)*100:.1f}%)")
