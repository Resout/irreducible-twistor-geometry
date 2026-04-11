#!/usr/bin/env python3
"""
INDEPENDENT REPRODUCTION: A₂ coupled Jacobian eigenvalues.
Uses the exact DS implementation from the codebase (spectral_gap_computation.py).
"""
import numpy as np

H = 3
FLOOR = 1.0 / H**3

def ds_combine(m, e, apply_floor=True):
    s, theta = m[:3], m[3]
    se, theta_e = e[:3], e[3]
    s_new = s * se + s * theta_e + theta * se
    theta_new = theta * theta_e
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    m_out = np.zeros(4)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    if apply_floor:
        born = m_out[3]**2 / np.sum(m_out**2)
        if born < FLOOR:
            m_out = enforce_floor(m_out)
    return m_out, K

def enforce_floor(m):
    s = m[:3].copy()
    lo, hi = m[3], 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        s_scale = (1.0 - mid) / np.sum(s) if np.sum(s) > 0 else 0
        s_trial = s * s_scale
        born = mid**2 / (np.sum(s_trial**2) + mid**2)
        if born < FLOOR:
            lo = mid
        else:
            hi = mid
    theta_new = (lo + hi) / 2
    s_scale = (1.0 - theta_new) / np.sum(s) if np.sum(s) > 0 else 0
    m_out = np.zeros(4)
    m_out[:3] = s * s_scale
    m_out[3] = theta_new
    return m_out

# ============================================================
# Step 1: Build the equilibrium evidence (from spectral_gap_computation.py)
# ============================================================
coeff = 3*np.sqrt(2) + 2 + np.sqrt(10/13)
s2_eq = 1.0 / coeff
s1_eq = 3*np.sqrt(2) * s2_eq
theta_eq = np.sqrt(10/13) * s2_eq
e_eq = np.array([s1_eq, s2_eq, s2_eq, theta_eq])

print("=" * 80)
print("A₂ COUPLED JACOBIAN EIGENVALUE REPRODUCTION")
print("=" * 80)
print(f"\nEvidence: ({e_eq[0]:.6f}, {e_eq[1]:.6f}, {e_eq[2]:.6f}, {e_eq[3]:.6f})")

# ============================================================
# Step 2: Find single-site fixed point
# ============================================================
m = np.array([0.4, 0.2, 0.2, 0.2])
for i in range(500):
    m_new, K = ds_combine(m, e_eq)
    if np.max(np.abs(m_new - m)) < 1e-15:
        break
    m = m_new
m_star = m_new
_, K_star = ds_combine(m_star, e_eq, apply_floor=False)
print(f"\nSingle-site fixed point:")
print(f"  m* = ({m_star[0]:.8f}, {m_star[1]:.8f}, {m_star[2]:.8f}, {m_star[3]:.8f})")
print(f"  K* = {K_star:.8f} (target: {7/30:.8f})")

# ============================================================
# Step 3: Build A₂ coupled system
# ============================================================
# A₂ Cartan matrix: [[2,-1],[-1,2]]
# Two sites, coupled at strength g = K* = 7/30
# Evidence for site i: e_base + g * A_ij * (m_j - uniform)

g = 7/30
uniform = np.array([0.25, 0.25, 0.25, 0.25])

def coupled_evidence(m_other, e_base, g_coupling):
    """Evidence for one site given the other site's mass function."""
    perturbation = g_coupling * (-1) * (m_other - uniform)  # A_ij = -1 for i≠j
    e = e_base + perturbation
    # Ensure positivity and L1=1
    e = np.maximum(e, 1e-10)
    e = e / np.sum(e)
    return e

def coupled_step(m1, m2, e_base, g_coupling):
    """One coupled update step for the A₂ system."""
    e1 = coupled_evidence(m2, e_base, g_coupling)
    e2 = coupled_evidence(m1, e_base, g_coupling)
    m1_new, K1 = ds_combine(m1, e1)
    m2_new, K2 = ds_combine(m2, e2)
    return m1_new, m2_new, K1, K2

print(f"\n{'='*80}")
print("A₂ COUPLED SYSTEM (g = K* = 7/30)")
print("=" * 80)

# Find coupled fixed point
m1 = np.array([0.5, 0.15, 0.15, 0.2])
m2 = np.array([0.15, 0.5, 0.15, 0.2])

for i in range(2000):
    m1_new, m2_new, K1, K2 = coupled_step(m1, m2, e_eq, g)
    diff = max(np.max(np.abs(m1_new - m1)), np.max(np.abs(m2_new - m2)))
    if i < 5 or i % 200 == 0:
        print(f"  Step {i}: diff={diff:.2e}, K1={K1:.6f}, K2={K2:.6f}")
    if diff < 1e-14:
        print(f"  Converged at step {i}")
        break
    m1, m2 = m1_new, m2_new

m1_star, m2_star = m1_new, m2_new
print(f"\n  Site 1: ({m1_star[0]:.8f}, {m1_star[1]:.8f}, {m1_star[2]:.8f}, {m1_star[3]:.8f})")
print(f"  Site 2: ({m2_star[0]:.8f}, {m2_star[1]:.8f}, {m2_star[2]:.8f}, {m2_star[3]:.8f})")

# ============================================================
# Step 4: Compute 8×8 coupled Jacobian
# ============================================================
print(f"\n{'='*80}")
print("8×8 COUPLED JACOBIAN")
print("=" * 80)

state = np.concatenate([m1_star, m2_star])  # 8-vector

def coupled_map(state_vec):
    m1 = state_vec[:4]
    m2 = state_vec[4:]
    m1_new, m2_new, _, _ = coupled_step(m1, m2, e_eq, g)
    return np.concatenate([m1_new, m2_new])

eps = 1e-8
f0 = coupled_map(state)
J8 = np.zeros((8, 8))
for j in range(8):
    s_pert = state.copy()
    s_pert[j] += eps
    f_pert = coupled_map(s_pert)
    J8[:, j] = (f_pert - f0) / eps

evals_8 = np.linalg.eigvals(J8)
evals_abs = np.sort(np.abs(evals_8))[::-1]

print("\n  8×8 Jacobian eigenvalues (by magnitude):")
for i, ev in enumerate(evals_abs):
    print(f"    |λ_{i}| = {ev:.6f}")

# The non-trivial eigenvalues (first 4, should match paper)
print(f"\n  Non-trivial eigenvalues (top 4):")
print(f"    λ₀ = {evals_abs[0]:.4f}  (expected: 0.5022)")
print(f"    λ₁ = {evals_abs[1]:.4f}  (expected: 0.4745)")
print(f"    λ₂ = {evals_abs[2]:.4f}  (expected: 0.3527)")
print(f"    λ₃ = {evals_abs[3]:.4f}  (expected: 0.3344)")

# Gaps
Delta_computed = -np.log(evals_abs[:4])
D0 = Delta_computed[0]
print(f"\n  Mass gaps Δ = -ln(λ):")
for i in range(4):
    print(f"    Δ_{i} = {Delta_computed[i]:.4f}")

print(f"\n  Mass ratios m/m(0++):")
for i in range(4):
    print(f"    Δ_{i}/Δ₀ = {Delta_computed[i]/D0:.4f}")

# Verify glueball spectrum
print(f"\n  Glueball ratios vs lattice QCD:")
print(f"    m(2++)/m(0++) = √2 = {np.sqrt(2):.4f} [exact theorem, not eigenvalue]")
print(f"    m(0-+)/m(0++) = Δ₂/Δ₀ = {Delta_computed[2]/D0:.4f} (lattice: 1.50)")
print(f"    m(0++*)/m(0++) = Δ₃/Δ₀ = {Delta_computed[3]/D0:.4f} (lattice: 1.56)")

# Eigenvectors
print(f"\n  Eigenvector analysis (site-1 vs site-2 symmetry):")
evecs = np.linalg.eig(J8)[1]
for i in range(4):
    idx = np.argsort(np.abs(np.linalg.eigvals(J8)))[::-1][i]
    vec = evecs[:, idx]
    # Check node symmetry: compare site-1 and site-2 parts
    v1, v2 = vec[:4], vec[4:]
    sym = np.linalg.norm(v1 + v2) / (np.linalg.norm(v1) + np.linalg.norm(v2) + 1e-30)
    asym = np.linalg.norm(v1 - v2) / (np.linalg.norm(v1) + np.linalg.norm(v2) + 1e-30)
    sym_type = "node-symmetric" if sym > asym else "node-antisymmetric"
    print(f"    Mode {i}: {sym_type} (sym={sym:.3f}, asym={asym:.3f})")

print(f"\n{'='*80}")
print("VERIFICATION COMPLETE")
print("=" * 80)
