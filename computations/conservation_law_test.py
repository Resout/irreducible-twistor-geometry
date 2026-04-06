"""
Conservation Law Test: Is K* = 7/30 a CONSEQUENCE of the fixed-point equations?

The question: at the S₂-symmetric equilibrium (s₂=s₃), the system has
6 unknowns (s₁, s₂, θ, e₁, e₂, φ) and potentially 6 constraints:
  1. L₁(m*) = 1
  2. L₁(e*) = 1
  3. Born(m*) = 1/27
  4-6. Fixed-point: m* = Φ(m*, e*) [3 component equations with floor]

If this system determines K* = 7/30 as OUTPUT (not input), then the
conservation law is a THEOREM of the fixed-point equations.

If the system is underdetermined and K* requires additional input,
then the conservation law adds independent information.

We test this by:
  (A) Solving the fixed-point system numerically WITHOUT imposing K*=7/30
  (B) Checking what K* emerges
  (C) Sweeping over evidence parameters to see if K*=7/30 is the unique
      self-consistent solution or one of many
"""

import numpy as np
from scipy.optimize import fsolve, minimize

H = 3
FLOOR = 1.0 / H**3  # 1/27

# ============================================================
# DS combination + floor (real positive masses)
# ============================================================
def ds_combine_raw(m, e):
    """DS combination WITHOUT floor. Returns (m_out, K)."""
    s, theta = m[:3], m[3]
    se, phi = e[:3], e[3]
    s_new = s * se + s * phi + theta * se
    theta_new = theta * phi
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    m_out = np.zeros(4)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    # L1 normalise
    total = np.sum(m_out)
    if total > 0:
        m_out = m_out / total
    return m_out, K

def enforce_floor(m):
    """Born floor enforcement for real positive masses.
    Uses the quadratic formula from popov_second_variation.py."""
    s, theta = m[:3], m[3]
    ssq = sum(si**2 for si in s)
    total_sq = ssq + theta**2
    born = theta**2 / total_sq if total_sq > 0 else 1
    if born >= FLOOR:
        return m.copy()
    ss = sum(s)
    if ss < 1e-15:
        return m.copy()
    r = ssq / ss**2
    # Solve (26-r)t² + 2rt - r = 0
    a_coeff = 26 - r
    b_coeff = 2 * r
    c_coeff = -r
    disc = b_coeff**2 - 4 * a_coeff * c_coeff
    t = (-b_coeff + np.sqrt(disc)) / (2 * a_coeff)
    alpha = (1 - t) / ss
    return np.array([s[0]*alpha, s[1]*alpha, s[2]*alpha, t])

def full_step(m, e):
    """One DS step + floor."""
    m_ds, K = ds_combine_raw(m, e)
    m_floored = enforce_floor(m_ds)
    return m_floored, K

def born(m):
    return m[3]**2 / np.sum(m**2) if np.sum(m**2) > 0 else 1

def K_conflict(m, e):
    return sum(m[i] * e[j] for i in range(3) for j in range(3) if i != j)

# ============================================================
# TEST A: Direct fixed-point solve (no K* constraint)
# ============================================================
# S₂-symmetric: s₂ = s₃, e₂ = e₃
# Unknowns: s₁, s₂, θ, e₁, e₂, φ (6 unknowns)
# Equations:
#   1. s₁ + 2s₂ + θ = 1             (L₁ state)
#   2. e₁ + 2e₂ + φ = 1             (L₁ evidence)
#   3. Born(m*) = 1/27               (Born floor active)
#   4-6. m* = Φ(m*, e*)             (fixed point, 3 components)

def residual_system(x):
    """
    x = [s1, s2, theta, e1, e2, phi]
    Returns 6 residuals that should all be zero at the fixed point.
    """
    s1, s2, theta, e1, e2, phi = x

    # Ensure positivity
    if any(v <= 0 for v in [s1, s2, theta, e1, e2, phi]):
        return [1e10] * 6

    m = np.array([s1, s2, s2, theta])
    e = np.array([e1, e2, e2, phi])

    # Constraint 1: L1(m) = 1
    res1 = s1 + 2*s2 + theta - 1.0

    # Constraint 2: L1(e) = 1
    res2 = e1 + 2*e2 + phi - 1.0

    # Constraint 3: Born(m*) = 1/27
    total_sq = s1**2 + 2*s2**2 + theta**2
    res3 = theta**2 / total_sq - FLOOR

    # Constraints 4-6: m* = Φ(m*, e*)
    m_out, K = full_step(m, e)
    # The output should equal the input
    res4 = m_out[0] - s1    # s1 component
    res5 = m_out[1] - s2    # s2 component
    res6 = m_out[3] - theta  # θ component

    return [res1, res2, res3, res4, res5, res6]

print("=" * 70)
print("TEST A: Solve fixed-point equations WITHOUT imposing K* = 7/30")
print("=" * 70)

# Try multiple initial guesses to be thorough
initial_guesses = [
    [0.75, 0.04, 0.15, 0.60, 0.12, 0.13],   # near known solution
    [0.50, 0.10, 0.20, 0.50, 0.10, 0.20],   # symmetric-ish
    [0.80, 0.02, 0.14, 0.70, 0.08, 0.10],   # concentrated
    [0.60, 0.08, 0.16, 0.40, 0.15, 0.15],   # moderate
    [0.85, 0.01, 0.12, 0.80, 0.05, 0.05],   # very concentrated
    [0.65, 0.06, 0.17, 0.55, 0.10, 0.15],   # another moderate
]

solutions = []
for i, x0 in enumerate(initial_guesses):
    try:
        sol = fsolve(residual_system, x0, full_output=True)
        x_sol, info, ier, msg = sol
        s1, s2, theta, e1, e2, phi = x_sol

        if ier == 1 and all(v > 0 for v in x_sol):
            m_sol = np.array([s1, s2, s2, theta])
            e_sol = np.array([e1, e2, e2, phi])
            K_val = K_conflict(m_sol, e_sol)

            residuals = residual_system(x_sol)
            max_res = max(abs(r) for r in residuals)

            if max_res < 1e-8:
                print(f"\nGuess {i+1}: CONVERGED (max residual = {max_res:.2e})")
                print(f"  m* = ({s1:.8f}, {s2:.8f}, {s2:.8f}, {theta:.8f})")
                print(f"  e* = ({e1:.8f}, {e2:.8f}, {e2:.8f}, {phi:.8f})")
                print(f"  K*  = {K_val:.10f}")
                print(f"  7/30 = {7/30:.10f}")
                print(f"  K* - 7/30 = {K_val - 7/30:.2e}")
                print(f"  Born(m*) = {born(m_sol):.10f} (target: {FLOOR:.10f})")

                # Verify it's actually a fixed point
                m_check, K_check = full_step(m_sol, e_sol)
                fp_err = np.max(np.abs(m_check - m_sol))
                print(f"  Fixed-point error: {fp_err:.2e}")

                solutions.append((x_sol, K_val, max_res))
            else:
                print(f"\nGuess {i+1}: converged but high residual = {max_res:.2e}")
        else:
            print(f"\nGuess {i+1}: did not converge ({msg.strip()})")
    except Exception as ex:
        print(f"\nGuess {i+1}: exception ({ex})")

# ============================================================
# TEST B: Is the solution unique? Sweep over evidence parameter
# ============================================================
print("\n" + "=" * 70)
print("TEST B: Uniqueness — sweep evidence concentration, find all fixed points")
print("=" * 70)

# For each value of e₁ (dominant evidence), solve for the rest
print(f"\n{'e1':>8s} {'K*':>12s} {'K*-7/30':>12s} {'theta':>10s} {'born':>10s} {'fp_err':>10s}")
print("-" * 70)

e1_values = np.linspace(0.40, 0.90, 51)
K_results = []

for e1_target in e1_values:
    # With e₁ fixed, we have 5 unknowns: s1, s2, θ, e2, φ
    # And 5 equations: L₁(m), L₁(e), Born(m), and 2 FP equations (s1, θ)
    # (s2 FP follows from S₂ symmetry + L₁)

    def residual_5(x):
        s1, s2, theta, e2, phi = x
        if any(v <= 0 for v in [s1, s2, theta, e2, phi]):
            return [1e10] * 5

        e1 = e1_target
        m = np.array([s1, s2, s2, theta])
        e = np.array([e1, e2, e2, phi])

        res1 = s1 + 2*s2 + theta - 1.0
        res2 = e1 + 2*e2 + phi - 1.0
        res3 = theta**2 / (s1**2 + 2*s2**2 + theta**2) - FLOOR

        m_out, K = full_step(m, e)
        res4 = m_out[0] - s1
        res5 = m_out[3] - theta

        return [res1, res2, res3, res4, res5]

    best = None
    for s1_0, s2_0, theta_0, e2_0, phi_0 in [
        (0.75, 0.04, 0.15, 0.12, 0.13),
        (0.80, 0.02, 0.14, 0.08, 0.10),
        (0.60, 0.08, 0.16, 0.15, 0.20),
    ]:
        try:
            sol = fsolve(residual_5, [s1_0, s2_0, theta_0, e2_0, phi_0],
                        full_output=True)
            x_s, info, ier, msg = sol
            if ier == 1 and all(v > 0 for v in x_s):
                res = residual_5(x_s)
                mr = max(abs(r) for r in res)
                if mr < 1e-8:
                    s1, s2, theta, e2, phi = x_s
                    m = np.array([s1, s2, s2, theta])
                    e = np.array([e1_target, e2, e2, phi])
                    K_val = K_conflict(m, e)
                    m_out, _ = full_step(m, e)
                    fp_err = np.max(np.abs(m_out - m))
                    if best is None or mr < best[1]:
                        best = (K_val, mr, theta, born(m), fp_err)
        except:
            pass

    if best is not None:
        K_val, mr, theta, b, fpe = best
        K_results.append((e1_target, K_val))
        marker = " <--- 7/30" if abs(K_val - 7/30) < 0.001 else ""
        print(f"{e1_target:8.4f} {K_val:12.8f} {K_val-7/30:12.2e} {theta:10.6f} {b:10.6f} {fpe:10.2e}{marker}")
    else:
        print(f"{e1_target:8.4f}    no solution found")

# ============================================================
# TEST C: Count degrees of freedom
# ============================================================
print("\n" + "=" * 70)
print("TEST C: Jacobian rank at solution — is the system fully determined?")
print("=" * 70)

if solutions:
    x_best = solutions[0][0]

    # Compute numerical Jacobian of the 6x6 system
    eps = 1e-8
    J = np.zeros((6, 6))
    f0 = residual_system(list(x_best))
    for j in range(6):
        x_pert = list(x_best)
        x_pert[j] += eps
        f_pert = residual_system(x_pert)
        for i in range(6):
            J[i, j] = (f_pert[i] - f0[i]) / eps

    sv = np.linalg.svd(J, compute_uv=False)
    rank = np.sum(sv > 1e-6)

    print(f"Jacobian singular values: {sv}")
    print(f"Numerical rank: {rank} (out of 6)")
    print(f"Condition number: {sv[0]/sv[-1]:.2e}" if sv[-1] > 0 else "Singular!")

    if rank == 6:
        print("\n>>> SYSTEM IS FULLY DETERMINED: K* = 7/30 is a THEOREM,")
        print("    not an independent axiom. The conservation law is a consequence.")
    elif rank == 5:
        print("\n>>> SYSTEM HAS 1 FREE PARAMETER: K* requires additional input.")
        print("    The conservation law adds independent information.")
    else:
        print(f"\n>>> SYSTEM HAS {6-rank} FREE PARAMETERS")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

if solutions:
    K_vals = [s[1] for s in solutions]
    K_mean = np.mean(K_vals)
    print(f"All converged solutions give K* = {K_mean:.10f}")
    print(f"7/30 = {7/30:.10f}")
    print(f"Deviation: {abs(K_mean - 7/30):.2e}")

if K_results:
    K_at_fp = [k for e, k in K_results]
    print(f"\nAcross e₁ sweep ({len(K_results)} solutions):")
    print(f"  K* range: [{min(K_at_fp):.8f}, {max(K_at_fp):.8f}]")
    print(f"  K* std:   {np.std(K_at_fp):.2e}")
    if np.std(K_at_fp) < 1e-6:
        print("  >>> K* is CONSTANT across all evidence concentrations")
        print(f"  >>> K* = {np.mean(K_at_fp):.10f} = 7/30? {abs(np.mean(K_at_fp) - 7/30) < 1e-6}")
    else:
        print("  >>> K* VARIES with evidence concentration")
        print("  >>> Conservation law adds independent information")
