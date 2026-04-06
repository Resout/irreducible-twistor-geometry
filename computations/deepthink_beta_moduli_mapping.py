#!/usr/bin/env python3
"""
SELF-CONTAINED COMPUTATION: Does the evidence distribution move along
the K*=7/30 moduli curve as the lattice coupling beta varies?

RUN THIS SCRIPT AS-IS. No modifications needed. Prints all results.

WHAT THIS COMPUTES:
  Part A: 4D SU(2) lattice simulation at 7 beta values
          100 configs each, L=4, plaquettes in different planes sharing links
          Extract conditional distribution C_beta, evidence (p_dom, phi)
          Run DS iteration, measure K*, lambda_0, Delta

  Part B: For each beta, map (p_dom, phi) onto the K*=7/30 moduli curve
          Check if position on curve changes with beta

  Part C: Analytical cross-check via strong-coupling expansion
          Compute <P_1 P_2>_connected as function of beta
          Predict how delta(beta) should vary

  Part D: If evidence moves with beta, compute Delta(beta) and check
          whether Delta(beta) * a(beta) gives finite physical mass

RUNTIME: ~60-90 minutes.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.special import iv as bessel_i  # modified Bessel function I_v(z)

# ============================================================================
# 4D SU(2) LATTICE
# ============================================================================
N_DIM = 4

def random_su2():
    x = np.random.randn(4)
    x /= np.linalg.norm(x)
    a, b, c, d = x
    return np.array([[a+1j*b, c+1j*d], [-c+1j*d, a-1j*b]])

def random_su2_near_id(eps=0.3):
    x = np.random.randn(4)
    x /= np.linalg.norm(x)
    x[0] = 1.0 / eps
    x /= np.linalg.norm(x)
    a, b, c, d = x
    return np.array([[a+1j*b, c+1j*d], [-c+1j*d, a-1j*b]])

class Lattice4D:
    def __init__(self, L):
        self.L = L
        self.links = np.zeros((N_DIM, L, L, L, L, 2, 2), dtype=complex)
        for mu in range(N_DIM):
            for x0 in range(L):
                for x1 in range(L):
                    for x2 in range(L):
                        for x3 in range(L):
                            self.links[mu, x0, x1, x2, x3] = random_su2()

    def shift(self, site, mu, d=1):
        s = list(site)
        s[mu] = (s[mu] + d) % self.L
        return tuple(s)

    def get(self, mu, site):
        return self.links[mu, site[0], site[1], site[2], site[3]]

    def put(self, mu, site, U):
        self.links[mu, site[0], site[1], site[2], site[3]] = U

    def plaq(self, site, mu, nu):
        x = site
        U1 = self.get(mu, x)
        U2 = self.get(nu, self.shift(x, mu))
        U3 = self.get(mu, self.shift(x, nu)).conj().T
        U4 = self.get(nu, x).conj().T
        return np.real(np.trace(U1 @ U2 @ U3 @ U4))

    def staple(self, site, mu):
        S = np.zeros((2, 2), dtype=complex)
        x = site
        x_mu = self.shift(x, mu)
        for nu in range(N_DIM):
            if nu == mu: continue
            x_nu = self.shift(x, nu)
            x_mu_mnu = self.shift(x_mu, nu, -1)
            x_mnu = self.shift(x, nu, -1)
            # Forward staple
            S += self.get(nu, x_mu) @ self.get(mu, x_nu).conj().T @ self.get(nu, x).conj().T
            # Backward staple
            S += self.get(nu, x_mu_mnu).conj().T @ self.get(mu, x_mnu).conj().T @ self.get(nu, x_mnu)
        return S

    def sweep(self, beta, eps=0.3):
        acc = 0
        tot = 0
        for mu in range(N_DIM):
            for x0 in range(self.L):
                for x1 in range(self.L):
                    for x2 in range(self.L):
                        for x3 in range(self.L):
                            site = (x0, x1, x2, x3)
                            old = self.get(mu, site).copy()
                            stap = self.staple(site, mu)
                            old_a = np.real(np.trace(old @ stap))
                            new = random_su2_near_id(eps) @ old
                            new_a = np.real(np.trace(new @ stap))
                            dS = beta * (new_a - old_a)
                            if dS > 0 or np.random.random() < np.exp(dS):
                                self.put(mu, site, new)
                                acc += 1
                            tot += 1
        return acc / tot

    def thermalize(self, beta, n=200):
        eps = min(0.5, 2.0 / max(beta, 0.1))
        rates = []
        for i in range(n):
            r = self.sweep(beta, eps)
            rates.append(r)
            # Adapt eps for ~50% acceptance
            if i > 10 and i % 20 == 0:
                avg_r = np.mean(rates[-20:])
                if avg_r > 0.55: eps *= 0.9
                elif avg_r < 0.45: eps *= 1.1
        return np.mean(rates[-20:]) if len(rates) >= 20 else np.mean(rates)

    def observables_shared_link(self):
        """Extract plaquettes in (0,1) and (0,2) planes at each site.
        These share the direction-0 link U_0(x)."""
        obs_01 = []
        obs_02 = []
        for x0 in range(self.L):
            for x1 in range(self.L):
                for x2 in range(self.L):
                    for x3 in range(self.L):
                        obs_01.append(self.plaq((x0,x1,x2,x3), 0, 1))
                        obs_02.append(self.plaq((x0,x1,x2,x3), 0, 2))
        return np.array(obs_01), np.array(obs_02)

    def observables_multi_pair(self):
        """Extract multiple correlated observable pairs.
        Pairs: (0,1)-(0,2), (0,1)-(0,3), (0,2)-(0,3), (1,2)-(1,3)
        Each pair shares one link direction."""
        pairs = [(0,1,0,2), (0,1,0,3), (0,2,0,3), (1,2,1,3)]
        all_obs1 = []
        all_obs2 = []
        for mu1, nu1, mu2, nu2 in pairs:
            for x0 in range(self.L):
                for x1 in range(self.L):
                    for x2 in range(self.L):
                        for x3 in range(self.L):
                            all_obs1.append(self.plaq((x0,x1,x2,x3), mu1, nu1))
                            all_obs2.append(self.plaq((x0,x1,x2,x3), mu2, nu2))
        return np.array(all_obs1), np.array(all_obs2)

    def observables_timeslice(self):
        """Time-slice averaged plaquettes — the paper's primary observable.
        O1(t) = (1/L^3) sum_{x1,x2,x3} P_{01}(t,x1,x2,x3)
        O2(t) = (1/L^3) sum_{x1,x2,x3} P_{02}(t,x1,x2,x3)
        These share direction-0 links. Spatial averaging amplifies signal by L^3.
        Returns L pairs per config (one per time slice)."""
        L = self.L
        obs1 = []
        obs2 = []
        for t in range(L):
            # Average P_{01} over spatial volume at time t
            sum_01 = 0
            sum_02 = 0
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        sum_01 += self.plaq((t, x1, x2, x3), 0, 1)
                        sum_02 += self.plaq((t, x1, x2, x3), 0, 2)
            obs1.append(sum_01 / L**3)
            obs2.append(sum_02 / L**3)

        # Also add (0,1)-(0,3) and (0,2)-(0,3) pairs for more statistics
        for t in range(L):
            sum_01 = 0
            sum_03 = 0
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        sum_01 += self.plaq((t, x1, x2, x3), 0, 1)
                        sum_03 += self.plaq((t, x1, x2, x3), 0, 3)
            obs1.append(sum_01 / L**3)
            obs2.append(sum_03 / L**3)

        for t in range(L):
            sum_02 = 0
            sum_03 = 0
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        sum_02 += self.plaq((t, x1, x2, x3), 0, 2)
                        sum_03 += self.plaq((t, x1, x2, x3), 0, 3)
            obs1.append(sum_02 / L**3)
            obs2.append(sum_03 / L**3)

        return np.array(obs1), np.array(obs2)

# ============================================================================
# DS FRAMEWORK
# ============================================================================
H = 3
FLOOR = 1.0 / H**3

def ds_combine(m, e):
    s, th = m[:3], m[3]
    se, ph = e[:3], e[3]
    sn = np.array([s[i]*se[i] + s[i]*ph + th*se[i] for i in range(3)])
    tn = th * ph
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)
    d = 1.0 - K
    if abs(d) < 1e-15: return m.copy(), K
    out = np.zeros(4)
    out[:3] = sn / d
    out[3] = tn / d
    out /= np.sum(out)
    return out, K

def enforce_floor(m):
    s, th = m[:3], m[3]
    ssq = np.sum(s**2)
    if th**2 / (ssq + th**2) >= FLOOR - 1e-15: return m.copy()
    ss = np.sum(s)
    if ss < 1e-15: return m.copy()
    r = ssq / ss**2
    disc = (2*r)**2 + 4*(26-r)*r
    t = (-2*r + np.sqrt(max(disc, 0))) / (2*(26-r))
    alpha = (1-t) / ss
    return np.array([s[0]*alpha, s[1]*alpha, s[2]*alpha, t])

def full_step(m, e):
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds), K

def find_fp(e, n=2000):
    m = np.array([0.7, 0.05, 0.05, 0.2])
    for _ in range(n):
        m2, K = full_step(m, e)
        if np.max(np.abs(m2 - m)) < 1e-14: return m, True
        m = m2
    return m, False

def measure_K(m, e):
    return sum(m[i]*e[j] for i in range(3) for j in range(3) if i != j)

def compute_eigs(m, e, eps=1e-8):
    m0, _ = full_step(m, e)
    J = np.zeros((3, 3))
    idx = [0, 1, 3]
    for j in range(3):
        mp = m.copy(); mp[idx[j]] += eps; mp[3 if idx[j] != 3 else 0] -= eps
        mo, _ = full_step(mp, e)
        for i in range(3):
            J[i, j] = (mo[idx[i]] - m0[idx[i]]) / eps
    return np.sort(np.abs(np.linalg.eigvals(J)))[::-1]

# ============================================================================
# CONDITIONAL DISTRIBUTION AND EVIDENCE EXTRACTION
# ============================================================================
def build_C(obs1, obs2):
    n = len(obs1)
    def tercile(vals):
        si = np.argsort(vals)
        b = np.zeros(len(vals), dtype=int)
        for rank, i in enumerate(si):
            b[i] = min(rank * H // n, H - 1)
        return b
    b1, b2 = tercile(obs1), tercile(obs2)
    C = np.zeros((H, H))
    for k in range(n):
        C[b1[k], b2[k]] += 1
    rs = C.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1
    return C / rs

def evidence_from_C_row(C, row_idx):
    """Convert one row of C into a DS evidence mass function.
    s_j = C[row_idx, j] * (1 - theta), theta = H(row)/log(H)."""
    row = C[row_idx]
    entropy = -np.sum(row * np.log(row + 1e-15)) / np.log(H)
    theta = entropy
    e = np.zeros(H + 1)
    for j in range(H):
        e[j] = row[j] * (1 - theta)
    e[H] = theta
    e /= np.sum(e)
    return e

def iterative_extraction(C, obs1, obs2, n_steps=50, n_runs=20):
    """Paper's iterative DS extraction procedure (Section 8).

    Starting from complete ignorance, at each step:
      1. Sample a cell (i,j) from C proportional to frequency
      2. Construct evidence mass from row i of C
      3. Combine running state with evidence via DS + floor
      4. Record K

    Returns K* (mean |K| over steps 35-50), the final running state,
    and the evidence that was used."""

    # Build joint frequency matrix for sampling
    # (approximation: use C itself as joint distribution)
    # Actually, we need the joint distribution P(i,j), not just conditional
    # Assume roughly uniform marginal for O1 (tercile binning ensures this)
    joint = C / H  # P(i,j) ≈ P(j|i) * P(i) = C[i,j] * (1/H)
    joint_flat = joint.flatten()
    joint_flat /= np.sum(joint_flat)

    all_K_stars = []
    all_final_m = []
    all_final_e = []

    for run in range(n_runs):
        # Start from complete ignorance
        m = np.array([0.0, 0.0, 0.0, 1.0])
        m[:3] = 1e-10  # tiny singletons to avoid division issues
        m /= np.sum(m)

        K_history = []

        for step in range(n_steps):
            # Sample cell (i,j) from joint distribution
            cell = np.random.choice(H*H, p=joint_flat)
            i_row = cell // H

            # Construct evidence from row i of C
            e = evidence_from_C_row(C, i_row)

            # Combine with running state
            m, K = full_step(m, e)
            K_history.append(abs(K))

        # Extract K* from steps 35-50
        if len(K_history) >= 50:
            K_star = np.mean(K_history[35:50])
        else:
            K_star = np.mean(K_history[-10:])

        all_K_stars.append(K_star)
        all_final_m.append(m.copy())
        all_final_e.append(e.copy())

    K_star_mean = np.mean(all_K_stars)
    K_star_std = np.std(all_K_stars)

    # The "evidence at convergence" is the average of the final evidence vectors
    avg_final_e = np.mean(all_final_e, axis=0)
    avg_final_e /= np.sum(avg_final_e)
    avg_final_m = np.mean(all_final_m, axis=0)
    avg_final_m /= np.sum(avg_final_m)

    phi = avg_final_e[H]
    s_sum = np.sum(avg_final_e[:H])
    p_dom = np.max(avg_final_e[:H]) / s_sum if s_sum > 0 else 1.0/H
    delta = np.trace(C) / H

    return K_star_mean, K_star_std, avg_final_m, avg_final_e, phi, p_dom, delta

# ============================================================================
# PART A: 4D LATTICE SIMULATIONS
# ============================================================================
print("=" * 100)
print("PART A: 4D SU(2) LATTICE SIMULATIONS AT MULTIPLE BETA VALUES")
print("=" * 100)

L = 4
beta_values = [1.5, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0]
n_configs = 80
n_therm = 250

np.random.seed(12345)

print(f"\nL = {L} (4D: {L}^4 = {L**4} sites)")
print(f"Beta values: {beta_values}")
print(f"Configs per beta: {n_configs}")
print(f"Thermalization: {n_therm} sweeps")
print(f"Observable pairs: 4 plane pairs sharing one link direction")
print(f"Samples per config: {L**4 * 4} = {L**4 * 4}")
print(f"Total samples per beta: {n_configs * L**4 * 4}")

lattice_results = []

for beta in beta_values:
    print(f"\n--- beta = {beta:.1f} ---")
    all_o1, all_o2 = [], []       # time-slice averaged (main observable)
    all_o1_raw, all_o2_raw = [], []  # site-by-site (for comparison)
    accept_rates = []

    for c in range(n_configs):
        lat = Lattice4D(L)
        ar = lat.thermalize(beta, n=n_therm)
        accept_rates.append(ar)
        # Time-slice averaged observables (paper's method — strong signal)
        o1, o2 = lat.observables_timeslice()
        all_o1.extend(o1)
        all_o2.extend(o2)
        # Site-by-site observables (for comparison — weak signal)
        o1r, o2r = lat.observables_shared_link()
        all_o1_raw.extend(o1r)
        all_o2_raw.extend(o2r)
        if (c+1) % 20 == 0:
            print(f"  config {c+1}/{n_configs}, accept_rate={ar:.3f}", flush=True)

    obs1 = np.array(all_o1)       # time-slice averaged
    obs2 = np.array(all_o2)
    obs1r = np.array(all_o1_raw)  # site-by-site
    obs2r = np.array(all_o2_raw)
    n_samples = len(obs1)
    n_samples_raw = len(obs1r)
    mean_plaq = np.mean(obs1r)
    corr_ts = np.corrcoef(obs1, obs2)[0, 1]
    corr_raw = np.corrcoef(obs1r, obs2r)[0, 1]
    mean_ar = np.mean(accept_rates)

    print(f"  Timeslice: {n_samples} pairs, corr={corr_ts:.6f}")
    print(f"  Site-raw:  {n_samples_raw} pairs, corr={corr_raw:.6f}")
    corr = corr_ts  # use timeslice for main analysis

    # Build conditional distribution from TIME-SLICE averaged observables
    C = build_C(obs1, obs2)
    # Also build from raw for comparison
    C_raw = build_C(obs1r, obs2r)

    # Paper's iterative DS extraction (50 steps, 20 independent runs)
    K_val, K_std, m_fp, e_vec, phi, p_dom, delta = iterative_extraction(
        C, obs1, obs2, n_steps=50, n_runs=20)

    # Also find the fixed point with the averaged evidence
    m_fp2, conv = find_fp(e_vec)
    K_val2 = measure_K(m_fp2, e_vec) if conv else K_val

    # Eigenvalues at the fixed point
    lam0, lam1, gap = 0, 0, 0
    if conv and 0 < K_val2 < 1:
        try:
            eigs = compute_eigs(m_fp2, e_vec)
            lam0 = eigs[0]
            lam1 = eigs[1] if len(eigs) > 1 else 0
            gap = -np.log(lam0) if 0 < lam0 < 1 else 0
        except:
            pass

    lattice_results.append({
        'beta': beta, 'L': L,
        'n_samples': n_samples, 'mean_plaq': mean_plaq,
        'corr': corr, 'accept_rate': mean_ar,
        'delta': delta, 'phi': phi, 'p_dom': p_dom,
        'K': K_val, 'K_std': K_std, 'K_fp': K_val2,
        'lam0': lam0, 'lam1': lam1, 'gap': gap,
        'conv': conv, 'C': C, 'e': e_vec, 'm': m_fp2,
    })

    print(f"  samples={n_samples}, <P>={mean_plaq:.4f}, corr={corr:.6f}, accept={mean_ar:.3f}")
    print(f"  C =")
    for row in C:
        print(f"    [{row[0]:.4f} {row[1]:.4f} {row[2]:.4f}]")
    delta_raw = np.trace(C_raw) / H
    print(f"  delta(timeslice)={delta:.6f}, delta(raw)={delta_raw:.6f}")
    print(f"  K*(iterative) = {K_val:.6f} +/- {K_std:.6f} (7/30={7/30:.6f})")
    print(f"  K*(fixed-pt)  = {K_val2:.6f}, conv={conv}")
    print(f"  phi={phi:.6f}, p_dom={p_dom:.6f}")
    print(f"  lam0={lam0:.6f}, Delta={gap:.4f}")


# ============================================================================
# PART B: MAPPING ONTO THE MODULI CURVE
# ============================================================================
print("\n\n" + "=" * 100)
print("PART B: MAPPING ONTO THE K*=7/30 MODULI CURVE")
print("=" * 100)

print(f"\n{'beta':>6s} {'delta':>8s} {'phi':>8s} {'p_dom':>8s} {'K*':>10s} "
      f"{'K*-7/30':>10s} {'lam0':>8s} {'Delta':>8s} {'<P>':>8s} {'corr':>8s}")
print("-" * 100)

for r in lattice_results:
    marker = " <--" if abs(r['K'] - 7/30) < 0.02 else ""
    print(f"{r['beta']:6.1f} {r['delta']:8.4f} {r['phi']:8.4f} {r['p_dom']:8.4f} "
          f"{r['K']:10.6f} {r.get('K_std',0):8.4f} {r['lam0']:8.4f} {r['gap']:8.4f} "
          f"{r['mean_plaq']:8.4f} {r['corr']:8.6f}{marker}")

# Trend analysis
print("\n\nTREND ANALYSIS:")
betas = np.array([r['beta'] for r in lattice_results])
phis = np.array([r['phi'] for r in lattice_results])
pdoms = np.array([r['p_dom'] for r in lattice_results])
deltas = np.array([r['delta'] for r in lattice_results])
Ks = np.array([r['K'] for r in lattice_results])
corrs = np.array([r['corr'] for r in lattice_results])
lam0s = np.array([r['lam0'] for r in lattice_results])
gaps = np.array([r['gap'] for r in lattice_results])
plaqs = np.array([r['mean_plaq'] for r in lattice_results])

if len(betas) > 2:
    for name, vals in [('delta', deltas), ('phi', phis), ('p_dom', pdoms),
                        ('K*', Ks), ('corr', corrs), ('lam0', lam0s), ('Delta', gaps)]:
        c = np.polyfit(betas, vals, 1)
        print(f"  {name:>8s}(beta) = {c[0]:+.6f}*beta + {c[1]:.6f}  "
              f"(slope {'significant' if abs(c[0]) > np.std(vals)/3 else 'flat'})")


# ============================================================================
# PART C: ANALYTICAL CROSS-CHECK (strong-coupling expansion)
# ============================================================================
print("\n\n" + "=" * 100)
print("PART C: ANALYTICAL CROSS-CHECK")
print("=" * 100)

print("""
For SU(2), the character expansion of the Boltzmann weight is:
  exp(beta/2 * Tr U) = sum_j (2j+1) * a_j(beta) * chi_j(U)
  where a_j(beta) = I_{2j+1}(beta) / (beta/2)  [modified Bessel function]

The plaquette expectation:
  <Tr U_P / 2> = I_2(beta) / I_1(beta)  [exact for fundamental rep]

The connected correlator of two plaquettes sharing one link:
  <P1*P2>_c / (<P1><P2>) ~ I_2(beta)^2 / I_1(beta)^2  [leading order]

This ratio controls the correlation strength and hence delta.
""")

print(f"{'beta':>6s} {'I1':>10s} {'I2':>10s} {'I2/I1':>10s} {'<P>_theory':>12s} {'<P>_lattice':>12s}")
print("-" * 65)

for i, beta in enumerate(beta_values):
    I1 = bessel_i(1, beta)
    I2 = bessel_i(2, beta)
    ratio = I2 / I1
    plaq_theory = ratio  # <Tr U / 2> for fundamental rep
    plaq_lattice = lattice_results[i]['mean_plaq'] / 2  # our plaq is Tr, not Tr/2
    print(f"{beta:6.1f} {I1:10.6f} {I2:10.6f} {ratio:10.6f} {plaq_theory:12.6f} {plaq_lattice:12.6f}")

# The I2/I1 ratio controls the nearest-neighbor correlation
# As beta increases: I2/I1 -> 1 (weak coupling, all links near identity)
# As beta decreases: I2/I1 -> 0 (strong coupling, random links)
print(f"\nI_2(beta)/I_1(beta) ranges from {bessel_i(2,beta_values[0])/bessel_i(1,beta_values[0]):.4f} "
      f"to {bessel_i(2,beta_values[-1])/bessel_i(1,beta_values[-1]):.4f}")
print(f"This ratio controls the structural correlation strength.")
print(f"If it changes with beta, the evidence distribution changes,")
print(f"and the position on the K*=7/30 curve changes.")


# ============================================================================
# PART D: PHYSICAL MASS GAP
# ============================================================================
print("\n\n" + "=" * 100)
print("PART D: PHYSICAL MASS GAP VS BETA")
print("=" * 100)

print("""
If Delta(beta) varies along the moduli curve:
  m_phys = Delta(beta) / a(beta)
  where a(beta) = lattice spacing from asymptotic freedom.

For SU(2): a(beta) ~ Lambda_L^{-1} * exp(-pi^2 * beta / 11)  [2-loop]
  (Lambda_L = Lambda parameter in lattice units)

If Delta(beta) * a(beta) = constant, then m_phys is finite and
the framework has a well-defined continuum limit.
""")

# Compute a(beta) relative to a reference (beta=2.3 as in paper)
beta_ref = 2.3
b0 = 11 / (24 * np.pi**2)  # 1-loop coefficient for SU(2)

print(f"{'beta':>6s} {'Delta':>8s} {'a/a_ref':>10s} {'Delta*a':>10s} {'ratio':>10s}")
print("-" * 50)

Delta_a_values = []
for r in lattice_results:
    if r['gap'] > 0:
        # Relative lattice spacing (1-loop asymptotic freedom)
        a_ratio = np.exp(-np.pi**2 * (r['beta'] - beta_ref) / 11)
        Delta_a = r['gap'] * a_ratio
        Delta_a_values.append(Delta_a)
        print(f"{r['beta']:6.1f} {r['gap']:8.4f} {a_ratio:10.4f} {Delta_a:10.4f} "
              f"{Delta_a/Delta_a_values[0]:10.4f}")

if Delta_a_values:
    cv = np.std(Delta_a_values) / np.mean(Delta_a_values)
    print(f"\nDelta*a: mean={np.mean(Delta_a_values):.4f}, std={np.std(Delta_a_values):.4f}, CV={cv:.4f}")
    if cv < 0.1:
        print(">>> Delta*a is approximately CONSTANT — continuum limit exists!")
    elif cv < 0.3:
        print(">>> Delta*a varies moderately — partial scaling")
    else:
        print(">>> Delta*a varies significantly — no simple scaling")


# ============================================================================
# PART E: SUMMARY
# ============================================================================
print("\n\n" + "=" * 100)
print("SUMMARY")
print("=" * 100)

print(f"""
QUESTION: Does the evidence distribution move along the K*=7/30
moduli curve as the lattice coupling beta varies?

LATTICE DATA ({len(lattice_results)} beta values, L={L}, {n_configs} configs each):
  delta range: [{min(deltas):.4f}, {max(deltas):.4f}]
  phi range:   [{min(phis):.4f}, {max(phis):.4f}]
  p_dom range: [{min(pdoms):.4f}, {max(pdoms):.4f}]
  K* range:    [{min(Ks):.4f}, {max(Ks):.4f}]
  corr range:  [{min(corrs):.6f}, {max(corrs):.6f}]

BESSEL FUNCTION ANALYSIS:
  I_2/I_1 ranges from {bessel_i(2,1.5)/bessel_i(1,1.5):.4f} (beta=1.5)
  to {bessel_i(2,4.0)/bessel_i(1,4.0):.4f} (beta=4.0)
  This ratio DOES change with beta, confirming that the structural
  correlation strength is beta-dependent.

INTERPRETATION:
  If phi and p_dom are constant across beta:
    -> Position on moduli curve is FIXED by gauge group
    -> No RG flow along the curve
    -> Path (B) is correct: the theory IS the continuum theory

  If phi and p_dom change with beta:
    -> Position on curve MOVES with beta
    -> The moduli curve IS the RG flow
    -> Delta runs with beta, potentially giving a continuum limit
""")

print("COMPUTATION COMPLETE.")
