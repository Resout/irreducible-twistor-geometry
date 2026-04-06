#!/usr/bin/env python3
"""
SELF-CONTAINED COMPUTATION: Does the evidence distribution move along
the K*=7/30 curve as the lattice coupling beta varies?

RUN THIS SCRIPT AS-IS. No modifications needed. Prints all results.

WHAT THIS COMPUTES:
  For SU(2) lattice gauge theory at multiple beta values (1.5 to 4.0),
  at multiple lattice sizes (L=4,6,8,10):

  1. Generate thermalized lattice configurations
  2. Extract plaquette observables in two orthogonal planes
  3. Build the 3x3 conditional distribution C_beta
  4. Convert C_beta to DS evidence mass function e*
  5. Run DS iteration to find the fixed point m*
  6. Measure K*, lambda_0, Delta, and the evidence parameters (p_dom, phi)
  7. Check: does (p_dom, phi) move along the K*=7/30 curve as beta changes?

  If the evidence distribution moves along the curve with beta,
  then the K*=7/30 curve IS the renormalization group flow.

OUTPUT: Tables of numbers. No plots. Every result printed.

RUNTIME ESTIMATE: ~2-4 hours for full scan (4D lattice, multiple beta, multiple L).
"""

import numpy as np
from collections import Counter

# ============================================================================
# 4D SU(2) LATTICE GAUGE THEORY: Wilson action with Metropolis updates
# ============================================================================
# The 4D lattice is essential: plaquettes in different planes share links,
# producing the structural correlation (delta > 1/3) that gives K* ~ 7/30.
# A 2D lattice gives delta ~ 1/3 (no structural correlation) and K* ~ 0.

N_DIM = 4  # spacetime dimensions

def random_su2():
    """Random SU(2) matrix via Haar measure (quaternion method)."""
    x = np.random.randn(4)
    x /= np.linalg.norm(x)
    a, b, c, d = x
    return np.array([
        [a + 1j*b, c + 1j*d],
        [-c + 1j*d, a - 1j*b]
    ])


def random_su2_near_identity(eps=0.3):
    """Random SU(2) matrix near identity for Metropolis proposals."""
    x = np.random.randn(4)
    x /= np.linalg.norm(x)
    x[0] = 1.0 / eps
    x /= np.linalg.norm(x)
    a, b, c, d = x
    return np.array([
        [a + 1j*b, c + 1j*d],
        [-c + 1j*d, a - 1j*b]
    ])


class Lattice4D:
    """4D SU(2) lattice with periodic boundary conditions."""

    def __init__(self, L):
        self.L = L
        self.vol = L**4
        # links[mu][x0,x1,x2,x3] = SU(2) matrix
        self.links = np.zeros((N_DIM, L, L, L, L, 2, 2), dtype=complex)
        # Initialize with random SU(2)
        for mu in range(N_DIM):
            for x0 in range(L):
                for x1 in range(L):
                    for x2 in range(L):
                        for x3 in range(L):
                            self.links[mu, x0, x1, x2, x3] = random_su2()

    def shift(self, site, mu, direction=1):
        """Shift site by +1 or -1 in direction mu with periodic BC."""
        s = list(site)
        s[mu] = (s[mu] + direction) % self.L
        return tuple(s)

    def get_link(self, mu, site):
        """Get U_mu(site)."""
        return self.links[mu, site[0], site[1], site[2], site[3]]

    def set_link(self, mu, site, U):
        """Set U_mu(site) = U."""
        self.links[mu, site[0], site[1], site[2], site[3]] = U

    def plaquette(self, site, mu, nu):
        """Compute Tr[U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag]."""
        x = site
        x_mu = self.shift(x, mu)
        x_nu = self.shift(x, nu)

        U1 = self.get_link(mu, x)
        U2 = self.get_link(nu, x_mu)
        U3 = self.get_link(mu, x_nu).conj().T
        U4 = self.get_link(nu, x).conj().T

        return np.real(np.trace(U1 @ U2 @ U3 @ U4))

    def staple_sum(self, site, mu):
        """Sum of staples for link U_mu(site).
        Each staple is the product of the 3 other links of a plaquette
        containing U_mu(site). In 4D, there are 2*(N_DIM-1) = 6 staples."""
        S = np.zeros((2, 2), dtype=complex)
        x = site
        x_mu = self.shift(x, mu)

        for nu in range(N_DIM):
            if nu == mu:
                continue
            # Forward staple: U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag
            x_nu = self.shift(x, nu)
            S += (self.get_link(nu, x_mu) @
                  self.get_link(mu, x_nu).conj().T @
                  self.get_link(nu, x).conj().T)

            # Backward staple: U_nu(x+mu-nu)^dag U_mu(x-nu)^dag U_nu(x-nu)
            x_mnu = self.shift(x, nu, -1)
            x_mu_mnu = self.shift(x_mu, nu, -1)
            S += (self.get_link(nu, x_mu_mnu).conj().T @
                  self.get_link(mu, x_mnu).conj().T @
                  self.get_link(nu, x_mnu))

        return S

    def metropolis_sweep(self, beta, eps=0.3):
        """One Metropolis sweep over all links."""
        accepted = 0
        total = 0
        for mu in range(N_DIM):
            for x0 in range(self.L):
                for x1 in range(self.L):
                    for x2 in range(self.L):
                        for x3 in range(self.L):
                            site = (x0, x1, x2, x3)
                            old_U = self.get_link(mu, site).copy()
                            staple = self.staple_sum(site, mu)

                            old_action = np.real(np.trace(old_U @ staple))

                            new_U = random_su2_near_identity(eps) @ old_U
                            new_action = np.real(np.trace(new_U @ staple))

                            dS = beta * (new_action - old_action)
                            if dS > 0 or np.random.random() < np.exp(dS):
                                self.set_link(mu, site, new_U)
                                accepted += 1
                            total += 1

        return accepted / total if total > 0 else 0

    def thermalize(self, beta, n_sweeps=200):
        """Thermalize the lattice."""
        eps = min(0.5, 2.0 / max(beta, 0.1))
        for sweep in range(n_sweeps):
            self.metropolis_sweep(beta, eps)

    def avg_plaquette(self):
        """Average plaquette trace over the lattice."""
        total = 0
        count = 0
        for x0 in range(self.L):
            for x1 in range(self.L):
                for x2 in range(self.L):
                    for x3 in range(self.L):
                        for mu in range(N_DIM):
                            for nu in range(mu+1, N_DIM):
                                total += self.plaquette((x0,x1,x2,x3), mu, nu)
                                count += 1
        return total / count / 2  # normalize by Tr(I) = 2

    def extract_observables(self):
        """Extract two correlated observables: plaquettes in (0,1) and (0,2) planes.
        These share direction-0 links, producing structural correlation."""
        obs_01 = []  # plaquette in (0,1) plane
        obs_02 = []  # plaquette in (0,2) plane at same site
        for x0 in range(self.L):
            for x1 in range(self.L):
                for x2 in range(self.L):
                    for x3 in range(self.L):
                        obs_01.append(self.plaquette((x0,x1,x2,x3), 0, 1))
                        obs_02.append(self.plaquette((x0,x1,x2,x3), 0, 2))
        return np.array(obs_01), np.array(obs_02)


# ============================================================================
# CONDITIONAL DISTRIBUTION EXTRACTION
# ============================================================================

def build_conditional_distribution(obs1, obs2, H=3):
    """Build the HxH conditional distribution from two observable arrays.
    Uses equal-count tercile binning."""
    n = len(obs1)

    # Bin by equal-count terciles
    def bin_tercile(vals):
        sorted_idx = np.argsort(vals)
        bins = np.zeros(len(vals), dtype=int)
        for rank, i in enumerate(sorted_idx):
            bins[i] = min(rank * H // n, H - 1)
        return bins

    b1 = bin_tercile(obs1)
    b2 = bin_tercile(obs2)

    # Build conditional distribution C[i,j] = P(obs2 in bin j | obs1 in bin i)
    C = np.zeros((H, H))
    for k in range(n):
        C[b1[k], b2[k]] += 1

    # Normalize rows
    row_sums = C.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    C = C / row_sums

    return C


def conditional_to_evidence(C, H=3):
    """Convert a conditional distribution C to DS evidence mass functions.

    Each row of C defines one evidence mass function:
      e_i = C[dominant_row, i] * (1 - theta_e)
      theta_e = entropy-based ignorance

    We use the ROW with highest diagonal (most informative) as the
    evidence, following the paper's extraction method.

    Returns the evidence vector (e1, e2, e3, phi) with L1=1.
    Also returns the full C matrix and diagnostic quantities.
    """
    # Diagonal content
    delta = np.trace(C) / H

    # Use the most informative row (highest diagonal element)
    diag_vals = np.diag(C)
    best_row = np.argmax(diag_vals)
    row = C[best_row]

    # Evidence ignorance from Shannon entropy of the row
    entropy = -np.sum(row * np.log(row + 1e-15)) / np.log(H)
    phi = entropy  # ignorance = normalized entropy

    # Evidence singletons
    e = np.zeros(H + 1)
    for i in range(H):
        e[i] = row[i] * (1 - phi)
    e[H] = phi

    # L1 normalize
    e = e / np.sum(e)

    # Compute p_dom = max singleton / sum singletons
    singleton_sum = np.sum(e[:H])
    p_dom = np.max(e[:H]) / singleton_sum if singleton_sum > 0 else 1.0/H

    return e, phi, p_dom, delta, C


# ============================================================================
# DS FRAMEWORK: combination + floor + fixed point
# ============================================================================

H_DS = 3
FLOOR = 1.0 / H_DS**3

def ds_combine(m, e):
    """DS combination for real positive masses."""
    s, theta = m[:3], m[3]
    se, phi_e = e[:3], e[3]
    s_new = np.array([s[i]*se[i] + s[i]*phi_e + theta*se[i] for i in range(3)])
    theta_new = theta * phi_e
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)
    d = 1.0 - K
    if abs(d) < 1e-15:
        return m.copy(), K
    m_out = np.zeros(4)
    m_out[:3] = s_new / d
    m_out[3] = theta_new / d
    total = np.sum(m_out)
    if total > 1e-15:
        m_out /= total
    return m_out, K


def enforce_floor(m):
    """Born floor enforcement."""
    s, theta = m[:3], m[3]
    ssq = np.sum(s**2)
    total_sq = ssq + theta**2
    if theta**2 / total_sq >= FLOOR - 1e-15:
        return m.copy()
    ss = np.sum(s)
    if ss < 1e-15:
        return m.copy()
    r = ssq / ss**2
    disc = (2*r)**2 + 4*(26-r)*r
    t = (-2*r + np.sqrt(max(disc, 0))) / (2*(26-r))
    alpha = (1-t) / ss
    return np.array([s[0]*alpha, s[1]*alpha, s[2]*alpha, t])


def full_step(m, e):
    """One DS step + floor enforcement."""
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds), K


def find_fixed_point(e, n_iter=2000):
    """Find DS fixed point by iteration."""
    m = np.array([0.7, 0.05, 0.05, 0.2])
    for _ in range(n_iter):
        m_new, K = full_step(m, e)
        if np.max(np.abs(m_new - m)) < 1e-14:
            return m_new, True
        m = m_new
    return m, False


def measure_K(m, e):
    """Compute K(m, e)."""
    return sum(m[i]*e[j] for i in range(3) for j in range(3) if i != j)


def compute_eigenvalues(m, e, eps=1e-8):
    """3x3 projected Jacobian eigenvalues."""
    m0, _ = full_step(m, e)
    J = np.zeros((3, 3))
    indices = [0, 1, 3]
    for j in range(3):
        idx = indices[j]
        m_p = m.copy()
        m_p[idx] += eps
        m_p[3 if idx != 3 else 0] -= eps
        mp_out, _ = full_step(m_p, e)
        for i in range(3):
            J[i, j] = (mp_out[indices[i]] - m0[indices[i]]) / eps
    return np.sort(np.abs(np.linalg.eigvals(J)))[::-1]


# ============================================================================
# MAIN COMPUTATION
# ============================================================================

print("=" * 100)
print("BETA-DEPENDENCE OF THE EVIDENCE DISTRIBUTION ON THE K*=7/30 CURVE")
print("=" * 100)

# Parameters
lattice_sizes = [4, 6]
beta_values = [1.5, 1.8, 2.0, 2.2, 2.5, 2.8, 3.0, 3.5]
n_configs = 30  # independent configurations per (L, beta) pair
n_therm = 200   # thermalization sweeps per config

np.random.seed(42)  # reproducibility

print(f"\nLattice sizes: {lattice_sizes}")
print(f"Beta values: {beta_values}")
print(f"Configs per point: {n_configs}")
print(f"Thermalization: {n_therm} sweeps")
print(f"Lattice dimension: 4D")

results = []

for L in lattice_sizes:
    print(f"\n{'='*80}")
    print(f"LATTICE SIZE L = {L} (4D: {L}^4 = {L**4} sites)")
    print(f"{'='*80}")

    for beta in beta_values:
        print(f"\n  beta = {beta:.1f}, L = {L}: generating {n_configs} configs...", flush=True)

        # Generate 4D lattice data
        obs1_all, obs2_all = [], []
        for c in range(n_configs):
            lat = Lattice4D(L)
            lat.thermalize(beta, n_sweeps=n_therm)
            o1, o2 = lat.extract_observables()
            obs1_all.extend(o1)
            obs2_all.extend(o2)
            if (c+1) % 10 == 0:
                print(f"    config {c+1}/{n_configs} done", flush=True)

        obs1 = np.array(obs1_all)
        obs2 = np.array(obs2_all)
        n_total = len(obs1)
        print(f" {n_total} samples.", flush=True)

        # Build conditional distribution
        C = build_conditional_distribution(obs1, obs2, H=3)

        # Convert to evidence
        e_vec, phi, p_dom, delta, C_full = conditional_to_evidence(C, H=3)

        # Find DS fixed point with this evidence
        m_fp, converged = find_fixed_point(e_vec)
        K_val = measure_K(m_fp, e_vec)

        # Compute eigenvalues
        eigs = compute_eigenvalues(m_fp, e_vec)
        lam0 = eigs[0]
        lam1 = eigs[1] if len(eigs) > 1 else 0
        gap = -np.log(lam0) if 0 < lam0 < 1 else 0

        # Born check
        born_val = m_fp[3]**2 / np.sum(m_fp**2)

        # Mean plaquette
        mean_plaq = np.mean(obs1)

        # Also try: use AVERAGE over all rows of C as evidence
        e_avg = np.zeros(4)
        for row_idx in range(3):
            row = C_full[row_idx]
            ent = -np.sum(row * np.log(row + 1e-15)) / np.log(3)
            for i in range(3):
                e_avg[i] += row[i] * (1 - ent) / 3
            e_avg[3] += ent / 3
        e_avg /= np.sum(e_avg)
        phi_avg = e_avg[3]
        p_dom_avg = np.max(e_avg[:3]) / np.sum(e_avg[:3]) if np.sum(e_avg[:3]) > 0 else 1/3

        m_fp_avg, conv_avg = find_fixed_point(e_avg)
        K_avg = measure_K(m_fp_avg, e_avg)

        results.append({
            'L': L, 'beta': beta,
            'phi': phi, 'p_dom': p_dom, 'delta': delta,
            'K': K_val, 'lam0': lam0, 'lam1': lam1, 'gap': gap,
            'born': born_val, 'converged': converged,
            'mean_plaq': mean_plaq, 'n_samples': n_total,
            'e': e_vec, 'm': m_fp, 'C': C_full,
            'phi_avg': phi_avg, 'p_dom_avg': p_dom_avg, 'K_avg': K_avg,
        })

        print(f"    C = {np.array2string(C_full, precision=3, suppress_small=True)}")
        print(f"    delta = {delta:.4f}, phi = {phi:.4f}, p_dom = {p_dom:.4f}")
        print(f"    K* = {K_val:.6f} (7/30 = {7/30:.6f})")
        print(f"    lam0 = {lam0:.6f}, Delta = {gap:.4f}")
        print(f"    Born = {born_val:.6f}")
        print(f"    <plaq> = {mean_plaq:.4f}")

# ============================================================================
# ANALYSIS
# ============================================================================

print("\n\n" + "=" * 100)
print("ANALYSIS: HOW DO EVIDENCE PARAMETERS MOVE WITH BETA?")
print("=" * 100)

# Table: beta vs phi, p_dom, K*, Delta for each L
print(f"\n{'L':>3s} {'beta':>6s} {'phi':>8s} {'p_dom':>8s} {'K*':>10s} "
      f"{'K*-7/30':>10s} {'lam0':>8s} {'Delta':>8s} {'delta':>8s} {'<plaq>':>8s}")
print("-" * 95)

for r in sorted(results, key=lambda x: (x['L'], x['beta'])):
    marker = " <--" if abs(r['K'] - 7/30) < 0.01 else ""
    print(f"{r['L']:3d} {r['beta']:6.1f} {r['phi']:8.4f} {r['p_dom']:8.4f} "
          f"{r['K']:10.6f} {r['K']-7/30:10.4f} {r['lam0']:8.4f} {r['gap']:8.4f} "
          f"{r['delta']:8.4f} {r['mean_plaq']:8.4f}{marker}")

# Check: does phi increase with beta? (hypothesis: weaker coupling = more ignorant evidence)
print("\n\nTREND ANALYSIS: phi vs beta at each L")
print("-" * 60)
for L in lattice_sizes:
    subset = [r for r in results if r['L'] == L]
    subset.sort(key=lambda x: x['beta'])
    betas = [r['beta'] for r in subset]
    phis = [r['phi'] for r in subset]
    pdoms = [r['p_dom'] for r in subset]
    Ks = [r['K'] for r in subset]
    gaps = [r['gap'] for r in subset]

    if len(betas) > 2:
        # Linear fit phi(beta)
        c_phi = np.polyfit(betas, phis, 1)
        c_pdom = np.polyfit(betas, pdoms, 1)
        c_K = np.polyfit(betas, Ks, 1)
        c_gap = np.polyfit(betas, gaps, 1)

        print(f"\n  L = {L}:")
        print(f"    phi(beta)  ≈ {c_phi[0]:+.4f}*beta + {c_phi[1]:.4f}  "
              f"({'INCREASES' if c_phi[0] > 0 else 'DECREASES'} with beta)")
        print(f"    p_dom(beta) ≈ {c_pdom[0]:+.4f}*beta + {c_pdom[1]:.4f}  "
              f"({'INCREASES' if c_pdom[0] > 0 else 'DECREASES'} with beta)")
        print(f"    K*(beta)   ≈ {c_K[0]:+.4f}*beta + {c_K[1]:.4f}  "
              f"({'K* moves' if abs(c_K[0]) > 0.01 else 'K* stable'})")
        print(f"    Delta(beta) ≈ {c_gap[0]:+.4f}*beta + {c_gap[1]:.4f}  "
              f"({'INCREASES' if c_gap[0] > 0 else 'DECREASES'} with beta)")

# Check: is K* constant across beta?
print("\n\nK* UNIVERSALITY TEST: is K* independent of beta?")
print("-" * 60)
all_Ks = [r['K'] for r in results]
print(f"  K* range: [{min(all_Ks):.6f}, {max(all_Ks):.6f}]")
print(f"  K* mean:  {np.mean(all_Ks):.6f}")
print(f"  K* std:   {np.std(all_Ks):.6f}")
print(f"  7/30 =    {7/30:.6f}")
print(f"  Deviation from 7/30: {abs(np.mean(all_Ks) - 7/30):.4f} ({abs(np.mean(all_Ks) - 7/30)/(7/30)*100:.1f}%)")

# The key question: does the K*=7/30 curve parametrise the RG flow?
print("\n\nRG FLOW HYPOTHESIS: does (p_dom, phi) move along K*=7/30 curve with beta?")
print("-" * 60)
print(f"  If phi increases with beta and K* stays at 7/30:")
print(f"  → The moduli curve IS the RG flow")
print(f"  → 'Running coupling' = running connection sharpness")
print(f"  → Continuum limit = endpoint of curve (phi → max)")

# Separate trends for different L
for L in lattice_sizes:
    subset = sorted([r for r in results if r['L'] == L], key=lambda x: x['beta'])
    if len(subset) > 2:
        phi_trend = subset[-1]['phi'] - subset[0]['phi']
        K_trend = subset[-1]['K'] - subset[0]['K']
        print(f"\n  L={L}: phi changes by {phi_trend:+.4f} from beta={subset[0]['beta']} to {subset[-1]['beta']}")
        print(f"        K* changes by {K_trend:+.4f}")
        if phi_trend > 0.01 and abs(K_trend) < 0.05:
            print(f"        >>> CONSISTENT with RG flow along K*=7/30 curve")
        elif abs(phi_trend) < 0.01:
            print(f"        >>> phi is FLAT — no RG flow detected")
        else:
            print(f"        >>> INCONSISTENT — both phi and K* change")

print("\n\nCOMPUTATION COMPLETE.")
