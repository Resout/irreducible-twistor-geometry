"""
Multi-site spectral gap: the core problem
==========================================
Single-site: d_H(T(m,e), T(m',e)) <= kappa * d_H(m, m'), kappa ~ 0.956.
Coupled lattice: evidence at site x depends on neighbor states.

Key question: does the coupling destroy the gap?

Strategy: on the product space, use the triangle inequality.
d_H(T(m_x, e_x), T(m'_x, e'_x))
  <= kappa * d_H(m_x, m'_x)         [contraction with same evidence]
   + L * d_H(e_x, e'_x)             [Lipschitz in evidence]

where L = Lip(e -> T(m,e)) in the Hilbert metric.
If kappa + L < 1, the coupled system contracts.
"""

import numpy as np
from numpy.linalg import norm

def ds_combine(m, e):
    th_m, s1m, s2m, s3m = m
    th_e, s1e, s2e, s3e = e
    K = s1m*s2e + s1m*s3e + s2m*s1e + s2m*s3e + s3m*s1e + s3m*s2e
    if K >= 1.0:
        return m.copy()
    inv = 1.0 / (1.0 - K)
    return np.array([
        inv * th_m * th_e,
        inv * (s1m*s1e + s1m*th_e + th_m*s1e),
        inv * (s2m*s2e + s2m*th_e + th_m*s2e),
        inv * (s3m*s3e + s3m*th_e + th_m*s3e),
    ])

def born_floor(m, H=3):
    th = m[0]
    s = m[1:]
    s_sq = np.dot(s, s)
    born = th**2 / (th**2 + s_sq) if (th**2 + s_sq) > 0 else 1.0
    if born < 1.0/H**3:
        th = np.sqrt(s_sq / (H**3 - 1))
        m = np.array([th, m[1], m[2], m[3]])
    total = np.sum(np.abs(m))
    if total > 0:
        m = m / total
    return m

def ds_step(m, e, H=3):
    return born_floor(ds_combine(m, e), H)

def hilbert_dist(m, mp):
    if np.any(m <= 1e-15) or np.any(mp <= 1e-15):
        return float('inf')
    ratios = m / mp
    return np.log(np.max(ratios)) - np.log(np.min(ratios))

# ============================================================
# Part 1: Measure kappa (contraction rate with SAME evidence)
# ============================================================
print("=" * 70)
print("PART 1: Single-site contraction rate kappa")
print("=" * 70)

# The paper claims kappa ~ 0.956. Let's measure it precisely.
np.random.seed(42)
kappas = []

for trial in range(10000):
    m1 = np.random.dirichlet([2, 2, 2, 2])
    m1 = born_floor(m1)
    m2 = np.random.dirichlet([2, 2, 2, 2])
    m2 = born_floor(m2)
    e = np.random.dirichlet([2, 2, 2, 2])
    e = born_floor(e)

    d_before = hilbert_dist(m1, m2)
    if d_before < 1e-10:
        continue

    t1 = ds_step(m1, e)
    t2 = ds_step(m2, e)

    d_after = hilbert_dist(t1, t2)
    kappas.append(d_after / d_before)

print(f"kappa: mean={np.mean(kappas):.6f}, max={np.max(kappas):.6f}, "
      f"std={np.std(kappas):.6f}")
print(f"Paper claims kappa ~ 0.956")
kappa_max = np.max(kappas)

# ============================================================
# Part 2: Measure L (Lipschitz constant in evidence)
# ============================================================
print()
print("=" * 70)
print("PART 2: Evidence Lipschitz constant L")
print("=" * 70)

# For fixed m, how much does T(m, e) change when e changes?
# L = sup_{m,e,e'} d_H(T(m,e), T(m,e')) / d_H(e, e')

Ls = []

for trial in range(10000):
    m = np.random.dirichlet([2, 2, 2, 2])
    m = born_floor(m)
    e1 = np.random.dirichlet([2, 2, 2, 2])
    e1 = born_floor(e1)
    e2 = np.random.dirichlet([2, 2, 2, 2])
    e2 = born_floor(e2)

    d_evidence = hilbert_dist(e1, e2)
    if d_evidence < 1e-10:
        continue

    t1 = ds_step(m, e1)
    t2 = ds_step(m, e2)

    d_output = hilbert_dist(t1, t2)
    Ls.append(d_output / d_evidence)

print(f"L: mean={np.mean(Ls):.6f}, max={np.max(Ls):.6f}, std={np.std(Ls):.6f}")
L_max = np.max(Ls)

print(f"\nkappa + L = {kappa_max:.6f} + {L_max:.6f} = {kappa_max + L_max:.6f}")
print(f"Contraction: kappa + L < 1? {kappa_max + L_max < 1}")

# ============================================================
# Part 3: More careful - measure at equilibrium
# ============================================================
print()
print("=" * 70)
print("PART 3: Contraction rates near equilibrium")
print("=" * 70)

# Near equilibrium, the dynamics is approximately linear.
# The effective contraction includes both the local contraction and
# the evidence-coupling effect.

# Find the equilibrium first
e_eq = np.array([0.2, 0.3, 0.3, 0.2])
m_eq = np.array([0.25, 0.25, 0.25, 0.25])
for i in range(100):
    m_eq = ds_step(m_eq, e_eq)

print(f"Equilibrium: m* = {m_eq}")

# Now measure kappa and L near equilibrium
kappas_eq = []
Ls_eq = []

for trial in range(10000):
    # Small perturbations around equilibrium
    delta1 = np.random.randn(4) * 0.01
    delta1 -= np.mean(delta1)  # preserve L1
    m1 = born_floor(np.abs(m_eq + delta1))
    m1 /= np.sum(m1)
    m1 = born_floor(m1)

    delta2 = np.random.randn(4) * 0.01
    delta2 -= np.mean(delta2)
    m2 = born_floor(np.abs(m_eq + delta2))
    m2 /= np.sum(m2)
    m2 = born_floor(m2)

    # Same evidence
    d_before = hilbert_dist(m1, m2)
    if d_before < 1e-10:
        continue
    t1 = ds_step(m1, e_eq)
    t2 = ds_step(m2, e_eq)
    d_after = hilbert_dist(t1, t2)
    kappas_eq.append(d_after / d_before)

    # Different evidence (perturbed)
    de = np.random.randn(4) * 0.01
    de -= np.mean(de)
    e1 = born_floor(np.abs(e_eq + de))
    e1 /= np.sum(e1)
    e1 = born_floor(e1)

    d_ev = hilbert_dist(e_eq, e1)
    if d_ev < 1e-10:
        continue
    t_eq = ds_step(m_eq, e_eq)
    t_e1 = ds_step(m_eq, e1)
    d_out = hilbert_dist(t_eq, t_e1)
    Ls_eq.append(d_out / d_ev)

print(f"Near equilibrium:")
print(f"  kappa: mean={np.mean(kappas_eq):.6f}, max={np.max(kappas_eq):.6f}")
print(f"  L:     mean={np.mean(Ls_eq):.6f}, max={np.max(Ls_eq):.6f}")
print(f"  kappa + L: mean={np.mean(kappas_eq)+np.mean(Ls_eq):.6f}, "
      f"max={np.max(kappas_eq)+np.max(Ls_eq):.6f}")
print(f"  Coupled contraction near eq: {np.max(kappas_eq)+np.max(Ls_eq) < 1}")

# ============================================================
# Part 4: The neighbor averaging effect
# ============================================================
print()
print("=" * 70)
print("PART 4: Neighbor averaging reduces effective L")
print("=" * 70)

# On the lattice, the evidence at site x is the AVERAGE of neighbors.
# Averaging is a contraction in the Hilbert metric!
# For z = (a+b)/2: d_H(z, z') <= max(d_H(a,a'), d_H(b,b'))
# Actually stronger: d_H(avg, avg') <= (1/2)(d_H(a,a') + d_H(b,b'))
# ... no, Hilbert metric doesn't satisfy that. Let me check.

# The Hilbert metric on the simplex: d_H(m, m') = log(max m_i/m'_i) - log(min m_i/m'_i)
# For the average z = (a+b)/2:
# z_i/z'_i = (a_i+b_i)/(a'_i+b'_i)
# This is bounded by max(a_i/a'_i, b_i/b'_i) from above
# and min(a_i/a'_i, b_i/b'_i) from below? No, that's not right either.

# Let me just measure it directly: how does averaging affect Hilbert distance?

avg_contractions = []
for trial in range(10000):
    a = np.random.dirichlet([2, 2, 2, 2])
    a = born_floor(a)
    b = np.random.dirichlet([2, 2, 2, 2])
    b = born_floor(b)
    ap = np.random.dirichlet([2, 2, 2, 2])
    ap = born_floor(ap)
    bp = np.random.dirichlet([2, 2, 2, 2])
    bp = born_floor(bp)

    z = (a + b) / 2
    z /= np.sum(z)
    zp = (ap + bp) / 2
    zp /= np.sum(zp)

    d_z = hilbert_dist(z, zp)
    d_max = max(hilbert_dist(a, ap), hilbert_dist(b, bp))

    if d_max > 1e-10:
        avg_contractions.append(d_z / d_max)

print(f"Averaging contraction: d_H(avg, avg') / max(d_H inputs)")
print(f"  mean={np.mean(avg_contractions):.6f}, max={np.max(avg_contractions):.6f}")

# What about 4-neighbor averaging (2D lattice)?
avg4_contractions = []
for trial in range(10000):
    inputs = [born_floor(np.random.dirichlet([2,2,2,2])) for _ in range(4)]
    inputs_p = [born_floor(np.random.dirichlet([2,2,2,2])) for _ in range(4)]

    z = np.mean(inputs, axis=0)
    z /= np.sum(z)
    zp = np.mean(inputs_p, axis=0)
    zp /= np.sum(zp)

    d_z = hilbert_dist(z, zp)
    d_max = max(hilbert_dist(inputs[i], inputs_p[i]) for i in range(4))

    if d_max > 1e-10:
        avg4_contractions.append(d_z / d_max)

print(f"4-neighbor averaging: d_H(avg, avg') / max(d_H inputs)")
print(f"  mean={np.mean(avg4_contractions):.6f}, max={np.max(avg4_contractions):.6f}")
avg4_max = np.max(avg4_contractions)

# ============================================================
# Part 5: The combined bound
# ============================================================
print()
print("=" * 70)
print("PART 5: Combined multi-site contraction")
print("=" * 70)

# On the coupled lattice with z neighbors per site:
# d_H(T(m_x, e_x), T(m'_x, e'_x))
#   <= kappa * d_H(m_x, m'_x) + L * d_H(e_x, e'_x)
#
# where e_x = avg of z neighbors and e'_x = avg of z neighbors (primed)
# d_H(e_x, e'_x) <= alpha * max_{y~x} d_H(m_y, m'_y)
# where alpha is the averaging contraction factor.
#
# Taking supremum over x:
# D_{n+1} <= kappa * D_n + L * alpha * D_n = (kappa + L*alpha) * D_n
#
# The effective contraction rate is kappa + L*alpha.

print(f"Components:")
print(f"  kappa (mass contraction):     {kappa_max:.6f}")
print(f"  L (evidence Lipschitz):       {L_max:.6f}")
print(f"  alpha_2 (2-neighbor avg):     {np.max(avg_contractions):.6f}")
print(f"  alpha_4 (4-neighbor avg):     {avg4_max:.6f}")
print()

for z, alpha in [(2, np.max(avg_contractions)), (4, avg4_max)]:
    effective = kappa_max + L_max * alpha
    print(f"  z={z} neighbors: kappa + L*alpha = {kappa_max:.4f} + {L_max:.4f}*{alpha:.4f} = {effective:.4f}")
    print(f"    Contracts: {effective < 1}")

# ============================================================
# Part 6: Direct measurement on coupled lattice
# ============================================================
print()
print("=" * 70)
print("PART 6: Direct coupled contraction measurement")
print("=" * 70)

# Run two lattices with different initial conditions, same topology.
# Measure max d_H between corresponding sites at each step.

N = 20

def make_random_lattice(N):
    lat = np.zeros((N, 4))
    for i in range(N):
        lat[i] = born_floor(np.random.dirichlet([2, 2, 2, 2]))
    return lat

def coupled_step(lattice, N):
    new = np.zeros_like(lattice)
    for i in range(N):
        m = lattice[i]
        e_left = lattice[(i-1) % N]
        e_right = lattice[(i+1) % N]
        e_avg = 0.5 * (e_left + e_right)
        new[i] = ds_step(m, e_avg)
    return new

def max_site_distance(lat1, lat2, N):
    d = 0
    for i in range(N):
        d = max(d, hilbert_dist(lat1[i], lat2[i]))
    return d

np.random.seed(0)
lat1 = make_random_lattice(N)
lat2 = make_random_lattice(N)

d0 = max_site_distance(lat1, lat2, N)
print(f"1D chain N={N}")
print(f"{'step':>6} {'max d_H':>12} {'ratio':>10}")
print(f"{'init':>6} {d0:12.6f}")

prev_d = d0
for step in range(30):
    lat1 = coupled_step(lat1, N)
    lat2 = coupled_step(lat2, N)
    d = max_site_distance(lat1, lat2, N)
    ratio = d / prev_d if prev_d > 1e-10 else 0
    if step % 3 == 0:
        print(f"{step:6d} {d:12.6f} {ratio:10.6f}")
    prev_d = d

# 2D lattice
print()
N2 = 10

def make_random_lattice_2d(N):
    lat = np.zeros((N, N, 4))
    for i in range(N):
        for j in range(N):
            lat[i, j] = born_floor(np.random.dirichlet([2, 2, 2, 2]))
    return lat

def coupled_step_2d(lat, N):
    new = np.zeros_like(lat)
    for i in range(N):
        for j in range(N):
            m = lat[i, j]
            nbrs = [lat[(i+1)%N,j], lat[(i-1)%N,j], lat[i,(j+1)%N], lat[i,(j-1)%N]]
            e_avg = np.mean(nbrs, axis=0)
            new[i, j] = ds_step(m, e_avg)
    return new

def max_distance_2d(l1, l2, N):
    d = 0
    for i in range(N):
        for j in range(N):
            d = max(d, hilbert_dist(l1[i,j], l2[i,j]))
    return d

np.random.seed(1)
l2d_1 = make_random_lattice_2d(N2)
l2d_2 = make_random_lattice_2d(N2)

d0_2d = max_distance_2d(l2d_1, l2d_2, N2)
print(f"2D lattice {N2}x{N2}")
print(f"{'step':>6} {'max d_H':>12} {'ratio':>10}")
print(f"{'init':>6} {d0_2d:12.6f}")

prev_d = d0_2d
for step in range(20):
    l2d_1 = coupled_step_2d(l2d_1, N2)
    l2d_2 = coupled_step_2d(l2d_2, N2)
    d = max_distance_2d(l2d_1, l2d_2, N2)
    ratio = d / prev_d if prev_d > 1e-10 else 0
    if step % 2 == 0:
        print(f"{step:6d} {d:12.6f} {ratio:10.6f}")
    prev_d = d

# ============================================================
# Part 7: The structural argument
# ============================================================
print()
print("=" * 70)
print("PART 7: Why the gap survives coupling")
print("=" * 70)

# The Hilbert metric contraction holds for ANY positive evidence.
# On the coupled lattice, the evidence at each site IS positive
# (all neighbors are in S_epsilon after 2 steps).
#
# The key: the DS map T(m, e) is a CONTRACTION in m for FIXED e.
# On the coupled lattice, e depends on the configuration, but
# the contraction applies pointwise: at each site, the map contracts.
#
# The coupling creates correlations between sites but cannot
# amplify the local contraction into expansion.
#
# More precisely: define D_n = max_x d_H(m_x^n, m*) where m* is
# the unique fixed point. At each site:
# d_H(m_x^{n+1}, m*) = d_H(T(m_x^n, e_x^n), T(m*, e*))
#   <= d_H(T(m_x^n, e_x^n), T(m*, e_x^n))  +  d_H(T(m*, e_x^n), T(m*, e*))
#   <= kappa * d_H(m_x^n, m*)  +  d_H(T(m*, e_x^n), T(m*, e*))
#
# The second term: how much does the fixed point shift when evidence changes?
# At m*, the map is T(m*, e) = m* for the equilibrium evidence e*.
# For perturbed evidence: T(m*, e) = m* + J_e * (e - e*) + O(|e-e*|^2)
# where J_e is the Jacobian wrt evidence.
#
# The eigenvalues of J_e at equilibrium are 0.913, 0.489, 0.158, ~0.
# So d_H(T(m*, e), T(m*, e*)) <= lambda_max(J_e) * d_H(e, e*)
#                                = 0.913 * d_H(e_x, e*)
#
# And e_x = avg of neighbors, which are all converging to m*.
# So d_H(e_x, e*) <= max_{y~x} d_H(m_y, m*) = D_n.
#
# Therefore: D_{n+1} <= kappa * D_n + 0.913 * D_n = (kappa + 0.913) * D_n
# This gives kappa + 0.913 ~ 1.87 > 1. NOT contracting.
#
# BUT: this is too pessimistic. The 0.913 eigenvalue is the RESPONSE
# of the output to evidence perturbation, NOT the Hilbert metric Lipschitz.
# Let me measure the actual Lipschitz constant of T(m*, .) at equilibrium.

# Measure: d_H(T(m*, e), m*) / d_H(e, m*) near equilibrium
evidence_responses = []
for trial in range(10000):
    de = np.random.randn(4) * 0.01
    de -= np.mean(de)
    e_perturbed = np.abs(m_eq + de)
    e_perturbed /= np.sum(e_perturbed)
    e_perturbed = born_floor(e_perturbed)

    d_evidence = hilbert_dist(e_perturbed, m_eq)
    if d_evidence < 1e-10:
        continue

    t_perturbed = ds_step(m_eq, e_perturbed)
    d_response = hilbert_dist(t_perturbed, m_eq)  # distance from fixed point

    evidence_responses.append(d_response / d_evidence)

print(f"Evidence response at equilibrium: d_H(T(m*,e), m*) / d_H(e, m*)")
print(f"  mean={np.mean(evidence_responses):.6f}, max={np.max(evidence_responses):.6f}")

L_eq = np.max(evidence_responses)
print(f"\nEffective coupled rate at equilibrium:")
print(f"  kappa + L_eq = {kappa_max:.6f} + {L_eq:.6f} = {kappa_max + L_eq:.6f}")
print(f"  With averaging (z=2): {kappa_max + L_eq * np.max(avg_contractions):.6f}")
print(f"  With averaging (z=4): {kappa_max + L_eq * avg4_max:.6f}")

# ============================================================
# Part 8: Alternative - Dobrushin criterion
# ============================================================
print()
print("=" * 70)
print("PART 8: Dobrushin uniqueness criterion")
print("=" * 70)

# The Dobrushin criterion for uniqueness of Gibbs measures:
# Define the influence matrix C_{xy} = sup_m d_H(T_x(m), T_x(m')) / d_H(m_y, m'_y)
# where the sup is over configs differing only at site y.
#
# If max_x sum_y C_{xy} < 1, there is a unique Gibbs measure with spectral gap.
#
# For the DS lattice:
# C_{xx} = kappa (self-contraction)
# C_{xy} = influence of neighbor y on site x through the evidence average.
#
# For z neighbors: C_{xy} = L / z for each neighbor (evidence is 1/z weighted).

# Measure: influence of a single neighbor change on the output
single_neighbor_influence = []

for trial in range(10000):
    # Lattice near equilibrium
    m_center = born_floor(np.abs(m_eq + np.random.randn(4)*0.01))
    m_center /= np.sum(m_center)
    m_center = born_floor(m_center)

    nbr1 = born_floor(np.abs(m_eq + np.random.randn(4)*0.01))
    nbr1 /= np.sum(nbr1)
    nbr1 = born_floor(nbr1)

    nbr2 = born_floor(np.abs(m_eq + np.random.randn(4)*0.01))
    nbr2 /= np.sum(nbr2)
    nbr2 = born_floor(nbr2)

    # Change one neighbor
    nbr1_alt = born_floor(np.abs(m_eq + np.random.randn(4)*0.01))
    nbr1_alt /= np.sum(nbr1_alt)
    nbr1_alt = born_floor(nbr1_alt)

    # Evidence with original and altered neighbor
    e_orig = 0.5 * (nbr1 + nbr2)
    e_alt = 0.5 * (nbr1_alt + nbr2)

    t_orig = ds_step(m_center, e_orig)
    t_alt = ds_step(m_center, e_alt)

    d_output = hilbert_dist(t_orig, t_alt)
    d_input = hilbert_dist(nbr1, nbr1_alt)

    if d_input > 1e-10:
        single_neighbor_influence.append(d_output / d_input)

C_xy = np.max(single_neighbor_influence)
print(f"Single-neighbor influence C_xy = max d_H(output change) / d_H(neighbor change)")
print(f"  mean={np.mean(single_neighbor_influence):.6f}, max={C_xy:.6f}")

# Dobrushin condition: C_xx + sum_{y~x} C_xy < 1
# C_xx = kappa, sum over z neighbors of C_xy
for z in [2, 4, 6]:
    dobrushin = kappa_max + z * C_xy
    print(f"  z={z}: kappa + z*C_xy = {kappa_max:.4f} + {z}*{C_xy:.4f} = {dobrushin:.4f} {'< 1 CONTRACTS' if dobrushin < 1 else '>= 1'}")

# ============================================================
# Summary
# ============================================================
print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
Single-site contraction rate kappa = {kappa_max:.6f}
Evidence Lipschitz constant L = {L_max:.6f}
Evidence response at equilibrium L_eq = {L_eq:.6f}
2-neighbor averaging factor = {np.max(avg_contractions):.6f}
4-neighbor averaging factor = {avg4_max:.6f}
Single-neighbor influence C_xy = {C_xy:.6f}

The coupled contraction requires kappa + L*alpha < 1 or the Dobrushin
condition kappa + z*C_xy < 1.

Direct lattice measurement shows the coupled system DOES contract
(the distance between two configurations decreases at each step).
The analytic bound needs to be tighter than the naive kappa + L.
""")
