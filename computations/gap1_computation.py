"""
Gap 1 computation: Quantitative F+ = F- verification.

4D lattice with DS dynamics as the update rule.
Hodge decomposition of plaquette curvature into self-dual and anti-self-dual.
The claim (Prop 19.10): at equilibrium, M is real => |F+| = |F-|.
"""
import numpy as np
from numpy.linalg import svd, det, inv, norm

# Pauli matrices
I2 = np.eye(2, dtype=complex)
s1_mat = np.array([[0, 1], [1, 0]], dtype=complex)
s2_mat = np.array([[0, -1j], [1j, 0]], dtype=complex)
s3_mat = np.array([[1, 0], [0, -1]], dtype=complex)

H = 3
BORN_FLOOR = 1.0 / H**3  # 1/27


def mass_to_M(mass):
    """su(2) embedding: M = (1/sqrt2)(theta*I + s1*sigma1 + s2*sigma2 + s3*sigma3)"""
    s1, s2, s3, th = mass[0], mass[1], mass[2], mass[3]
    return (th * I2 + s1 * s1_mat + s2 * s2_mat + s3 * s3_mat) / np.sqrt(2)


def enforce_born_floor(mass, floor=BORN_FLOOR):
    m2 = np.abs(mass)**2
    S = np.sum(m2)
    if S < 1e-15:
        return mass
    if m2[3] / S >= floor:
        return mass
    lo, hi = 0.0, 1.0
    for _ in range(30):
        mid = (lo + hi) / 2
        trial = mass.copy()
        trial[:3] *= mid
        trial[3] = 1.0 - np.sum(trial[:3])
        t2 = np.abs(trial)**2
        St = np.sum(t2)
        bt = t2[3] / St if St > 0 else 1.0
        if bt < floor:
            hi = mid
        else:
            lo = mid
    result = mass.copy()
    result[:3] *= lo
    result[3] = 1.0 - np.sum(result[:3])
    return result


def ds_combine(m1, m2, floor=BORN_FLOOR):
    s1, th1 = m1[:H], m1[H]
    s2, th2 = m2[:H], m2[H]
    s_new = s1 * s2 + s1 * th2 + th1 * s2
    th_new = th1 * th2
    K = np.sum(s1) * np.sum(s2) - np.sum(s1 * s2)
    n = 1.0 - K
    if abs(n) < 1e-15:
        n = 1e-15
    combined = np.append(s_new, th_new) / n
    combined = enforce_born_floor(combined, floor)
    return combined, K


# Fixed-point mass function
m_star = np.array([0.5959, 0.1405, 0.1405, 0.1232])
S_star = np.sum(m_star**2)
born_star = m_star**2 / S_star
M_star = mass_to_M(m_star)

print("Fixed-point mass:", m_star)
print("L1 sum:", np.sum(m_star))
print("Born probs:", born_star)
print("Born(theta):", born_star[3], " target:", BORN_FLOOR)
print("det(M*):", det(M_star))
print()

# Kink computation (equation 25)
theta = m_star[3]
c = np.sqrt(np.sum(m_star[:3]**2) / 26)
abs_theta = abs(theta)

dbar_Phi = np.zeros((4, 4))
for i in range(3):
    dbar_Phi[3, i] = theta * m_star[i] / (52 * abs_theta * c)
dbar_Phi[3, 3] = -c * theta**2 / (2 * abs_theta**3)

# Kink as 2x2 matrix (theta-bar direction)
kink_dir = dbar_Phi[:, 3]
delta_M = (kink_dir[3] * I2 + kink_dir[0] * s1_mat +
           kink_dir[1] * s2_mat + kink_dir[2] * s3_mat) / np.sqrt(2)
U_k, S_k, Vh_k = svd(delta_M)
print("Kink singular values:", S_k)
print("Rank:", np.sum(S_k > 1e-10))
S_ent_dist = S_k**2 / np.sum(S_k**2)
S_ent = -np.sum(S_ent_dist * np.log2(S_ent_dist + 1e-30))
print("Entanglement entropy:", S_ent, "bits")
print()

# 4D lattice computation
L = 4
np.random.seed(42)

# Initialize link masses as perturbations of fixed point
link_masses = np.zeros((L, L, L, L, 4, 4))
for idx in np.ndindex(L, L, L, L, 4):
    pert = np.random.randn(4) * 0.05
    m_link = np.abs(m_star + pert)
    m_link = m_link / np.sum(m_link)
    link_masses[idx] = m_link


def get_M_at(x, d):
    return mass_to_M(link_masses[x[0] % L, x[1] % L, x[2] % L, x[3] % L, d])


def thermalize(n_steps):
    global link_masses
    for _ in range(n_steps):
        new = link_masses.copy()
        for idx in np.ndindex(L, L, L, L, 4):
            x0, x1, x2, x3, d = idx
            nx = [x0, x1, x2, x3]
            nx[d] = (nx[d] + 1) % L
            neighbor = link_masses[nx[0], nx[1], nx[2], nx[3], d]
            combined, K = ds_combine(link_masses[idx], neighbor)
            new[idx] = combined
        link_masses = new


def compute_chirality_ratio():
    """Compute |F+|/|F-| across all sites using proper 4D Hodge decomposition."""
    ratios = []
    for site in np.ndindex(L, L, L, L):
        Fc = {}
        for mu in range(4):
            for nu in range(mu + 1, 4):
                x = list(site)
                M1 = get_M_at(x, mu)
                x_mu = list(site)
                x_mu[mu] = (x_mu[mu] + 1) % L
                M2 = get_M_at(x_mu, nu)
                x_nu = list(site)
                x_nu[nu] = (x_nu[nu] + 1) % L
                M3 = get_M_at(x_nu, mu)
                M4 = get_M_at(x, nu)
                U_p = M1 @ M2 @ inv(M3) @ inv(M4)
                Fc[(mu, nu)] = U_p - I2

        F01 = Fc[(0, 1)]
        F02 = Fc[(0, 2)]
        F03 = Fc[(0, 3)]
        F12 = Fc[(1, 2)]
        F13 = Fc[(1, 3)]
        F23 = Fc[(2, 3)]

        # Hodge star in 4D Euclidean:
        # *F01 = F23, *F02 = -F13, *F03 = F12
        # Self-dual: F+ = (F + *F)/2
        # Anti-self-dual: F- = (F - *F)/2
        Fp_sq = (norm(F01 + F23)**2 + norm(F02 - F13)**2 + norm(F03 + F12)**2) / 4
        Fm_sq = (norm(F01 - F23)**2 + norm(F02 + F13)**2 + norm(F03 - F12)**2) / 4

        if Fm_sq > 1e-20 and Fp_sq > 1e-20:
            ratios.append(np.sqrt(Fp_sq / Fm_sq))

    return np.array(ratios)


# Run thermalization and measure
print("Thermalizing and measuring |F+|/|F-| with proper Hodge decomposition...")
print()

for total_steps in [5, 10, 20, 40, 80]:
    if total_steps == 5:
        thermalize(5)
    else:
        prev = [5, 10, 20, 40, 80]
        idx = prev.index(total_steps)
        thermalize(total_steps - prev[idx - 1])

    ratios = compute_chirality_ratio()
    if len(ratios) > 0:
        print(f"After {total_steps:3d} DS steps: "
              f"|F+|/|F-| = {ratios.mean():.6f} +/- {ratios.std():.6f}  "
              f"(n={len(ratios)}, min={ratios.min():.4f}, max={ratios.max():.4f})")
    else:
        print(f"After {total_steps:3d} DS steps: no valid plaquettes")

# Check if masses are converging to real
print()
sample = link_masses[0, 0, 0, 0, 0]
print("Sample mass after thermalization:", sample)
print("Imaginary content:", np.max(np.abs(np.imag(sample))) if np.iscomplexobj(sample) else "all real")

# Extract K* from the last round
Ks = []
for idx in np.ndindex(L, L, L, L, 4):
    x0, x1, x2, x3, d = idx
    nx = [x0, x1, x2, x3]
    nx[d] = (nx[d] + 1) % L
    _, K = ds_combine(link_masses[idx], link_masses[nx[0], nx[1], nx[2], nx[3], d])
    Ks.append(abs(K))
Ks = np.array(Ks)
print(f"\nConflict K: mean={Ks.mean():.6f}, std={Ks.std():.6f}  (predicted K*=7/30={7/30:.6f})")
