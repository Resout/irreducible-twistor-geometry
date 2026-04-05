"""
Wilson loop area law — version 2.

The area law is about the EXPECTATION <W(C)> over configurations,
not the value at a single (uniform) configuration.

At uniform equilibrium: W = Tr(I)/2 = 1 for all loops.
The area law comes from fluctuations around equilibrium.

Strategy: sample mass function configurations from B,
compute Wilson loops, average.
"""

import numpy as np
from numpy.linalg import norm, det

def ds_combine(m, e):
    th_m, s1m, s2m, s3m = m
    th_e, s1e, s2e, s3e = e
    K = s1m*s2e + s1m*s3e + s2m*s1e + s2m*s3e + s3m*s1e + s3m*s2e
    if abs(K) >= 1.0:
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
    total_sq = th**2 + s_sq
    if total_sq < 1e-30:
        return np.array([1.0, 0.0, 0.0, 0.0])
    born = th**2 / total_sq
    if born < 1.0/H**3:
        th_new = np.sqrt(s_sq / (H**3 - 1))
        m = np.array([th_new, m[1], m[2], m[3]])
    total = np.sum(np.abs(m))
    if total > 0:
        m = m / total
    return m

def pauli_matrix(m):
    th, s1, s2, s3 = m
    return np.array([
        [th + s3, s1 - 1j*s2],
        [s1 + 1j*s2, th - s3]
    ], dtype=complex) / np.sqrt(2)

def wilson_loop_2d(lat, x, y, R, T):
    """Wilson loop of size R×T at position (x,y) on a 2D lattice of mass functions."""
    N = lat.shape[0]
    W = np.eye(2, dtype=complex)

    # Build transition functions around the contour
    # Right leg
    for dx in range(R):
        M1 = pauli_matrix(lat[(x+dx)%N, y])
        M2 = pauli_matrix(lat[(x+dx+1)%N, y])
        d = det(M1)
        if abs(d) > 1e-10:
            W = W @ np.linalg.solve(M1, M2)

    # Up leg
    for dy in range(T):
        M1 = pauli_matrix(lat[(x+R)%N, (y+dy)%N])
        M2 = pauli_matrix(lat[(x+R)%N, (y+dy+1)%N])
        d = det(M1)
        if abs(d) > 1e-10:
            W = W @ np.linalg.solve(M1, M2)

    # Left leg
    for dx in range(R):
        M1 = pauli_matrix(lat[(x+R-dx)%N, (y+T)%N])
        M2 = pauli_matrix(lat[(x+R-dx-1)%N, (y+T)%N])
        d = det(M1)
        if abs(d) > 1e-10:
            W = W @ np.linalg.solve(M1, M2)

    # Down leg
    for dy in range(T):
        M1 = pauli_matrix(lat[x, (y+T-dy)%N])
        M2 = pauli_matrix(lat[x, (y+T-dy-1)%N])
        d = det(M1)
        if abs(d) > 1e-10:
            W = W @ np.linalg.solve(M1, M2)

    return np.real(np.trace(W)) / 2

# ============================================================
# Sample configurations from B and measure Wilson loops
# ============================================================
print("=" * 70)
print("WILSON LOOP AREA LAW — CONFIGURATION AVERAGE")
print("=" * 70)

N = 12  # lattice size
n_configs = 200  # number of configurations to average over

# For each configuration: random mass functions on B (not thermalized)
# These represent fluctuations around equilibrium.

results = {}
for R in range(1, 6):
    for T in range(R, 6):
        W_values = []
        for config in range(n_configs):
            np.random.seed(config * 1000 + R * 100 + T)
            # Sample a random configuration on B
            lat = np.zeros((N, N, 4))
            for i in range(N):
                for j in range(N):
                    lat[i,j] = born_floor(np.random.dirichlet([2,2,2,2]))

            # Measure Wilson loop at a random position
            x, y = np.random.randint(0, N, 2)
            W = wilson_loop_2d(lat, x, y, R, T)
            W_values.append(W)

        W_mean = np.mean(W_values)
        W_abs_mean = np.mean(np.abs(W_values))
        area = R * T
        perim = 2 * (R + T)
        results[(R, T)] = {
            'area': area, 'perim': perim,
            'W_mean': W_mean, 'W_abs': W_abs_mean,
        }

print(f"\n{'R×T':>6} {'Area':>5} {'Perim':>5} {'<W>':>12} {'<|W|>':>12} "
      f"{'-ln<|W|>':>10} {'-ln/A':>8}")

for (R, T), d in sorted(results.items()):
    if d['W_abs'] > 1e-10:
        logW = -np.log(d['W_abs'])
        logW_per_area = logW / d['area']
    else:
        logW = float('inf')
        logW_per_area = float('inf')

    print(f"{R}×{T:>2} {d['area']:5d} {d['perim']:5d} {d['W_mean']:12.6f} "
          f"{d['W_abs']:12.6f} {logW:10.4f} {logW_per_area:8.4f}")

# ============================================================
# Fit area law vs perimeter law
# ============================================================
print("\n" + "=" * 70)
print("AREA LAW VS PERIMETER LAW FIT")
print("=" * 70)

areas = []
perimeters = []
log_Ws = []

for (R, T), d in sorted(results.items()):
    if d['W_abs'] > 1e-10:
        areas.append(d['area'])
        perimeters.append(d['perim'])
        log_Ws.append(-np.log(d['W_abs']))

areas = np.array(areas)
perimeters = np.array(perimeters)
log_Ws = np.array(log_Ws)

# Area fit: -ln<|W|> = σ·A + c
A_area = np.column_stack([areas, np.ones(len(areas))])
coeffs_area = np.linalg.lstsq(A_area, log_Ws, rcond=None)[0]
pred_area = A_area @ coeffs_area
SS_res_area = np.sum((log_Ws - pred_area)**2)
SS_tot = np.sum((log_Ws - np.mean(log_Ws))**2)
R2_area = 1 - SS_res_area / SS_tot if SS_tot > 0 else 0

# Perimeter fit: -ln<|W|> = μ·P + c
A_perim = np.column_stack([perimeters, np.ones(len(perimeters))])
coeffs_perim = np.linalg.lstsq(A_perim, log_Ws, rcond=None)[0]
pred_perim = A_perim @ coeffs_perim
SS_res_perim = np.sum((log_Ws - pred_perim)**2)
R2_perim = 1 - SS_res_perim / SS_tot if SS_tot > 0 else 0

# Combined: -ln<|W|> = σ·A + μ·P + c
A_both = np.column_stack([areas, perimeters, np.ones(len(areas))])
coeffs_both = np.linalg.lstsq(A_both, log_Ws, rcond=None)[0]
pred_both = A_both @ coeffs_both
SS_res_both = np.sum((log_Ws - pred_both)**2)
R2_both = 1 - SS_res_both / SS_tot if SS_tot > 0 else 0

print(f"Area law:      -ln<|W|> = {coeffs_area[0]:.6f}·A + {coeffs_area[1]:.6f}")
print(f"               R² = {R2_area:.4f}")
print(f"Perimeter law: -ln<|W|> = {coeffs_perim[0]:.6f}·P + {coeffs_perim[1]:.6f}")
print(f"               R² = {R2_perim:.4f}")
print(f"Combined:      -ln<|W|> = {coeffs_both[0]:.6f}·A + {coeffs_both[1]:.6f}·P + {coeffs_both[2]:.6f}")
print(f"               R² = {R2_both:.4f}")

print(f"\nString tension σ = {coeffs_area[0]:.6f}")
print(f"Predicted: σ = -ln(1-K*) = {-np.log(23/30):.6f}")

if R2_area > R2_perim:
    print(f"\n*** AREA LAW WINS: R²(area)={R2_area:.4f} > R²(perim)={R2_perim:.4f} ***")
else:
    print(f"\n*** PERIMETER LAW WINS: R²(perim)={R2_perim:.4f} > R²(area)={R2_area:.4f} ***")
