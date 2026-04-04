"""Compute the linearised eigenvalue of the det=0 instability."""
import numpy as np

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


# Base point on det=0
b = 0.01 - 0.3j
denom = -2 + 4*b
a = -(1 - 4*b + 2*b**2) / denom
theta = 1 - a - 2*b
m0 = np.array([a, b, b, theta], dtype=complex)
print(f"Base point: |det| = {abs(det_M(m0)):.2e}")

# Transverse direction
grad_det = np.array([-m0[0], -m0[1], -m0[2], m0[3]])
grad_det = grad_det / np.linalg.norm(grad_det)

# Amplification ratio at multiple epsilon scales
print("\nAmplification ratio (det_out / det_in) at different scales:")
for eps in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]:
    m_pert = m0 + eps * grad_det
    m_pert = m_pert / np.sum(m_pert)
    det_in = abs(det_M(m_pert))
    m_out = F_map(m_pert)
    det_out = abs(det_M(m_out))
    if det_in > 1e-30:
        ratio = det_out / det_in
        print(f"  eps={eps:.0e}  |det_in|={det_in:.4e}  |det_out|={det_out:.4e}  ratio={ratio:.6f}")

# Amplification at many different det=0 points
print("\nAmplification ratio across different det=0 states:")
np.random.seed(42)
ratios = []
for trial in range(200):
    b_t = (0.05 * np.random.randn()) + (0.05 * np.random.randn() - 0.2) * 1j
    denom_t = -2 + 4 * b_t
    if abs(denom_t) < 0.01: continue
    a_t = -(1 - 4*b_t + 2*b_t**2) / denom_t
    theta_t = 1 - a_t - 2*b_t
    m_t = np.array([a_t, b_t, b_t, theta_t], dtype=complex)
    if abs(np.sum(m_t) - 1) > 1e-8: continue
    if abs(det_M(m_t)) > 1e-6: continue

    g = np.array([-m_t[0], -m_t[1], -m_t[2], m_t[3]])
    g = g / np.linalg.norm(g)
    m_p = m_t + 1e-8 * g
    m_p = m_p / np.sum(m_p)

    det_in = abs(det_M(m_p))
    m_out = F_map(m_p)
    det_out = abs(det_M(m_out))

    if det_in > 1e-20:
        ratios.append(det_out / det_in)

ratios = np.array(ratios)
print(f"  n = {len(ratios)}")
print(f"  mean = {ratios.mean():.4f}")
print(f"  std  = {ratios.std():.4f}")
print(f"  min  = {ratios.min():.4f}")
print(f"  max  = {ratios.max():.4f}")
print(f"  ALL > 1: {(ratios > 1).all()}")

# Check if the ratio is a derived constant
print(f"\nChecking if ratio is a derived constant:")
r = ratios.mean()
print(f"  ratio = {r:.4f}")
print(f"  1/(1-K*)^2 = {1/(1-7/30)**2:.4f}")
print(f"  (H^2+1) = {10}")
print(f"  H^3 = {27}")
print(f"  H^2+1 / (1-K*)^2 = {10/(23/30)**2:.4f}")
print(f"  1/eta = {27/4:.4f}")
print(f"  H/(H-1) = {3/2:.4f}")

# Also try asymmetric det=0 states (s2 != s3)
print("\nAsymmetric det=0 states:")
ratios2 = []
for trial in range(200):
    s1_t = 0.3 + 0.1*np.random.randn() + (0.1*np.random.randn())*1j
    s2_t = 0.1 + 0.05*np.random.randn() + (0.05*np.random.randn())*1j
    s3_t = 0.15 + 0.05*np.random.randn() + (0.05*np.random.randn())*1j
    theta_t = 1 - s1_t - s2_t - s3_t
    m_t = np.array([s1_t, s2_t, s3_t, theta_t], dtype=complex)

    # Project onto det=0: adjust theta so theta^2 = sum s_i^2
    ss = np.sum(m_t[:3]**2)
    theta_new = np.sqrt(ss)  # one of the two square roots
    m_t[3] = theta_new
    m_t = m_t / np.sum(m_t)  # re-normalize L1

    if abs(det_M(m_t)) > 1e-4: continue  # projection didn't work well

    g = np.array([-m_t[0], -m_t[1], -m_t[2], m_t[3]])
    g = g / np.linalg.norm(g)
    m_p = m_t + 1e-8 * g
    m_p = m_p / np.sum(m_p)

    det_in = abs(det_M(m_p))
    m_out = F_map(m_p)
    det_out = abs(det_M(m_out))

    if det_in > 1e-20:
        ratios2.append(det_out / det_in)

r2 = np.array(ratios2)
if len(r2) > 0:
    print(f"  n = {len(r2)}")
    print(f"  mean = {r2.mean():.4f}")
    print(f"  std  = {r2.std():.4f}")
    print(f"  ALL > 1: {(r2 > 1).all()}")
