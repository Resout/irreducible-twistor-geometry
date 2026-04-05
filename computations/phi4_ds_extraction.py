#!/usr/bin/env python3
"""
φ⁴ scalar field DS extraction test — Prediction 2 (critical falsification).

The framework predicts:
  - Gauge theory (non-abelian): δ > 1/3 for ALL pairs → K* = 7/30 (mass gap)
  - φ⁴ (no gauge invariance): δ > 1/3 for correlated pairs, δ ≈ 1/3 for
    independent pairs → K* = 7/30 only for correlated, K* = 0 for independent

The test: measure δ at various separations in φ⁴ and test whether each is
significantly above 1/3 using bootstrap confidence intervals.

If ALL pairs in φ⁴ show δ significantly above 1/3, the framework's gauge/non-gauge
distinction is undermined.
"""
import numpy as np

H = 3


# ============================================================
# 4D φ⁴ lattice simulation (Metropolis)
# ============================================================
class Phi4Lattice:
    def __init__(self, L, kappa, lam):
        self.L = L
        self.kappa = kappa
        self.lam = lam
        self.phi = np.random.randn(L, L, L, L) * 0.5

    def local_action(self, x, y, z, t, phi_val):
        L = self.L
        nn_sum = (
            self.phi[(x+1)%L, y, z, t] + self.phi[(x-1)%L, y, z, t] +
            self.phi[x, (y+1)%L, z, t] + self.phi[x, (y-1)%L, z, t] +
            self.phi[x, y, (z+1)%L, t] + self.phi[x, y, (z-1)%L, t] +
            self.phi[x, y, z, (t+1)%L] + self.phi[x, y, z, (t-1)%L]
        )
        return -self.kappa * phi_val * nn_sum + phi_val**2 + self.lam * (phi_val**2 - 1)**2

    def sweep(self, delta=1.5):
        L = self.L
        accepted = 0
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        old_val = self.phi[x, y, z, t]
                        new_val = old_val + np.random.uniform(-delta, delta)
                        dS = (self.local_action(x, y, z, t, new_val) -
                              self.local_action(x, y, z, t, old_val))
                        if dS < 0 or np.random.random() < np.exp(-dS):
                            self.phi[x, y, z, t] = new_val
                            accepted += 1
        return accepted / L**4

    def thermalize(self, n_sweeps=200, delta=1.5):
        for _ in range(n_sweeps):
            self.sweep(delta)

    def time_slice_average(self, t):
        return np.mean(self.phi[:, :, :, t])

    def spatial_subvolume_average(self, t, half):
        L = self.L
        mid = L // 2
        if half == 0:
            return np.mean(self.phi[:mid, :, :, t])
        else:
            return np.mean(self.phi[mid:, :, :, t])


# ============================================================
# DS extraction: tercile binning → δ
# ============================================================
def bin_tercile(vals):
    sorted_idx = np.argsort(vals)
    bins = np.zeros(len(vals), dtype=int)
    n = len(vals)
    for rank, idx in enumerate(sorted_idx):
        bins[idx] = min(rank * H // n, H - 1)
    return bins


def compute_delta(obs1, obs2):
    b1 = bin_tercile(obs1)
    b2 = bin_tercile(obs2)
    C = np.zeros((H, H))
    for i in range(len(b1)):
        C[b1[i], b2[i]] += 1
    row_sums = C.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    C = C / row_sums
    return np.trace(C) / H, C


def bootstrap_delta(obs1, obs2, n_boot=2000):
    """Bootstrap CI for δ."""
    n = len(obs1)
    deltas = np.zeros(n_boot)
    for b in range(n_boot):
        idx = np.random.randint(0, n, size=n)
        d, _ = compute_delta(obs1[idx], obs2[idx])
        deltas[b] = d
    return deltas


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("φ⁴ SCALAR FIELD — DS EXTRACTION TEST (Prediction 2)")
    print("=" * 70)

    L = 8
    kappa = 0.15
    lam = 1.0
    n_configs = 800  # more stats than before
    n_therm = 300
    n_skip = 5

    print(f"\nParameters: L={L}, κ={kappa}, λ={lam}")
    print(f"Configs: {n_configs}, thermalization: {n_therm}, skip: {n_skip}")

    print("\nThermalizing...")
    lat = Phi4Lattice(L, kappa, lam)
    lat.thermalize(n_therm)

    # Collect all time-slice averages for every config
    all_slices = np.zeros((n_configs, L))
    sub0 = np.zeros(n_configs)
    sub1 = np.zeros(n_configs)

    print(f"Generating {n_configs} configurations...")
    for cfg in range(n_configs):
        for _ in range(n_skip):
            lat.sweep()
        for t in range(L):
            all_slices[cfg, t] = lat.time_slice_average(t)
        sub0[cfg] = lat.spatial_subvolume_average(0, 0)
        sub1[cfg] = lat.spatial_subvolume_average(0, 1)
        if (cfg + 1) % 200 == 0:
            print(f"  {cfg+1}/{n_configs}")

    # ============================================================
    print("\n" + "=" * 70)
    print("SEPARATION SCAN with bootstrap significance")
    print("=" * 70)
    print(f"  Null hypothesis: δ = 1/H = {1/H:.4f}")
    print(f"  Structural correlation: δ significantly > 1/H")
    print()
    print(f"  {'Δt':>4}  {'δ':>7}  {'95% CI':>16}  {'p(δ>1/3)':>9}  {'Pearson r':>10}  {'Verdict':>12}")
    print(f"  {'---':>4}  {'---':>7}  {'---':>16}  {'---':>9}  {'---':>10}  {'---':>12}")

    ref = all_slices[:, 0]
    results = []

    for dt in range(L):
        target = all_slices[:, dt]
        delta_obs, C = compute_delta(ref, target)
        boot = bootstrap_delta(ref, target, n_boot=3000)
        ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])
        p_above = np.mean(boot > 1/H)
        r = np.corrcoef(ref, target)[0, 1]

        if dt == 0:
            verdict = "SELF"
        elif ci_lo > 1/H + 0.005:
            verdict = "STRUCTURAL"
        elif ci_hi < 1/H + 0.005:
            verdict = "INDEPENDENT"
        else:
            verdict = "MARGINAL"

        results.append((dt, delta_obs, ci_lo, ci_hi, p_above, r, verdict))
        label = " ← adj" if dt == 1 else " ← max" if dt == L//2 else ""
        print(f"  {dt:4d}  {delta_obs:7.4f}  [{ci_lo:.4f}, {ci_hi:.4f}]  {p_above:9.3f}  {r:+10.4f}  {verdict:>12}{label}")

    # Spatial sub-volumes
    print()
    delta_sub, C_sub = compute_delta(sub0, sub1)
    boot_sub = bootstrap_delta(sub0, sub1, n_boot=3000)
    ci_lo_s, ci_hi_s = np.percentile(boot_sub, [2.5, 97.5])
    p_above_s = np.mean(boot_sub > 1/H)

    if ci_lo_s > 1/H + 0.005:
        v_sub = "STRUCTURAL"
    elif ci_hi_s < 1/H + 0.005:
        v_sub = "INDEPENDENT"
    else:
        v_sub = "MARGINAL"

    r_sub = np.corrcoef(sub0, sub1)[0, 1]
    print(f"  {'sub':>4}  {delta_sub:7.4f}  [{ci_lo_s:.4f}, {ci_hi_s:.4f}]  {p_above_s:9.3f}  {r_sub:+10.4f}  {v_sub:>12}")

    # ============================================================
    # Compare to gauge theory reference values
    print("\n" + "=" * 70)
    print("COMPARISON TO GAUGE THEORY")
    print("=" * 70)
    print(f"  SU(2) cross-rep δ (L=4,6,8):  0.343 ± 0.001  [ALWAYS > 1/3]")
    print(f"  U(1) coprime δ (L=4,6,8):      0.335 ± 0.004  [approaching 1/3]")
    print()

    n_structural = sum(1 for _, _, _, _, _, _, v in results if v == "STRUCTURAL" and _ > 0)
    n_independent = sum(1 for _, _, _, _, _, _, v in results if v == "INDEPENDENT")
    n_marginal = sum(1 for _, _, _, _, _, _, v in results if v == "MARGINAL")

    print(f"  φ⁴ results (Δt > 0): {n_structural} structural, {n_independent} independent, {n_marginal} marginal")

    # ============================================================
    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)

    if n_structural == L - 1:  # all non-self pairs structural
        print("  FALSIFICATION WARNING: ALL pairs show δ significantly > 1/3")
        print("  φ⁴ behaves like gauge theory — the gauge/non-gauge distinction fails.")
    elif n_structural == 0:
        print("  CONSISTENT WITH PREDICTION:")
        print("  No pairs in φ⁴ show structural correlation (δ significantly > 1/3).")
        print("  φ⁴ has no gauge constraint → no universal structural dependence.")
        print("  Compare: SU(2) shows δ = 0.343 systematically above 1/3 at all L.")
    elif n_independent > 0:
        print("  CONSISTENT WITH PREDICTION:")
        print(f"  {n_structural} correlated pairs (short range), {n_independent} independent pairs.")
        print("  φ⁴ has some correlation (kinetic coupling) but NOT universal.")
        print("  Gauge theory prediction: ALL pairs structural. φ⁴: only some.")
    else:
        print(f"  INCONCLUSIVE: {n_structural} structural, {n_marginal} marginal.")
        print("  Needs more statistics or larger lattice.")


if __name__ == '__main__':
    np.random.seed(42)
    main()
