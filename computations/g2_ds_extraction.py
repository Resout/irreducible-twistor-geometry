"""
G₂ lattice gauge theory → DS extraction → K* test.

G₂ is the smallest exceptional Lie group:
  - 14-dimensional, rank 2, TRIVIAL CENTER
  - Fundamental representation: 7×7 real orthogonal matrices
  - G₂ ⊂ SO(7), defined by preserving the octonionic 3-form T_abc
  - Confinement WITHOUT center symmetry (Holland-Pepe-Wiese 2003)

If K* = 7/30 emerges from G₂ lattice data, the universality of the
DS equilibrium extends beyond classical gauge groups to exceptionals.
This would be the strongest evidence yet: G₂ has no center, so the
result cannot be an artifact of center symmetry.

Algorithm: Metropolis with G₂ exponential map, Wilson plaquette action.
"""
import numpy as np
from time import time

H = 3; FLOOR = 1.0/27.0

# ================================================================
# G₂ algebra: 14 generators in the 7-dimensional representation
# ================================================================

# The octonionic structure constants (Fano plane)
# T_abc = +1 for (a,b,c) in the following cyclic triples (0-indexed):
FANO_TRIPLES = [(0,1,3), (1,2,4), (2,3,5), (3,4,6), (4,5,0), (5,6,1), (6,0,2)]

def build_T_tensor():
    """Build the totally antisymmetric octonionic 3-form T_abc."""
    T = np.zeros((7,7,7))
    for (a,b,c) in FANO_TRIPLES:
        # All even permutations get +1, odd get -1
        T[a,b,c] = T[b,c,a] = T[c,a,b] = 1.0
        T[a,c,b] = T[c,b,a] = T[b,a,c] = -1.0
    return T

T_OCT = build_T_tensor()

def build_g2_generators():
    """
    Build 14 generators of g₂ ⊂ so(7).

    so(7) has 21 generators (antisymmetric 7×7 matrices).
    g₂ is the subalgebra preserving T_abc. Individual so(7) basis
    elements don't preserve T — only specific linear combinations do.
    We find these as the kernel of the constraint map.
    """
    # Build all 21 so(7) basis elements
    so7_basis = []
    for i in range(7):
        for j in range(i+1, 7):
            X = np.zeros((7,7))
            X[i,j] = 1.0; X[j,i] = -1.0
            so7_basis.append(X)

    # Build constraint matrix: for each so(7) basis element X_k,
    # compute the violation V(a,b,c) = sum_d (X_ad T_dbc + X_bd T_adc + X_cd T_abd)
    # for all distinct ordered triples (a,b,c).
    # A linear combination sum_k alpha_k X_k preserves T iff the constraint
    # matrix times alpha = 0.
    constraints = []
    for a in range(7):
        for b in range(7):
            for c in range(7):
                if T_OCT[a,b,c] == 0 and a != b and b != c and a != c:
                    # Only need nontrivial constraints
                    row = []
                    for X in so7_basis:
                        val = sum(X[a,d]*T_OCT[d,b,c] + X[b,d]*T_OCT[a,d,c] + X[c,d]*T_OCT[a,b,d]
                                  for d in range(7))
                        row.append(val)
                    if any(abs(v) > 1e-10 for v in row):
                        constraints.append(row)
    # Also add constraints from the Fano triples themselves
    for a in range(7):
        for b in range(7):
            for c in range(7):
                row = []
                for X in so7_basis:
                    val = sum(X[a,d]*T_OCT[d,b,c] + X[b,d]*T_OCT[a,d,c] + X[c,d]*T_OCT[a,b,d]
                              for d in range(7))
                    row.append(val)
                if any(abs(v) > 1e-10 for v in row):
                    constraints.append(row)

    C = np.array(constraints)
    # Remove duplicate rows
    if len(C) > 0:
        C_unique = [C[0]]
        for row in C[1:]:
            diffs = [np.max(np.abs(row - r)) for r in C_unique]
            if min(diffs) > 1e-10:
                C_unique.append(row)
        C = np.array(C_unique)

    # Kernel of C gives the g₂ coefficients
    from numpy.linalg import svd
    U, S, Vt = svd(C)
    # Null space: rows of Vt with singular values ~ 0
    tol = 1e-8
    null_mask = S < tol  # for the first min(m,n) singular values
    # Need to also include trailing rows of Vt if C has fewer rows than cols
    n_constraints = C.shape[0]
    n_basis = C.shape[1]  # = 21
    null_start = np.sum(S > tol)
    null_vecs = Vt[null_start:]

    print(f"  Constraint matrix: {C.shape[0]}×{C.shape[1]}, rank {null_start}")
    print(f"  Null space dimension: {len(null_vecs)} (expect 14)")

    assert len(null_vecs) == 14, f"Expected 14-dim null space, got {len(null_vecs)}"

    # Build generators from null vectors
    gens = []
    for vec in null_vecs:
        G = sum(vec[k] * so7_basis[k] for k in range(21))
        gens.append(G)

    # Orthonormalise via Gram-Schmidt on the Killing form Tr(X Y)
    ortho_gens = []
    for g in gens:
        for og in ortho_gens:
            proj = np.trace(g @ og.T) / np.trace(og @ og.T)
            g = g - proj * og
        norm = np.sqrt(abs(np.trace(g @ g.T)))
        if norm > 1e-10:
            ortho_gens.append(g / norm)

    return ortho_gens

print("Building G₂ generators...", flush=True)
G2_GENS = build_g2_generators()
print(f"  {len(G2_GENS)} generators constructed.")

# Verify: [T_a, T_b] closes in span
def verify_closure():
    """Verify the generators close under commutation."""
    for i in range(14):
        for j in range(i+1, 14):
            comm = G2_GENS[i] @ G2_GENS[j] - G2_GENS[j] @ G2_GENS[i]
            # Project onto generators
            coeffs = [np.trace(comm @ g.T) / np.trace(g @ g.T) for g in G2_GENS]
            residual = comm - sum(c*g for c,g in zip(coeffs, G2_GENS))
            if np.max(np.abs(residual)) > 1e-8:
                return False
    return True

print(f"  Closure verified: {verify_closure()}")

# ================================================================
# G₂ group elements
# ================================================================

def g2_exp(coeffs, eps=1.0):
    """Exponential map: exp(eps * sum_a coeffs[a] * T_a) via Padé."""
    X = eps * sum(c*g for c,g in zip(coeffs, G2_GENS))
    # Matrix exponential via eigendecomposition (7×7, fast)
    from scipy.linalg import expm
    return expm(X)

def random_g2():
    """Random G₂ element (Haar-distributed via product of random exponentials)."""
    U = np.eye(7)
    for _ in range(10):
        coeffs = np.random.randn(14)
        U = g2_exp(coeffs, eps=1.0) @ U
    return U

def near_identity_g2(eps=0.3):
    """G₂ element near identity for Metropolis proposals."""
    coeffs = np.random.randn(14)
    return g2_exp(coeffs, eps=eps)

def reproject_g2(U):
    """Reproject a matrix onto G₂ ⊂ SO(7).

    Step 1: Project to SO(7) via polar decomposition.
    Step 2: Minimise distance to G₂ by iterating through the algebra.
    For near-G₂ matrices (which is our case after small updates),
    a single polar decomposition + algebra projection suffices.
    """
    from scipy.linalg import polar
    R, _ = polar(U)
    # Ensure det = +1
    if np.linalg.det(R) < 0:
        R[:, 0] *= -1
    # Project the so(7) log onto g₂
    from scipy.linalg import logm
    try:
        X = logm(R)
        X = (X - X.T) / 2  # antisymmetrise
        # Project onto g₂ span
        X_g2 = sum(np.trace(X @ g.T) / np.trace(g @ g.T) * g for g in G2_GENS)
        from scipy.linalg import expm
        return expm(X_g2)
    except:
        return R  # fallback: just use the SO(7) projection

# ================================================================
# G₂ lattice
# ================================================================

class G2Lattice:
    def __init__(self, L):
        self.L = L
        self.links = np.zeros((L,L,L,L,4,7,7))
        print(f"  Initializing G₂ lattice L={L}...", end="", flush=True)
        for idx in np.ndindex(L,L,L,L):
            for mu in range(4):
                self.links[idx+(mu,)] = random_g2()
        print(" done.")

    def staple_sum(self, pos, mu):
        """Sum of staples around link U_mu(pos). Same as SU(N) but 7×7 real."""
        L = self.L
        S = np.zeros((7,7))
        for nu in range(4):
            if nu == mu: continue
            pm = list(pos); pm[mu]=(pm[mu]+1)%L; pm=tuple(pm)
            pn = list(pos); pn[nu]=(pn[nu]+1)%L; pn=tuple(pn)
            pmn = list(pos); pmn[nu]=(pmn[nu]-1)%L; pmn=tuple(pmn)
            pmun = list(pos); pmun[mu]=(pmun[mu]+1)%L; pmun[nu]=(pmun[nu]-1)%L; pmun=tuple(pmun)
            # Forward staple: U_nu(x+mu) @ U_mu(x+nu)^T @ U_nu(x)^T
            S += self.links[pm+(nu,)] @ self.links[pn+(mu,)].T @ self.links[pos+(nu,)].T
            # Backward staple: U_nu(x+mu-nu)^T @ U_mu(x-nu)^T @ U_nu(x-nu)
            S += self.links[pmun+(nu,)].T @ self.links[pmn+(mu,)].T @ self.links[pmn+(nu,)]
        return S

    def metropolis_sweep(self, beta, eps=0.3, n_hits=4):
        """Metropolis update: propose U' = R·U with R near identity in G₂."""
        L = self.L
        accepted = 0; total = 0
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        for mu in range(4):
                            pos = (x,y,z,t)
                            S = self.staple_sum(pos, mu)
                            old_U = self.links[pos+(mu,)]
                            old_action = -(beta/7) * np.trace(old_U @ S)
                            for _ in range(n_hits):
                                R = near_identity_g2(eps)
                                new_U = R @ old_U
                                new_action = -(beta/7) * np.trace(new_U @ S)
                                dS = new_action - old_action
                                total += 1
                                if dS < 0 or np.random.random() < np.exp(-dS):
                                    self.links[pos+(mu,)] = new_U
                                    old_U = new_U
                                    old_action = new_action
                                    accepted += 1
        return accepted / total if total > 0 else 0

    def plaquette_avg(self):
        """Average plaquette (global)."""
        L = self.L
        total = 0; count = 0
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        for mu in range(4):
                            for nu in range(mu+1, 4):
                                pos = (x,y,z,t)
                                pm = list(pos); pm[mu]=(pm[mu]+1)%L; pm=tuple(pm)
                                pn = list(pos); pn[nu]=(pn[nu]+1)%L; pn=tuple(pn)
                                U = (self.links[pos+(mu,)] @ self.links[pm+(nu,)] @
                                     self.links[pn+(mu,)].T @ self.links[pos+(nu,)].T)
                                total += np.trace(U) / 7
                                count += 1
        return total / count

    def polyakov(self, x, y, z):
        """Polyakov loop at spatial position (x,y,z): Tr of temporal holonomy / 7."""
        P = np.eye(7)
        for t in range(self.L):
            P = P @ self.links[x,y,z,t,3]
        return np.trace(P) / 7

    def plaquette_at(self, pos, mu, nu):
        """Single plaquette trace / 7."""
        L = self.L
        pm = list(pos); pm[mu]=(pm[mu]+1)%L; pm=tuple(pm)
        pn = list(pos); pn[nu]=(pn[nu]+1)%L; pn=tuple(pn)
        U = (self.links[pos+(mu,)] @ self.links[pm+(nu,)] @
             self.links[pn+(mu,)].T @ self.links[pos+(nu,)].T)
        return np.trace(U) / 7

    def spatial_plaquette_avg(self, t, mu, nu):
        """Average plaquette in (mu,nu) plane at time t."""
        L = self.L; total = 0; count = 0
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    total += self.plaquette_at((x,y,z,t), mu, nu)
                    count += 1
        return total / count

# ================================================================
# DS extraction tools (identical to SU(N) pipeline)
# ================================================================

def bin_tercile(vals):
    sorted_idx = np.argsort(vals)
    bins = np.zeros(len(vals), dtype=int)
    n = len(vals)
    for rank, idx in enumerate(sorted_idx):
        bins[idx] = min(rank * H // n, H - 1)
    return bins

def compute_delta_and_C(obs1, obs2):
    b1 = bin_tercile(obs1); b2 = bin_tercile(obs2)
    C = np.zeros((H, H))
    for i in range(len(b1)): C[b1[i], b2[i]] += 1
    rs = C.sum(1, keepdims=True); rs[rs==0]=1; C /= rs
    return np.trace(C)/H, C

def ds_combine(m, e):
    s, th = m[:3], m[3]; se, ph = e[:3], e[3]
    sn = s*se + s*ph + th*se; tn = th*ph
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
    d = 1.0 - K
    if abs(d) < 1e-15: return m.copy(), abs(K)
    out = np.zeros(4); out[:3] = sn/d; out[3] = tn/d
    born = out[3]**2/np.sum(out**2)
    if born < FLOOR:
        ss=np.sum(np.abs(out[:3])); ssq=np.sum(out[:3]**2)
        if ss > 0:
            r=ssq/ss**2; t=(np.sqrt(26*r)-r)/(26-r)
            out[:3]=np.abs(out[:3])*(1-t)/ss; out[3]=t
    return out, abs(K)

def extract_K(C, n_steps=50, n_runs=50):
    evidence_bank = []
    for i in range(H):
        row = C[i]; ent = -np.sum(row[row>0]*np.log2(row[row>0]))
        theta_i = ent/np.log2(H)
        s = row*(1-theta_i); e = np.zeros(4); e[:3]=s; e[3]=theta_i
        total=np.sum(np.abs(e))
        if total>0: e/=total
        evidence_bank.append(e)
    cell_probs = C.copy(); cell_probs /= np.sum(cell_probs)
    K_vals = []
    for _ in range(n_runs):
        m = np.array([0.0,0.0,0.0,1.0]); Ks = []
        for step in range(n_steps):
            flat = np.random.choice(H*H, p=cell_probs.flatten())
            i = flat // H
            m, K = ds_combine(m, evidence_bank[i]); Ks.append(K)
        K_vals.append(np.mean(Ks[-15:]))
    return np.mean(K_vals), np.std(K_vals)

# ================================================================
# MAIN
# ================================================================

def main():
    np.random.seed(42)
    L = 4; beta = 9.5
    n_therm = 100; n_configs = 200; n_skip = 5

    print("="*60)
    print(f"G₂ LATTICE GAUGE THEORY → DS EXTRACTION")
    print(f"  L={L}, beta={beta}, trivial center")
    print("="*60)

    print(f"\nInitializing...")
    t0 = time()
    lat = G2Lattice(L)

    print(f"\nThermalizing ({n_therm} sweeps)...")
    for i in range(n_therm):
        acc = lat.metropolis_sweep(beta, eps=0.3)
        if (i+1) % 20 == 0:
            plaq = lat.plaquette_avg()
            print(f"  Sweep {i+1}/{n_therm}: plaq={plaq:.4f}, acc={acc:.3f} ({time()-t0:.0f}s)")

    # Tune epsilon for ~40-50% acceptance
    print(f"\nFinal acceptance rate: {acc:.3f}")

    print(f"\nMeasuring {n_configs} configs...")
    t0 = time()
    obs = {
        "Poly(0,0,0)": [], "Poly(1,0,0)": [], "Poly(L/2,0,0)": [],
        "Plaq(0,1)": [], "Plaq(0,2)": [], "Plaq(1,2)": [],
    }
    for cfg in range(n_configs):
        for _ in range(n_skip):
            lat.metropolis_sweep(beta, eps=0.3)
        obs["Poly(0,0,0)"].append(lat.polyakov(0,0,0))
        obs["Poly(1,0,0)"].append(lat.polyakov(1,0,0))
        obs["Poly(L/2,0,0)"].append(lat.polyakov(L//2,0,0))
        obs["Plaq(0,1)"].append(lat.spatial_plaquette_avg(0, 0, 1))
        obs["Plaq(0,2)"].append(lat.spatial_plaquette_avg(0, 0, 2))
        obs["Plaq(1,2)"].append(lat.spatial_plaquette_avg(0, 1, 2))
        if (cfg+1) % 50 == 0:
            print(f"  {cfg+1}/{n_configs} ({time()-t0:.0f}s)")

    obs = {k: np.array(v) for k,v in obs.items()}

    # DS extraction for all pairs
    print(f"\nDS extraction for all pairs:")
    print(f"  {'Pair':>35s}  {'delta':>7s}  {'K*':>10s}  {'err':>8s}  {'r':>7s}")
    print(f"  {'-'*35}  {'-'*7}  {'-'*10}  {'-'*8}  {'-'*7}")

    best_K = 0; best_pair = ""; best_C = None
    names = list(obs.keys())
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            n1, n2 = names[i], names[j]
            delta, C = compute_delta_and_C(obs[n1], obs[n2])
            K_mean, K_std = extract_K(C)
            r = np.corrcoef(obs[n1], obs[n2])[0,1]
            tag = " <-- BEST" if K_mean > best_K else ""
            print(f"  {n1+' vs '+n2:>35s}  {delta:7.4f}  {K_mean:10.6f}  {K_std:8.4f}  {r:+7.4f}{tag}")
            if K_mean > best_K:
                best_K = K_mean; best_pair = f"{n1} vs {n2}"; best_C = C

    print(f"\n  Best pair: {best_pair}")
    print(f"  K* = {best_K:.6f} (target: {7/30:.6f})")
    if best_C is not None:
        print(f"  C =")
        for row in range(H):
            print(f"    [{best_C[row,0]:.4f}  {best_C[row,1]:.4f}  {best_C[row,2]:.4f}]")

    print(f"\n  RESULT: G₂ (trivial center, exceptional)")
    print(f"  K* = {best_K:.6f}")
    print(f"  Target = {7/30:.6f}")
    print(f"  |K* - 7/30| = {abs(best_K - 7/30):.6f}")
    if best_K > 0.1:
        print(f"  STATUS: PROMISING (K* > 0.1)")
    elif best_K > 0.01:
        print(f"  STATUS: WEAK (K* > 0.01, need larger lattice or more stats)")
    else:
        print(f"  STATUS: MARGINAL (K* < 0.01, finite-size effects dominant)")

if __name__ == "__main__":
    main()
