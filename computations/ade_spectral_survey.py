#!/usr/bin/env python3
"""
ADE Spectral Survey
===================
Computes the complete eigenvalue tower for all ADE Dynkin diagrams,
using the DS transfer operator at H=3, K*=7/30, Born floor=1/27.

Questions I'm trying to answer:
  1. What does E8 look like? (h(E8)=30 is structurally derived in this framework)
  2. Is the 1.082 state in A2 a pattern? Does it appear in all groups?
  3. Does the folding invariance hold for the FULL tower, not just rho?
  4. Are there more exact ratio results (like sqrt(2)) hiding in higher groups?
  5. What is the A_n -> A_inf spectral behaviour?

Written during freeform autonomous session, 2026-04-08.
"""

import numpy as np
from scipy.optimize import fsolve

# ============================================================
# DS PRIMITIVES
# ============================================================

def ds_combine(m, e):
    s = m[:3]; theta = m[3]; ev = e[:3]; phi = e[3]
    s_pre = s*ev + s*phi + theta*ev
    theta_pre = theta*phi
    K = 1.0 - np.sum(s_pre) - theta_pre
    if abs(1.0-K) < 1e-15: return m.copy(), K
    out = np.zeros(4)
    out[:3] = s_pre/(1.0-K); out[3] = theta_pre/(1.0-K)
    return out, K

def born_prob(m):
    L2 = np.sum(m**2)
    return m[3]**2/L2 if L2 > 1e-30 else 0.0

def enforce_floor(m, fv=1.0/27.0):
    if born_prob(m) >= fv - 1e-14: return m.copy()
    S = np.sum(m[:3])
    if S < 1e-15: return m.copy()
    Sq = np.sum(m[:3]**2)
    A = 26.0*S**2-Sq; B = 2.0*Sq; C = -Sq
    disc = B**2-4*A*C
    if disc < 0: return m.copy()
    cands = [t for t in [(-B+np.sqrt(disc))/(2*A), (-B-np.sqrt(disc))/(2*A)] if 0<t<1]
    if not cands: return m.copy()
    t = min(cands, key=lambda x: abs(x-m[3]))
    out = np.zeros(4); out[:3] = m[:3]*(1-t)/S; out[3] = t
    return out

def ds_step(m, e):
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds), K

# ============================================================
# SINGLE-SITE EQUILIBRIUM
# ============================================================

def find_single_site_eq():
    def eqs(p):
        s1,th,w1,phi = p
        s2=(1-s1-th)/2; w2=(1-w1-phi)/2
        m=np.array([s1,s2,s2,th]); e=np.array([w1,w2,w2,phi])
        eq1=th**2/(s1**2+2*s2**2+th**2)-1/27
        eq2=phi**2/(w1**2+2*w2**2+phi**2)-1/27
        eq3=2*s1*w2+2*s2*w1+2*s2*w2-7/30
        eq4=ds_step(m,e)[0][0]-s1
        return [eq1,eq2,eq3,eq4]
    sol = fsolve(eqs,[0.787,0.155,0.631,0.129])
    s1,th,w1,phi=sol; s2=(1-s1-th)/2; w2=(1-w1-phi)/2
    return np.array([s1,s2,s2,th]), np.array([w1,w2,w2,phi])

m_star, e_star = find_single_site_eq()
uniform = np.array([0.25]*4)
g = 7.0/30.0

# ============================================================
# ADE DYNKIN DIAGRAMS (all simply-laced)
# ============================================================

def make_An(n):
    C = 2*np.eye(n, dtype=int)
    for i in range(n-1): C[i,i+1]=C[i+1,i]=-1
    return C

def make_Dn(n):
    """D_n: chain 0-1-..-(n-3), then fork (n-3)-(n-2) and (n-3)-(n-1)"""
    assert n>=4
    C = 2*np.eye(n, dtype=int)
    for i in range(n-3): C[i,i+1]=C[i+1,i]=-1  # chain up to fork
    C[n-3,n-2]=C[n-2,n-3]=-1  # fork branch 1
    C[n-3,n-1]=C[n-1,n-3]=-1  # fork branch 2
    return C

def make_E6():
    """E6: chain 0-1-2-3-4, branch from node 2 to node 5"""
    C = 2*np.eye(6, dtype=int)
    for i in range(4): C[i,i+1]=C[i+1,i]=-1
    C[2,5]=C[5,2]=-1
    return C

def make_E7():
    """E7: chain 0-1-2-3-4-5, branch from node 3 to node 6"""
    C = 2*np.eye(7, dtype=int)
    for i in range(5): C[i,i+1]=C[i+1,i]=-1
    C[3,6]=C[6,3]=-1
    return C

def make_E8():
    """E8: chain 0-1-2-3-4-5-6, branch from node 4 to node 7"""
    C = 2*np.eye(8, dtype=int)
    for i in range(6): C[i,i+1]=C[i+1,i]=-1
    C[4,7]=C[7,4]=-1
    return C

# ============================================================
# GENERAL COUPLED SYSTEM
# ============================================================

def coupled_step(state, cartan, e0=None, g=7.0/30.0):
    if e0 is None: e0 = e_star
    n = len(cartan)
    masses = [state[4*i:4*(i+1)] for i in range(n)]
    new_masses = []
    for i in range(n):
        e = e0.copy()
        for j in range(n):
            off = cartan[i,j]
            if off != 0 and j != i:
                e += g * off * (masses[j] - uniform)
        e = np.maximum(e, 1e-10)
        e /= np.sum(e)
        nm, _ = ds_step(masses[i], e)
        new_masses.append(nm)
    return np.concatenate(new_masses)

def find_eq(cartan, max_iter=200000, tol=1e-13):
    n = len(cartan)
    state = np.tile(m_star, n)
    for it in range(max_iter):
        ns = coupled_step(state, cartan)
        d = np.linalg.norm(ns - state)
        if d < tol:
            return state, it, d
        state = ns
    return state, max_iter, d

def jacobian(state, cartan, eps=1e-9):
    n = len(state)
    J = np.zeros((n,n))
    for j in range(n):
        sp=state.copy(); sm=state.copy()
        sp[j]+=eps; sm[j]-=eps
        J[:,j]=(coupled_step(sp,cartan)-coupled_step(sm,cartan))/(2*eps)
    return J

def get_spectrum(cartan, name="?"):
    n = len(cartan)
    state, niters, resid = find_eq(cartan)
    J = jacobian(state, cartan)
    evals = np.sort(np.abs(np.linalg.eigvals(J)))[::-1]
    sig = [e for e in evals if e > 1e-6]
    rho = sig[0] if sig else 0.0
    D0 = -np.log(rho) if rho > 1e-15 else np.inf
    ratios = [-np.log(e)/D0 if e > 1e-15 else np.inf for e in sig]
    return {
        'name': name, 'rank': n, 'state': state,
        'niters': niters, 'resid': resid,
        'evals': sig, 'rho': rho, 'D0': D0, 'ratios': ratios
    }

# ============================================================
# RUN SURVEY
# ============================================================

GROUPS = [
    ('A1', make_An(1)),
    ('A2', make_An(2)),
    ('A3', make_An(3)),
    ('A4', make_An(4)),
    ('A5', make_An(5)),
    ('A6', make_An(6)),
    ('D4', make_Dn(4)),
    ('D5', make_Dn(5)),
    ('D6', make_Dn(6)),
    ('D7', make_Dn(7)),
    ('E6', make_E6()),
    ('E7', make_E7()),
    ('E8', make_E8()),
]

print("=" * 70)
print("ADE SPECTRAL SURVEY")
print("DS transfer operator | H=3, K*=7/30, Born floor=1/27")
print("=" * 70)

results = {}
for name, C in GROUPS:
    print(f"  {name:4s} (rank {len(C)})...", end='', flush=True)
    r = get_spectrum(C, name)
    results[name] = r
    print(f" ρ={r['rho']:.5f}  Δ₀={r['D0']:.4f}  [{r['niters']} iters, resid={r['resid']:.1e}]")

# ============================================================
# DETAILED OUTPUT
# ============================================================

print()
print("=" * 70)
print("SPECTRAL RADIUS TABLE")
print("=" * 70)
print(f"{'Group':6s}  {'rank':4s}  {'rho':8s}  {'Delta0':8s}  {'n_evals':7s}  {'notes'}")
print("-" * 70)
for name, r in results.items():
    rr = results[name]
    notes = []
    if abs(rr['D0'] - 1.2626) < 0.01: notes.append("~A1_single")
    print(f"{name:6s}  {rr['rank']:4d}  {rr['rho']:.6f}  {rr['D0']:.4f}    {len(rr['evals']):3d}      {'  '.join(notes)}")

print()
print("=" * 70)
print("FULL EIGENVALUE TOWERS")
print("=" * 70)

for name, r in results.items():
    rr = results[name]
    print(f"\n{name} (rank {rr['rank']}, {len(rr['evals'])} non-trivial eigenvalues):")
    print(f"  ρ = {rr['rho']:.6f}   Δ₀ = {rr['D0']:.5f}")
    # Print in pairs for readability
    for i in range(0, len(rr['evals']), 4):
        chunk_e = rr['evals'][i:i+4]
        chunk_r = rr['ratios'][i:i+4]
        ev_str = "  ".join(f"|λ|={e:.4f} (Δ/Δ₀={r:.4f})" for e,r in zip(chunk_e, chunk_r))
        print(f"  [{i:2d}..{i+len(chunk_e)-1:2d}] {ev_str}")

# ============================================================
# RATIO ANALYSIS ACROSS GROUPS
# ============================================================

print()
print("=" * 70)
print("SECOND EIGENVALUE RATIO (the '1.082 pattern')")
print("=" * 70)
print("In A2, the second eigenvalue has ratio 1.082 (unidentified state)")
print("Does this appear in all groups? What is the pattern?")
print()
for name, r in results.items():
    rr = results[name]
    ratios = rr['ratios']
    if len(ratios) >= 2:
        print(f"  {name:4s}: ratio[1]/ratio[0] = {ratios[1]:.4f}", end='')
        if len(ratios) >= 4:
            print(f"   ratio[2]={ratios[2]:.4f}  ratio[3]={ratios[3]:.4f}", end='')
        print()

# ============================================================
# FOLDING INVARIANCE CHECK (Full Tower)
# ============================================================

print()
print("=" * 70)
print("FOLDING INVARIANCE: Does rho(fold(G)) = rho(G) for the FULL TOWER?")
print("=" * 70)

# Known ADE folding pairs (within simply-laced groups):
# A_{2n} folds to D_{n+1} ... but D is also simply-laced
# A4 folds to ... hmm. The standard ADE folding:
# The non-simply-laced come from folding simply-laced. For simply-laced only:
# D_4 -> A_2 (fold by Z_3 outer automorphism -> G_2, but G_2 is not ADE)
# Actually within ADE, the folding pairs are:
# A_{2n-1} folded by Z_2 -> C_n (not ADE)
# D_{n+1} folded by Z_2 -> B_n (not ADE)
# E_6 folded by Z_2 -> F_4 (not ADE)
# D_4 folded by Z_3 -> G_2 (not ADE)
# So within ADE, folding goes to non-ADE groups.
# But the PAPER says "folding invariance" applies to the ADE folding PAIRS
# for the spectral radius. Let me check the pairs mentioned in the paper.

# From memory: the paper claims rho(fold(G)) = rho(G) for ALL 10 folding pairs.
# The ADE folding pairs include:
# (A_2, A_2) <- trivial
# (D_4, A_2+A_1) <- fold D_4 by S_3
# (E_6, A_2) <- fold E_6 by Z_3
# (E_7, A_3) <- fold by Z_2
# (E_8, D_4) <- fold by Z_2

# Since we have all these groups, let's check spectral radius pairs:
print()
print("Checking spectral radius of related groups:")
folding_checks = [
    ("E8", "D4", "E8 folds to D4 (by Z_2 outer automorphism)"),
    ("E7", "A3", "E7 folds to ... (related)"),
    ("E6", "A2", "E6 folds to A2 (by Z_3) -> G_2, but rho comparison"),
    ("D4", "A2", "D4 folds to G_2 (Z_3), rho comparison"),
    ("D5", "A3", "D5 vs A3 spectral radii"),
    ("D6", "A5", "D6 vs A5"),
]
for g1, g2, note in folding_checks:
    r1 = results.get(g1, {}).get('rho', None)
    r2 = results.get(g2, {}).get('rho', None)
    if r1 and r2:
        diff = abs(r1-r2)
        print(f"  {g1:4s} ρ={r1:.5f}  vs  {g2:4s} ρ={r2:.5f}  diff={diff:.2e}  -- {note}")

# ============================================================
# EXACT RATIO SEARCH
# ============================================================

print()
print("=" * 70)
print("EXACT RATIO SEARCH")
print("Search for algebraically clean ratios in eigenvalue towers")
print("Comparing to: sqrt(2), 3/2, 4/3, 5/4, phi=(1+sqrt(5))/2, etc.")
print("=" * 70)

candidates = {
    'sqrt(2)': np.sqrt(2),
    '3/2':     1.5,
    '4/3':     4/3,
    '5/4':     1.25,
    'phi':     (1+np.sqrt(5))/2,
    '1+1/H':   1+1/3,
    '1+K*':    1+7/30,
    '1+Born':  1+1/27,
    '1+1/30':  1+1/30,
    'sqrt(3)': np.sqrt(3),
    '7/6':     7/6,
    '8/7':     8/7,
}

print()
for name, r in results.items():
    rr = results[name]
    ratios = rr['ratios']
    hits = []
    for i, ratio in enumerate(ratios):
        if ratio > 50: continue
        for cname, cval in candidates.items():
            if abs(ratio - cval) < 0.005:
                hits.append(f"  ratio[{i}]={ratio:.4f} ≈ {cname}={cval:.4f} (diff={abs(ratio-cval):.4f})")
    if hits:
        print(f"{name}:")
        for h in hits:
            print(h)

# ============================================================
# E8 SPECIAL ANALYSIS
# ============================================================

print()
print("=" * 70)
print("E8 SPECIAL ANALYSIS")
print("h(E8)=30=H(H^2+1) is structurally derived in this framework.")
print("Does the spectral structure reflect this?")
print("=" * 70)

e8 = results['E8']
print(f"\nE8 spectral radius: ρ = {e8['rho']:.8f}")
print(f"E8 mass gap:        Δ = {e8['D0']:.8f}")
print(f"1/h(E8) = 1/30    = {1/30:.8f}")
print(f"K*/h(E8)          = {7/30/30:.8f}")
print(f"e^(-1/30)         = {np.exp(-1/30):.8f}")
print(f"e^(-K*)           = {np.exp(-7/30):.8f}")
print(f"1 - K*            = {1 - 7/30:.8f}")

print(f"\nE8 eigenvalue tower ({len(e8['evals'])} states):")
for i, (ev, ratio) in enumerate(zip(e8['evals'], e8['ratios'])):
    mark = ""
    for cname, cval in candidates.items():
        if abs(ratio - cval) < 0.01:
            mark = f"  ← {cname}"
    print(f"  [{i:2d}] |λ|={ev:.6f}  Δ={-np.log(ev):.5f}  Δ/Δ₀={ratio:.5f}{mark}")

print(f"\nDoes Δ₀(E8) relate to h(E8)=30?")
D0_E8 = e8['D0']
print(f"  Δ₀ = {D0_E8:.6f}")
print(f"  Δ₀ × h = {D0_E8 * 30:.6f}   (if Δ₀ = 1/h, this should be 1)")
print(f"  Δ₀ × K* = {D0_E8 * 7/30:.6f}")
print(f"  e^(-Δ₀ × h) = {np.exp(-D0_E8 * 30):.6f}")

# ============================================================
# A_n LARGE-n LIMIT
# ============================================================

print()
print("=" * 70)
print("A_n LARGE-n LIMIT")
print("Paper: Δ_A(inf) ≈ 0.2498, NOT exactly 1/4 (diff = 2.4e-4)")
print("Investigating convergence...")
print("=" * 70)

an_results = [(name, r) for name, r in results.items() if name.startswith('A')]
an_results.sort(key=lambda x: int(x[0][1:]))
print(f"\n{'n':4s}  {'rho':8s}  {'Delta0':8s}  {'1/4 - Delta0':12s}")
for name, r in an_results:
    n = int(name[1:])
    D0 = r['D0']
    diff = 0.25 - D0
    print(f"  {n:2d}  {r['rho']:.6f}  {D0:.6f}  {diff:+.6f}")

print(f"\nConjectured limit: Δ_A(∞) = ?")
print(f"1/4 = {0.25:.6f}")
print(f"K*/h_A(∞) = ??? (Coxeter number of A_inf = ∞)")
# Note: h(A_n) = n+1. As n→∞, h→∞. So K*/h→0. Not the right formula.
# Try: does Δ_A(∞) = 1/4 - something·K*?
if len(an_results) >= 3:
    D0_large = [r['D0'] for _, r in an_results[-3:]]
    # Fit exponential convergence: D0(n) = D0_inf + a*exp(-b*n)
    ns = [int(name[1:]) for name, _ in an_results[-3:]]
    # Simple extrapolation
    print(f"\nLast three A_n Delta0 values: {[f'{d:.6f}' for d in D0_large]}")
    if len(D0_large) >= 2:
        trend = D0_large[-1] - D0_large[-2]
        print(f"Delta per step: {trend:.6f}")
        extrapolated = D0_large[-1] + trend * 10
        print(f"Rough extrapolation (10 more steps): {extrapolated:.6f}")

# ============================================================
# THE 1.082 STATE: WHAT IS IT ACROSS GROUPS?
# ============================================================

print()
print("=" * 70)
print("THE 1.082 STATE: ACROSS ALL GROUPS")
print("In A2, ratio[1]=1.082 is the node-antisymmetric angular mode (s2-s3)")
print("It has no lattice QCD match. Is it universal? Is it a pattern?")
print("=" * 70)

print()
print("Second-state ratio across all groups:")
for name, r in results.items():
    rr = results[name]
    if len(rr['ratios']) >= 2:
        r2 = rr['ratios'][1]
        print(f"  {name:4s} rank={rr['rank']}: ratio[1] = {r2:.4f}", end='')
        # Is this ratio approaching something?
        diff_from_1 = r2 - 1.0
        print(f"  (gap from ground: {diff_from_1:.4f})")
    else:
        print(f"  {name:4s} rank={rr['rank']}: only 1 eigenvalue")

print()
print("Pattern: does the second ratio grow/shrink with rank?")
print("If it approaches 1.000 as rank grows -> soft Goldstone-like mode")
print("If it approaches sqrt(2) -> same mechanism as 2++ ratio")
print("If it's constant -> universal feature of the framework")

# ============================================================
# KOOPMAN TWO-PARTICLE STATES (across groups)
# ============================================================

print()
print("=" * 70)
print("KOOPMAN TWO-PARTICLE STATES")
print("Products of eigenvalue pairs = two-glueball bound state energies?")
print("=" * 70)

for name, r in results.items():
    rr = results[name]
    evals = rr['evals'][:4]  # top 4
    D0 = rr['D0']
    if len(evals) < 2: continue
    print(f"\n{name}: two-particle states (lowest 3):")
    products = []
    for i in range(len(evals)):
        for j in range(i, len(evals)):
            prod = evals[i]*evals[j]
            if prod > 1e-6:
                delta_sum = -np.log(prod)
                ratio = delta_sum/D0
                products.append((ratio, i, j, prod))
    products.sort()
    for ratio, i, j, prod in products[:3]:
        print(f"  ({i},{j}): |λ|={prod:.4f}  Δ/Δ₀={ratio:.4f}")

# ============================================================
# EQUILIBRIUM STRUCTURE ACROSS GROUPS
# ============================================================

print()
print("=" * 70)
print("EQUILIBRIUM STRUCTURE: How does m* shift with coupling?")
print("=" * 70)

print(f"\nSingle-site equilibrium: m* = {m_star}")
print(f"  s1={m_star[0]:.4f}  s2=s3={m_star[1]:.4f}  theta={m_star[3]:.4f}")
print(f"  Ratio s1/s2 = {m_star[0]/m_star[1]:.2f}")
print(f"  Ratio s1/theta = {m_star[0]/m_star[3]:.2f}")
print()

for name, r in results.items():
    rr = results[name]
    n = rr['rank']
    state = rr['state']
    # Get equilibrium of first node
    m0 = state[:4]
    print(f"{name:4s} node 0: s1={m0[0]:.4f}  s2={m0[1]:.4f}  theta={m0[3]:.4f}  "
          f"s1/s2={m0[0]/m0[1]:.2f}  Born={born_prob(m0):.4f}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("Key results from the ADE survey:")
print()
print("1. Spectral radius ordering:")
an_rho = [(int(n[1:]), results[n]['rho']) for n in results if n.startswith('A')]
an_rho.sort()
for n, rho in an_rho:
    print(f"   A{n}: rho = {rho:.4f}  Delta = {-np.log(rho):.4f}")

print()
for name in ['D4','D5','D6','D7','E6','E7','E8']:
    if name in results:
        rr = results[name]
        print(f"   {name}: rho = {rr['rho']:.4f}  Delta = {rr['D0']:.4f}")

print()
print("2. E8 special values:")
e8r = results.get('E8', {})
if e8r:
    print(f"   rho(E8) = {e8r['rho']:.6f}")
    print(f"   Delta(E8) = {e8r['D0']:.6f}")
    print(f"   h(E8) = 30,  Delta(E8) * 30 = {e8r['D0']*30:.4f}")

print()
print("3. The 1.082-type state:")
a2_ratio1 = results.get('A2', {}).get('ratios', [1,1])[1] if results.get('A2') else None
if a2_ratio1:
    print(f"   A2: {a2_ratio1:.4f}")
for name in ['A3','A4','A5','D4','D5','E6','E7','E8']:
    if name in results and len(results[name]['ratios']) >= 2:
        print(f"   {name}: {results[name]['ratios'][1]:.4f}")
