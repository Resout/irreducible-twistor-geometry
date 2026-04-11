#!/usr/bin/env python3
"""
Investigate the exact structure of Λ from framework quantities alone.

Don't anchor to measured Λ. Instead: what dimensionless number does
the framework produce at the fourth rung? Work from the inside out.

The three known quantities:
  Floor = 1/H³
  Born_DS = determined by equilibrium
  Instanton = det(I-J)^{-1/2} × e^{-H³/K*}

The question: what is the fourth rung, expressed purely in terms
of H, K*, and the equilibrium?
"""

from mpmath import mp, mpf, sqrt, ln, exp, fabs, nstr, pi, log, power

mp.dps = 80
H = 3
K = mpf(7)/30
FLOOR = mpf(1)/H**3
S = mpf(H)**3 / K  # = 810/7

# Equilibrium values (high precision)
# These are determined by 4 equations in 4 unknowns
# Let's get them fresh

def ds_combine(m, e):
    s = m[:3]; theta = m[3]; ev = e[:3]; phi = e[3]
    s_pre = [s[i]*ev[i] + s[i]*phi + theta*ev[i] for i in range(3)]
    theta_pre = theta * phi
    total_pre = sum(s_pre) + theta_pre
    Kv = mpf(1) - total_pre
    if fabs(mpf(1) - Kv) < mpf(10)**(-70): return list(m), Kv
    denom = mpf(1) - Kv
    return [sp / denom for sp in s_pre] + [theta_pre / denom], Kv

def born_prob(m):
    L2sq = sum(x**2 for x in m)
    if L2sq < mpf(10)**(-60): return mpf(0)
    return m[3]**2 / L2sq

def enforce_floor(m):
    b = born_prob(m)
    if b >= FLOOR - mpf(10)**(-60): return list(m)
    Sv = sum(m[:3])
    if Sv < mpf(10)**(-60): return list(m)
    Sq = sum(x**2 for x in m[:3])
    A_c = mpf(26)*Sv**2 - Sq; B_c = mpf(2)*Sq; C_c = -Sq
    disc = B_c**2 - 4*A_c*C_c
    if disc < 0: return list(m)
    t1 = (-B_c + sqrt(disc))/(2*A_c)
    t2 = (-B_c - sqrt(disc))/(2*A_c)
    cands = [t for t in [t1,t2] if mpf(0)<t<mpf(1)]
    if not cands: return list(m)
    t = min(cands, key=lambda x: fabs(x-m[3]))
    alpha = (mpf(1)-t)/Sv
    return [m[i]*alpha for i in range(3)] + [t]

def ds_step(m, e):
    m_ds, Kv = ds_combine(m, e)
    return enforce_floor(m_ds), Kv

def find_eq():
    from scipy.optimize import fsolve
    from mpmath import findroot
    def eq_r(p):
        s1,th,w1,ph = p
        s2=(1-s1-th)/2; w2=(1-w1-ph)/2
        m=[s1,s2,s2,th]; e=[w1,w2,w2,ph]
        L2m=s1**2+2*s2**2+th**2; L2e=w1**2+2*w2**2+ph**2
        Kv=s1*w2+s1*w2+s2*w1+s2*w2+s2*w1+s2*w2
        m_mp=[mpf(x) for x in m]; e_mp=[mpf(x) for x in e]
        mo,_=ds_step(m_mp,e_mp)
        return[th**2/L2m-1/27,ph**2/L2e-1/27,Kv-7/30,float(mo[0])-s1]
    sol=fsolve(eq_r,[0.787,0.155,0.631,0.129])
    def eq_mp(*p):
        s1,th,w1,ph=p
        s2=(mpf(1)-s1-th)/2; w2=(mpf(1)-w1-ph)/2
        m=[s1,s2,s2,th]; e=[w1,w2,w2,ph]
        eq1=th**2/(s1**2+2*s2**2+th**2)-FLOOR
        eq2=ph**2/(w1**2+2*w2**2+ph**2)-FLOOR
        Kv=s1*w2+s1*w2+s2*w1+s2*w2+s2*w1+s2*w2
        mo,_=ds_step(m,e)
        return[eq1,eq2,Kv-K,mo[0]-s1]
    x=findroot(eq_mp,[mpf(str(v)) for v in sol])
    s1,th,w1,ph=x
    s2=(mpf(1)-s1-th)/2; w2=(mpf(1)-w1-ph)/2
    return [s1,s2,s2,th],[w1,w2,w2,ph]

print("Computing equilibrium...")
m_star, e_star = find_eq()
s1,s2,s3,theta = m_star
e1,e2,e3,phi = e_star

# DS output
m_ds, K_val = ds_combine(m_star, e_star)
Born_DS = born_prob(m_ds)
theta_ds = m_ds[3]

# Single-site Jacobian eigenvalues
eps = mpf(10)**(-35)
J4 = [[mpf(0)]*4 for _ in range(4)]
for j in range(4):
    sp=list(m_star); sm=list(m_star)
    sp[j]+=eps; sm[j]-=eps
    fp,_=ds_step(sp,e_star); fm,_=ds_step(sm,e_star)
    for i in range(4):
        J4[i][j]=(fp[i]-fm[i])/(2*eps)

from mpmath import matrix, eig, chop
J4m = matrix(4,4)
for i in range(4):
    for j in range(4):
        J4m[i,j] = J4[i][j]
evals,_ = eig(J4m)
evals = sorted([chop(e) for e in evals], key=lambda x: -float(fabs(x)))
lam0 = fabs(evals[0])  # ~0.2829
lam1 = fabs(evals[1])  # ~0.2813

# Instanton
det_IJ = (1-lam0)*(1-lam1)
prefactor = det_IJ**(-mpf(1)/2)
instanton = prefactor * exp(-S)

print(f"\nEquilibrium:")
print(f"  θ* = {nstr(theta, 20)}")
print(f"  φ* = {nstr(phi, 20)}")
print(f"  Born_DS = {nstr(Born_DS, 20)}")
print(f"  λ₀ = {nstr(lam0, 20)}")
print(f"  λ₁ = {nstr(lam1, 20)}")
print(f"  S = {nstr(S, 20)}")
print(f"  instanton = {nstr(instanton, 20)}")

# ============================================================
# APPROACH 1: What exact expressions give clean Λ?
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 1: EXACT EXPRESSIONS FROM FRAMEWORK QUANTITIES")
print("=" * 70)

# All the building blocks
blocks = {
    "S": S,
    "K*": K,
    "1-K*": 1-K,
    "H": mpf(H),
    "H-1": mpf(H-1),
    "H+1": mpf(H+1),
    "H²": mpf(H**2),
    "H³": mpf(H**3),
    "H²+1": mpf(H**2+1),
    "H²-1": mpf(H**2-1),
    "H(H+1)": mpf(H*(H+1)),
    "5!": mpf(120),
    "θ*": theta,
    "φ*": phi,
    "Born_DS": Born_DS,
    "λ₀": lam0,
    "λ₁": lam1,
    "e^{-S}": exp(-S),
    "inst": instanton,
}

# The Λ formula: inst × Born_DS^n with n ≈ 23.55
# What exactly IS n?

n_exact = (ln(mpf("1.15e-123")) - ln(instanton)) / ln(Born_DS)
print(f"\n  Exact n (from Λ_obs = 1.15e-123): {nstr(n_exact, 20)}")

# Test: is n related to S or other quantities?
print(f"\n  n vs framework:")
print(f"    47/2 = {23.5}  err = {float(fabs(n_exact-23.5)/n_exact*100):.3f}%")
print(f"    S/5 = {float(S/5):.4f}  err = {float(fabs(n_exact-S/5)/n_exact*100):.3f}%")
print(f"    (S-K*)/5 = {float((S-K)/5):.4f}")
print(f"    H²(H-1)+H/2 = {H**2*(H-1)+H/2:.1f}")

# ============================================================
# APPROACH 2: Build Λ from STRUCTURE, not measurement
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 2: THE FOURTH RUNG FROM STRUCTURE")
print("=" * 70)

# The ladder in ln-space:
# rung 0: ln(floor) = -ln(H³) = -3ln(H)
# rung 1: ln(Born_DS) = some function of equilibrium
# rung 2: ln(inst) ≈ -S + small correction = -H³/K* + ...
# rung 3: ln(Λ) = ???

# The gaps between rungs:
g01 = ln(Born_DS) - ln(FLOOR)    # floor → depth
g12 = ln(instanton) - ln(Born_DS)  # depth → instanton
g23_obs = ln(mpf("1.15e-123")) - ln(instanton)  # instanton → Λ (observed)

print(f"\n  Gap 0→1 (floor→depth): {nstr(g01, 15)}")
print(f"  Gap 1→2 (depth→inst):  {nstr(g12, 15)}")
print(f"  Gap 2→3 (inst→Λ obs):  {nstr(g23_obs, 15)}")

# Ratios of gaps
r12 = g12/g01
r23 = g23_obs/g12
r_total = g23_obs/g01

print(f"\n  Gap ratios:")
print(f"    g12/g01 = {nstr(r12, 15)}")
print(f"    g23/g12 = {nstr(r23, 15)}")
print(f"    g23/g01 = {nstr(r_total, 15)}")

# Test: g23/g12 vs framework numbers
print(f"\n  g23/g12 = {nstr(r23, 10)} vs:")
r23_tests = {
    "3/2": mpf(3)/2,
    "H/(H-1)": mpf(H)/(H-1),
    "S/(S-H³+H)": S/(S-H**3+H),
    "(H²+1)/(H²-1)": mpf(H**2+1)/(H**2-1),
    "47/H(H+1)": mpf(47)/(H*(H+1)),
    "(H+1)/H": mpf(H+1)/H,
    "1+K*": 1+K,
    "1+1/H": 1+mpf(1)/H,
    "S/(H³+H)": S/(H**3+H),
    "(1-K*)/K* × K*...": (1-K)/K * K,  # = 1-K = 23/30
}
for name, val in r23_tests.items():
    err = float(fabs(r23-val)/r23*100)
    if err < 5:
        print(f"    {name:30s} = {nstr(val, 10):>12s}  err={err:.3f}%")

# Test: g12/g01 vs framework
print(f"\n  g12/g01 = {nstr(r12, 10)} vs:")
r12_tests = {
    "(H+1)²": mpf((H+1)**2),
    "H(H+2)": mpf(H*(H+2)),
    "S/H": S/H,
    "S/ln(H³)": S/ln(mpf(H**3)),
    "H⁴/K*/(something)": mpf(1),
}
for name, val in r12_tests.items():
    err = float(fabs(r12-val)/fabs(r12)*100)
    if err < 10:
        print(f"    {name:30s} = {nstr(val, 10):>12s}  err={err:.3f}%")

# ============================================================
# APPROACH 3: Λ from the instanton action at multiple depths
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 3: MULTIPLE INSTANTON ACTIONS")
print("=" * 70)

# The hierarchy tower uses e^{-S/d} at different depths d.
# The instanton is e^{-S/1} = e^{-S} (d=1, trivial sector).
# What if Λ involves the PRODUCT of multiple sectors?

# Product of all hierarchy entries:
depths = [1, 5, 10, 16, 24, 26]
product_all = sum(S/d for d in depths)
print(f"\n  Sum of S/d for d in {depths}:")
print(f"  = S × Σ(1/d) = {nstr(S, 6)} × {float(sum(1/d for d in depths)):.6f} = {nstr(product_all, 10)}")
print(f"  e^{{-sum}} = 10^{float(-product_all/ln(10)):.2f}")

# What about d=1 and d=26 together?
sum_1_26 = S/1 + S/26
print(f"\n  S/1 + S/26 = S(1+1/26) = S × 27/26 = {nstr(sum_1_26, 10)}")
print(f"  e^{{-(S+S/26)}} = 10^{float(-sum_1_26/ln(10)):.2f}")

# S + S/H = S(H+1)/H
sum_S_SH = S * (H+1)/H
print(f"\n  S(H+1)/H = {nstr(sum_S_SH, 10)}")
print(f"  e^{{-S(H+1)/H}} = 10^{float(-sum_S_SH/ln(10)):.2f}")

# S × (1 + 1/d) for various d
for d in [2, 3, 5, 10, 26]:
    total_action = S * (1 + mpf(1)/d)
    l10 = float(-total_action/ln(10))
    print(f"  S(1+1/{d}) = S×{float(1+1/d):.4f} = {float(total_action):.2f} → 10^{l10:.2f}")

# ============================================================
# APPROACH 4: Λ = e^{-S} × e^{-S/d_eff} = e^{-S(1+1/d)}
# with d chosen to give exact Λ
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 4: DOUBLE INSTANTON e^{-S(1+1/d)}")
print("=" * 70)

# If Λ = e^{-S(1+1/d)}, what d gives the observed value?
# ln(Λ) = -S(1+1/d)
# -283.08 = -115.71(1+1/d)
# 1+1/d = 283.08/115.71 = 2.446
# 1/d = 1.446
# d = 0.6916

# But what if we include the prefactor?
# ln(Λ) = ln(prefactor) - S(1+1/d)
# -283.08 = 0.331 - S(1+1/d)
# S(1+1/d) = 283.41
# 1+1/d = 283.41/115.71 = 2.4494
# 1/d = 1.4494
# d = 0.6899

ln_Lambda_obs = ln(mpf("1.15e-123"))
ln_pref = ln(prefactor)

ratio_needed = (-ln_Lambda_obs + ln_pref) / S
d_eff = 1 / (ratio_needed - 1)

print(f"\n  If Λ = prefactor × e^{{-S(1+1/d)}}:")
print(f"  ratio = {nstr(ratio_needed, 15)}")
print(f"  d = {nstr(d_eff, 15)}")

# Is d_eff a framework number?
d_tests = {
    "K*": K,
    "K*/(1-K*)": K/(1-K),
    "1/√2": mpf(1)/sqrt(2),
    "ln(H)/H": ln(mpf(H))/H,
    "(H-1)/H": mpf(H-1)/H,
    "K*×H": K*H,
    "7/10": mpf(7)/10,
    "K*×3": K*3,
    "sin²θ_W = 3/8": mpf(3)/8,
}
for name, val in d_tests.items():
    err = float(fabs(d_eff-val)/d_eff*100)
    if err < 15:
        print(f"    d ≈ {name:25s} = {nstr(val, 10):>12s}  err={err:.3f}%")

# ============================================================
# APPROACH 5: Pure algebraic — what if Λ = e^{-5!/2}?
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 5: PURE ALGEBRAIC CANDIDATES")
print("=" * 70)

# S = 810/7, and S + 1/K* = 120 = 5!
# So 5!/2 = 60, S+S/d, etc.

candidates = {
    "e^{-5!/2} = e^{-60}": exp(-mpf(60)),
    "e^{-5!} = e^{-120}": exp(-mpf(120)),
    "e^{-5!} × prefactor": prefactor * exp(-mpf(120)),
    "e^{-2S} × prefactor²": prefactor**2 * exp(-2*S),
    "e^{-2S} × prefactor": prefactor * exp(-2*S),
    "e^{-2S+ln(H)} = H×e^{-2S}": mpf(H) * exp(-2*S),
    "e^{-S} × e^{-5!+S} = e^{-5!}": exp(-mpf(120)),
    "inst × e^{-5!+S}": instanton * exp(-(120-S)),
    "inst × e^{-S×K*/(1-K*)}": instanton * exp(-S*K/(1-K)),
    "inst² / floor": instanton**2 / FLOOR,
    "inst² × H³": instanton**2 * H**3,
    "inst² × (H³-1)": instanton**2 * (H**3-1),
    "inst^{5/2} / (H-1)": instanton**(mpf(5)/2) / (H-1),
    "inst^{5/2} × H^{9/2} / (H-1)": instanton**(mpf(5)/2) * mpf(H)**(mpf(9)/2) / (H-1),
    "inst × Born_DS^{47/2}": instanton * Born_DS**(mpf(47)/2),
    "inst × Born_DS^{S/5}": instanton * Born_DS**(S/5),
    "inst × floor^{47/2}": instanton * FLOOR**(mpf(47)/2),
    "e^{-S} × (θ*φ*)^{S/H}": exp(-S) * (theta*phi)**(S/H),
    "e^{-S} × θ*^S": exp(-S) * theta**S,
    "inst × θ*^{S/(H-1)}": instanton * theta**(S/(H-1)),
    "inst × (θ*φ*/(1-K*))^S": instanton * (theta*phi/(1-K))**S,
    "inst × Born_DS^{(H³-1)/K*}": instanton * Born_DS**(mpf(H**3-1)/K),
}

print(f"\n  {'Formula':50s} {'log10':>10s} {'vs -122.9':>10s}")
print(f"  {'-'*50} {'-'*10} {'-'*10}")

results = []
for name, val in candidates.items():
    if float(fabs(val)) > 0:
        l10 = float(ln(fabs(val))/ln(10))
        err = abs(l10 - (-122.935))
        results.append((err, name, l10))

results.sort()
for err, name, l10 in results[:15]:
    print(f"  {name:50s} {l10:>10.2f} {err:>10.2f}")

# ============================================================
# APPROACH 6: The θ connection
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 6: θ* AND φ* AS THE BRIDGE")
print("=" * 70)

# θ_DS = θ*φ*/(1-K*) is the vacuum depth's θ component
# The Born_DS = θ_DS² / ||m_DS||²
# But what if Λ involves θ* or φ* directly?

print(f"\n  θ* = {nstr(theta, 30)}")
print(f"  φ* = {nstr(phi, 30)}")
print(f"  θ*φ* = {nstr(theta*phi, 30)}")
print(f"  θ*φ*/(1-K*) = θ_DS = {nstr(theta*phi/(1-K), 20)}")

# θ* raised to various powers
theta_tests = {
    "θ*^S = θ*^{810/7}": theta**S,
    "θ*^{5!}": theta**(mpf(120)),
    "θ*^{H⁴}": theta**(mpf(H**4)),
    "(θ*φ*)^{S/2}": (theta*phi)**(S/2),
    "(θ*φ*/(1-K*))^{S}": (theta*phi/(1-K))**S,
    "θ*^S × φ*^S": theta**S * phi**S,
    "θ*^{S+S/(H³-1)}": theta**(S+S/(H**3-1)),
}

print(f"\n  {'Expression':40s} {'log10':>12s}")
for name, val in theta_tests.items():
    if float(fabs(val)) > 0:
        l10 = float(ln(fabs(val))/ln(10))
        print(f"  {name:40s} {l10:>12.2f}")

# ============================================================
# KEY INSIGHT TEST: inst × Born_DS^n where n = S/|ln(Born_DS)| × K
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 7: THE DEPTH AS A MODULAR EXPONENT")
print("=" * 70)

# What if n = S × something / ln(1/Born_DS)?
# n ≈ 23.55
# ln(1/Born_DS) = 7.12
# S = 115.71
# S/ln(1/Born_DS) = 16.25
# n × ln(1/Born_DS) / S = 23.55 × 7.12 / 115.71 = 167.7/115.71 = 1.449

# So n = 1.449 × S / ln(1/Born_DS)
# And 1.449 ≈ 1 + 1/d with d ≈ 2.23 ≈ ln(10)/ln(3)?

depth_factor = float(n_exact * ln(1/Born_DS) / S)
print(f"\n  n × ln(1/Born_DS) / S = {depth_factor:.6f}")
print(f"  = 1 + {depth_factor-1:.6f}")

df_tests = {
    "3/2": 1.5,
    "(H+1)/H = 4/3": (H+1)/H,
    "1+K*": float(1+K),
    "1+1/(H-1)": float(1+1/(H-1)),
    "ln(27)/ln(10)": float(ln(27)/ln(10)),
    "1/K*×(1-K*)/H": float((1-K)/(K*H)),
    "H/(H-1)": float(H/(H-1)),
    "√2": float(sqrt(2)),
}
for name, val in df_tests.items():
    err = abs(depth_factor - float(val))/depth_factor*100
    if err < 10:
        print(f"    ≈ {name:30s} = {float(val):.6f}  err={err:.3f}%")


if __name__ == "__main__":
    pass
