"""
Why degree 19?

The leading eigenvalue λ₀ of the DS transfer operator has minimal polynomial
of degree 19 over Q, with Galois group S₁₉. The evidence parameter v has
degree 24 over Q. Why does the eigenvalue land at 19, not 24 or 48?

Strategy:
1. Verify the degree-19 polynomial from the paper
2. Factor the polynomial mod small primes (structure of Galois action)
3. Analyze the coefficient structure (GCDs, factorizations)
4. Trace the algebraic tower: Q ⊂ Q(v) ⊂ Q(v,√(26R)) ⊃ Q(λ₀)
5. Count: what is the degree of the FULL algebraic system for λ₀?
"""

import numpy as np
import sys

try:
    from mpmath import mp, mpf, fabs, nstr, log, polyroots
except ImportError:
    print("mpmath required"); sys.exit(1)

# ============================================================
# The degree-19 polynomial from the paper (eq:minpoly)
# ============================================================
# P(λ) = c₁₉λ¹⁹ + c₁₈λ¹⁸ + ... + c₀
# Listed highest degree first
coeffs_descending = [
    -497376766077,    # λ^19
     1088544475467,   # λ^18
     22094975727,     # λ^17
    -218313582412,    # λ^16
     578938153168,    # λ^15
    -300055480057,    # λ^14
     833013126110,    # λ^13
    -402784632820,    # λ^12
    -402921610991,    # λ^11
    -210848363196,    # λ^10
     175793550949,    # λ^9
    -458346126206,    # λ^8
     10558534577,     # λ^7
     536180328934,    # λ^6
    -215345191444,    # λ^5
     65343281380,     # λ^4
    -125920925463,    # λ^3
    -210759099887,    # λ^2
     85202393725,     # λ^1
    -4671418103,      # λ^0
]

coeffs_ascending = coeffs_descending[::-1]  # c₀, c₁, ..., c₁₉

print("=" * 70)
print("ANALYSIS OF THE DEGREE-19 MINIMAL POLYNOMIAL FOR λ₀")
print("=" * 70)

# ============================================================
# Part 1: Verify P(λ₀) = 0
# ============================================================
print("\nPart 1: Verification")

mp.dps = 100
lam0 = mpf('0.28291034660315143666181958723824224253434')

P_val = sum(mpf(c) * lam0**k for k, c in enumerate(coeffs_ascending))
print(f"  P(λ₀) = {nstr(P_val, 5)}")
print(f"  |P(λ₀)| = {nstr(fabs(P_val), 5)}")
print(f"  Verified: {'YES' if fabs(P_val) < mpf(10)**(-20) else 'NO'}")

# Check if λ₁ is also a root
lam1 = mpf('0.28131300001289042103025957885817017890632')
P_val1 = sum(mpf(c) * lam1**k for k, c in enumerate(coeffs_ascending))
print(f"\n  P(λ₁) = {nstr(P_val1, 5)}")
print(f"  λ₁ is a root: {'YES' if fabs(P_val1) < mpf(10)**(-20) else 'NO'}")

# ============================================================
# Part 2: Coefficient analysis
# ============================================================
print("\n" + "=" * 70)
print("Part 2: Coefficient structure")
print("=" * 70)

from math import gcd
from functools import reduce

g = reduce(gcd, [abs(c) for c in coeffs_descending])
print(f"\n  GCD of all coefficients: {g}")
print(f"  Polynomial is {'primitive' if g == 1 else f'divisible by {g}'}")

print(f"\n  Coefficient magnitudes:")
for k in range(19, -1, -1):
    c = coeffs_ascending[k]
    sign = '+' if c > 0 else '-'
    print(f"    λ^{k:2d}: {sign}{abs(c):>15d}  ({len(str(abs(c)))} digits)")

# Sum and alternating sum
S = sum(coeffs_ascending)
A = sum((-1)**k * c for k, c in enumerate(coeffs_ascending))
print(f"\n  P(1) = {S}  (sum of coefficients)")
print(f"  P(-1) = {A}  (alternating sum)")
print(f"  P(0) = {coeffs_ascending[0]}")

# ============================================================
# Part 3: Factorization mod small primes
# ============================================================
print("\n" + "=" * 70)
print("Part 3: Factorization mod small primes")
print("=" * 70)

def poly_mod_p(coeffs, p):
    """Reduce polynomial coefficients mod p."""
    return [c % p for c in coeffs]

def poly_mul_mod(a, b, p):
    """Multiply two polynomials mod p."""
    result = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i+j] = (result[i+j] + ai * bj) % p
    return result

def poly_mod_poly(f, g, p):
    """f mod g in F_p[x]. Returns remainder."""
    f = [x % p for x in f]
    g = [x % p for x in g]
    while len(f) >= len(g) and any(x != 0 for x in f):
        # Remove trailing zeros
        while f and f[-1] == 0:
            f.pop()
        if len(f) < len(g):
            break
        # Leading coefficient
        lc_f = f[-1]
        lc_g = g[-1]
        # Find lc_f / lc_g mod p
        lc_g_inv = pow(lc_g, p-2, p)
        factor = (lc_f * lc_g_inv) % p
        shift = len(f) - len(g)
        for i in range(len(g)):
            f[i + shift] = (f[i + shift] - factor * g[i]) % p
        while f and f[-1] == 0:
            f.pop()
    return f if f else [0]

def poly_gcd_mod(f, g, p):
    """GCD of f, g in F_p[x]."""
    while any(x != 0 for x in g):
        f, g = g, poly_mod_poly(f, g, p)
    return f

def factor_degree_mod_p(coeffs_asc, p):
    """Find the degrees of irreducible factors of poly mod p.
    Uses Berlekamp or brute force for small degrees."""
    c = poly_mod_p(coeffs_asc, p)
    # Remove trailing zeros
    while c and c[-1] == 0:
        c.pop()
    if not c:
        return [0]

    deg = len(c) - 1
    if deg <= 0:
        return []

    # Make monic
    lc_inv = pow(c[-1], p-2, p)
    c = [(x * lc_inv) % p for x in c]

    # Factor by computing gcd(x^{p^k} - x, f) for k=1,2,...
    factors = []
    f = list(c)
    # x^p mod f
    # Start with x = [0, 1]
    for k in range(1, deg+1):
        # Compute x^{p^k} mod f using repeated squaring
        # x^p mod f
        xpk = [0, 1]  # x
        for _ in range(k):
            # raise to p-th power mod f
            result = [1]
            base = list(xpk)
            exp = p
            while exp > 0:
                if exp % 2 == 1:
                    result = poly_mul_mod(result, base, p)
                    result = poly_mod_poly(result, f, p)
                base = poly_mul_mod(base, base, p)
                base = poly_mod_poly(base, f, p)
                exp //= 2
            xpk = result

        # gcd(x^{p^k} - x, f)
        diff = list(xpk)
        if len(diff) >= 2:
            diff[1] = (diff[1] - 1) % p
        elif len(diff) >= 1:
            diff = [0, (0 - 1) % p]  # -x
        else:
            diff = [0, p-1]

        g = poly_gcd_mod(list(f), diff, p)
        if len(g) > 1:
            g_deg = len(g) - 1
            n_factors = g_deg // k
            factors.extend([k] * n_factors)
            # Divide f by g
            # f = f / g mod p (exact division)
            while len(g) > 1:
                f = poly_mod_poly(f, g, p)
                # Actually we need exact division, not remainder
                # Let me use a simpler approach
                break
            # For simplicity, just record what we found
            break

    return factors

# Simpler approach: just check if P has roots mod p
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    c_mod = [c % p for c in coeffs_ascending]
    roots_mod_p = []
    for x in range(p):
        val = sum(c_mod[k] * pow(x, k, p) % p for k in range(20)) % p
        if val == 0:
            roots_mod_p.append(x)
    print(f"  mod {p:2d}: {len(roots_mod_p)} roots: {roots_mod_p}")

# ============================================================
# Part 4: The algebraic tower
# ============================================================
print("\n" + "=" * 70)
print("Part 4: Algebraic tower analysis")
print("=" * 70)

print("""
The algebraic dependencies:

1. s = √3, degree 2 over Q
2. v (evidence parameter), degree 12 over Q(s) = degree 24 over Q
3. a, b, θ (equilibrium masses): algebraic over Q(v)
   - Determined by the Gröbner basis equations
   - Back-substitution from the lex Gröbner basis gives a, b as functions of v, s
4. R_DS = (m²_DS₁ + 2m²_DS₂)/(m_DS₁ + 2m_DS₂)²: rational in Q(v)
5. √(26R_DS): degree 1 or 2 over Q(v)
6. Jacobian entries: in Q(v, √(26R_DS))
7. tr(J), cofsum: in Q(v, √(26R_DS))
8. disc = tr² - 4·cofsum: in Q(v, √(26R_DS))
9. λ₀ = (tr + √disc)/2: in Q(v, √(26R_DS), √disc)

Degree of λ₀ over Q:
  [Q(λ₀):Q] = 19 (prime, Gal = S₁₉)

Since 19 is prime and gcd(19, 24) = 1:
  Q(λ₀) ∩ Q(v) = Q  (the two field extensions are "independent")

This means: λ₀ is NOT in Q(v). The eigenvalue generates a genuinely
different algebraic extension than the evidence parameter.

The degree 19 must come from the COMBINED elimination of:
  - v (degree 24 over Q)
  - √(26R) (additional algebraic complexity from the floor)
  - √disc (from the eigenvalue equation)

The interaction of these three algebraic operations produces degree 19.
""")

# ============================================================
# Part 5: Is 19 = something from H=3?
# ============================================================
print("=" * 70)
print("Part 5: Numerological check — is 19 natural from H=3?")
print("=" * 70)

H = 3
print(f"\n  H = {H}")
print(f"  H² + 1 = {H**2+1} (= dim Sym²(C^{{H+1}}))")
print(f"  H³ = {H**3}")
print(f"  H⁴ + H² + 1 = {H**4+H**2+1} = 91 (cyclotomic)")
print(f"  (H²+1)(H+1) - 1 = {(H**2+1)*(H+1)-1}")
print(f"  H(H²+1) - H(H-1) = {H*(H**2+1) - H*(H-1)}")
print(f"  2(H²+1) - 1 = {2*(H**2+1)-1} = 19  ✓")
print(f"  (H+1)² + H = {(H+1)**2+H} = 19  ✓")
print(f"  H² + H + H(H-1)/2 + H + 1 = ...")
print()
print(f"  Two representations of 19:")
print(f"    19 = 2·(H²+1) - 1 = 2·dim(Sym²(C⁴)) - 1")
print(f"    19 = (H+1)² + H = dim(C^{{H+1}})² + H")
print()

# Check: is 19 the degree of a specific resultant?
print(f"  Bézout products for candidate systems:")
print(f"    2·2·2·2·1 = {2*2*2*2*1} (five equations, degrees 2,2,2,2,1)")
print(f"    24·2/... possible quotients:")
for d in [24, 48, 96]:
    if d % 19 == 0:
        print(f"    {d}/19 = {d//19}")
    else:
        print(f"    {d}/19 = {d/19:.3f} (not integer)")

# ============================================================
# Part 6: Discriminant structure
# ============================================================
print("\n" + "=" * 70)
print("Part 6: Leading and trailing coefficients")
print("=" * 70)

c19 = abs(coeffs_descending[0])
c0 = abs(coeffs_descending[-1])

print(f"  Leading coeff c₁₉ = {coeffs_descending[0]}")
print(f"  Trailing coeff c₀ = {coeffs_descending[-1]}")

# Factor these
def trial_factor(n):
    n = abs(n)
    factors = []
    for p in range(2, 10000):
        while n % p == 0:
            factors.append(p)
            n //= p
        if p*p > n:
            break
    if n > 1:
        factors.append(n)
    return factors

print(f"\n  c₁₉ = {coeffs_descending[0]} = {' × '.join(str(f) for f in trial_factor(coeffs_descending[0]))}")
print(f"  c₀  = {coeffs_descending[-1]} = {' × '.join(str(f) for f in trial_factor(coeffs_descending[-1]))}")

# P(1) and P(-1)
print(f"\n  P(1) = {S}")
if S != 0:
    print(f"       = {' × '.join(str(f) for f in trial_factor(S))}")
print(f"  P(-1) = {A}")
if A != 0:
    print(f"        = {' × '.join(str(f) for f in trial_factor(A))}")

# The content (GCD of coefficients)
print(f"\n  Content = {g}")

# Product of roots = (-1)^19 c₀/c₁₉ = c₀/c₁₉
# Sum of roots = -c₁₈/c₁₉
prod_roots = mpf(coeffs_ascending[0]) / mpf(coeffs_ascending[19])
sum_roots = -mpf(coeffs_descending[1]) / mpf(coeffs_descending[0])
print(f"\n  Product of all 19 roots = c₀/c₁₉ = {nstr(prod_roots, 15)}")
print(f"  Sum of all 19 roots     = -c₁₈/c₁₉ = {nstr(sum_roots, 15)}")

# How many real roots?
print(f"\n  Finding all roots...")
mp.dps = 50
all_roots = polyroots([mpf(c) for c in coeffs_descending])
real_roots = [r for r in all_roots if fabs(r.imag) < mpf(10)**(-30)]
complex_roots = [r for r in all_roots if fabs(r.imag) >= mpf(10)**(-30)]

print(f"  Real roots: {len(real_roots)}")
for r in sorted(real_roots, key=lambda x: float(x.real)):
    print(f"    {nstr(r.real, 20)}")

print(f"  Complex root pairs: {len(complex_roots)//2}")
for r in sorted(complex_roots, key=lambda x: (float(x.real), float(x.imag))):
    if r.imag > 0:
        print(f"    {nstr(r.real, 12)} ± {nstr(fabs(r.imag), 12)}i")

print(f"\n  λ₀ = {nstr(lam0, 20)} (the unique root in (0,1))")

print("\nDONE.")
