"""
Full symbolic derivation of the spectral gap.

The six constraint equations determine a 6×6 Jacobian.
Its characteristic polynomial, evaluated at the solution,
has Δ as the modulus of a complex root.

Strategy:
1. Build symbolic Jacobian
2. Get characteristic polynomial (symbolic in variables + λ)
3. Eliminate variables using the six equations (Gröbner basis or sequential substitution)
4. Get a polynomial in λ alone with coefficients in K* and H
"""

import sympy as sp
from sympy import sqrt, Rational, symbols, solve, Matrix, simplify, factor
from sympy import expand, collect, Poly, groebner, resultant

s, w, th, a, b, phi, lam = symbols('s w theta a b phi lambda')

K_val = Rational(7, 30)

# ============================================================
# Step 1: The six equations and the Jacobian
# ============================================================
print("Step 1: Building Jacobian...")

eqs = [
    s + 2*w + th - 1,
    a + 2*b + phi - 1,
    26*th**2 - s**2 - 2*w**2,
    26*phi**2 - a**2 - 2*b**2,
    2*s*b + 2*w*(a + b) - K_val,
    b*s*(w + th) - a*w*(s + th),
]

F = Matrix(eqs)
J = F.jacobian([s, w, th, a, b, phi])
char_matrix = J - lam * sp.eye(6)
print("Computing characteristic polynomial...")
char_poly = char_matrix.det()
char_poly = expand(char_poly)
print(f"Characteristic polynomial computed. Degree in λ: {sp.degree(char_poly, lam)}")

# Collect by powers of λ
char_collected = Poly(char_poly, lam)
print(f"Coefficients (as expressions in s,w,θ,a,b,φ):")
for power, coeff in sorted(char_collected.as_dict().items(), reverse=True):
    deg = power[0]
    coeff_simplified = simplify(coeff)
    print(f"  λ^{deg}: {coeff_simplified}")

# ============================================================
# Step 2: Reduce using the linear equations first
# ============================================================
print("\nStep 2: Eliminating θ and φ using eqs 1,2...")

# θ = 1 - s - 2w, φ = 1 - a - 2b
subs_linear = {th: 1 - s - 2*w, phi: 1 - a - 2*b}
char_reduced = char_poly.subs(subs_linear)
char_reduced = expand(char_reduced)
print(f"After eliminating θ,φ: polynomial in s,w,a,b,λ")

# ============================================================
# Step 3: Use Born floor equations to eliminate w and b
# ============================================================
print("\nStep 3: Using Born floor equations...")

# eq3 with θ = 1-s-2w: 26(1-s-2w)² = s² + 2w²
# eq4 with φ = 1-a-2b: 26(1-a-2b)² = a² + 2b²

eq3_reduced = expand(26*(1 - s - 2*w)**2 - s**2 - 2*w**2)
eq4_reduced = expand(26*(1 - a - 2*b)**2 - a**2 - 2*b**2)
print(f"eq3: {eq3_reduced} = 0")
print(f"eq4: {eq4_reduced} = 0")

# Solve eq3 for w in terms of s
print("\nSolving eq3 for w...")
w_sols = solve(eq3_reduced, w)
print(f"w(s) = {w_sols}")

# Solve eq4 for b in terms of a
print("Solving eq4 for b...")
b_sols = solve(eq4_reduced, b)
print(f"b(a) = {b_sols}")

# Both give two branches. We need the physical ones (small w, small b).
# From the discriminant: 104(r²+2) where r=s/w, the solutions are
# w = (26s + 52 ± √(26s²+52)) / (25s² + 104s + 102)... let's just pick.

# The physical branches (verified numerically):
# w: the one with minus sign (smaller w)
# b: the one with minus sign (smaller b)

# Let me check which is which
import numpy as np
for i, ws in enumerate(w_sols):
    val = float(ws.subs(s, Rational(787, 1000)))
    print(f"  w_{i}(0.787) = {val:.6f}")

w_phys = w_sols[0] if float(w_sols[0].subs(s, Rational(787,1000))) < 0.1 else w_sols[1]

for i, bs in enumerate(b_sols):
    val = float(bs.subs(a, Rational(631, 1000)))
    print(f"  b_{i}(0.631) = {val:.6f}")

b_phys = b_sols[0] if float(b_sols[0].subs(a, Rational(631,1000))) < 0.2 else b_sols[1]

print(f"\nPhysical branches:")
print(f"  w(s) = {w_phys}")
print(f"  b(a) = {b_phys}")

# ============================================================
# Step 4: Substitute into the coupling equations
# ============================================================
print("\nStep 4: Reducing to (s, a) using Born floor...")

# eq5 reduced: 2s·b(a) + 2w(s)·(a + b(a)) = 7/30
eq5_reduced = 2*s*b_phys + 2*w_phys*(a + b_phys) - K_val
eq5_reduced = simplify(eq5_reduced)
print(f"eq5(s,a) = {eq5_reduced}")

# eq6 reduced: b(a)·s·(w(s)+θ(s)) - a·w(s)·(s+θ(s))
th_of_s = 1 - s - 2*w_phys
eq6_reduced = b_phys*s*(w_phys + th_of_s) - a*w_phys*(s + th_of_s)
eq6_reduced = simplify(eq6_reduced)
print(f"eq6(s,a) = {eq6_reduced}")

# ============================================================
# Step 5: Substitute w(s) and b(a) into the char polynomial
# ============================================================
print("\nStep 5: Substituting w(s) and b(a) into characteristic polynomial...")

subs_floor = {w: w_phys, th: 1 - s - 2*w_phys, phi: 1 - a - 2*b_phys}
# Also substitute b → b(a)
subs_floor[b] = b_phys

char_in_sa = char_poly.subs(subs_linear)  # first eliminate θ, φ
print("  Eliminated θ, φ...")

# Now substitute w(s) and b(a)
char_in_sa = char_in_sa.subs(w, w_phys)
print("  Substituted w(s)...")
char_in_sa = char_in_sa.subs(b, b_phys)
print("  Substituted b(a)...")

# This is now a function of s, a, and λ (with square roots from w(s) and b(a))
# To clear the square roots, we might need to rationalize.

# Let's first try evaluating at the numerical solution to verify
s_num = Rational(786898446158, 10**12)
a_num = Rational(631200887903, 10**12)

print("\nEvaluating at numerical solution to verify...")
char_numerical = char_in_sa.subs(s, s_num).subs(a, a_num)
# This should be a polynomial in λ
char_numerical = simplify(char_numerical)
print("Simplified.")

# Get coefficients
char_poly_lam = Poly(expand(char_numerical), lam)
print(f"Coefficients of λ:")
for power, coeff in sorted(char_poly_lam.as_dict().items(), reverse=True):
    deg = power[0]
    print(f"  λ^{deg}: {float(coeff):.8f}")

# Roots
import numpy as np
coeffs_list = [float(c) for c in char_poly_lam.all_coeffs()]
roots = np.roots(coeffs_list)
print(f"\nRoots:")
for r in sorted(roots, key=lambda x: -abs(x)):
    print(f"  {r}")
print(f"\nModuli: {sorted(np.abs(roots), reverse=True)}")

# ============================================================
# Step 6: Try to get the polynomial in λ with exact K*, H coefficients
# ============================================================
print("\n" + "=" * 60)
print("Step 6: Exact elimination")
print("=" * 60)

# The hard part: eliminate s and a from the system
# {char_poly(s,a,λ)=0, eq5(s,a)=0, eq6(s,a)=0}
# to get a polynomial in λ alone.
#
# This requires resultants or Gröbner bases, and the square roots
# from w(s) and b(a) complicate things.
#
# Alternative: work with the squared quantities to avoid radicals.
# Define u = s², v = a². Then the Born floor equations become
# polynomial in u and w² (or use the quadratic solutions directly).

# Actually, let me try a different approach: use the NUMERICAL
# characteristic polynomial but recognize its coefficients as
# rational functions of K* = 7/30.

print("Characteristic polynomial coefficients at the solution:")
print("Looking for structure in terms of K* = 7/30...")

K = 7/30
H = 3

for i, c in enumerate(coeffs_list):
    deg = 6 - i
    # Try to express c in terms of K and H
    # c = p(K, H) / q(K, H) for small integer polynomials
    print(f"\n  λ^{deg}: c = {c:.10f}")
    # Check simple forms
    candidates = {}
    for n1 in range(-10, 11):
        for d1 in range(1, 31):
            for n2 in range(-5, 6):
                val = n1/d1 + n2*K
                if abs(val - c) < 0.001:
                    candidates[f"{n1}/{d1} + {n2}·K"] = val
                val2 = n1/d1 * K**abs(n2)
                if abs(val2 - c) < 0.001:
                    candidates[f"{n1}/{d1}·K^{n2}"] = val2
    if candidates:
        best = min(candidates.items(), key=lambda x: abs(x[1]-c))
        print(f"    Best match: {best[0]} = {best[1]:.10f} (diff = {best[1]-c:.2e})")
    else:
        print(f"    No simple match found")
