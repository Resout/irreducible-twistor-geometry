"""
Can the crystal detect number-theoretic structure beyond periodicity?

The θ-fingerprint detects invertibility (5.44x separation). What
OTHER mathematical properties does it detect?

Test: prime vs composite, multiplicative vs additive, quadratic
residues vs non-residues. If the crystal can detect these, the
solver gains a STRUCTURAL ORACLE for mathematical problems.

Also: the crystal as a Legendre symbol detector.
For p prime, the Legendre symbol (n/p) = {0, 1, -1} has exactly
3 values — perfectly matched to H=3. Can the crystal learn it?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
from solver.algebra import H, MASS_DIM, schmidt_number
from solver.crystals import Entangler, classify_relationship, slot_measure
from solver.functional import build_function_crystal, _detect_output_period
from solver import compute

torch.set_grad_enabled(False)


def theta_fingerprint(joint):
    fp = torch.stack([
        joint[0, 3], joint[1, 3], joint[2, 3],
        joint[3, 0], joint[3, 1], joint[3, 2],
    ])
    norm = fp.abs().pow(2).sum().sqrt()
    if norm > 1e-10:
        fp = fp / norm
    return fp

def theta_asymmetry(fp):
    return (fp[:3].abs() - fp[3:].abs()).abs().sum().item()

def theta_distance(fp1, fp2):
    return (fp1 - fp2).abs().pow(2).sum().sqrt().item()


domain = range(1, 301)

# Build canonical fingerprint for reference
fc_id = build_function_crystal("identity", lambda x: int(x) % 3, domain)
fp_id = theta_fingerprint(fc_id.joint)


# ═══════════════════════════════════════════════════════════════
#  Test 1: Legendre symbol as a 3-valued function
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 1: THE LEGENDRE SYMBOL — a natural H=3 function")
print("(n/p) ∈ {-1, 0, 1} — exactly 3 values for H=3!")
print("=" * 70)

def legendre(n, p):
    """Legendre symbol (n/p). Returns -1, 0, or 1."""
    if n % p == 0:
        return 0
    result = pow(n, (p - 1) // 2, p)
    return result if result <= 1 else result - p

# Test for several primes
for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
    f = lambda x, _p=p: legendre(int(x), _p)
    fc = build_function_crystal(f"(n/{p})", f, domain)
    fp = theta_fingerprint(fc.joint)
    d_id = theta_distance(fp, fp_id)
    asymm = theta_asymmetry(fp)
    outputs = [f(x) for x in range(1, 51)]
    period = _detect_output_period(outputs)
    n_values = len(set(outputs))

    print(f"  (n/{p:2d}): Schmidt={fc.schmidt:.3f}, rel={fc.relationship:13s}, "
          f"asymm={asymm:.4f}, d(id)={d_id:.4f}, distinct={n_values}, period={period}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: Primality detector
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 2: PRIMALITY — is_prime(n)")
print("Binary function (2 values in 3 bins)")
print("=" * 70)

# is_prime maps to {0, 1} — only 2 distinct values
f_prime = lambda n: 1 if compute.is_prime(max(2, int(n))) else 0
fc_prime = build_function_crystal("is_prime", f_prime, domain)
fp_prime = theta_fingerprint(fc_prime.joint)
print(f"  is_prime: Schmidt={fc_prime.schmidt:.3f}, rel={fc_prime.relationship}, "
      f"asymm={theta_asymmetry(fp_prime):.4f}")

# Primes mod 3: residues of primes
f_pmod3 = lambda n: max(2, int(n)) % 3 if compute.is_prime(max(2, int(n))) else -1
# Better: count primes up to n
def prime_counting(n):
    n = max(1, int(n))
    return sum(1 for i in range(2, n + 1) if compute.is_prime(i))

fc_pi = build_function_crystal("π(n)", prime_counting, range(1, 101))
print(f"  π(n):     Schmidt={fc_pi.schmidt:.3f}, rel={fc_pi.relationship}, "
      f"asymm={theta_asymmetry(theta_fingerprint(fc_pi.joint)):.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 3: Multiplicative vs additive functions
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 3: MULTIPLICATIVE vs ADDITIVE")
print("Multiplicative: f(ab) = f(a)f(b) for gcd(a,b)=1")
print("Additive: f(ab) = f(a) + f(b) for gcd(a,b)=1")
print("=" * 70)

# Multiplicative functions
mult_crystals = []
for name, f in [
    ("φ(n)", lambda n: compute.euler_totient(max(1, int(n)))),
    ("σ(n)", lambda n: compute.divisor_sum(max(1, int(n)))),
    ("d(n)", lambda n: compute.divisor_count(max(1, int(n)))),
    ("n mod 3", lambda n: int(n) % 3),
    ("μ(n)²", lambda n: 0 if any(v > 1 for v in compute.factorize(max(2, int(n))).values()) else 1),
]:
    fc = build_function_crystal(name, f, domain)
    fp = theta_fingerprint(fc.joint)
    mult_crystals.append({
        "name": name, "crystal": fc, "theta_fp": fp,
        "asymm": theta_asymmetry(fp),
        "d_id": theta_distance(fp, fp_id),
    })
    print(f"  {name:<12s}: Schmidt={fc.schmidt:.3f}, rel={fc.relationship:13s}, "
          f"asymm={theta_asymmetry(fp):.4f}, d(id)={theta_distance(fp, fp_id):.4f}")

# Additive functions
print()
add_crystals = []
for name, f in [
    ("ω(n)", lambda n: len(compute.factorize(max(2, int(n))))),
    ("Ω(n)", lambda n: sum(compute.factorize(max(2, int(n))).values())),
    ("digit_sum", lambda n: compute.digit_sum(max(1, int(n)))),
    ("log₂(n)", lambda n: int(n).bit_length() - 1),
]:
    fc = build_function_crystal(name, f, domain)
    fp = theta_fingerprint(fc.joint)
    add_crystals.append({
        "name": name, "crystal": fc, "theta_fp": fp,
        "asymm": theta_asymmetry(fp),
        "d_id": theta_distance(fp, fp_id),
    })
    print(f"  {name:<12s}: Schmidt={fc.schmidt:.3f}, rel={fc.relationship:13s}, "
          f"asymm={theta_asymmetry(fp):.4f}, d(id)={theta_distance(fp, fp_id):.4f}")

# Compare
avg_mult_asymm = sum(c["asymm"] for c in mult_crystals) / len(mult_crystals)
avg_add_asymm = sum(c["asymm"] for c in add_crystals) / len(add_crystals)
avg_mult_d = sum(c["d_id"] for c in mult_crystals) / len(mult_crystals)
avg_add_d = sum(c["d_id"] for c in add_crystals) / len(add_crystals)
print(f"\n  Multiplicative: avg asymm={avg_mult_asymm:.4f}, avg d(id)={avg_mult_d:.4f}")
print(f"  Additive:       avg asymm={avg_add_asymm:.4f}, avg d(id)={avg_add_d:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 4: Quadratic residues — period (p-1)/2 structure
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 4: QUADRATIC RESIDUES")
print("x² mod p has period (p-1)/2 — multiplicative structure")
print("=" * 70)

for p in [5, 7, 11, 13, 17, 19, 23]:
    f = lambda x, _p=p: (int(x) ** 2) % _p
    fc = build_function_crystal(f"x² mod {p}", f, domain)
    fp = theta_fingerprint(fc.joint)

    # The period of x² mod p
    outputs = [f(x) for x in range(1, p + 1)]
    n_distinct = len(set(outputs))
    # Number of QRs = (p-1)/2 + 1 (counting 0)

    print(f"  x² mod {p:2d}: Schmidt={fc.schmidt:.3f}, rel={fc.relationship:13s}, "
          f"asymm={theta_asymmetry(fp):.4f}, distinct={n_distinct}, "
          f"(p-1)/2={((p-1)//2)}")


# ═══════════════════════════════════════════════════════════════
#  Test 5: The crystal as Legendre symbol computer
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 5: CAN THE CRYSTAL COMPUTE THE LEGENDRE SYMBOL?")
print("Build crystal from (n/p) on small domain, query at large n")
print("=" * 70)

# Pick p=7: Legendre symbol has values {-1, 0, 1}
# QRs mod 7: 1, 2, 4 are QR (residues of 1², 3², 2²)
# NQRs: 3, 5, 6
p = 7
f_leg7 = lambda x: legendre(int(x), p)
fc_leg7 = build_function_crystal("(n/7)", f_leg7, domain)

print(f"\n  Crystal (n/7): Schmidt={fc_leg7.schmidt:.3f}")
print(f"  Relationship: {fc_leg7.relationship}")

# Query at specific values
print(f"\n  {'n':>5s} {'(n/7)':>6s} {'crystal_bin':>11s} {'conf':>6s}")
print("  " + "-" * 32)

for n in [1, 2, 3, 4, 5, 6, 8, 10, 14, 15, 49, 50, 100]:
    exact = legendre(n, p)
    result = fc_leg7.query(float(n))
    bin_label = ["QR(1)", "zero(0)", "NQR(-1)"][result["dominant_bin"]]
    print(f"  {n:>5d} {exact:>6d} {bin_label:>11s} {result['confidence']:>6.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 6: Cross-fingerprint of Legendre symbols
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 6: DO LEGENDRE SYMBOLS FOR DIFFERENT PRIMES CLUSTER?")
print("=" * 70)

leg_fps = {}
for p in [5, 7, 11, 13, 17, 19, 23, 29]:
    f = lambda x, _p=p: legendre(int(x), _p)
    fc = build_function_crystal(f"(n/{p})", f, domain)
    leg_fps[p] = theta_fingerprint(fc.joint)

# Pairwise distances
print(f"\n  Pairwise θ-distances between Legendre crystals:")
primes = sorted(leg_fps.keys())
print(f"  {'':>5s}", end="")
for p in primes:
    print(f" {p:>5d}", end="")
print()

for p1 in primes:
    print(f"  {p1:>5d}", end="")
    for p2 in primes:
        d = theta_distance(leg_fps[p1], leg_fps[p2])
        print(f" {d:>5.3f}", end="")
    print()

# Average pairwise distance
dists = []
for i, p1 in enumerate(primes):
    for p2 in primes[i+1:]:
        dists.append(theta_distance(leg_fps[p1], leg_fps[p2]))
avg_d = sum(dists) / len(dists)
print(f"\n  Average pairwise distance: {avg_d:.4f}")
print(f"  Min: {min(dists):.4f}, Max: {max(dists):.4f}")

# Comparison: distance from Legendre crystals to other function types
fc_mod = build_function_crystal("x mod 3", lambda x: int(x) % 3, domain)
fp_mod = theta_fingerprint(fc_mod.joint)
fc_phi = build_function_crystal("φ(n)", lambda x: compute.euler_totient(max(1, int(x))), domain)
fp_phi = theta_fingerprint(fc_phi.joint)

print(f"\n  Legendre vs x mod 3 (avg): {sum(theta_distance(fp, fp_mod) for fp in leg_fps.values()) / len(leg_fps):.4f}")
print(f"  Legendre vs φ(n) (avg):    {sum(theta_distance(fp, fp_phi) for fp in leg_fps.values()) / len(leg_fps):.4f}")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
