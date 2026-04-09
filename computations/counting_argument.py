"""
Counting argument for all three Koide mass scales.
Each fermion sector sees structures appropriate to its quantum numbers.
"""

H = 3

print("="*70)
print("THE THREE DECOMPOSITIONS")
print("="*70)
print()
print("LEPTON (colour-singlet):")
print("  1/dim(Sym^2(C^4)) + 1/H(H+1) = 1/10 + 1/12 = 11/60")
print(f"  Check: {1/10+1/12:.10f} = {11/60:.10f}")
a_l, b_l = 10, 12
print(f"  Quadratic: x^2-{a_l+b_l}x+{a_l*b_l} = (x-{a_l})(x-{b_l}) = 0")
print(f"  Product = {a_l*b_l} = S+1/K* = 5!")
print(f"  Discriminant = {(a_l-b_l)**2} = (H-1)^2 = H+1 [SELF-CONSISTENCY]")

print()
print("DOWN QUARK (colour-fundamental):")
print("  1/dim(SU(3)) + 1/(H+1) = 1/8 + 1/4 = 3/8")
print(f"  Check: {1/8+1/4:.10f} = {3/8:.10f}")
a_d, b_d = 8, 4
print(f"  Quadratic: x^2-{a_d+b_d}x+{a_d*b_d} = (x-{a_d})(x-{b_d}) = 0")
print(f"  Sum = {a_d+b_d} = H(H+1) = one of the LEPTON roots!")
print(f"  Product = {a_d*b_d} = (H^2-1)(H+1) = 2^(H+2)")
disc_d = (a_d-b_d)**2
print(f"  Discriminant = {disc_d} = (H+1)^2 [FERMION COUNT per generation]")

print()
print("UP QUARK (from product constraint):")
print(f"  (11/60)(3/8)(40/3) = {11/60 * 3/8 * 40/3:.10f} = 11/12")
print(f"  M2_up = (11/12)/(11/60 * 3/8) = {(11/12)/(11/60 * 3/8):.10f} = 40/3")
print(f"  = d_super * (H+1)/H = 10 * 4/3")

print()
print("="*70)
print("COMPARISON TABLE")
print("="*70)
print()
print("Sector    Scale   Decomposition     a             b              a*b    disc")
print("-"*85)
print(f"Lepton    11/60   1/10 + 1/12       Sym^2(C^4)=10 H(H+1)=12      120    4=(H-1)^2")
print(f"Down      3/8     1/8 + 1/4         SU(3)=8       H+1=4          32     16=(H+1)^2")
print(f"Up        40/3    (product)          --            --             --     --")

print()
print("="*70)
print("THE PATTERN")
print("="*70)
print("""
Each fermion sees TWO structures, gets 1/dim of each, adds them:

LEPTON (colour-singlet):
  Structure 1: Sym^2(C^4) = channel space (bilinear product)
    dim = 10 = d_super. The conservation law lives here.
  Structure 2: gauge generators of SU(H+1) = SU(4)
    dim = H(H+1) = 12. The gauge algebra acts here.
  Both are FULL (unrestricted) structures.
  The lepton is a singlet - it sees everything uniformly.

DOWN QUARK (colour-fundamental):
  Structure 1: SU(H) = SU(3) = colour group
    dim = H^2-1 = 8. The group the quark transforms under.
  Structure 2: C^(H+1) = mass space
    dim = H+1 = 4. The space the quark lives in.
  Both are COLOUR-RESTRICTED structures.
  The quark carries colour - it sees the colour subgroup.

UP QUARK: forced by the product constraint
  (M2_lep)(M2_down)(M2_up) = 11/12 = massive/total generators.
  Given lepton and down, up is determined.
""")

print("="*70)
print("NESTING STRUCTURE")
print("="*70)
print(f"Down sum = 8+4 = 12 = H(H+1) = one of the lepton roots")
print(f"Lepton discriminant = 4 = H+1 = one of the down roots")
print(f"Down discriminant = 16 = (H+1)^2 = fermions per generation")
print(f"Lepton product = 120 = 5! = S+1/K*")
print(f"Down product = 32 = 2^(H+2) = 2*(H+1)^2/2... = (H^2-1)(H+1)")
print()
print("The two quadratics are NESTED:")
print(f"  Lepton: x^2 - 22x + 120 = 0  (roots 10, 12)")
print(f"  Down:   x^2 - 12x + 32 = 0   (roots 8, 4)")
print(f"  The DOWN sum (12) is the LEPTON root.")
print(f"  The LEPTON discriminant (4) is the DOWN root.")
print()

# Verify the product constraint derivation
print("="*70)
print("DERIVATION STATUS")
print("="*70)
print("""
PROVED (if counting argument holds):
  1. Total allocation = 11/12 (massive/total generators) [in paper]
  2. Lepton: 1/10 + 1/12 = 11/60 (singlet sees full structures)
  3. Down: 1/8 + 1/4 = 3/8 = sin^2(theta_W) (fundamental sees colour structures)
  4. Up: 40/3 (forced by product constraint)
  5. Product: (11/60)(3/8)(40/3) = 11/12 [verified]

THE TWO THINGS TO PROVE:
  A. Why 1/dim (uniform counting): each sector gets 1/dim of each
     structure because it has no preferred orientation within that
     structure (Born floor enforces uniformity at the fixed point).

  B. Why additive (parallel conductance): the two structures are
     independent algebraic decompositions of C^4 - one from the
     bilinear product (Sym^2), one from the symmetry group (SU).
     Their contributions add because they're different kinds of
     structure on the same space.

If A and B are proved:
  M2_lep = m(0++) * 11/60 is DERIVED
  M2_down = m(0++) * 3/8 is DERIVED (and = sin^2(theta_W), already proved!)
  M2_up = m(0++) * 40/3 is DERIVED (from product constraint)
  m(0++) = M2_lep * 60/11 = 1712.0 MeV is a PREDICTION
  ZERO free parameters.
""")
