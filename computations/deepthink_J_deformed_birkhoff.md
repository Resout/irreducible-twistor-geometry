# Task: J-Deformed Birkhoff Factorisation on CP¹

## The problem

At the K*=7/30 equilibrium of the DS framework, the transition function M̃(Z) = M(π(Z)) is a pullback from S⁴ to CP³. It does not depend on the fibre coordinate ζ. With respect to the standard complex structure J₀, the Birkhoff factorisation is trivial and the Ward connection is flat.

But the Born floor deforms the almost complex structure from J₀ to J, with J ≠ J₀. Mason's theorem says F⁺ ∝ N_J ≠ 0. The physical Ward connection uses the J-deformed ∂̄-operator, not the standard one.

Compute the Ward-reconstructed curvature using the **deformed** Birkhoff factorisation.

## What you have

The equilibrium at 100-digit precision:

- m* = [0.786898446158, 0.029282259968, 0.029282259968, 0.154537033907]
- e* from make_evidence(p_dom) with K(m*,e*) = 7/30 exactly
- M* = (θI + s₁σ₁ + s₂σ₂ + s₃σ₃)/√2, det(M*) = -0.2985
- The anti-holomorphic Jacobian ∂̄Φ at m* has rank 1, ‖∂̄Φ‖ = 0.7305
- The floor correction Jacobian J_corr = (J_floor - I)·J_DS is known analytically

## What you need to compute

The deformed ∂̄-operator on the bundle restricted to a twistor line L_x ≅ CP¹ is:

∂̄_J = ∂̄_{J₀} + A_J

where A_J is a matrix-valued (0,1)-form on CP¹ encoding the Born floor's deformation of the complex structure. The task is:

### Step 1: Construct A_J on the twistor line

The deformation J - J₀ is determined by the Born floor's anti-holomorphic content. On the twistor line L_x parameterised by ζ, the deformation acts on the bundle through:

A_J(ζ) = M(x)⁻¹ · [∂̄_J M̃ - ∂̄_{J₀} M̃](ζ)

Since M̃ = M(π(Z)) is ∂̄_{J₀}-holomorphic (ζ-independent), we have ∂̄_{J₀} M̃ = 0. So:

A_J(ζ) = M(x)⁻¹ · ∂̄_J M̃(ζ)

The term ∂̄_J M̃ encodes how the J-deformation makes the pullback transition function non-holomorphic. It involves the Nijenhuis tensor N_J contracted with the fibre tangent direction at each ζ.

Concretely: the DS map Φ = Floor∘DS has anti-holomorphic Jacobian ∂̄Φ (the 4×4 matrix computed in Stage 2 of the previous script). The deformation A_J at fibre position ζ is the contraction of ∂̄Φ with the fibre tangent vector at ζ, embedded into the bundle via the Pauli map.

On the standard twistor line through x = 0 with π_{A'} = (1, ζ): the fibre tangent at ζ is the direction ∂/∂ζ̄ in CP³. This direction, in the mass function coordinates (s₁,s₂,s₃,θ), depends on ζ through the incidence relation. The anti-holomorphic derivative of the incidence gives the ζ-dependent fibre direction.

Use the incidence relation coefficients from ward_penrose_integral.py:

dz̄₁/dx^μ = A1[μ] + C1[μ]/ζ
dz̄₂/dx^μ = A2[μ] + C2[μ]/ζ

with A1[0] = A1[3] = -i/√2, C1[1] = -i/√2, C1[2] = 1/√2, A2[1] = -i/√2, A2[2] = 1/√2, C2[0] = -i/√2, C2[3] = i/√2 (all others zero).

The deformation A_J(ζ) on L₀ is:

A_J(ζ) = Σ_μ B[μ] · c_μ(ζ)

where B[μ] = M_ds⁻¹ · Σ_j (σ_j/√2) · J_corr[j,μ] are the connection basis matrices and c_μ(ζ) = (A1[μ]+A2[μ]) + (C1[μ]+C2[μ])/ζ are the incidence coefficients.

This gives A_J as a Laurent polynomial in ζ: A_J(ζ) = ρ₀ + ρ₋₁/ζ.

### Step 2: Solve the J-deformed Birkhoff equation

Find h₊(ζ) holomorphic inside |ζ| = 1 satisfying:

∂̄_{J₀} h₊ + A_J · h₊ = 0

Since A_J is a (0,1)-form on CP¹, this is an ODE in ζ̄. On the equator |ζ| = 1 where ζ̄ = 1/ζ, the equation becomes:

∂h₊/∂(1/ζ) = -A_J(ζ) · h₊(ζ)

For ‖A_J‖ < 1, solve perturbatively:

h₊ = I + h₊⁽¹⁾ + h₊⁽²⁾ + ...

where h₊⁽¹⁾ = -P₊[A_J] (Cauchy projection of A_J onto holomorphic modes), h₊⁽²⁾ = P₊[A_J · P₋[A_J]] (second-order correction from the non-commutative structure), etc.

Compute h₊ to at least second order.

### Step 3: Extract the Ward connection

A_μ(x) = h₊(ζ)⁻¹ · ∂_μ h₊(ζ)

The spacetime derivative ∂_μ acts through the incidence relation (changing which twistor line we're on). At x = 0, ∂_μ h₊ involves the ζ-derivative of h₊ weighted by the incidence coefficients.

### Step 4: Compute F

F_μν = ∂_μ A_ν - ∂_ν A_μ + [A_μ, A_ν]

Extract |F⁺|² = |F_{0'0'}|² (the Penrose residue component).

### Step 5: Report

- Is F_Ward = 0 or ≠ 0 at the equilibrium?
- If ≠ 0: what is |F⁺|²? What fraction is non-abelian?
- Does the result connect to the commutator content 8332/625?

## The key insight this tests

A pullback bundle M̃ = M(π(Z)) has no ζ-dependence in the J₀ sense. But the J-deformation rotates what "holomorphic" means. In the J-rotated frame, M̃ acquires effective ζ̄-dependence through A_J. The Birkhoff factorisation with respect to ∂̄_J is non-trivial. This is the content of Mason's theorem: non-integrable J produces F⁺ ≠ 0 even for pullback bundles.

If this computation gives F ≠ 0: the Ward reconstruction works and Gap 7 is closed.
If F = 0 even with the J-deformation: something deeper is wrong with the geometric bridge.

## Precision

100 digits is sufficient. The inputs (m*, e*, J_corr, B[μ], incidence coefficients) are all available from the previous scripts. The Cauchy projection is a discrete Fourier transform on CP¹. The perturbative expansion converges because ‖A_J‖ < 1.
