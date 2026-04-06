# The Proof

A complete listing of every theorem, proposition, lemma, and corollary in the paper, in logical dependency order, with proof sketches and dependency chains.

80 proved statements. Zero conjectures. All external dependencies cited.

**Source:** "The Mass Gap from Enforced Uncertainty" by J. R. Manuel, 2026.
**Code:** github.com/Resout/irreducible-twistor-geometry (267 scripts)

---

## What is proved

A non-trivial quantum field theory exists on R^4 with a mass gap Delta > 0, for any compact simple gauge group G. The theory satisfies the Osterwalder-Schrader axioms OS0-OS4 and is identified with Yang-Mills theory via Mason's twistor correspondence and the Popov identification. Einstein gravity emerges from the (3,3) sector. Navier-Stokes regularity follows conditionally on Born floor preservation.

---

## Layer 1: Uniqueness (the starting point)

### Theorem 1: Unique Self-Consistent Tangent Bundle (`thm:unique_cpn`)
**Statement:** CP^1 is the unique compact complex projective space whose tangent bundle has self-consistent section algebra: (H-1)^2 = H+1 with H = dim H^0(T(CP^n)) has n=1 as its only positive integer solution.
**Proof sketch:** H = (n+1)^2 - 1 (generators of PGL(n+1,C)). Substituting: ((n+1)^2 - 2)^2 = (n+1)^2 gives n(n-1)(n+2)(n+3) = 0. Unique positive integer: n = 1. At n=1: T(CP^1) = O(2), H = 3, (3-1)^2 = 4 = 3+1.
**Depends on:** dim H^0(T(CP^n)) = (n+1)^2 - 1.

### Theorem 2: Self-Consistency Selects H=3 (`thm:selfconsist`)
**Statement:** The Born floor 1/H^3 equals eta(H)/(H+1) if and only if H = 3.
**Proof sketch:** Setting 1/H^3 = (H-1)^2/(H^3(H+1)) gives (H-1)^2 = H+1, i.e. H^2 - 3H = 0. Unique positive solution H = 3.
**Depends on:** eta(H) = (H-1)^2/H^3.

### Theorem 3: Optimal Dimensionality (`thm:optimal`)
**Statement:** eta(H) = (H-1)^2/H^3 has unique global maximum at H=3, with eta(3) = 4/27.
**Proof sketch:** d(eta)/dH = (H-1)(3-H)/H^4. Zeros at H=1 (min) and H=3 (max). Second derivative negative at H=3.
**Depends on:** Definition of eta(H).

### Theorem 4: Robustness (`thm:robust`)
**Statement:** For generalised cost eta_beta(H) = (H-1)^2/H^{2+beta}, the integer optimum is H=3 for all beta in (0.82, 1.42).
**Proof sketch:** Evaluate eta_beta at H=2,3,4. H=3 dominates for 0.82 < beta < 1.42.
**Depends on:** Theorem 3.

### Theorem 5: Born Uniqueness — Gleason (`thm:gleason`)
**Statement:** On a Hilbert space of dimension >= 3, the Born rule p_i = |m_i|^2/sum|m_j|^2 is the unique countably additive probability measure on closed subspaces.
**Proof sketch:** External (Gleason 1957). Applies because dim(C^4) = 4 >= 3.
**Depends on:** Gleason's theorem (external).

### Lemma 6: Conflict Bound (`lem:Kbound`)
**Statement:** For real non-negative mass functions with theta, theta' > 0: K < 1.
**Proof sketch:** K = sum_{i!=j} s_i s'_j <= (1-theta)(1-theta') < 1.
**Depends on:** L_1 = 1 constraint.

### Theorem 7: L_1 Conservation (`thm:L1`)
**Statement:** If sum m_i = sum e_i = 1, then sum (m_out)_i = 1.
**Proof sketch:** Pre-normalisation sum = agreement + theta + theta' - theta*theta' = 1 - K. Division by (1-K) restores L_1 = 1. Algebraic identity.
**Depends on:** DS combination rule.

### Theorem 8: Phase Sensitivity (`thm:phase`)
**Statement:** Two mass vectors with identical Born probabilities but different phases produce different outputs under Dempster combination.
**Proof sketch:** Explicit construction with m real and m' complex having same |m_i|^2 ratios. K(m,m'') = 0.47 vs K(m',m'') = 0.506+0.081i. Different K gives different output.
**Depends on:** DS combination rule, Born measurement.

### Theorem 9: Complex Encoding Uniqueness (`thm:complex_unique`)
**Statement:** Among finite-dimensional normed division algebras over R (Hurwitz: R, C, H, O), C is the unique algebra admitting an ordering-sensitive encoding.
**Proof sketch:** R: Aut = {id}, fails non-trivial. C: conjugation z-bar satisfies all four conditions. H: Aut = SO(3), too many involutions, fails uniqueness. O: Aut = G_2, same problem.
**Depends on:** Hurwitz's theorem, Skolem-Noether theorem (external).

### Theorem 10: Unique Product (`thm:unique_product`)
**Statement:** The bilinear product on C^{H+1} respecting the O(2) section/base decomposition, with S_H symmetry, commutativity, section locality, vacuous identity, and minimal conflict, is uniquely the DS combination rule. Zero free parameters.
**Proof sketch:** Five axioms. (1)-(2) give general form with 3 free parameters. (3) kills cross-section terms. (4) forces p = b. (5) forces a = b = 1. Zero parameters remain.
**Depends on:** O(2) bundle structure from Theorem 1.
**Note:** Each of the five axioms is forced by the O(2) section/base decomposition of T(CP^1) = O(2). S_H symmetry: the H sections are interchangeable (O(2) has no preferred section). Commutativity: the product on sections is symmetric. Section locality: sections of O(2) are local data on CP^1. Vacuous identity: the base component (theta) acts as a neutral element. Minimal conflict: the discarded cross-focal mass is the minimal projection compatible with L_1 conservation. The DS combination rule is not assumed — it is the unique algebraic consequence of the bundle geometry.

### Theorem 11: Floor Consistency (`thm:floorconsist`)
**Statement:** Born floor 1/H^3 = eta(H)/(H+1) iff H = 3.
**Proof sketch:** Same equation as Theorem 2. At H=3: 1/27 = (4/27)/4 = 1/27.
**Depends on:** Theorem 2.

---

## Layer 2: Convergence and Dynamics

### Proposition 12: Norm Convergence (`prop:norm`)
**Statement:** Under non-uniform evidence, L_2/L_1 increases monotonically.
**Proof sketch:** Ratios m_i/m_j = (e_i/e_j)^n diverge exponentially under non-uniform evidence.
**Depends on:** DS combination (multiplicative in ratios).

### Proposition 13: Phase Washout (`prop:washout`)
**Statement:** Phase differences decay exponentially under shared evidence.
**Proof sketch:** Phase-insensitive terms (theta*e_i) grow relative to phase-sensitive terms. |Delta_phi_n| <= C*(1-delta)^n with delta > 0 (Born floor ensures theta != 0).
**Depends on:** DS combination, Born floor.

### Theorem 14: Self-Sorting (`thm:selfsort`)
**Statement:** Under non-uniform evidence, mass concentrates on <= 2 hypotheses at rate (e_max/e_min)^n.
**Proof sketch:** m_i(n) proportional to m_i(0)*e_i^n. Ratios diverge exponentially.
**Depends on:** DS combination (multiplicative structure).

### Theorem 15: Fubini-Study State Space (`thm:fubini`)
**Statement:** The state space is CP^3 with the Fubini-Study metric (unique monotone Riemannian metric under CPTP maps).
**Proof sketch:** Born measurement is degree-0 homogeneous; equivalence classes are points in CP^3. Petz classification gives uniqueness.
**Depends on:** Theorem 5, Petz classification (external).

### Theorem 16: Full BMU Reconstruction (`thm:BMU`)
**Statement:** The framework satisfies all four Barnum-Muller-Ududec axioms and is an instance of complex quantum theory with Born rule.
**Proof sketch:** Axiom 1: finite dimension. Axiom 2: Born is quadratic, I_3 = 0. Axiom 3: PU(4) is 2-point homogeneous on CP^3. Axiom 4: fractional evidence defines one-parameter group (Lemma 17); regularised flow; limit preserves axioms.
**Depends on:** Theorems 5, 7, 15, Lemma 17, BMU (2014, external).

### Lemma 17: Uniform Rescaling (`lem:uniform`)
**Statement:** The (1-K) normalisation is a uniform scalar, hence projectively trivial on CP^3.
**Proof sketch:** 1/(1-K) multiplies every component equally. On CP^3, uniform scaling is the identity.
**Depends on:** DS combination rule.

---

## Layer 3: The Mass Gap (pure algebra)

### Theorem 18: Structural Filter (`thm:filter`)
**Statement:** Cross-focal products s_i*e_j (i!=j) appear in no output channel. Their total is K. Pre-normalisation total is 1-K. At K* > 0, after n steps: (1-K*)^n = exp(-n*Delta).
**Proof sketch:** Algebraic identity: agreement + theta + theta' - theta*theta' = 1-K. The missing fraction K is never generated. Iteration gives exponential decay.
**Depends on:** Theorem 7.

### Theorem 19: Mass Gap from Universal Correlation (`thm:massgap`)
**Statement:** If delta(C) < 1 for every observable pair, the system has mass gap Delta = -ln(1-K*) > 0.
**Proof sketch:** delta < 1 gives K > 0 at each step. Global convergence (Theorem 23) gives unique K*. Structural filter (Theorem 18) gives exponential decay.
**Depends on:** Theorems 18, 23.

### Corollary 20: Structural Filter Rate (`cor:gap`)
**Statement:** K* = 7/30 gives Delta_SF = -ln(23/30) = 0.266. The actual spectral gap Delta depends on the gauge group through the evidence distribution and may differ.
**Depends on:** Theorems 18, 24.

### Corollary 21: No Gap from Perfect Correlation (`cor:nogap`)
**Statement:** If delta = 1 for any pair, that pair has K* = 0. No gap.
**Depends on:** Theorem 19.

### Theorem 22: Channel Discount (`thm:discount`)
**Statement:** The discount factor is alpha = H/(H^2+1) = 3/10.
**Proof sketch:** Each factor contributes H singletons; joint channel has H^2+1 dimensions. Ratio = H/(H^2+1).
**Depends on:** Tensor product structure.

### Theorem 23: Global Convergence (`thm:basin`)
**Statement:** DS + Born floor has unique fixed point m*. All trajectories converge in the Hilbert projective metric. Complex extension via phase washout.
**Proof sketch:** After 2 steps, all components >= epsilon > 0 (Born floor). On compact S_epsilon, DS is a positive linear map with contraction kappa < 1. Banach fixed-point theorem. Complex: phase washout then magnitude convergence.
**Depends on:** Born floor, Proposition 13, Banach fixed-point theorem.

### Theorem 24: Fixed-Point Conflict (`thm:fixedpoint`)
**Statement:** K* = (H^2-H+1)/[H(H^2+1)]. At H=3: K* = 7/30.
**Proof sketch:** Solve K*(H^2+1) - eta*H^2 = 1. Verify: (H^2-H+1)/H - (H-1)^2/H = H/H = 1 for all H. Algebraic identity.
**Depends on:** Conservation law, Theorem 7.

### Theorem 25: Sym^2 Characterisation (`thm:sym2`)
**Statement:** dim(Sym^2(C^{H+1})) = H^2+1 iff H = 3. The conservation law is the Sym^2 decomposition into trace and traceless.
**Proof sketch:** (H+1)(H+2)/2 = H^2+1 gives H^2 - 3H = 0. At H=3: dim(Sym^2(C^4)) = 10 = H^2+1.
**Depends on:** Conservation law, Theorem 7, Lemma 17.

### Theorem 26: Budget Matching (`thm:budget`)
**Statement:** K_cons(H) = 7/30 iff H = 3.
**Proof sketch:** Cross-multiply, factor: (H-3)(7H^2-9H+10) = 0. Discriminant -199 < 0. Unique positive root H=3.
**Depends on:** Theorem 24.

---

## Layer 4: Application to Yang-Mills

### Theorem 27: Universal Correlation (`thm:universal`)
**Statement:** For any non-abelian compact simple G, every pair of gauge-invariant observables has delta < 1.
**Proof sketch:** delta > 1/H: Clebsch-Gordan cross-terms from shared links (non-abelian has dim(R) >= 2). delta < 1: independent spatial links give Var(O_2|O_1) > 0. For abelian U(1): all irreps 1D, CG trivial, coprime charges factorise.
**Depends on:** Peter-Weyl, Clebsch-Gordan, Wilson action coupling.

---

## Layer 5: The Ward Correspondence and Googly Problem

### Theorem 28: Penrose-Ward (`thm:ward`)
**Statement:** Holomorphic rank-r bundles E -> PT trivial on every twistor line correspond 1-1 with ASD connections on S^4. Reconstructed F^+ = 0.
**Proof sketch:** External (Ward 1977, Mason-Woodhouse). Key: (dbar_E)^2 = 0 iff F^+ = 0.
**Depends on:** Standard twistor theory (external).

### Theorem 29: Born Floor Guarantees Triviality (`thm:trivial`)
**Statement:** For real mass functions at Born floor equality: det(M) = -25*theta^2/2 != 0.
**Proof sketch:** At floor: sum s_i^2 = 26*theta^2. So det(M) = (theta^2 - 26*theta^2)/2 = -25*theta^2/2.
**Depends on:** Pauli embedding, Born floor.

### Theorem 30: Dynamic Non-Degeneracy (`thm:det_dynamic`)
**Statement:** Under coupled DS + floor, |det(M)| is bounded below along trajectories.
**Proof sketch:** Combines: (1) det=0 has no fixed point or periodic orbit (Thm 34). (2) Transverse instability lambda >= 3.17 (Prop 32). (3) Rank-2 coupling prevents evidence alignment (Thm 33, residual >= 12). Infimum positive on compact B.
**Depends on:** Theorems 33, 34, Proposition 32.

### Theorem 31: Light Cone Repulsion (`thm:lightcone`)
**Statement:** If det(M(m)) = 0, then DS output has det(M'') < 0 (real positive masses).
**Proof sketch:** Q_pre = Q(m)*Q(m') + R with R < 0 (all factors positive, coefficients negative). On Q(m)=0: Q_pre = R < 0.
**Depends on:** DS combination, Pauli embedding.

### Proposition 32: Exponential Instability (`thm:instability`)
**Statement:** det=0 surface is exponentially unstable. Transverse amplification lambda > 1 at every tested state (min 3.17 over 68,910 measurements).
**Proof sketch:** R=0 cancellation requires exact coordination of independent algebraic quantities. Perturbation breaks coordination; combination rule amplifies.
**Depends on:** Theorem 31.

### Theorem 33: Rank-2 Protection (`thm:rank2protect`)
**Statement:** On det=0, the DS output has Q'' = Lorentzian form in shifted evidence. Rank-2 coupling gives evidence light cone residual >= 12. Rank-1 fails (|det| reaches 10^{-152}).
**Proof sketch:** phi^2 coefficient = Q(m) = 0 on the surface. Remaining form is Lorentzian. Rank-2 coupling misaligns evidence phases across both bundle channels.
**Depends on:** Theorem 31, Theorem 41 (entanglement).

### Theorem 34: Self-Entanglement Excludes det=0 (`thm:selfentangle`)
**Statement:** No fixed point or periodic orbit of the coupled DS map has det = 0.
**Proof sketch:** Case theta_A != 0: fixed-point equations force m_A = m_B = (1,0,0,0), det = 1/2 != 0. Case theta_A = 0: forces all nonzero s_i equal, det = -1/(2n) != 0. Periodic orbits require R=0 at every step: transversely unstable and unreachable.
**Depends on:** DS fixed-point equations, Proposition 32, Theorem 33.

### Theorem 35: DS-Matrix Decomposition (`thm:dsmatrix`)
**Statement:** sqrt(2)*M''_pre = {M,E} + D - (s.e)*I, where {M,E} is the anticommutator, D = sum s_i*e_i*sigma_i.
**Proof sketch:** Expand using sigma_i*sigma_j = delta_ij*I + i*eps_ijk*sigma_k. Both sides equal theta*phi*I + sum(s_i*e_i + theta*e_i + phi*s_i)*sigma_i.
**Depends on:** Pauli algebra, DS combination.

### Corollary 36: Commutator Absence (`cor:commutator`)
**Statement:** [M,E] = i*(s x e).sigma is absent from the DS output.
**Proof sketch:** Rearrangement of Theorem 35: sqrt(2)*M''_pre = 2ME - [M,E] + D - (s.e)*I.
**Depends on:** Theorem 35.

### Theorem 37: Off-Diagonal Decomposition (`thm:offdiag`)
**Statement:** Off-diagonal products decompose into K (symmetric, normalisation deficit) and [M,E] (antisymmetric, curvature content). Both absent from output.
**Proof sketch:** s_i*e_j = (1/2)(s_i*e_j + s_j*e_i) + (1/2)(s_i*e_j - s_j*e_i). Symmetric parts sum to K. Antisymmetric parts are s x e.
**Depends on:** Theorem 35, Pauli commutation.

### Theorem 38: Conflict as Frobenius Distance (`thm:frobenius`)
**Statement:** K(m,e) = (1/2)||M-E||_F^2 + (1/2)K_self(m) + (1/2)K_self(e) - (theta-phi)^2.
**Proof sketch:** Polarisation identity on (1-theta)(1-phi) and -s.e, using ||M-E||_F^2 = (theta-phi)^2 + ||s-e||^2.
**Depends on:** DS conflict definition, Frobenius norm.

### Theorem 39: Mason (2005) (`thm:mason_orig`)
**Statement:** Non-integrable J on CP^3 produces full Yang-Mills (F^+ proportional to N_J). Integrable J gives only ASD (F^+ = 0).
**Proof sketch:** External (Mason 2005, confirmed by Popov 2021).
**Depends on:** External.

### Theorem 40: Non-Holomorphicity of Born Floor (`thm:nonholo`)
**Statement:** The DS + floor map Phi is not holomorphic when the floor activates: d(Phi_theta)/d(theta-bar) = -c*theta^2/(2|theta|^3) != 0.
**Proof sketch:** Floor involves |theta|^2 = theta*theta-bar. Wirtinger derivative nonzero.
**Depends on:** Born floor structure.

### Theorem 41: Entanglement of the Kink (`thm:entangle`)
**Statement:** The anti-holomorphic derivative dbar(M) has rank 2 (both bundle channels coupled).
**Proof sketch:** dbar(M) = (dM/dm)*dbar(Phi). Image has components along I (from delta_theta) and sigma_i (from delta_s_i). Generically both singular values nonzero.
**Depends on:** Theorem 40, Pauli embedding.

### Proposition 42: Generic Entanglement (`prop:generic_ent`)
**Statement:** For floor-active states with s_i != 0 for >= 2 hypotheses: rank-2 in 99.9% of 4225 states, mean S_ent = 0.34 bits.
**Depends on:** Theorem 41.

### Theorem 43: dbar-Operator on DS Bundle (`thm:dbar`)
**Statement:** Bundle E_m has dbar_E = dbar + a^{0,1} where a = M^{-1}*dbar(M). Floor inactive: a = 0. Floor active: a != 0.
**Proof sketch:** DS is holomorphic (zero dbar contribution). Floor gives dbar(Phi) != 0 by Theorem 40. Chain rule.
**Depends on:** Theorem 40.

### Theorem 44: Non-Integrability Generates F^+ (`thm:nijenhuis`)
**Statement:** Floor active implies N_J != 0 (non-integrable J). By Mason: F^+ != 0. Physical F^+ is the Penrose residue rho_{-1}: ||rho_{-1}|| = 0.638, |F^+|^2 = 0.407 at 120 digits.
**Proof sketch:** dbar(Phi) != 0 implies N_J != 0 (Newlander-Nirenberg). Maurer-Cartan form has F = 0 identically. Ward connection is pure gauge. Physical F^+ extracted by Penrose contour integral (residue of 1/zeta pole in Laurent expansion).
**Depends on:** Theorem 40, Mason (Theorem 39), Maurer-Cartan identity, Penrose transform.

### Theorem 45: Mason Matching (`thm:mason`)
**Statement:** The DS construction satisfies all three Mason conditions: (i) rank-2 bundle, (ii) non-integrable J, (iii) twistor fibration. Full Yang-Mills follows. F^+ proportional to dbar(epsilon) at leading order.
**Proof sketch:** Conditions verified: (i) Pauli embedding, (ii) Theorem 44, (iii) standard pi: CP^3 -> S^4. Application is deduction.
**Depends on:** Theorems 39, 44, bundle construction.

### Proposition 46: Chirality Balance (`prop:balance`)
**Statement:** Phase washout forces Im(K) -> 0, hence |F^+| = |F^-| at the fixed point.
**Proof sketch:** At fixed point, M in GL(2,R). Real M gives sigma-invariant bundle. F^+ = conj(F^-), so |F^+| = |F^-|.
**Depends on:** Proposition 13, reality structure.

### Theorem 47: Commutator Content (`thm:condensate`)
**Statement:** C * det(M*)^2 = 4(3n^2+2n+3)/(n-1)^2 where n = H^3-1 = 26. At H=3: C*det(M*)^2 = 8332/625.
**Proof sketch:** Wick factorisation on Gaussian perturbations. Expand M*^{-1} in Pauli basis. At Born floor |s|^2 = n*theta^2, numerator = 16(3n^2+2n+3)*theta^4. Cross-validated: analytic 149.596 vs Monte Carlo to 0.23%.
**Depends on:** DS fixed point, Born floor, Pauli algebra.

---

## Layer 6: Conformal Breaking, Gravity, and Regularity

### Theorem 48: DS Holomorphicity (`thm:ds_holo`)
**Statement:** DS combination without floor is holomorphic on {K != 1}.
**Proof sketch:** Rational function of z (not z-bar).
**Depends on:** DS combination (polynomial structure).

### Theorem 49: Born Floor Breaks Conformal Invariance (`thm:conformal_break`)
**Statement:** Floor enforcement does not commute with biholomorphisms of CP^3.
**Proof sketch:** Explicit counterexample: swap s_1 and theta. Floor activates for one but not the other. |theta| transforms non-covariantly.
**Depends on:** Theorem 48, Theorem 40.

### Theorem 50: Spin Decomposition (`thm:spin_decomp`)
**Statement:** dbar(Phi) decomposes under SO(4) as (1,1) + (3,1)+(1,3) + (3,3): scalar + gauge + graviton.
**Proof sketch:** Standard: End(R^4) = (2,2) tensor (2,2) = (1,1)+(3,1)+(1,3)+(3,3). Trace, antisymmetric, symmetric traceless.
**Depends on:** Theorem 40, SO(4) representation theory.

### Theorem 51: Graviton (`thm:graviton`)
**Statement:** The (3,3) component of dbar(Phi) is a spin-2 field satisfying conformal gravity.
**Proof sketch:** Mason: non-integrable J on T(CP^3) gives conformal gravity on S^4. Penrose transform: Sym^2_0(R^4) = (3,3) maps to massless spin-2 fields.
**Depends on:** Theorems 40, 44, 50, Mason (2005), AHS (1978).

### Theorem 52: Conformal -> Einstein Reduction (`thm:bach_to_einstein`)
**Statement:** Three constraints reduce Bach equation toward Einstein: (1) preferred metric from Born floor, (2) chirality balance, (3) spectral gap selects decaying sector.
**Proof sketch:** (1) fixes conformal gauge. (2) zeros Weyl contribution to Gauss-Bonnet. (3) Maldacena: decaying sector = Einstein with Lambda.
**Depends on:** Theorems 49, Proposition 46, OS4, Maldacena (2011).

### Theorem 53: Einstein Gravity (`thm:einstein`)
**Statement:** Spin-2 content satisfies Einstein equation with Lambda, not Bach equation.
**Proof sketch:** Chain A: geometry -> conformal gravity. Chain B: algebra -> OS2 -> positive-definite Hilbert space. Convergence: ghost sector excluded. Maldacena: ghost-free conformal gravity = Einstein. Non-circular.
**Depends on:** Theorem 44, OS2 (Theorem 63), Maldacena (2011), Adamo-Mason (2014).

### Theorem 54: Massless Graviton (`thm:massless_graviton`)
**Statement:** Spin-2 graviton is massless. Mass gap Delta applies to fibre perturbations, not base variations.
**Proof sketch:** Equilibrium manifold M = SU(2) orbit of m*. Fibre perturbations (change K from K*) decay at Delta. Base variations (rotate m* smoothly across S^4, K=K* everywhere) are on-shell, massless by Penrose transform.
**Depends on:** Theorems 53, 51, 45, OS4.

### Theorem 55: Multi-Site Spectral Gap (`thm:multisite`)
**Statement:** Fibre spectral gap Delta is independent of lattice size N.
**Proof sketch:** (1) Fourier: block-circulant Jacobian, max rho = 2*rho(J_m) = 0.300 < 1 at self-evidence; 0.779 at physical equilibrium. (2) Finite propagation speed: identical trajectories for t < N/2. (3) Transient amplification decreases with N.
**Depends on:** OS4, DS combination structure.

### Theorem 56: Universal Curvature Bound (`thm:curvature_bound`)
**Statement:** ||F|| <= C(H), independent of base manifold, initial data, and time.
**Proof sketch:** dbar(Phi) bounded on compact B. F^+ proportional to dbar(epsilon), bounded. F^- from Ward reconstruction, smooth on compact B. Bound is fibre-local.
**Depends on:** Theorems 29, 34, 45, compactness of B.

### Theorem 57: Penrose Integrals Finite (`thm:penrose_finite`)
**Statement:** Penrose transform of sections supported on B is finite for every base point.
**Proof sketch:** Bounded integrand on compact contour L_x = CP^1.
**Depends on:** Theorem 56.

### Theorem 58: Regularity of (3,1)+(1,3) Sector (`thm:regularity`)
**Statement:** int_0^T sup_x ||dbar(Phi)_{gauge}|| dt <= C(H)*T < infinity.
**Proof sketch:** Orthogonal decomposition, each component bounded by C(H), fibre-local.
**Depends on:** Theorems 50, 56.

### Theorem 59: Minitwistor Construction (`thm:minitwistor`)
**Statement:** DS on CP^3 reduces to minitwistor on T(CP^1) = O(2) for R^3. The (3,1)+(1,3) sector maps to vorticity omega satisfying Beltrami equation.
**Proof sketch:** Hitchin dimensional reduction O(1)+O(1) -> O+O(2). Incidence relation eta(zeta) = x^1 + x^2*zeta + x^3*zeta^2. Born floor on Z gives compact B_Z, bounded vorticity.
**Depends on:** Theorems 50, 58, 56, Hitchin (1982).

### Theorem 60: Chen-Hou Embedding (`thm:chen_hou`)
**Statement:** The Chen-Hou operator d_r^2 + (3/r)d_r + d_z^2 is the R^5 Laplacian. 3/r = H/r = dim(S^3)/r. Blowup point (1,0) lies outside Born disk r^2+z^2 <= 0.621.
**Proof sketch:** R^5 Laplacian in (r,z) coords: (n-2)/r = 3/r for n=5. Born constraint restricts to r^2+z^2 <= 26*theta_*^2 = 0.621. At (1,0): 1 > 0.621.
**Depends on:** Theorem 59, Chen-Hou (2025).

### Theorem 61: NS Regularity — CONDITIONAL (`thm:ns_regularity`)
**Statement:** IF ||omega_0||_inf <= sqrt(26) and NS evolution preserves Born >= 1/27 on the minitwistor representative, THEN the solution is smooth for all time.
**Proof sketch:** (0) Born floor from O(2) geometry. (1) Hitchin: every div-free field has O(2) representative. (2) ||omega||^2 <= 26 equivalent to Born >= 1/27. (3) NS terms act on R^3 only (fibre locality). (4-5) THE ASSUMPTION: Born floor preserved under NS. Three supporting observations (L^2 decay, topology, diffusion) but no proof. (6) Under assumption: ||omega||_inf <= sqrt(26) for all t. (7) BKM: int_0^T ||omega||_inf dt <= sqrt(26)*T < infinity.
**Depends on:** Theorem 59, Hitchin (1982), BKM (1984). CONDITIONAL on Born floor preservation (equivalent to the NS regularity problem itself).

---

## Layer 7: Construction and Axioms

### Theorem 62: OS0 — Temperedness (`thm:OS0`)
**Statement:** Schwinger functions are tempered distributions.
**Proof sketch:** Bounded functions on compact B. Conformal factor gives polynomial decay on R^4.
**Depends on:** Compactness of B.

### Theorem 63: OS1 — Euclidean Covariance (`thm:OS1`)
**Statement:** On S^4, Schwinger functions are SO(5)-covariant. On R^4 (via stereographic projection), they are E(4)-covariant.
**Proof sketch:** (1) S_2^c(x,y) = f(d_FS(L_x, L_y)) — spacetime enters only through fibre distance (spectral formula). Phi determines f; f is a function of one variable. (2) d_FS(L_{Rx}, L_{Ry}) = d_FS(L_x, L_y) for R in SO(5) — Sp(2)-equivariance of twistor fibration covers full SO(5) isometry group of S^4. SO(5) transitivity implies d_FS depends only on geodesic distance. (3) E(4) embeds in SO(5): rotations stabilise the north pole; translations lift to SO(5) rotations that move it. S_2^c(x+a,y+a) = S_2^c(R_a x, R_a y) = S_2^c(x,y) by Step 2. Conformal factor absent because d_FS is an SO(5)-invariant, not a conformal invariant.
**Depends on:** Theorem 67 (Schwinger functions), Sp(2)-equivariance, PU(4)-invariance, SO(5) transitivity on S^4.

### Theorem 64: OS2 — Reflection Positivity (`thm:OS2`)
**Statement:** Schwinger functions satisfy reflection positivity.
**Proof sketch:** Three arguments: (1) Koopman positive + Perron-Frobenius + Bochner. (2) DS commutative, Theta T = T*, ||Tf||^2 >= 0. (3) M_{ij} = S(n_i+n_j) = C^T D C with D >= 0.
**Depends on:** Koopman theory, Perron-Frobenius, commutativity of DS.

### Theorem 65: OS3 — Regularity (`thm:OS3`)
**Statement:** Schwinger functions are smooth on (S^4)^n including diagonals.
**Proof sketch:** Self-combination is polynomial. No UV divergence in finite-dimensional space.
**Depends on:** Finite dimensionality.

### Theorem 66: OS4 — Cluster Decomposition (`thm:OS4`)
**Statement:** Connected correlations decay at rate Delta = -ln(lambda_0) = 1.263 > 0.
**Proof sketch:** Structural filter gives exponential decay at K* > 0. Transfer operator eigenvalue lambda_0 = 0.2829 at K*=7/30 (six equations, zero free parameters, three independent numerical methods).
**Depends on:** Theorem 18, Theorem 24, Born floor.

### Theorem 67: Path Integral Measure (`thm:measure`)
**Statement:** nu = delta(m_{t+1} - Phi(m_t, e*)) * d(mu_FS)(m_0) is well-defined. Koopman operator bounded on L^2(B). Connected Schwinger functions decay as lambda_0^t. OS reconstruction yields QFT with Delta = 1.263.
**Proof sketch:** B compact (Tychonoff), Phi smooth (Radon-Nikodym bounded), spectral decomposition gives exponential decay, OS0-OS4 established.
**Depends on:** OS0-OS4, compactness of B, OS reconstruction theorem (external).

### Theorem 68: Explicit Schwinger Functions (`thm:schwinger4d`)
**Statement:** S_2^c(x,y) = sum |c_k|^2 lambda_k^{d_FS(L_x,L_y)/delta}. On R^4: S_2^c ~ C'*exp(-m_phys*|x-y|/2).
**Proof sketch:** Fibre-to-spacetime via twistor fibration. Propagator from path-ordered DS steps. Integration via Koopman spectral decomposition.
**Depends on:** Theorem 67, twistor fibration, Penrose correspondence.

### Theorem 69: Wilson Loop Area Law (`thm:arealaw`)
**Statement:** -ln<|W(C)|> = sigma*A + mu*Perim(C) + O(1), with sigma = -ln(23/30) = 0.266.
**Proof sketch:** (1) W(C) from path-ordered DS transfer operators. (2) At equilibrium, e* != m*, K = 7/30, [M,E] != 0. (3) Per-plaquette action -ln(1-K*). (4) Cumulant expansion + OS4 decorrelation. (5) No screening (no charged matter). Confinement chain: non-abelian -> delta < 1 -> K* > 0 -> area law.
**Depends on:** Theorems 24, 27, 18, OS4.

---

## Layer 8: The Identification

### Theorem 70: Mason Spectral Identification (`thm:ward_spectral`)
**Statement:** Mason's conditions hold at the equilibrium (||dbar(Phi)|| = 1.12, floor fires every step). DS dynamics IS Yang-Mills dynamics. DS spectral gap IS Yang-Mills mass gap. Identification is independent of F_Ward = 0 at equilibrium.
**Proof sketch:** Three Mason conditions verified at m* (not only during transients). Spectral continuity (Proposition 73): the identification holds throughout B.
**Depends on:** Theorem 39 (Mason), Theorem 44 (nijenhuis), Lemma 71 (rank-1).

### Theorem 71: Popov Identification (`thm:popov_identification`)

**This is the second most important result in the paper, after the unique product theorem (Theorem 10).** It is what makes "DS gap = YM gap" an algebraic identity rather than an analogy.

**Statement:** The Popov Hessian restricted to the DS fibre IS the DS transfer operator. dPhi|_{m*} = (d^2S/da^2)|_{a*, T_{m*}}.
**Proof sketch:** (1) Equilibrium form a* is rank-1 (all a*_alpha proportional, from rank-1 dbar(Phi) — Lemma 72). (2) [a*, delta_a] = 0 for all delta_a in T_{m*} (proportional matrices commute — this is the key step). (3) Popov Hessian L(delta_a) = dbar_{J'} delta_a + [a*, delta_a] reduces to dbar_{J'} delta_a = dPhi|_{m*}. (4) Eigenvalues match to 120 digits.
**Depends on:** Lemma 72, Theorem 44, Popov (2021).

**Why this matters:** The Popov action S[a] = (1/2pi i) int Omega wedge tr(a wedge dbar_{J'} a + 2/3 a^3) is the twistor-space action whose equations of motion ARE full Yang-Mills (Mason 2005, Popov 2021). Its second variation (Hessian) at the equilibrium determines the quantum spectrum. The rank-1 property of the Born floor — the fact that the floor adjusts exactly one real quantity — kills the commutator correction [a*, delta_a] identically, reducing the infinite-dimensional Popov Hessian on the DS fibre to the finite-dimensional DS transfer operator. Without rank-1, the commutator would be generically nonzero, the two operators would differ, and proving spectral equivalence would be a separate (harder) problem. The Born floor's minimality — the smallest possible non-integrable deformation — is what makes the identification exact.

### Lemma 72: Rank-1 Anti-Holomorphic Jacobian (`rem:rank1`)
**Statement:** At any floor-active point: rank(dbar(Phi)) = 1.
**Proof sketch:** DS is holomorphic (dbar(G) = 0), so dbar(Phi) = dbar(F) . G where F is the floor map. Floor adjusts one real quantity (|theta|); singleton rescaling proportional to theta correction. Every row of dbar(F) is a scalar multiple of the row vector (dt/dz-bar). Column times row = rank 1. Verified: sigma_2/sigma_1 = 10^{-201} at 300 digits.
**Depends on:** Theorem 48 (DS holomorphic), chain rule.

### Proposition 73: Fibre-Varying Bound (`prop:fibre_bound`)
**Statement:** The mass gap infimum over all CP^1-homogeneity modes k >= 0 is achieved at k=0 (fibre-constant).
**Proof sketch:** DS map commutes with CP^1 Fourier decomposition (fibre locality). Lambda_0 = 0.2829 for every k. Higher k = higher spin = heavier glueballs in confining theory.
**Depends on:** Lemma 72, Theorem 71, Theorem 69, Penrose transform.

### Proposition 74: Continuous Spectral Family (`prop:continuous_spectrum`)
**Statement:** Eigenvalues lambda_k(m) are continuous on B. |lambda_k(m)| < 1 for all m in B, all k >= 1.
**Proof sketch:** Phi smooth on B. Eigenvalues continuous in matrix entries. Along trajectories m_t -> m*, lambda_k(m_t) -> lambda_k(m*). Max |lambda_0| = 0.964 < 1 over 500 random samples.
**Depends on:** Smoothness of Phi, Theorem 23 (global convergence).

### Theorem 75: Universal Popov Identification (`thm:popov_universal`)
**Statement:** The Popov identification holds at every point of the K*=7/30 moduli curve M. rank(dbar(Phi)) = 1, [a*, delta_a] = 0, and Popov Hessian = DS transfer operator at every self-consistent equilibrium. Consequently, DS spectral gap = Yang-Mills mass gap for every compact simple non-abelian G.
**Proof sketch:** Lemma 72 proves rank-1 at ALL floor-active points (not equilibrium-specific). Rank-1 forces proportional a*_alpha, which commute with any perturbation. With [a*, delta_a] = 0, Popov Hessian reduces to transfer operator at every point on M. Verified: sigma_2/sigma_1 < 10^{-8} and ||[a*, .]||/||dPhi|| < 10^{-8} at 196 points spanning M.
**Depends on:** Lemma 72 (rank-1), Theorem 71 (Popov identification), Theorem 27 (universal correlation).

**Why this matters:** Theorem 71 proved the identification at one specific equilibrium (p_dom = 0.932). Theorem 75 extends it to the entire moduli curve — every self-consistent equilibrium at K*=7/30, regardless of which gauge group produces it. The mass gap existence is universal; the mass gap value depends on the gauge group through the evidence structure.

### Remark: The Moduli Curve (`rem:moduli`)
The K*=7/30 moduli curve M is a 1D family of self-consistent equilibria parametrised by the evidence ignorance phi in (0, 0.61). The conservation law (Sym^2(C^4) decomposition) provides the self-consistency constraint that selects K*=7/30 from a 2D family of fixed points. The gauge group selects a specific point on M.

**Universal on M** (CV < 0.5%): K*=7/30, sigma=-ln(23/30), C*det(M*)^2=8332/625, rank(dbar(Phi))=1, tr(dbar(Phi)) approx -0.491.

**Gauge-specific** (varies along M): lambda_0 (0.108 to 0.923), Delta (2.23 to 0.08), Penrose residue, SO(4) fractions, eigenvalue splitting.

**Key relationship:** lambda_0 approx 1.33*phi + 0.075 (linear in evidence ignorance). Mass gap existence is universal; mass gap value depends on how specific the gauge constraint is (smaller phi = larger gap).

**Verification:** 196 points computed (deepthink_moduli_curve.py), independently confirmed.

### Theorem 76: K*=7/30 is Independent of the Coupling (`thm:K_coupling_independent`)
**Statement:** K*=7/30 is independent of the lattice coupling beta, the gauge group G, and the evidence distribution e*. The position on the moduli curve M depends on beta and G; the value K* does not.
**Proof sketch:** K* = (H^2-H+1)/(H(H^2+1)) is derived from the conservation law, which is a consequence of the Sym^2(C^{H+1}) decomposition. The inputs are: (1) V = S tensor S = C^4 from the axioms of the self-consistent relational structure, (2) the decomposition V = Sym^2(S) + Lambda^2(S) = C^3 + C giving H=3, (3) eta = 4/27, (4) L_1 = 1, (5) channel counts from dim(Sym^2(C^4)) = 10. None reference beta, G, or e*. The verification K*(H^2+1) - eta*H^2 = 1 is an algebraic identity for all H. Beta enters through C_beta, which determines the position on M (the evidence), but every point on M has K*=7/30 because K* is determined by the state space decomposition upstream of the evidence.
**Depends on:** Theorem 24 (fixed-point conflict), Theorem 25 (Sym^2 characterisation), self-consistent relational structure axioms.

**Why this matters:** This theorem separates what is universal (K*=7/30, string tension, conservation law) from what is scale-dependent (the position on the moduli curve, the spectral gap Delta, the evidence structure). The coupling beta changes the evidence through C_beta, moving the system along M, but cannot change K* because K* is an algebraic invariant of V = C^4 = S tensor S.

---

## Layer 9: Universality (for any G)

### Theorem 78: Folding Invariance of the Mass Gap (`thm:folding_invariance`)
**Statement:** Let Γ̃ be a simply-laced Dynkin diagram with automorphism σ of order k, and Γ = Γ̃/σ the folded diagram. Then ρ(Γ) = ρ(Γ̃), hence Δ(Γ) = Δ(Γ̃).
**Proof sketch:** By Cartan matrix symmetry, the coupled equilibrium is σ-invariant. The Jacobian commutes with σ, so decomposes into eigenspaces. The folded algebra = symmetric projection. The spectral radius lies in the symmetric sector (verified for all 10 testable pairs to 10⁻⁹). The folding map is: C_n → A_{2n-1}, B_n → D_{n+1}, G₂ → D₄, F₄ → E₆.
**Depends on:** Uniqueness of coupled equilibrium, Dynkin diagram automorphisms.

### Theorem 79: ADE Classification of Mass Gaps (`thm:ade_classification`)
**Statement:** For any compact simple G, the mass gap Δ(G) at K*=7/30 depends only on the simply-laced covering type (ADE) of its Dynkin diagram. Δ > 0 for every compact simple G. The complete classification:
- A-series: SU(2) Δ=1.263, SU(3) Δ=0.689, ..., A_∞ Δ≈0.250
- D-series: SO(8)/SO(7)/G₂ Δ=0.118, D_∞ Δ≈0.115
- E-series: E₆/F₄ Δ=0.087 (tightest), E₇ Δ=0.098, E₈ Δ=0.095
**Proof sketch:** By folding invariance, reduce to ADE. Verify ρ < 1 for all 36 groups in the Killing-Cartan classification. The A-series is analysed by Fourier decomposition of the bulk Jacobian. The D-series bottleneck is the 3-valent fork node. The E-series has three discrete values.
**Depends on:** Theorem 78 (folding), coupled DS dynamics, Chevalley-Serre presentation.

### Remark: Topology, Not Spectrum (`rem:topology`)
A₅ and D₄ have identical adjacency spectral radius (√3) but Δ_{A₅}=0.297 vs Δ_{D₄}=0.118. The coupling changes the equilibrium itself: the 3-valent fork node shifts 53% from the single-site m*, with J_self eigenvalues [0.844, 0.820] vs chain [0.556, 0.546]. Graph topology (fork vs chain) controls the gap, not spectral properties.

**Gauntlet verification:** 36 groups tested at g = K* = 7/30. All pass (ρ < 1).
Safety margins (g_crit/K*): SU(3): 2.79×. F₄: 6.06×. E₈: 1.97×.

---

## Complete Dependency Chain

```
n(n-1)(n+2)(n+3) = 0, n=1
        |
  T(CP^1) = O(2), H = dim H^0 = 3
        |
  (H-1)^2 = H+1 at H=3 (four independent routes)
        |
  C^4, CP^3, Born floor 1/27, DS combination (unique product)
        |
  K* = 7/30 (conservation law = algebraic identity for all H)
        |
  DS + floor on CP^3
        |
  +---> dbar(Phi) != 0 (floor non-holomorphic) ---> N_J != 0
  |         |
  |     Mason (external): non-integrable J ---> full Yang-Mills
  |         |
  |     Rank-1 dbar(Phi) ---> [a*, delta_a] = 0 ---> Popov = DS
  |         |
  |     DS gap = Popov gap = YM gap (value depends on G)
  |     Universal Popov: holds at EVERY point on K*=7/30 curve
  |     K*=7/30 independent of beta (algebraic invariant of Sym^2(C^4))
  |
  +---> OS0: temperedness (compact B)
  +---> OS1: covariance (Sp(2)-equivariance of fibration)
  +---> OS2: reflection positivity (Koopman + commutativity)
  +---> OS3: regularity (finite-dimensional, polynomial)
  +---> OS4: cluster decomposition (lambda_0 = 0.2829 < 1)
  +---> Non-triviality (kurtosis 9 sigma, G_4 14.5 sigma)
        |
  OS reconstruction ---> Wightman QFT with mass gap
        |
  For any G: delta < 1 (universal correlation) + rho < 1 (36 groups)
        |
  Folding invariance: rho(fold(G)) = rho(G) -> ADE classification
        |
  Mass gap depends only on ADE type (A_n, D_n, E_6, E_7, E_8)
        |
  +---> (3,3) sector: graviton (conformal -> Einstein via OS2)
  +---> Minitwistor: NS regularity (conditional on floor preservation)
```

## External Dependencies

1. **Mason (2005)**: Non-integrable J -> full YM. Confirmed by Popov (2021).
2. **Osterwalder-Schrader (1973/75)**: OS axioms -> Wightman QFT.
3. **Gleason (1957)**: Born rule uniqueness.
4. **Maldacena (2011), Adamo-Mason (2014)**: Ghost-free conformal gravity = Einstein.
5. **Hitchin (1982)**: Minitwistor correspondence.
6. **BKM (1984)**: Bounded vorticity -> NS regularity (conditional result only).
7. **BMU (2014)**: Quantum reconstruction from four axioms.
8. **Hurwitz**: Classification of normed division algebras.

## What is NOT Proved

1. Asymptotic freedom / running coupling (K* = 7/30 is fixed).
2. OS axioms for curvature-derived Schwinger functions (proved for Born observables only).
3. Coupling g = K* from first principles (sufficient, not derived).
4. NS regularity unconditionally (conditional on Born floor preservation).
5. UV coincidence with conventional Yang-Mills at all scales.

## Computational Verification

267 scripts. Key cross-checks: K* to machine precision, lambda_0 to 500 digits, Delta to 50 sig fig, rho_{-1} to 120 digits, [R_0,R_1] traceless to 10^{-123}, tr(rho^2) = 3.70 (cup product), 36 gauge groups all rho < 1 (independently confirmed), kurtosis 9 sigma, G_4 14.5 sigma.

All code public: github.com/Resout/irreducible-twistor-geometry

## Author

J. R. Manuel, 2026
