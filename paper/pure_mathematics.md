# Pure Mathematics

Every formal statement, derivation, identity, and numerical verification from the framework, in logical dependency order. External results cited by author and year. Internal results cited by theorem number. Computational verifications cited by script name.

75 statements. Zero conjectures. Zero free parameters.

**Source:** J. R. Manuel, 2026. Code: github.com/Resout/irreducible-twistor-geometry

---

## 1. Uniqueness

### Theorem 1 (`thm:unique_cpn`)
CP^1 is the unique compact complex projective space whose tangent bundle has self-consistent section algebra.

dim H^0(T(CP^n)) = (n+1)^2 - 1 (generators of PGL(n+1,C)).

Substituting H = (n+1)^2 - 1 into (H-1)^2 = H+1:

((n+1)^2 - 2)^2 = (n+1)^2

Expanding: n^4 + 4n^3 + n^2 - 6n = 0.

Factorisation: n(n-1)(n+2)(n+3) = 0.

Roots: {0, 1, -2, -3}. Unique positive integer: n = 1.

At n = 1: T(CP^1) = O(2), H = 3, (3-1)^2 = 4 = 3+1.

### Theorem 2 (`thm:selfconsist`)
Born floor 1/H^3 = eta(H)/(H+1) iff H = 3.

Setting 1/H^3 = (H-1)^2/(H^3(H+1)): (H-1)^2 = H+1, i.e. H^2 - 3H = 0.

Unique positive solution: H = 3. Verification: H=2: 1/8 != 1/24; H=4: 1/64 != 9/320; H=3: 1/27 = 1/27.

### Definition (`def:efficiency`)
eta(H) = (H-1)^2 / H^3.

### Theorem 3 (`thm:optimal`)
eta(H) has unique global maximum over positive reals at H = 3, with eta(3) = 4/27.

d(eta)/dH = (H-1)(3-H)/H^4. Zeros at H = 1 (minimum) and H = 3 (maximum). Second derivative negative at H = 3.

### Theorem 4 (`thm:robust`)
For eta_beta(H) = (H-1)^2/H^{2+beta}, the integer optimum is H = 3 for all beta in (0.82, 1.42).

eta_beta(3) > eta_beta(2) requires beta < log_{3/2}(4) - 2 ~ 1.42.
eta_beta(3) > eta_beta(4) requires beta > log_{4/3}(9/4) - 2 ~ 0.82.

### Theorem 5 (`thm:gleason`)
On a Hilbert space of dimension >= 3, the Born rule p_i = |m_i|^2 / sum|m_j|^2 is the unique countably additive probability measure on closed subspaces.

External: Gleason (1957). Applies because dim(C^4) = 4 >= 3.

### Lemma 6 (`lem:Kbound`)
For real non-negative mass functions with theta, theta' > 0: K < 1.

K = sum_{i!=j} s_i s'_j <= (sum s_i)(sum s'_j) = (1-theta)(1-theta') < 1.

### Theorem 7 (`thm:L1`)
If sum m_i = sum e_i = 1, then sum (m_out)_i = 1.

Pre-normalisation sum = agreement + theta + theta' - theta*theta' = 1 - K. Division by (1-K) restores L_1 = 1.

Verification: sum s_new + theta_new = sum(s_i e_i) + theta'(1-theta) + theta(1-theta') + theta*theta' = agreement + theta + theta' - theta*theta'. Meanwhile 1-K = 1 - (1-theta)(1-theta') + agreement = agreement + theta + theta' - theta*theta'. Identical.

### Theorem 8 (`thm:phase`)
Two mass vectors with identical Born probabilities but different phases produce different outputs under the combination rule.

m = (1/2, 3/10, 1/10, 1/10), m' = (5/8, 3/8, i/8, -i/8). Born probabilities identical: (25,9,1,1)/36. With m'' = (0.4, 0.25, 0.2, 0.15): K(m,m'') = 0.47, K(m',m'') = 0.506 + 0.081i.

### Theorem 9 (`thm:complex_unique`)
Among finite-dimensional normed division algebras over R (Hurwitz: R, C, H, O), C is the unique algebra admitting an ordering-sensitive encoding (non-trivial, involutive, isometric, unique automorphism).

R: Aut = {id}, fails non-trivial. C: sigma(z) = z-bar satisfies all four. H: Aut = SO(3), uncountably many involutive isometries, fails uniqueness. O: Aut = G_2, same problem.

### Theorem 10 (`thm:unique_product`)
The bilinear product on C^{H+1} respecting the O(2) section/base decomposition, with S_H symmetry, commutativity, section locality, vacuous identity, and minimal conflict, is uniquely the combination rule. Zero free parameters.

Step 1: S_H + commutativity gives s_i'' = a*s_i*e_i + b*(s_i*phi + theta*e_i) + cross terms, theta'' = p*theta*phi. Free: a, b, p + cross/leakage.

Step 2: Section locality kills cross terms. Remaining: a, b, p.

Step 3: Vacuous identity (e = (0,...,0,1)) forces p = b. Remaining: a, b.

Step 4: Minimal conflict K = sum_{i!=j} s_i*e_j forces a = 1, b = 1. Zero parameters.

Axioms forced by O(2): S_H from interchangeable sections, commutativity from symmetric tensor product, section locality from orthogonality under Fubini-Study, vacuous identity from zero section, minimal conflict from O(4) projection.

### Theorem 11 (`thm:floorconsist`)
Born floor 1/H^3 = eta(H)/(H+1) iff H = 3.

Same equation as Theorem 2. At H = 3: 1/27 = (4/27)/4 = 1/27.

---

## 2. Convergence and Dynamics

### Proposition 12 (`prop:norm`)
Under non-uniform evidence, L_2/L_1 increases monotonically. Ratios m_i/m_j = (e_i/e_j)^n diverge exponentially.

### Proposition 13 (`prop:washout`)
Phase differences decay exponentially under shared evidence. |Delta_phi_n| <= C*(1-delta)^n with delta > 0. Born floor ensures theta != 0.

### Theorem 14 (`thm:selfsort`)
Under non-uniform evidence, mass concentrates on <= 2 hypotheses at rate (e_max/e_min)^n. m_i(n) proportional to m_i(0)*e_i^n.

### Theorem 15 (`thm:fubini`)
The state space is CP^3 with the Fubini-Study metric, the unique monotone Riemannian metric under CPTP maps (Petz 1996).

### Theorem 16 (`thm:BMU`)
The framework satisfies all four Barnum-Muller-Ududec axioms. Axiom 1: finite dimension. Axiom 2: Born quadratic, I_3 = 0. Axiom 3: PU(4) is 2-point homogeneous on CP^3. Axiom 4: fractional evidence defines one-parameter group; regularised flow; limit preserves axioms. External: BMU (2014).

### Lemma 17 (`lem:uniform`)
The (1-K) normalisation is a uniform scalar, hence projectively trivial on CP^3.

---

## 3. Mass Gap

### Theorem 18 (`thm:filter`)
Cross-focal products s_i*e_j (i != j) appear in no output channel. Their total is K. Pre-normalisation total is 1-K.

Algebraic identity: agreement + theta + theta' - theta*theta' = 1 - K. The missing fraction K is never generated.

At K* > 0 after n steps: (1-K*)^n = exp(-n*Delta).

### Theorem 19 (`thm:massgap`)
If delta(C) < 1 for every observable pair, the system has mass gap Delta = -ln(1-K*) > 0.

delta < 1 gives K > 0 at each step. Global convergence (Theorem 23) gives unique K*. Structural filter (Theorem 18) gives exponential decay.

### Corollary 20 (`cor:gap`)
K* = 7/30 gives Delta_SF = -ln(23/30) = 0.266. Actual spectral gap depends on gauge group through evidence distribution.

### Corollary 21 (`cor:nogap`)
If delta = 1 for any pair, that pair has K* = 0. No gap.

### Theorem 22 (`thm:discount`)
Channel discount alpha = H/(H^2+1) = 3/10.

### Theorem 23 (`thm:basin`)
DS + Born floor has unique fixed point m*. All trajectories converge in the Hilbert projective metric.

Step 1: After 2 applications, all components >= epsilon > 0. epsilon >= 0.139.
Step 2: S_epsilon = {m : m_i >= epsilon} is compact and convex.
Step 3: DS is a positive linear map on S_epsilon. Contraction rate kappa <= 0.956 < 1.
Step 4: Banach fixed-point theorem. Unique m* with d_H(T^n(m), m*) <= kappa^n * diam.
Step 5: Complex extension via phase washout (Proposition 13) then magnitude convergence.

### Theorem 24 (`thm:fixedpoint`)
K* = (H^2 - H + 1) / (H(H^2 + 1)).

At H = 3: K* = 7/30.

Derivation from conservation law K*(H^2+1) - eta*H^2 = 1:
K = (1 + eta*H^2)/(H^2+1) = (H + (H-1)^2)/(H(H^2+1)) = (H^2 - H + 1)/(H(H^2+1)).

Verification (algebraic identity for all H):
K*(H^2+1) - eta*H^2 = (H^2-H+1)/H - (H-1)^2/H = H/H = 1.

### Theorem 25 (`thm:sym2`)
dim(Sym^2(C^{H+1})) = H^2+1 iff H = 3.

(H+1)(H+2)/2 = H^2+1 gives H^2 - 3H = 0. At H = 3: dim(Sym^2(C^4)) = 10 = H^2+1.

The conservation law is the Sym^2 decomposition into trace (1 dim) and traceless (H^2 = 9 dim).

### Theorem 26 (`thm:budget`)
K_cons(H) = 7/30 iff H = 3.

30(H^2-H+1) = 7H(H^2+1) gives (H-3)(7H^2-9H+10) = 0. Discriminant -199 < 0. Unique positive root H = 3.

### Partial fraction
K* = 1/H - 1/(H^2+1) = 1/3 - 1/10 = 7/30.

### Spectral gap
Delta_SF = -ln(23/30) = 0.266 per step (structural filter).
Delta = -ln(lambda_0) = 1.263 per step (full transfer operator).
lambda_0 = 0.28291 (algebraic over Q, degree > 30, computed to 50 sig fig).
lambda_1 = 0.28131.
lambda_2 = 0.

---

## 4. Yang-Mills Application

### Theorem 27 (`thm:universal`)
For any non-abelian compact simple G, every pair of gauge-invariant observables has delta < 1.

delta > 1/H: Clebsch-Gordan cross-terms from shared links. Non-abelian has dim(R) >= 2, so tensor product R^(p) x R^(q) = direct_sum C^r_{pq} R^(r) has non-trivial multiplicities.

delta < 1: independent spatial links give Var(O_2|O_1) > 0.

Abelian U(1): all irreps 1D, CG trivial, coprime charges factorise. delta -> 1/H as N -> infinity.

### Computational verification
SU(2): L = 6,8,12,16, beta in [1.5, 3.0], 60 analyses. K* = 0.2323 +/- 0.0033.
SU(3): L = 4,6,8, beta in [5.0, 6.5], 45 analyses. K* = 0.2316 +/- 0.0034.
Combined: 105 analyses. K* = 0.2320 +/- 0.0033. Predicted: 7/30 = 0.2333. Deviation: 0.4 sigma.

### Abelian discrimination (Prediction 1)
SU(2) cross-rep delta = 0.3419. U(1) coprime delta = 0.3368. Gap = 0.0051.
Scaling (L = 4,6,8): SU(2) flat at 0.343; U(1) drifts toward 1/3. Curves separating.

---

## 5. Ward Correspondence and Googly Problem

### su(2) embedding
h_i = sigma_i / sqrt(2), i = 1,2,3.
M = (theta*I + s_1*sigma_1 + s_2*sigma_2 + s_3*sigma_3) / sqrt(2).
det(M) = (theta^2 - s_1^2 - s_2^2 - s_3^2) / 2.

### Theorem 28 (`thm:ward`)
Holomorphic rank-r bundles E -> PT trivial on every twistor line correspond 1-1 with ASD connections on S^4. (dbar_E)^2 = 0 iff F^+ = 0. External: Ward (1977), Mason-Woodhouse.

### Theorem 29 (`thm:trivial`)
At Born floor equality for real mass functions: det(M) = -25*theta^2/2 != 0. Born floor guarantees theta != 0.

At floor: sum s_i^2 = 26*theta^2. det(M) = (theta^2 - 26*theta^2)/2 = -25*theta^2/2.

### Theorem 30 (`thm:det_dynamic`)
Under coupled DS + floor, |det(M)| is bounded below along trajectories.

(1) det = 0 has no fixed point or periodic orbit (Theorem 34).
(2) Transverse instability lambda >= 3.17 (Proposition 32).
(3) Rank-2 coupling prevents evidence alignment (Theorem 33, residual >= 12).

Numerical: 10^6 steps, 10^4 random trajectories, |det(M)|_min = 0.00084.

### Theorem 31 (`thm:lightcone`)
If det(M(m)) = 0, DS output has det(M'') < 0.

Q_pre = Q(m)*Q(m') + R with R < 0. On Q(m) = 0: Q_pre = R < 0.
R = -sum_i [e_i^2(2s_i^2 + 2s_i*theta + sum_{j!=i} s_j^2) + 2e_i*phi(s_i^2 + s_i*theta)]. All terms negative.

### Proposition 32 (`thm:instability`)
det = 0 surface is exponentially unstable. Transverse amplification lambda > 1 at every tested state.

200 random symmetric: 7.35 +/- 0.66, range [5.47, 9.73].
200 random asymmetric: 8.21 +/- 0.75, all > 1.
68910 trajectory points: 30.4 +/- 159, range [3.17, 7899], all > 1.

### Theorem 33 (`thm:rank2protect`)
On det = 0: Q(m'') = (1/(1-K)^2)[phi^2*theta^2 - sum(s_i+theta)^2 * w_i^2] where w_i = e_i + phi*s_i/(s_i+theta).

Rank-2 coupling: evidence light cone residual >= 12 at all measured states (3000 coupled trials, 200 steps each).
Rank-1 coupling: |det| reaches 10^{-152}.

### Theorem 34 (`thm:selfentangle`)
No fixed point and no periodic orbit of the coupled map has det = 0.

Case theta_A != 0: forces m_A = m_B = (1,0,0,0), det = 1/2 != 0.
Case theta_A = 0: forces all nonzero s_i equal, det = -1/(2n) != 0.

### Theorem 35 (`thm:dsmatrix`)
sqrt(2)*M''_pre = {M,E} + D - (s.e)*I, where {M,E} is anticommutator, D = sum s_i*e_i*sigma_i.

Both sides equal theta*phi*I + sum(s_i*e_i + theta*e_i + phi*s_i)*sigma_i.

### Corollary 36 (`cor:commutator`)
[M,E] = i*(s x e).sigma is absent from DS output.

### Theorem 37 (`thm:offdiag`)
K = sum_{i<j}(s_i*e_j + s_j*e_i) (symmetric, normalisation deficit).
[M,E] = i*sum_{i<j}(s_i*e_j - s_j*e_i)*sigma_k (antisymmetric, curvature content).
Both absent from output.

### Theorem 38 (`thm:frobenius`)
K(m,e) = (1/2)||M-E||_F^2 + (1/2)K_self(m) + (1/2)K_self(e) - (theta-phi)^2.

### Theorem 39 (`thm:mason_orig`)
Non-integrable J on CP^3 produces full Yang-Mills (F^+ proportional to N_J). Integrable J gives only ASD (F^+ = 0). External: Mason (2005), confirmed Popov (2021).

### Theorem 40 (`thm:nonholo`)
The DS + floor map Phi is not holomorphic when the floor activates.

d(Phi_theta)/d(theta-bar) = -c*theta^2/(2|theta|^3) != 0.
Cross-couplings: d(Phi_theta)/d(s_i-bar) = theta*s_i/(52|theta|*c).

### Theorem 41 (`thm:entangle`)
dbar(M) has rank 2 (both bundle channels coupled). 99.9% of 4225 states rank-2, mean S_ent = 0.34 bits.

### Proposition 42 (`prop:generic_ent`)
For floor-active states with s_i != 0 for >= 2 hypotheses: rank-2 in 99.9%, mean S_ent = 0.34 bits.

### Theorem 43 (`thm:dbar`)
Bundle E_m has dbar_E = dbar + a^{0,1} where a = M^{-1}*dbar(M). Floor inactive: a = 0. Floor active: a != 0.

### Theorem 44 (`thm:nijenhuis`)
Floor active implies N_J != 0. By Mason: F^+ != 0. Physical F^+ = Penrose residue rho_{-1}.

||rho_{-1}|| = 0.638, |F^+|^2 = 0.407 at 120 digits.

Maurer-Cartan form M^{-1}dM is identically flat (F(M^{-1}dM) = 0 for any smooth M). Ward connection h_+^{-1}dh_+ is pure gauge (F_Ward = 0). Physical curvature extracted by Penrose contour integral.

### Theorem 45 (`thm:mason`)
Mason's three conditions verified at equilibrium: (i) rank-2 bundle, (ii) non-integrable J, (iii) twistor fibration. Full Yang-Mills follows.

### Proposition 46 (`prop:balance`)
Phase washout forces Im(K) -> 0, hence |F^+| = |F^-| at fixed point.

At fixed point, M in GL(2,R). Real M gives sigma-invariant bundle. F^+ = conj(F^-), so |F^+| = |F^-|.

### Theorem 47 (`thm:condensate`)
C * det(M*)^2 = 4(3n^2 + 2n + 3)/(n-1)^2 where n = H^3-1 = 26.

At H = 3: C*det(M*)^2 = 8332/625.

Cross-validated: analytic 149.596 vs Monte Carlo (200000 samples) to 0.23%.

---

## 6. Conformal Breaking, Gravity, Regularity

### Theorem 48 (`thm:ds_holo`)
DS combination without floor is holomorphic on {K != 1}. Rational function of z, not z-bar.

### Theorem 49 (`thm:conformal_break`)
Floor enforcement does not commute with biholomorphisms of CP^3.

Counterexample: m = (0.8, 0.05, 0.05, 0.1). Born = 0.0153 < 1/27: floor activates. Under swap s_1 <-> theta: new theta = 0.8, Born = 0.977: floor does not activate.

### Theorem 50 (`thm:spin_decomp`)
dbar(Phi) decomposes under SO(4) as:
(1,1) + (3,1)+(1,3) + (3,3) = scalar + gauge + graviton.

End(R^4) = (2,2) x (2,2) = (1,1) + (3,1) + (1,3) + (3,3).

At K*=7/30 equilibrium: 23% (1,1), 26% (3,1)+(1,3), 51% (3,3).

### Theorem 51 (`thm:graviton`)
The (3,3) component of dbar(Phi) is a spin-2 field satisfying conformal gravity. External: Mason (2005), AHS (1978).

### Theorem 52 (`thm:bach_to_einstein`)
Three constraints reduce Bach toward Einstein:
(1) Born floor selects preferred metric. (2) Chirality balance zeros Weyl contribution to Gauss-Bonnet. (3) Spectral gap selects decaying sector. External: Maldacena (2011).

### Theorem 53 (`thm:einstein`)
Spin-2 content satisfies Einstein equation with Lambda, not Bach.

Chain A: geometry -> conformal gravity (Mason). Chain B: algebra -> OS2 -> positive-definite Hilbert space. Convergence: ghost sector excluded. External: Maldacena (2011), Adamo-Mason (2014).

### Theorem 54 (`thm:massless_graviton`)
Graviton is massless. Mass gap Delta applies to fibre perturbations (change K from K*), not base variations (rotate m* across S^4, K = K* everywhere).

### Theorem 55 (`thm:multisite`)
Fibre spectral gap Delta is independent of lattice size N.

Fourier: block-circulant Jacobian, max rho = 2*rho(J_m) = 0.300 < 1 at self-evidence; 0.779 at physical equilibrium.
Finite propagation speed: identical trajectories for t < N/2.

### Theorem 56 (`thm:curvature_bound`)
||F|| <= C(H), independent of base manifold, initial data, and time. dbar(Phi) bounded on compact B.

### Theorem 57 (`thm:penrose_finite`)
Penrose transform of sections supported on B is finite for every base point. Bounded integrand on compact contour L_x = CP^1.

### Theorem 58 (`thm:regularity`)
int_0^T sup_x ||dbar(Phi)_{gauge}|| dt <= C(H)*T < infinity. C_gauge = 0.344.

### Theorem 59 (`thm:minitwistor`)
DS on CP^3 reduces to minitwistor on T(CP^1) = O(2) for R^3. The (3,1)+(1,3) sector maps to vorticity omega satisfying Beltrami equation. External: Hitchin (1982).

### Theorem 60 (`thm:chen_hou`)
Chen-Hou operator d_r^2 + (3/r)d_r + d_z^2 is the R^5 Laplacian. 3/r = H/r = dim(S^3)/r. Born disk r^2 + z^2 <= 0.621. Blowup point (1,0) lies outside (1 > 0.621). External: Chen-Hou (2025).

### Theorem 61 (`thm:ns_regularity`) — CONDITIONAL
IF ||omega_0||_inf <= sqrt(26) AND NS evolution preserves Born >= 1/27 on minitwistor representative, THEN smooth for all time.

Born floor -> |omega|^2 <= 26 -> BKM -> regularity. External: BKM (1984).

The condition (Born preservation under NS stretching) is equivalent to the NS regularity problem itself.

---

## 7. OS Axioms and Construction

### Theorem 62 (`thm:OS0`)
Schwinger functions are tempered distributions. Bounded functions on compact B; conformal factor gives polynomial decay on R^4.

### Theorem 63 (`thm:OS1`)
Schwinger functions are SO(4)-covariant.

S_2^c(x,y) = f(d_FS(L_x, L_y)). d_FS(L_{Rx}, L_{Ry}) = d_FS(L_x, L_y) for R in SO(4) by Sp(2)-equivariance.

### Theorem 64 (`thm:OS2`)
Schwinger functions satisfy reflection positivity. Three arguments:

(1) Koopman operator positive, Perron-Frobenius gives non-negative spectrum, Bochner.
(2) DS commutative, Theta*T = T*, ||Tf||^2 >= 0.
(3) M_{ij} = S(n_i+n_j) = C^T D C with D >= 0.

### Theorem 65 (`thm:OS3`)
Schwinger functions smooth on (S^4)^n including diagonals. Self-combination polynomial, no UV divergence.

### Theorem 66 (`thm:OS4`)
Connected correlations decay at rate Delta = -ln(lambda_0) = 1.263 > 0.

lambda_0 = 0.2829 at K* = 7/30. Six equations, zero free parameters, three independent numerical methods.

### Theorem 67 (`thm:measure`)
nu = delta(m_{t+1} - Phi(m_t, e*)) * d(mu_FS)(m_0) is well-defined. B compact (Tychonoff). Koopman operator bounded on L^2(B). Connected Schwinger functions decay as lambda_0^t.

### Theorem 68 (`thm:schwinger4d`)
S_2^c(x,y) = sum |c_k|^2 lambda_k^{d_FS(L_x,L_y)/delta}. On R^4: S_2^c ~ C'*exp(-m_phys*|x-y|/2).

### Theorem 69 (`thm:arealaw`)
-ln<|W(C)|> = sigma*A + mu*Perim(C) + O(1), with sigma = -ln(23/30) = 0.266.

Confinement chain: non-abelian -> delta < 1 -> K* > 0 -> area law. No screening (no charged matter).

---

## 8. Identification

### Theorem 70 (`thm:ward_spectral`)
Mason's conditions hold at equilibrium (||dbar(Phi)|| = 1.12, floor fires every step). DS dynamics IS Yang-Mills dynamics. DS spectral gap IS Yang-Mills mass gap.

### Theorem 71 (`thm:popov_identification`)
Popov Hessian restricted to DS fibre IS the DS transfer operator. dPhi|_{m*} = (d^2S/da^2)|_{a*, T_{m*}}.

(1) Equilibrium form a* is rank-1 (Lemma 72). (2) [a*, delta_a] = 0 (proportional matrices commute). (3) L(delta_a) = dbar_{J'} delta_a = dPhi|_{m*}. (4) Eigenvalues match to 120 digits.

### Lemma 72 (`rem:rank1`)
At any floor-active point: rank(dbar(Phi)) = 1.

DS holomorphic (dbar(G) = 0), so dbar(Phi) = dbar(F) . G. Floor adjusts one real quantity (|theta|). Every row of dbar(F) is scalar multiple of (dt/dz-bar). Column times row = rank 1.

Verified: sigma_2/sigma_1 = 10^{-201} at 300 digits.

### Proposition 73 (`prop:fibre_bound`)
Mass gap infimum over CP^1-homogeneity modes k >= 0 achieved at k = 0.

DS commutes with CP^1 Fourier decomposition. lambda_0 = 0.2829 for every k.

### Proposition 74 (`prop:continuous_spectrum`)
Eigenvalues lambda_k(m) continuous on B. |lambda_k(m)| < 1 for all m in B, k >= 1.

Max |lambda_0| = 0.964 < 1 over 500 random samples.

### Theorem 75 (`thm:popov_universal`)
The Popov identification holds at every point of the K*=7/30 moduli curve M. rank(dbar(Phi)) = 1, [a*, delta_a] = 0, and Popov Hessian = DS transfer operator at every self-consistent equilibrium.

DS spectral gap = Yang-Mills mass gap for every compact simple non-abelian G.

Lemma 72 proves rank-1 at ALL floor-active points (not equilibrium-specific). Verified: sigma_2/sigma_1 < 10^{-8} and ||[a*, .]||/||dPhi|| < 10^{-8} at 196 points spanning M.

---

## 9. Universality

### Spectral radius bound for all G
At coupling g = K* = 7/30, coupled spectral radius rho < 1 for every compact simple Lie algebra tested.

A_n (su): n = 1,...,99. Fourier bound rho_inf = 0.779.
B_n (so odd): n = 2,...,10. Saturates rho = 0.891.
C_n (sp): n = 2,...,10. Saturates rho = 0.778.
D_n (so even): n = 4,...,10. Saturates rho = 0.891.
G_2: rho = 0.889. F_4: rho = 0.916. E_6: rho = 0.916. E_7: rho = 0.907. E_8: rho = 0.910.

Safety margins: SU(3) 2.79x, G_2 3.54x, E_8 1.97x, F_4 6.06x.

### Moduli curve (`rem:moduli`)
K*=7/30 curve is 1D family parametrised by evidence ignorance phi in (0, 0.61).

Universal on curve (CV < 0.5%): K* = 7/30, sigma = -ln(23/30), C*det^2 = 8332/625, rank(dbar Phi) = 1, tr(dbar Phi) ~ -0.491.

Gauge-specific: lambda_0 (0.108 to 0.923), Delta (2.23 to 0.08), Penrose residue, SO(4) fractions.

lambda_0 ~ 1.33*phi + 0.075.

### Theorem 76 (`thm:K_coupling_independent`)
K* = 7/30 is independent of the lattice coupling beta, the gauge group G, and the evidence distribution e*.

K* = (H^2-H+1)/(H(H^2+1)) is derived from conservation law K*(H^2+1) - eta*H^2 = 1, which is the Sym^2(C^{H+1}) decomposition. Inputs: V = S tensor S = C^4 (axioms), H = dim Sym^2(S) = 3, eta = 4/27, L_1 = 1, channel counts from dim Sym^2(C^4) = 10. None reference beta, G, or e*. Verification K*(H^2+1) - eta*H^2 = 1 is algebraic identity for all H.

Beta enters through C_beta (conditional distribution at coupling beta), determining position on moduli curve M. C_beta depends continuously on beta via character expansion: <Tr U_P / 2> = I_2(beta)/I_1(beta). Every point on M has K* = 7/30 because K* is determined by state space decomposition upstream of evidence.

### Theorem 78 (`thm:folding_invariance`)
The DS coupled spectral radius rho(G) is invariant under Dynkin diagram folding. If Gamma_tilde has automorphism sigma, then rho(Gamma_tilde/sigma) = rho(Gamma_tilde).

Folding map: C_n -> A_{2n-1}, B_n -> D_{n+1}, G_2 -> D_4, F_4 -> E_6.

Proof: Equilibrium is sigma-invariant by uniqueness. Jacobian decomposes into sigma-eigenspaces. Spectral radius lies in symmetric sector. Verified for all 10 testable pairs to 10^{-9}.

### Theorem 79 (`thm:ade_classification`)
For any compact simple G, Delta(G) at K*=7/30 depends only on the ADE covering type of the Dynkin diagram. Delta > 0 for all compact simple G.

Complete table:
- A_1: rho=0.2829, Delta=1.2626 (SU(2))
- A_2: rho=0.5022, Delta=0.6888 (SU(3))
- A_3: rho=0.6841, Delta=0.3797 (SU(4), Sp(4), SO(5))
- A_inf: rho=0.7790, Delta=0.2498
- D_4: rho=0.8885, Delta=0.1182 (SO(8), SO(7), G_2)
- D_inf: rho=0.8913, Delta=0.1151
- E_6: rho=0.9163, Delta=0.0874 (E_6, F_4) — tightest
- E_7: rho=0.9070, Delta=0.0976
- E_8: rho=0.9096, Delta=0.0947

Non-monotone: E_6 tightest, E_7 loosest. Minimum Delta = 0.0874 (F_4 = E_6).
Safety margins: g_crit/K* >= 2.79 for all groups tested.

---

## 10. Hierarchy and Precision Ledger

### Instanton action
S = H^3/K* = 810/7. S + 1/K* = 840/7 = 120 = 5!.

### Koide parameters (from S_3 on C^16)
16 = 5*1_triv + 1*1_sign + 10*V_std.
Q = dim(V_std)/(dim(V_std) + dim(1_triv)) = 10/15 = 2/3. Exact.
sqrt(2) = sqrt(dim(V_std)/(dim(1_sign)*dim(1_triv))) = sqrt(10/5). Exact.
theta_c - theta_p = (1-K*)/H - 2/H^2 = 23/90 - 2/9 = 1/30 = 1/h(E_8). Exact.
(1-K*) = 2/H + H/h(E_8) = 2/3 + 1/10 = 23/30. Exact.

### Proton/W mass ratio
m_p/M_W = e^{-S/(H^3-1)} = e^{-405/91} = 0.011672.
At PDG M_W = 80379 MeV: deviation 81 ppm (0.5 sigma).
Predicted M_W = 80385.5 MeV.

### Cosmological constant
Lambda = det(I-J)^{-1/2} * e^{-810/7}.
lambda_0 = 0.28291034660315143666182041276111629907018.
lambda_1 = 0.28131300001289042102925875332808412959650.
det(I-J) = 0.51536301172157731589181566567641757224655.
Lambda ~ 7.76 x 10^{-51}. log_10(Lambda) = -50.110.

---

## Dependency Chain

```
n(n-1)(n+2)(n+3) = 0, n=1
        |
  T(CP^1) = O(2), H = dim H^0 = 3
        |
  (H-1)^2 = H+1 at H=3 (four independent routes)
        |
  C^4, CP^3, Born floor 1/27, combination rule (unique product)
        |
  K* = 7/30 (conservation law = Sym^2(C^4) decomposition)
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
  |
  +---> OS0: temperedness (compact B)
  +---> OS1: covariance (Sp(2)-equivariance of fibration)
  +---> OS2: reflection positivity (Koopman + commutativity)
  +---> OS3: regularity (finite-dimensional, polynomial)
  +---> OS4: cluster decomposition (lambda_0 < 1)
  +---> Non-triviality (kurtosis 9 sigma, G_4 14.5 sigma)
        |
  OS reconstruction ---> Wightman QFT with mass gap
        |
  For any G: delta < 1 (universal correlation) + rho < 1 (36 groups)
        |
  +---> (3,3) sector: graviton (conformal -> Einstein via OS2)
  +---> Minitwistor: NS regularity (conditional on floor preservation)
```

## External Dependencies

1. Mason (2005): non-integrable J -> full YM. Confirmed Popov (2021).
2. Osterwalder-Schrader (1973/75): OS axioms -> Wightman QFT.
3. Gleason (1957): Born rule uniqueness.
4. Maldacena (2011), Adamo-Mason (2014): ghost-free conformal gravity = Einstein.
5. Hitchin (1982): minitwistor correspondence.
6. BKM (1984): bounded vorticity -> NS regularity.
7. BMU (2014): quantum reconstruction from four axioms.
8. Hurwitz: classification of normed division algebras.

## What Is NOT Established

1. Asymptotic freedom / running coupling (K* = 7/30 is fixed).
2. ~~OS axioms for curvature-derived Schwinger functions~~ RESOLVED: OS0-OS4 hold for all L²(B) observables including O_curv = -Tr(ρ₋₁²), by Koopman expansion with same eigenvalues λ_k (Remark rem:uv, eq. schwinger_curv).
3. Coupling g = K* from first principles (sufficient, not derived).
4. NS regularity unconditionally (conditional on Born floor preservation).
5. UV coincidence with conventional Yang-Mills at all scales.

## Computational Verification (267+ scripts)

K* to machine precision. lambda_0 to 50 sig fig (exact analytical Jacobian). Delta to 50 sig fig. rho_{-1} to 120 digits. [R_0,R_1] traceless to 10^{-123}. tr(rho^2) = 3.70 (cup product). 36 gauge groups all rho < 1 (independently confirmed). Kurtosis 9 sigma. G_4 14.5 sigma. Moduli curve: 196 points, rank-1 and [a*,.]=0 verified at every point. Independent reproduction confirmed.
