# Irreducible Twistor Geometry

The framework of self-evident irreducible geometry forcing many, if not all, queries into a common and complete language. Found from a single premise, that thoughts are a shape.

## What This Is

257 verification scripts for the framework described in:

**"The Mass Gap from Enforced Uncertainty: Evidence Combination, Twistor Geometry, and Yang-Mills Theory"**
*J. R. Manuel, 2026*

A single algebraic identity — the Dempster-Shafer combination rule at H=3 hypotheses — forces a framework that produces:

- **Quantum theory** (BMU reconstruction on CP³)
- **Yang-Mills mass gap** (Δ = 1.2626, from K* = 7/30)
- **Einstein gravity** (conformal → Einstein via OS2)
- **Koide formula** (three lepton masses to 0.01%, zero free parameters)
- **Hierarchy tower** (cosmological constant to fine structure constant from S = 810/7)
- **Regularity of descended spin-1 field** (vorticity bounded by Born floor on minitwistor; the descended equation is a 4-component reaction-diffusion system, not Navier-Stokes)
- **Penrose residue bridge** (DS spectral gap = YM mass gap via non-integrable twistor geometry)

All from one equation. Zero free parameters. 21 exact algebraic identities (Class A). The state space CP³ is independently Penrose's twistor space for R⁴ — this is not assumed, it is forced by H=3.

## Constants

Everything derives from H = 3:

| Quantity | Value | Status |
|---|---|---|
| H | 3 | Forced (4 independent routes) |
| K* | 7/30 | Exact (conservation law) |
| Born floor | 1/27 | Exact (self-consistency) |
| Δ (spectral gap) | 1.2626 | Transcendental, = −ln(λ₀), 500-digit verified |
| λ₀ | 0.28291… | Algebraic over Q, degree >30 (500-digit verified) |
| S (instanton action) | 810/7 | Exact |
| S + 1/K* | 120 = 5! | Exact |

## The Identification: DS = Yang-Mills

The Born floor on CP³ induces a non-integrable almost complex structure J'. Mason's theorem (2005) guarantees that the resulting twistor geometry produces full Yang-Mills equations (both chiralities). The physical curvature F⁺ lives in the Penrose residue — the 1/ζ pole on twistor lines that the standard Ward construction discards. The Nijenhuis tensor N (from the Born floor) and the Maurer-Cartan form a = M⁻¹dM (from spatial variation) produce curvature only through their product N·a: the entanglement of entanglement.

Key results (confirmed by independent computation):
- ||ρ₋₁|| = 0.638 at uniform equilibrium (nonzero vacuum condensate)
- Response matrices [R₀, R₁] are genuinely su(2) (traceless to 10⁻¹²³)
- Correlator decay Δ_eff = 1.2626 matches transfer operator spectral gap to 6 figures
- The DS spectral gap IS the YM mass gap, bridged by the Penrose transform

## Structure

```
paper/          LaTeX source + PDF
computations/   257 standalone verification scripts (flat, no interdependencies)
```

Each `.py` file in `computations/` is a standalone verification of a specific claim. Key scripts:

**Core framework:**
- `spectral_gap_computation.py` — Δ = 1.2626 from transfer operator eigenvalue
- `conservation_law_derivation.py` — K* = 7/30 from Sym²(C⁴) decomposition
- `ds_gravity_k7_30.py` — exact K*=7/30 equilibrium, ∂̄Φ computation, SO(4) decomposition
- `abelian_bulletproof.py` — SU(2) vs U(1) discrimination at L=4,6,8,10,12
- `nontriviality_4point.py` — connected 4-point function at K*=7/30 equilibrium

**Yang-Mills identification:**
- `ward_penrose_integral.py` — Penrose contour integral, F⁺ from Nijenhuis residue
- `ward_highprecision.py` — 100-digit precision ratio verification
- `ward_reconstruction.py` — Ward reconstruction at equilibrium
- `ward_second_order_deepthink.py` — O(ε²) Birkhoff factorisation package

**Gauge theory:**
- `sun_from_su2.py` — N−1 root SU(2) subalgebras generate su(N)
- `sun_coupled_gap.py` — coupled mass gap verified for SU(N), N=2 through 6
- `wilson_area_law.py` — Wilson loop area law (confinement)

**Predictions:**
- `phi4_ds_extraction.py` — φ⁴ scalar field test (Prediction 2)
- `koide_from_crystal.py` — Koide formula from S₃ representation theory

**Integrity:**
- `verify_paper.py` — structural integrity check (environments, refs, citations, numerical claims)

Some scripts import from a `solver` module not included here. The scripts document what was computed; the results are in the paper.

## Paper

Available on Zenodo. The paper has 55 theorems, 0 conjectures, and compiles to ~90 pages. It documents its own errors and their corrections — including the Maurer-Cartan identity (M⁻¹dM is identically flat, not a gauge connection) and the degree-19 PSLQ mirage (artefact at 500 digits).

From Zenodo:

> The Mass Gap from Enforced Uncertainty: Evidence Combination, Twistor Geometry, and Yang-Mills Theory. Also happily sparring with whatever other question seems to be naturally constructing itself upon the path as we navigate it. On a personal note I'm starting to think quite a lot of the unsolved or unresolved mathematical and physics/computational geometric/theoretical framework forks in the modern scientific fields as they present them are actually very nearly just one, singular disagreement. I'll be continuing at the leisure of a laborer until it seems finished.

You are very welcome to use the computational techniques and discoveries here for any purpose. I humbly ask that if you do, you give me a reference — it'd make mum happy. And my sister embarrassed. And give the old man transient cognitive accolades to speak of at dinner.

Me no do math when shape do math good. Learn read shape fun. Math hard. Shape simple.

Man I really hope it's not only correct, but it's correct and Penrose sees this in time — that legend had it figured out way back when!

## Author

J. R. Manuel
