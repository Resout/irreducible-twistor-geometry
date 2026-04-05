# Irreducible Twistor Geometry

The framework of self evident irreducible geometry forcing many, if not all, queries into a common and complete language. Found from a single premise, that thoughts are a shape.

## What This Is

230 verification scripts for the framework described in:

**"The Mass Gap from Enforced Uncertainty: Evidence Combination, Twistor Geometry, and Yang-Mills Theory"**
*J. R. Manuel, 2026*

A single algebraic identity — the Dempster-Shafer combination rule at H=3 hypotheses — forces a framework that produces:

- **Quantum theory** (BMU reconstruction on CP³)
- **Yang-Mills mass gap** (Δ = 1.263, from K* = 7/30)
- **Einstein gravity** (conformal → Einstein via OS2)
- **Koide formula** (three lepton masses to 0.01%, zero free parameters)
- **Hierarchy tower** (cosmological constant to fine structure constant from S = 810/7)
- **Regularity of descended spin-1 field** (vorticity bounded by Born floor on minitwistor; the descended equation is a 4-component reaction-diffusion system, not Navier-Stokes)

All from one equation. Zero free parameters. 20 exact algebraic identities (Class A). The state space CP³ is independently Penrose's twistor space for R⁴ — this is not assumed, it is forced by H=3.

## Constants

Everything derives from H = 3:

| Quantity | Value | Status |
|---|---|---|
| H | 3 | Forced (4 independent routes) |
| K* | 7/30 | Exact (conservation law) |
| Born floor | 1/27 | Exact (self-consistency) |
| Δ (spectral gap) | 1.2626 | Transcendental, = -ln(λ₀), 50 sig fig |
| λ₀ | 0.28291 | Algebraic, degree 19 over Q, Gal=S₁₉ |
| S (instanton action) | 810/7 | Exact |
| S + 1/K* | 120 = 5! | Exact |

## Structure

```
paper/          LaTeX source + PDF
computations/   230 standalone verification scripts (flat, no interdependencies)
```

Each `.py` file in `computations/` is a standalone verification of a specific claim. Key scripts:

- `spectral_gap_computation.py` — Δ = 1.263 from transfer operator eigenvalue
- `conservation_law_derivation.py` — K* = 7/30 from Sym²(C⁴) decomposition
- `ds_gravity_k7_30.py` — exact K*=7/30 equilibrium, ∂̄Φ computation, SO(4) decomposition
- `abelian_bulletproof.py` — SU(2) vs U(1) discrimination at L=4,6,8,10,12
- `sun_from_su2.py` — N-1 root SU(2) subalgebras generate su(N)
- `sun_coupled_gap.py` — coupled mass gap verified for SU(N), N=2 through 6
- `nontriviality_4point.py` — connected 4-point function at K*=7/30 equilibrium
- `phi4_ds_extraction.py` — φ⁴ scalar field test (Prediction 2)
- `koide_from_crystal.py` — Koide formula from S₃ representation theory
- `verify_paper.py` — structural integrity check (environments, refs, citations, numerical claims)

Some scripts import from a `solver` module not included here. The scripts document what was computed; the results are in the paper.

## Paper

Available on Zenodo. The paper has 53 theorems, 0 conjectures, and compiles to ~90 pages.

From Zenodo:

The Mass Gap from Enforced Uncertainty: Evidence Combination, Twistor Geometry, and Yang-Mills Theory. Also happily sparring with whatever other question seems to be naturally constructing itself upon the path as we navigate it. On a personal note I'm starting to think quite a lot of the unsolved or unresolved mathematical and physics/computational geometric/theoretical framework forks in the modern scientific fields as they present them are actually very nearly just one, singular disagreement. I'll be continuing at the leisure of a laborer until it seems finished.

Disclaimer, you are very, very welcome to use the computational techniques and discoveries unfolded and demonstrated/described here in this/these paper(s) for any of the borderline infinite useful utilities it could be sicced on. I very humbly ask that if you do, then you might also give me a reference, it'd make mum happy I'm sure. And my sister embarrassed. And give the old man transient cognitive accolaydes to speak of at dinner. And teach - Cameron/Kranzy/Jodah/McCrae less dissapointed.

Also, yes I am aware I'm toying with almost every degree of indescriminitely cutting edge theory and hyperspecified ultra endpoint en-masse iterated contributory evolution of mathematical excess with the cadence of such tantamount to an unwashed troglodyte - from an institutional perspective that is.

It's what makes it fun.

Me no do math when shape do math good. learn read shape fun. math hard. shape simple.

All in good-ish faith!

Man I really hope it's not only correct, but it's correct and Penrose sees this in time, that legend had it figured out way back when!

J. R. Manuel


## Author

J. R. Manuel
