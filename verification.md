# Exact verification report

Checks passed: **2900/2900**.
Failed checks: **0**.

| Category | Checks |
|---|---:|
| source matrices | 52 |
| det sigma table | 52 |
| PBW overlaps | 104 |
| monomial basis changes | 70 |
| normal elements | 12 |
| central elements at h=1 | 32 |
| normalization of G | 4 |
| rational coordinate matrices | 2 |
| degenerate presentations | 25 |
| negative controls | 10 |
| PI table examples | 10 |
| root-of-unity arithmetic | 7 |
| additional PBW identities | 61 |
| manuscript discrepancy controls | 1 |
| localization normality | 48 |
| anti-isomorphism identities | 12 |
| central rational certificates | 102 |
| non-PI affine matrices | 14 |
| central power certificates | 9 |
| normal-generator matrices | 4 |
| Laurent identities | 490 |
| plane charts | 28 |
| formal parameter charts | 42 |
| symbolic torus splitting | 12 |
| Laurent inverse controls | 19 |
| crossed-product closure | 38 |
| torus action matrices | 19 |
| general algebraic certificates | 47 |
| plane lemma identities | 20 |
| finite-module certificates | 7 |
| general Dickson certificate | 2 |
| C rank-nine module | 6 |
| C module reduction regressions | 55 |
| nonprimitive orbit invariants | 8 |
| primitive localization obstructions | 6 |
| conditional orbit regressions | 24 |
| DG center lattice certificates | 5 |
| unbounded algebraic proof objects | 40 |
| unbounded integer ideal certificates | 5 |
| all-family degree certificates | 84 |
| center-action parity branches | 9 |
| fourth-root degree branches | 6 |
| K/L parity certificates | 6 |
| C bound certificates | 3 |
| weight matrix regressions | 6 |
| weight spectrum regressions | 6 |
| DG central generator regressions | 128 |
| Dickson samples n=3..8 | 12 |
| DG degree samples N<=10 | 412 |
| L degree samples N<=8 | 720 |
| independent Smith controls | 4 |

The parameter equations are reduced symbolically; no floating point is used.
For the 26 presentations, all four cubic ambiguities are resolved in the admissible parameter field.
The 26 Sigma matrices and independently transcribed M matrices agree after tensor reshuffling.
The overall factor h missing in the printed M for J is restored explicitly.
All 26 det(sigma) matrices agree with the manuscript table.
All 14 monomial basis changes are checked at symbolic parameters.

## Scope

All 21 closed degree models are implemented with exact torus matrices, center lattices and crossed-product action images, without a root-order cutoff.
The derived degree is compared independently with the manuscript formula. General recurrence and module identities have separate algebraic certificates.
Finite regression inputs remain explicitly labeled. Abstract division-ring, ideal-theoretic and independence proofs are documented mathematical dependencies, not proof-assistant formalizations.
The PI table evaluator still instantiates Table tab:main; it is not presented as an independent proof.
The exact degree of C and remaining primitivity cases stay open, as in the manuscript.

## Manuscript correction found by the completed checks

In the L torus ordered as w,y1,y1^-1*y2,u, the first upper exponent must be -(2*eta+phi).
The manuscript prints +(2*eta+phi). The relation w*y1=(h^2*f)^(-1)*y1*w verifies the correction.
Its Pfaffian and alternating Smith invariants are unchanged, so the stated degree formula remains unchanged.
See docs/MANUSCRIPT_CORRECTIONS.md for the exact replacement.

## Negative control for the printed Z relations

Replacing x1*y1 by x1*y2 in the two indicated positions produces nonzero cubic overlap residuals.
The four explicit residuals are retained in verification.json.

Python 3.12.14; SymPy 1.14.0.
