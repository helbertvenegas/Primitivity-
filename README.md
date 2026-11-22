# Double Ore extensions (14641): exact computational companion

Version **0.3.0** accompanies **Polynomial identities and primitivity in double
extension regular algebras of type (14641)**, by James Gómez and Helbert Venegas.

The 26 defining presentations come from Section 4 of J. J. Zhang and J. Zhang,
*Double extension regular algebras of type (14641)*, J. Algebra 322 (2009), 373--409.
[DOI](https://doi.org/10.1016/j.jalgebra.2009.03.041) ·
[arXiv](https://arxiv.org/abs/0712.2550).

## Start here

- **Colab/Jupyter:** open `notebooks/Double_Ore_14641_PI_verification.ipynb`.
  It is self-contained and needs no checkout, private files, Drive access or GPU.
- **Defining matrices:** `docs/MATRICES.md` lists P,Q,Sigma,M and det(sigma) for A--Z.
- **Degree matrices:** `docs/TORUS_ACTION_MATRICES.md` lists the torus commutation
  matrices and crossed-product actions for all 21 closed degree models.
- **Full results:** `reports/verification.json` contains every residual/certificate;
  `reports/verification.md` summarizes it and `reports/coverage.json` maps the
  checks to all 108 manuscript labels.
- **Independent checker:** `independent/check.py` uses only the Python standard library.
  It checks the source data, different Laurent normal forms, generic identities,
  and unimodular Smith/kernel witnesses.
- **Unbounded arguments:** `docs/UNBOUNDED_PROOFS.md` gives the general weight
  induction, all parity classes, and D/G/K/L integer ideal derivations.
- **General certificates:** `docs/SYMBOLIC_CERTIFICATES.md` explains the algebraic
  arguments connecting the exact computations to the general formulas.
- **Correction found:** see `docs/MANUSCRIPT_CORRECTIONS.md` for the L matrix sign.
- **Español:** see `docs/USO_RAPIDO_ES.md`.

```bash
python -m pip install -r requirements.txt
python -m ore14641.verify
python -m ore14641.certificates
python -I independent/check.py reports/certificates.json reports/independent_verification.json
python -m unittest discover -s tests -v
```

Tested with Python 3.12 and SymPy 1.14.0. The verifier fails if a check fails or
one of the required computational blocks has no labeled certificate.

## Scope

| Part | Implemented verification |
|---|---|
| Presentations | All 26 P,Q,Sigma,M data, parameter equations, determinants, four cubic PBW ambiguities per family |
| Monomial criterion | All 14 basis changes, admissible plane charts, inverse identities, rescaling invariants and plane/module calculations |
| No-PI cases | Affine matrices derived from Sigma, eigenvalue obstructions, and degenerate presentations |
| General Dickson argument | Symbolic eigenvectors, trace/determinant identities, and a root-power certificate valid for the general order argument |
| Degree models | All 21 torus matrices; exact crossed-product multiplication and conjugation matrices; tensor-product splittings |
| Degrees at roots of unity | Arbitrary common root orders; Hermite center lattice, independent Smith index, center-action image, comparison with closed formulas |
| F,I,Q,R,V | Theta and Z1,Z2 identities, squares, plane normalizers, torus and action tables |
| C | Exact nine-element module reduction, Hilbert-series identity, subtorus degree, exhaustive degree-ratio constraints |
| Primitivity ingredients | Rational central witnesses, central orbit polynomials, localization normality, outer-action parity obstructions, Omega/Xi and weight matrices |
| Uniform proofs | General spectral projector/interval induction, integral ideal equivalences, complete center parametrizations and exhaustive parity classes |
| Traceability | Hypotheses, checked claims and certificate IDs for required blocks; all 108 manuscript labels; computational checks, mathematical dependencies and open questions distinguished |

Finite regression inputs are labeled separately. The APIs are not restricted
to the orders used in those regressions. In particular, the suite also checks
examples at common order 3988 without enumerating all residues.

The mathematical division-ring, independence, simplicity, ideal-theory and
Nullstellensatz arguments remain explicit proof dependencies. This is a
computational companion, not a proof-assistant formalization. The exact degree
of C and the remaining primitivity questions are left open as in the manuscript.

## Exact parameter interface

```python
from sympy import Rational, sqrt
from ore14641 import f, h
from ore14641.degrees import degree_from_parameters

z = degree_from_parameters('Z', {f: 3, h: -Rational(1, 2)})
w = degree_from_parameters('W', {f: 3, h: -Rational(1, 2)})
assert z['algebra_pi_degree'] == 6
assert w['algebra_pi_degree'] == 12

# All three possibilities are returned; no exact value is invented.
c = degree_from_parameters('C', {h: 1/sqrt(3)})
assert c['exact_degree_is_open']
print(c['possible_pi_degrees'])
```

Default algebraic branches for omitted p or q are recorded in the returned
certificate. Use exact integers, Rational, radicals or algebraic numbers.
Floating point phases and parameters are rejected.

For large root orders, supply integer phases directly:

```python
from ore14641.words import alpha, beta
from ore14641.degrees import degree_certificate

result = degree_certificate('O', 3988, {alpha: 4, beta: 8})
assert result['center_action_order'] == 2
assert result['algebra_pi_degree'] == 2*997**2
```

The common cyclic group may be larger than the subgroup generated by the
parameters. For the full scalar charts it is convenient to include i.

## Important conventions and source corrections

The input/output orders are

\[
(y_1x_1,y_1x_2,y_2x_1,y_2x_2)^T
=\Sigma(x_1y_1,x_2y_1,x_1y_2,x_2y_2)^T.
\]

If R interchanges positions 2 and 3, M=R Sigma R. The ordinary scalar determinant
of Sigma is different from the determinant automorphism det(sigma).
The twist h is retained; N has no independent h; G is stored in the original
coordinates and its shift is checked separately. Square-root monomial charts
are not applied at f=0.

The printed Z expanded relations and the missing overall h in the printed M
for J are documented in `docs/SOURCE_NOTES.md`. The completed checks also find
that the manuscript's first L commutation exponent should have the opposite
sign. The closed degree formula is unaffected; its Pfaffian and alternating
Smith invariants are unchanged. The source manuscript is not silently edited.

## Files, rebuilding and citation

- `ore14641/`: defining data, exact polynomial and Laurent engines, degree and
  primitivity certificates, coverage ledger and verifier.
- `tests/`: independent regression and mutation/negative controls.
- `notebooks/`: standalone notebook, including code and saved outputs.
- `reports/`: full certificates, machine-readable torus models, provenance and hashes.
- `docs/`: both matrix catalogues, proof bridges, correction note and Spanish guide.

```bash
python scripts/build_notebook.py
python -m pip install -r requirements-notebook.txt
python scripts/execute_notebook.py
```

The runner uses nbclient and a real ipykernel, storing Jupyter output messages.
The default uses a separate kernel process. Where socket binding is restricted,
`python scripts/execute_notebook.py --kernel inprocess` uses the official
IPython in-process kernel in a fresh Python invocation and an empty working
directory. The saved execution record states the exact mode. An interactive
Colab runtime requires a Google sign-in and is not claimed as executed.
The GitHub workflow runs the primary suite, independent checker, mutation tests
and the notebook with a separate kernel process.

The supplied journal PDF and manuscript are not redistributed. Their SHA-256
fingerprints identify the snapshots checked. Citation metadata are supplied
in CITATION.cff. No public repository URL or DOI is invented; add them after
publication. No license has been assigned to this prepared draft.
