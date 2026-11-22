# How the computations support the general arguments

This document distinguishes exact symbolic ingredients, arbitrary-order
algorithms and finite regression inputs. It does not turn a list of passing
examples into a proof for all parameters. Labels refer to the supplied
manuscript snapshot.

## Defining data and admissible charts

All 26 Sigma and independently transcribed M matrices are checked, with the
parameter equations reduced exactly. The four degree-three ambiguities of
each terminating quadratic presentation resolve. Each of the fourteen basis
changes is checked in the original algebra, including the defining relation
of its balanced or coordinate plane.

Thus the Laurent word models use verified relations, rather than assuming
their charts solely because they are listed in a table. The inverse calculus
checks `lem:A1`, `eq:A2`, `eq:A3` and `lem:A3`, and has independent inverse and
mutation controls. Formal scalar charts are compared with the original
parameters before being used for degree computations.

## General Dickson certificate

For M(t)=[[pt,c12],[c21,p^-1*t]], set

\[
\alpha=-\frac{c_{12}c_{21}p^2}{(p^2-1)^2},\quad
t=w+\alpha/w,\quad
\lambda_1=pw+\alpha/(pw),\quad
\lambda_2=w/p+\alpha p/w.
\]

The program checks M V=V diag(lambda1,lambda2), with columns
(c12,lambda_i-pt), as well as the trace and determinant identities. The
determinant of V is explicitly nonzero in the stated rational parameter chart.
The c12=0 chart is covered separately by an explicit triangular eigenbasis,
including the diagonal case c12=c21=0.

The final nth-power step is checked with formal variables P=p^n, W=w^n and
A=alpha^n. Both

\[
PW+A/(PW)-(W+A/W),\qquad W/P+AP/W-(W+A/W)
\]

vanish at P=1. Together with the defining Dickson identity, this is the
algebraic content of the argument for every order n; it has no numerical
order cutoff. The extra direct checks at n=3,...,8 are regressions.

The q-integer induction used for powers of y_i is supported by the exact
geometric recurrence. `dg_central_elements(name,n)` then constructs the
displayed polynomial central generators for any specified primitive order
n>=3, at h=1. Their centrality is additionally checked at several orders.
The finite-order twist theorem transfers the PI argument to torsion h.

## Torus matrices, actions and crossed-product factors

The Laurent models derive the six commutation scalars of each torus, every
normalizer's integral exponent matrix and its scalar multipliers. They also
check all products of the chosen crossed-product cosets and the lattice index
that distinguishes those cosets. The polynomial checks for Theta,Z1,Z2 and
their squares connect the F,I,Q,R,V models to the original presentations.

For each PI input, let A be its integer alternating matrix and N a common
root order. The center lattice Gamma is computed as follows:

\[
L=\mathbb Z^d+(A^T/N)\mathbb Z^d,
\qquad \Gamma=L^*.
\]

If H is a column Hermite basis of N Z^d+A^T Z^d, then Gamma has basis
N H^(-T). The program checks integrality, A Gamma=0 modulo N and its determinant.
Its index is checked independently by the ordinary integer Smith diagonal.
All phase inputs are exact integers; composite orders and singular pairings
are supported. No enumeration of the entire residue space is needed.

Normalizer images on central monomials include the full quantum multiplication
cocycle, including negative powers. Their distinct signatures on a lattice
basis give the action image order on the center. For a four-coset model the
order is two or four, handled by `prop:crossed`; for the F/I models all eight
actions are distinct. These mathematical crossed-product results are stated
dependencies, not assumptions that every coset always contributes a factor.

`degree_certificate` computes the algebra degree from this model and then
compares it with the separately entered closed formula. It accepts arbitrary
common orders. `degree_from_parameters` converts exact algebraic roots to a
single cyclotomic field and obtains exact discrete logarithms, without
numerical arguments or floating point tolerances.

## Why the monomial formulas are general

The symbolic unimodular splittings in `thm:monodeg` are checked before any
root specialization. For B,E,J,S,T,U, the first rank-two factor has parameter
-1 and the second has the displayed scalar h^4,4h^4 or16h^4. Every nontrivial
center action contains the indicated inversion, giving action order four.

For M,N,O,P,W,Z, the pair parameters are u^(+/-2) and v^2. Thus the torus degree
is ord(u^2)*ord(v^2). With m=ord(u^2), one has u^m=1 exactly when ord(u) is odd.
The action scalars derived by the word calculus give image order two for O,W,Z
exactly when both u and v have odd order, and four otherwise. In M,N,P the
remaining signs require simultaneously u^m=1 and (-1)^m=1, which is impossible.
This is the parity argument behind the formulas; it is separate from the
finite branch examples. K and L additionally have index-two source-coordinate models and generic
integer ideal equivalences. The complete general derivations are in
`UNBOUNDED_PROOFS.md`.

For D/G the generic Pfaffian is checked to be 2k(k+k'). After taking the actual
generated cyclic group, the gcd of the entries is coprime to N, yielding the
displayed degree N^2/gcd(N,2k(k+k')). L's corrected exponent has the same
Pfaffian and gcd as the printed one; its degree formula is unchanged.

## C: a general rank-nine reduction

For a=y1^2*y2 and b=y1^3-y2^3, the program checks their commutation and the
recurrence identities, including a*y1=p*y1*a and

\[
y_1^5=b y_1^2+a y_2^2.
\]

The module reducer uses these recurrences without a degree cutoff. Its basis is

\[
1,y_1,y_2,y_1^2,y_1y_2,y_2^2,y_1^3,y_1y_2^2,y_1^4.
\]

The sum of its degree shifts is (1+t+t^2)^2, and the program checks
(1+t+t^2)^2/(1-t^3)^2=1/(1-t)^2. Generation follows from the checked reducing
identities, and the Hilbert-series equality gives freeness. Direct reductions
at bounded degrees are extra regression tests. The C torus calculation and
the exhaustive positive-integer solutions of 9=[L:K][L:F]ell^2 recover the
possible degree ratios 1,3,9. The exact choice remains open.

## Primitivity ingredients

Rational central witnesses are certified by a common normalizing automorphism
of numerator and denominator, not by informal cancellation of noncommuting
fractions. Orbit polynomial coefficients are checked to be invariant under
all torus generators and all crossed-product cosets. Conditional orbit APIs
keep the other parameter symbolic when only u or v is torsion.

For the primitive monomial branch, the program derives either a nonidentity
exponent action or an odd/even exponent obstruction to an inner monomial
action. Normality of Nx,Ny is checked in the PBW algebra. The mathematical
simplicity, prime-ideal localization, local-closedness and Nullstellensatz
steps remain those in the article.

For D/G, Omega is checked to commute with all Ore generators, Xi's four
normalizing relations are checked, and the generic trace, eigenvalue and
nonzero-Omega obstruction identities in `lem:hw` and `lem:lower` are verified.
The right-action matrix API accepts arbitrary n. Its closed tridiagonal form
and spectra have additional finite regression checks. The induction now has generic complementary spectral projectors, both
nonvanishing-branch obstructions, and an affine interval cover certificate for
every degree. See `UNBOUNDED_PROOFS.md`. The concluding ideal-theoretic deductions
remain the mathematical arguments in the manuscript, using these ingredients.

## Coverage and remaining open statements

`reports/coverage.json` inventories every manuscript label and records the
proof dependencies and open questions explicitly. The required computational
block check fails if a block has no certificate. The original mathematical
attributions are retained. This package does not claim a proof-assistant
formalization or resolve questions that the manuscript leaves open.
