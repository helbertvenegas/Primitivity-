# Source notes and corrections

The matrix source is the supplied **published** Zhang--Zhang article,
J. Algebra 322 (2009), 373--409, Section 4, rather than a secondary list of
expanded relations. The accompanying manuscript source is identified by
SHA-256 in `reports/source_provenance.json`. The source PDFs and manuscript
are not included in this distribution.

| Families | Subcase | Printed pages | PDF pages in the supplied file |
|---|---|---:|---:|
| A | 4.1.4 | 395 | 23 |
| B,C | 4.2.2 | 397 | 25 |
| D,E,F | 4.2.3 | 397 | 25 |
| G | 4.2.4 | 398 | 26 |
| H | 4.3.1 | 398 | 26 |
| I,J,K,L | 4.3.2 | 398 | 26 |
| M,N,O,P,Q,R,S,T,U | 4.3.3 | 398--399 | 26--27 |
| V,W,X,Y | 4.3.4 | 399 | 27 |
| Z | 4.4.4 | 400 | 28 |

## Family Z: matrix versus expanded relations

On printed p. 400, the Sigma matrix gives

\[
y_1x_1=h(x_1y_1+x_2y_2),\qquad
y_2x_2=h(fx_1y_1-x_2y_2).
\]

The expanded list below that matrix instead uses x1*y2 in place of x1*y1
in both positions (and suppresses the common h, as in the expanded A list).
The implementation follows the Sigma matrix and restores the twist parameter.

For an explicit negative control, the verifier makes exactly those two changes
to Sigma. All four cubic overlap residuals become nonzero. For example, the
left-minus-right residual for y1*x2*x1 is

\[
h^2\big((f-1)x_1^2y_2+f x_1x_2y_1-(f+2)x_1x_2y_2\big).
\]

This cannot vanish on the admissible Z locus because f!=0 and h!=0. The full
four certificates are retained under category `negative controls` in the JSON
report. This is a consistency check of the two presentations, not a change of
the family specified by the original Sigma matrix.

## Family J: missing overall factor in M

On printed p. 398 the display reads Sigma=h Sigma0 but M is printed without
the corresponding h. By the definition of M in (E3.7.1), M=R Sigma R, so the
factor must also occur in M. The independently transcribed M0 is stored with
its required h multiplier. The source note and code comment expose this
normalization explicitly.

## Family G: original and shifted coordinates

The source Sigma includes p*x1*y1 in the relation for y1*x2 and -x1*y2 in
that for y2*x2. The manuscript uses e2=x2+(p-1)^(-1)*x1. The registry retains
the original source matrix; four additional PBW identities check that this
shift gives precisely the manuscript's normalized relations. Thus source and
normalized matrices are not inadvertently interchanged.

## Parameter endpoints

D allows p=+/-1, although those members are non-PI. M,O,P,W allow f=0;
their mixed matrices remain invertible, but their square-root monomial basis
changes may become singular. N has f^2!=g^2 and no independent h. The engine
checks the admissible degenerate presentations directly. In particular,
O(0,h) and D(-1,h) have identical P,Q,Sigma data.

The PI criteria and expected determinant automorphism matrices are separately
transcribed from the supplied manuscript. Expected data are not generated
by the functions that are checked against them.
