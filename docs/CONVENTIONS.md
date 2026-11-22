# Matrix and rewriting conventions

Write P=(p12,p11), Q=(q12,q11). The plane relations are

\[
x_2x_1=q_{12}x_1x_2+q_{11}x_1^2,\qquad
y_2y_1=p_{12}y_1y_2+p_{11}y_1^2.
\]

The mixed relation is

\[
\begin{pmatrix}y_1x_1\\y_1x_2\\y_2x_1\\y_2x_2\end{pmatrix}
=\Sigma\begin{pmatrix}x_1y_1\\x_2y_1\\x_1y_2\\x_2y_2\end{pmatrix}.
\]

M uses orders (y1x1,y2x1,y1x2,y2x2) and (x1y1,x1y2,x2y1,x2y2).
Thus M=R Sigma R, where R interchanges positions 2 and 3. The M arrays in
`families.py` were entered independently from the published display. For J,
the omitted overall h is restored; see `SOURCE_NOTES.md`.

## The determinant automorphism

Let Sigma have 2 by 2 blocks Aij, with sigma_ij(X)=Aij X for X=(x1,x2)^T.
Composition of these maps reverses their displayed coefficient matrix order:

\[
(\sigma_{kl}\circ\sigma_{ij})(X)=A_{ij}A_{kl}X.
\]

The defining expression

\[
\det\sigma=-p_{11}\sigma_{12}\sigma_{11}
              +\sigma_{22}\sigma_{11}
              -p_{12}\sigma_{12}\sigma_{21}
\]

therefore has the coefficient matrix

\[
D=-p_{11}A_{11}A_{12}+A_{11}A_{22}-p_{12}A_{21}A_{12}.
\]

The code also checks that D preserves the Q-plane relation. The ordinary scalar
determinant det(Sigma), used to check invertibility of the mixed matrix, is
a separate quantity. Confusing these objects or reversing composition order
can give incorrect answers even if a few scalar examples happen to agree.

## Exact parameter arithmetic

- For B,E,F, reduce coefficients modulo p^2+1.
- For C, reduce modulo p^2+p+1.
- For I,J, reduce modulo q^2+1.
- For square-root monomial bases, substitute f=s^2 and require s!=0.
- All other parameters remain symbolic, with the nonvanishing restrictions in
  the family record. No root-of-unity assumption is imposed during PBW checks.
- Assigned algebraic parameters must satisfy their equations exactly; unchecked
  symbolic replacements of a constrained p or q are rejected.

Each rational expression is reduced modulo the parameter equation, including
inverting its denominator when necessary and possible. Floating point
coefficients are rejected. The parameter equations cover both branches; e.g.
using p^2+1 rather than just p=i checks both p=i and p=-i.

## PBW and termination

Words are tuples in the alphabet x1<x2<y1<y2. Each of the six quadratic rules
strictly decreases degree-lexicographic order. Their right sides contain only
ordered pairs, so the only possible inclusion-free overlap words are

\[
y_1x_2x_1,\quad y_2x_2x_1,\quad y_2y_1x_1,\quad y_2y_1x_2.
\]

The program first reduces the left pair and the right pair independently in
each word, then compares the resulting normal forms. All 104 generic residuals
vanish. By termination and the Diamond Lemma, this is a finite confluence
check for the presentations on the admissible parameter locus; it is stronger
than testing a collection of random words of bounded degree.

## Monomial basis matrices

Rows of E and Z are the new x-basis and y-basis coefficients. The mixed matrix
in the new bases is

\[
\Sigma'=(Z\otimes E)\Sigma(Z\otimes E)^{-1}.
\]

It is compared with the sparse matrix for z_a e_b=lambda_ab e_sigma(b) z_tau(a).
Type 0 has no transposition; type I has both; type II transposes only the
y-basis; type III transposes only the x-basis. The four lambda scalars retain
their full h factors (except N, which has no independent h).

## Smith normal form

For v_i v_j=zeta^a_ij v_j v_i, let A=(a_ij) be integer alternating and zeta
primitive of order N. If d_i are the ordinary integer Smith diagonal entries,
including zeros, then

\[
|\operatorname{im}(A\bmod N)|=\prod_i\frac{N}{\gcd(N,d_i)},\qquad
\operatorname{PIdeg}(T)=\sqrt{|\operatorname{im}(A\bmod N)|}.
\]

This works for composite N and singular matrices, without a field-rank shortcut.
The code computes the degree of the torus supplied. For L, the algebra degree
uses the additional crossed-product factor 2 established in the manuscript.
For D and G it uses the birational torus identification from the manuscript.

SymPy API: [smith_normal_form](https://docs.sympy.org/latest/modules/matrices/normalforms.html).
