# Certificates for unbounded degrees and root orders

These are algebraic certificates for the general arguments, not a claim of
formalization in a proof assistant. The independent checker verifies identities,
integral equivalences, Laurent words and exhaustive parity classes. The
structural algebra theorems used to interpret them are listed at the end.
All expressions below are in characteristic zero with the admissibility
conditions of the original family retained.

## 1. Weight induction for every n

Write a=bc, beta1=c(p^2-1)/p^2, beta2=-b(p^2-1), and
lambda(z)=z+a/z. For eta=p^nu*w, right multiplication by t on
(my1,my2) has the matrix

\[
B(\eta)=\begin{pmatrix}p\lambda(\eta)&-p\beta_1\\
-p^{-1}\beta_2&p^{-1}\lambda(\eta)\end{pmatrix}.
\]

Set lambda_plus=lambda(p*eta) and lambda_minus=lambda(eta/p). The certificates
retain the two rational matrices

\[
\Pi_+=\frac{B-\lambda_-I}{\lambda_+-\lambda_-},\qquad
\Pi_-=\frac{\lambda_+I-B}{\lambda_+-\lambda_-}.
\]

Both engines check Pi_plus+Pi_minus=I, Pi_plus^2=Pi_plus,
Pi_minus^2=Pi_minus, Pi_plus*Pi_minus=0 and
B=lambda_plus*Pi_plus+lambda_minus*Pi_minus. Their denominator is nonzero:

\[
\lambda_+-\lambda_-=(p-p^{-1})(\eta-bc/\eta).
\]

Here p is not a root of unity, bc!=0 and w is transcendental. If the lower
component vanishes, the first relation forces my2=(b/eta)*my1. The exact
identity lambda(eta)*(b/eta)-b-c*(b/eta)^2=0 then implies m*Omega=0.
Regularity of Omega in the flat extension rules this out. Replacing eta by
bc/eta interchanges the two branches, proving the same nonvanishing assertion
for the upper component. These injectivity and flatness facts are explicit
mathematical dependencies, not numerical tests.

The distinctness certificate is

\[
\lambda(Pw)-\lambda(Qw)=(P-Q)\left(w-\frac{bc}{PQw}\right).
\]

For distinct powers P,Q of p it cannot vanish. At degree n the index is
nu=-n+2j, 0<=j<=n. The next two branches have indices

\[
\nu-1=-(n+1)+2j,\qquad \nu+1=-(n+1)+2(j+1).
\]

Their index intervals are [0,n] and [1,n+1], with union [0,n+1]. This is a
single induction scheme, not a loop up to a chosen n. The base is degree zero
and the PBW dimension at degree n+1 is n+2. Nonvanishing of both branches and
this dimension show that every eigenvalue occurs on exactly one line.
Centrality and regularity of Omega identify the zero-weight line in degree
2k with F(w)*Omega^k.

For the ideal consequence, the integer measure |nu| decreases by one at each
chosen branch, so it reaches zero after finitely many steps. Clearing the
localization denominators and applying the noetherian prime-product argument
gives the ideal statement in `prop:weights`. These last ideal-theoretic steps
remain the proof in the manuscript; the certificate does not replace them.

## 2. A general rank-two pair and center action

If ab=z*ba and z has exact finite order d, a monomial a^i*b^j is central
exactly when d divides i and j. Thus the center index is d^2 and the pair
has PI degree d. Two commuting such pairs have center lattice
 diag(m,m,l,l) and degree m*l. The word checker verifies each unimodular
splitting before applying this rule; no root order is substituted.

For B,E,J,S,T,U the pair degrees are 2 and ord(z), where z is H^4 or
-alpha^4 in the formal chart. All three nontrivial cosets have a nonidentity
integral exponent action. Such an action cannot be the identity on any
full-rank central lattice. Coset closure and the group law then give a faithful
center action of order four, and hence the factor 8 in the degree formula.

For M,N,O,P,W,Z, use u=alpha/beta, except Z where u=-alpha/beta, and
v=alpha*beta. The pair parameters are u^(+/-2) and v^2, so their degrees
are m=ord(u^2), l=ord(v^2). Exactly one nonidentity coset has identity exponent
matrix. The checker obtains its scalar characters from the signed Laurent words
and tests them on the four central generators of exponents m,m,l,l.

The following are exhaustive classes of exact orders, with r an integer in the
range making the order positive. They are not representative sample orders.

| Exact order d | m=ord(z^2) | Parity of m | z^m |
|---|---|---|---|
| 2r+1 | 2r+1 | odd | 1 |
| 4r+2 | 2r+1 | odd | -1 |
| 4r | 2r | even | -1 |

The rule follows from m=d/gcd(d,2), and the unique nontrivial square root of
one is -1. It applies to u and v independently. All nine combinations are
checked using the actual scalar characters. The exceptional action is trivial
exactly when both orders are odd for O,W,Z; it is never trivial for M,N,P.
Consequently the center-action factor is two in the stated O/W/Z branch and
four in the other branches. This proves the parity formulas for all orders.

## 3. F, I, Q, R and V

For F/I the only nontrivial torus pair has parameter H^4. The two remaining
central generators are r and P. Restricting the eight source-derived actions
to r and P gives eight distinct sign/inversion signatures, independently of
ord(H). Hence the degree is 8*ord(H^4).

For Q/V, let m=ord(h^2). In the order n,Theta,rho,Y, the complete central
congruences are

\[
\begin{array}{ll}
Q:&a\equiv0\ (2),\ c-d\equiv0\ (2),\ a-d\equiv0\ (m),\ a+b\equiv0\ (m),\\
V:&b\equiv0\ (2),\ c-d\equiv0\ (2),\ b+d\equiv0\ (m),\ a+b\equiv0\ (m).
\end{array}
\]

The exported basis matrices are

\[
\Gamma_Q=\begin{pmatrix}2&0&0&0\\-2&m&0&0\\2&0&m&2\\2&0&m&0\end{pmatrix},
\quad
\Gamma_V=\begin{pmatrix}-2&m&0&0\\2&0&0&0\\-2&0&m&2\\-2&0&m&0\end{pmatrix}.
\]

If C is the matrix of the four congruence forms and
D=diag(2,2,m,m), the checker verifies D^(-1)*C*Gamma is an **integer
unimodular** matrix independent of m. Therefore these columns parametrize
all solutions, rather than only a sublattice. Their index is 4m^2. The second
coset has a nonidentity exponent action on the full-rank center lattice,
giving degree 2*(2m)=4m.

For R, in the order r,X,Z1,W, centrality is m|b and m|(c+2d).
A complete parametrization has columns
(1,0,0,0), (0,m,0,0), (0,0,m,0), (0,0,-2,1).
The same unimodular parametrization check proves the index m^2. The four
actions on the central variable r are r,-r,-r^(-1),r^(-1), so the degree is 4m.

## 4. D and G: a primitive entry ideal

The alternating four by four matrix has Pfaffian 2k(k+k'). The certificates
include explicit integer matrices proving equality of the ideal generated by
its entries with (M,k,k') for D (N=2M), and with (k,k') for G. Generation
of the scalar cyclic group implies that this entry gcd s1 is coprime to N.
The alternating Smith invariants satisfy s1|s2 and s1*s2=|Pfaffian|. Since
s1 is a unit modulo N, gcd(N,s2)=gcd(N,2k(k+k')), including Pfaffian zero.
The torus degree is therefore N^2/gcd(N,2k(k+k')).

## 5. K and L: complete general Smith and parity reductions

For K use an additional index-two torus in the ORIGINAL coordinates, ordered
as (y1,y2,x1,x2^2). Its upper exponent entries are

\[
(M,\eta,2\eta+\varphi,\eta,2\eta+\varphi,-2\theta),\quad N=2M.
\]

Its word relations are checked against the source and its Pfaffian is
-N*theta. Explicit ideal equivalences give
s1=gcd(M,eta,phi,2theta). In particular s1|N and
 gcd(s1,theta)=gcd(M,eta,phi,theta)=1 by generation. The second alternating
Smith invariant is N*theta/s1, so the torus degree is N. The remaining
x2 coset has a nonidentity exponent action, swapping the two y exponents,
on the full-rank center. The algebra degree is 2N. This also agrees with
the independently calculated four-coset model used by the degree API.

For L use the index-two torus in the manuscript, ordered as
(w,y1,y1^(-1)*y2,u). The source-word checker verifies the corrected entries

\[
(-2\eta-\varphi,0,2\theta,M,2\eta+2\varphi,0).
\]

The Pfaffian is N*theta. Integer ideal equivalences give
s1=gcd(2eta,phi,2theta,M). Since s1|M, the torus degree is

\[
\frac{N}{d},\qquad d=\gcd(\theta,M,\varphi,2\eta).
\]

The base ideal (theta,M,phi,eta) is the unit ideal. An explicit inclusion
of twice this base ideal in (theta,M,phi,2eta) proves d|2. Thus d is 1 or 2,
with d=2 exactly when theta,M,phi are all even. Also
N/N0=gcd(M,theta,phi), whose parity is even in precisely the same case.
The x1 coset has a nonidentity exponent action, giving algebra degree 2N/d.
This is the formula in the manuscript, for every generated cyclic group.
The first sign correction leaves both this entry ideal and the Pfaffian
unchanged.

## 6. What the second checker independently verifies

`python -I independent/check.py reports/certificates.json` uses the standard
library alone. It rechecks:

- the 26 source reshufflings, determinant-automorphism matrices and all cubic PBW ambiguities;
- the 14 monomial basis changes and their invertibility;
- signed Laurent torus words and every action image, using coordinate monomials or central squares plus alternating words;
- the generic rational identities, integer ideal witnesses, Q/V/R complete center parametrizations, induction index cover, and exhaustive parity partition;
- U*A*V=D with integral unimodular U,V, Smith divisibility, complete modular kernels, center-action signatures including negative powers, and degrees in each exported finite regression certificate.

For the last item, a Smith-derived kernel basis K=V*diag(N/gcd(N,Dii)) is
compared with the exported Gamma by Gamma=K*W, where W is an integer
unimodular matrix. This proves completeness, not just central membership.
The second checker contains no imports from SymPy or ore14641. Mutation tests
alter source matrices, projector identities, Smith matrices, kernel changes,
Laurent action scalars, parity classes and lattice parametrizations; each
alteration must be rejected.

The source definitions and mathematical bridge theorems remain trusted inputs.
An independent arithmetic implementation does not independently recover
Zhang--Zhang's classification from scratch. Structural dependencies are PBW
independence/domain, flatness, the center theorem for quantum tori, alternating
integer Smith theory, the finite crossed-product theorem, and the
ideal-theoretic/Nullstellensatz arguments explicitly recorded in the coverage
contracts. Open questions in the manuscript remain open.
