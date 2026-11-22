# Torus commutation and action matrices

All columns are generator images. If an action has matrix M and scalars c_i,
the i-th generator is sent to c_i times the ordered monomial with exponent column M_i.
These matrices are derived by exact Laurent word reduction, not copied from the closed degree formula.

## Family B

Formal scalar chart: $H=h$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & -1 & -1 & 1\\-1 & 1 & 1 & -1\\-1 & 1 & 1 & - \frac{1}{H^{4}}\\1 & -1 & - H^{4} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}-1 & 0 & 2 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  - p, \  - H^{2} p, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 2\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - p, \  -1, \  1, \  - \frac{p}{H^{2}}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}-1 & 0 & 2 & 0\\0 & -1 & 0 & 2\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - p, \  p, \  H^{2} p, \  - \frac{p}{H^{2}}\right)$$

## Family E

Formal scalar chart: $\alpha=h \left(1 + i\right)$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & -1 & 1 & 1\\-1 & 1 & 1 & -1\\1 & 1 & 1 & - \frac{1}{\alpha^{4}}\\1 & -1 & - \alpha^{4} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( i, \  - p, \  i \alpha^{2}, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 2\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  i, \  1, \  \frac{i p}{\alpha^{2}}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 2\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( i, \  - i p, \  i \alpha^{2}, \  \frac{i p}{\alpha^{2}}\right)$$

## Family J

Formal scalar chart: $\alpha=h \left(1 + i\right)$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & -1 & -1 & 1\\-1 & 1 & 1 & 1\\-1 & 1 & 1 & - \frac{1}{\alpha^{4}}\\1 & 1 & - \alpha^{4} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}-1 & 0 & 2 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - i, \  1, \  - i \alpha^{2} q, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - q, \  - i, \  1, \  - \frac{i}{\alpha^{2}}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}-1 & 0 & 2 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - i q, \  i, \  i \alpha^{2} q, \  - \frac{i}{\alpha^{2}}\right)$$

## Family K

Formal scalar chart: $H=h$, $S=s$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & -1 & q^{2} & \frac{1}{S^{2}}\\-1 & 1 & 1 & 1\\\frac{1}{q^{2}} & 1 & 1 & \frac{1}{H^{4}}\\S^{2} & 1 & H^{4} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( S, \  1, \  H^{2}, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( \frac{1}{q}, \  1, \  1, \  \frac{1}{H^{2}}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( \frac{S}{q}, \  1, \  H^{2}, \  \frac{1}{H^{2}}\right)$$

## Family L

Formal scalar chart: $H=h$, $F=f$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & 1 & q^{2} & F^{2}\\1 & 1 & 1 & 1\\\frac{1}{q^{2}} & 1 & 1 & \frac{1}{F^{4} H^{4}}\\\frac{1}{F^{2}} & 1 & F^{4} H^{4} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( \frac{1}{F}, \  -1, \  F^{2} H^{2}, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 2\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( \frac{1}{q}, \  1, \  1, \  - \frac{1}{F^{2} H^{2}}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 2\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( \frac{1}{F q}, \  -1, \  F^{2} H^{2}, \  - \frac{1}{F^{2} H^{2}}\right)$$

## Family M

Formal scalar chart: $\alpha=h \left(s + 1\right)$, $\beta=h \left(s - 1\right)$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & \frac{\alpha^{2}}{\beta^{2}} & 1 & 1\\\frac{\beta^{2}}{\alpha^{2}} & 1 & 1 & 1\\1 & 1 & 1 & \frac{1}{\alpha^{2} \beta^{2}}\\1 & 1 & \alpha^{2} \beta^{2} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - \frac{\alpha}{\beta}, \  1, \  - \alpha \beta, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  \frac{\beta}{\alpha}, \  1, \  \frac{1}{\alpha \beta}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - \frac{\alpha}{\beta}, \  \frac{\alpha}{\beta}, \  - \alpha \beta, \  \frac{1}{\alpha \beta}\right)$$

## Family N

Formal scalar chart: $\alpha=f - g$, $\beta=f + g$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & \frac{\beta^{2}}{\alpha^{2}} & 1 & 1\\\frac{\alpha^{2}}{\beta^{2}} & 1 & \frac{\beta^{2}}{\alpha^{2}} & 1\\1 & \frac{\alpha^{2}}{\beta^{2}} & 1 & \frac{1}{\alpha^{2} \beta^{2}}\\1 & 1 & \alpha^{2} \beta^{2} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}-1 & 0 & 2 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( \frac{\beta}{\alpha}, \  1, \  - \alpha^{2}, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  - \frac{\alpha}{\beta}, \  1, \  \frac{1}{\alpha \beta}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}-1 & 0 & 2 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - \frac{\beta}{\alpha}, \  - \frac{\beta}{\alpha}, \  - \alpha^{2}, \  \frac{1}{\alpha \beta}\right)$$

## Family O

Formal scalar chart: $\alpha=h \left(s + 1\right)$, $\beta=h \left(1 - s\right)$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & \frac{\beta^{2}}{\alpha^{2}} & 1 & 1\\\frac{\alpha^{2}}{\beta^{2}} & 1 & 1 & 1\\1 & 1 & 1 & \frac{1}{\alpha^{2} \beta^{2}}\\1 & 1 & \alpha^{2} \beta^{2} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( \frac{\beta}{\alpha}, \  1, \  \alpha \beta, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  \frac{\alpha}{\beta}, \  1, \  \frac{1}{\alpha \beta}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( \frac{\beta}{\alpha}, \  \frac{\beta}{\alpha}, \  \alpha \beta, \  \frac{1}{\alpha \beta}\right)$$

## Family P

Formal scalar chart: $\alpha=h \left(s + 1\right)$, $\beta=h \left(1 - s\right)$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & \frac{\alpha^{2}}{\beta^{2}} & 1 & \frac{\beta^{2}}{\alpha^{2}}\\\frac{\beta^{2}}{\alpha^{2}} & 1 & 1 & 1\\1 & 1 & 1 & \frac{1}{\alpha^{2} \beta^{2}}\\\frac{\alpha^{2}}{\beta^{2}} & 1 & \alpha^{2} \beta^{2} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( \frac{\alpha}{\beta}, \  -1, \  \alpha \beta, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 2\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  \frac{\beta}{\alpha}, \  1, \  - \frac{1}{\beta^{2}}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 2\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( \frac{\alpha}{\beta}, \  - \frac{\beta}{\alpha}, \  \alpha \beta, \  - \frac{1}{\beta^{2}}\right)$$

## Family S

Formal scalar chart: $H=2 h$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & -1 & 1 & 1\\-1 & 1 & 1 & 1\\1 & 1 & 1 & \frac{1}{H^{4}}\\1 & 1 & H^{4} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  1, \  H^{2}, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  -1, \  1, \  \frac{1}{H^{2}}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  -1, \  H^{2}, \  \frac{1}{H^{2}}\right)$$

## Family T

Formal scalar chart: $H=2 h$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & -1 & 1 & -1\\-1 & 1 & 1 & 1\\1 & 1 & 1 & \frac{1}{H^{4}}\\-1 & 1 & H^{4} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - i, \  1, \  H^{2}, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  -1, \  1, \  - \frac{1}{H^{2}}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - i, \  -1, \  H^{2}, \  - \frac{1}{H^{2}}\right)$$

## Family U

Formal scalar chart: $H=2 h$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & -1 & 1 & 1\\-1 & 1 & -1 & 1\\1 & -1 & 1 & \frac{1}{H^{4}}\\1 & 1 & H^{4} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  1, \  - H^{2}, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  - i, \  1, \  \frac{1}{H^{2}}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  i, \  - H^{2}, \  \frac{1}{H^{2}}\right)$$

## Family W

Formal scalar chart: $\alpha=h \left(s - i\right)$, $\beta=h \left(s + i\right)$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & \frac{\beta^{2}}{\alpha^{2}} & 1 & 1\\\frac{\alpha^{2}}{\beta^{2}} & 1 & \frac{\beta^{2}}{\alpha^{2}} & 1\\1 & \frac{\alpha^{2}}{\beta^{2}} & 1 & \frac{1}{\alpha^{2} \beta^{2}}\\1 & 1 & \alpha^{2} \beta^{2} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}-1 & 0 & 2 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - \frac{\beta}{\alpha}, \  1, \  - \alpha^{2}, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  \frac{\alpha}{\beta}, \  1, \  \frac{1}{\alpha \beta}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}-1 & 0 & 2 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - \frac{\beta}{\alpha}, \  \frac{\beta}{\alpha}, \  - \alpha^{2}, \  \frac{1}{\alpha \beta}\right)$$

## Family Z

Formal scalar chart: $\alpha=h \left(s - i\right)$, $\beta=h \left(s + i\right)$.

Torus generator order: `r`, `rho`, `X`, `Y`.

$$C=\left[\begin{matrix}1 & \frac{\alpha^{2}}{\beta^{2}} & 1 & \frac{\beta^{2}}{\alpha^{2}}\\\frac{\beta^{2}}{\alpha^{2}} & 1 & 1 & 1\\1 & 1 & 1 & \frac{1}{\alpha^{2} \beta^{2}}\\\frac{\alpha^{2}}{\beta^{2}} & 1 & \alpha^{2} \beta^{2} & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - \frac{\alpha}{\beta}, \  1, \  \alpha \beta, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 2\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  - \frac{\beta}{\alpha}, \  1, \  - \frac{1}{\beta^{2}}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & -1 & 0 & 2\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - \frac{\alpha}{\beta}, \  - \frac{\beta}{\alpha}, \  \alpha \beta, \  - \frac{1}{\beta^{2}}\right)$$

## Family F

Formal scalar chart: $H=- h \left(1 - p\right)$.

Torus generator order: `r`, `X`, `Y`, `P`.

$$C=\left[\begin{matrix}1 & 1 & 1 & 1\\1 & 1 & \frac{1}{H^{4}} & 1\\1 & H^{4} & 1 & 1\\1 & 1 & 1 & 1\end{matrix}\right]$$

Coset `000`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `001`:

$$M=\left[\begin{matrix}-1 & 2 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  -1, \  1, \  1\right)$$

Coset `010`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  H^{2}, \  1, \  1\right)$$

Coset `011`:

$$M=\left[\begin{matrix}-1 & 2 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  - H^{2}, \  1, \  1\right)$$

Coset `100`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  1, \  \frac{1}{H^{2}}, \  -1\right)$$

Coset `101`:

$$M=\left[\begin{matrix}-1 & 2 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  -1, \  \frac{1}{H^{2}}, \  -1\right)$$

Coset `110`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  H^{2}, \  \frac{1}{H^{2}}, \  -1\right)$$

Coset `111`:

$$M=\left[\begin{matrix}-1 & 2 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  - H^{2}, \  \frac{1}{H^{2}}, \  -1\right)$$

## Family I

Formal scalar chart: $H=\frac{1}{h \left(1 - q\right)}$.

Torus generator order: `r`, `X`, `Y`, `P`.

$$C=\left[\begin{matrix}1 & 1 & 1 & 1\\1 & 1 & \frac{1}{H^{4}} & 1\\1 & H^{4} & 1 & 1\\1 & 1 & 1 & 1\end{matrix}\right]$$

Coset `000`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `001`:

$$M=\left[\begin{matrix}-1 & 2 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  -1, \  1, \  1\right)$$

Coset `010`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  H^{2}, \  1, \  1\right)$$

Coset `011`:

$$M=\left[\begin{matrix}-1 & 2 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  - H^{2}, \  1, \  1\right)$$

Coset `100`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  1, \  \frac{1}{H^{2}}, \  -1\right)$$

Coset `101`:

$$M=\left[\begin{matrix}-1 & 2 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  -1, \  \frac{1}{H^{2}}, \  -1\right)$$

Coset `110`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  H^{2}, \  \frac{1}{H^{2}}, \  -1\right)$$

Coset `111`:

$$M=\left[\begin{matrix}-1 & 2 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  - H^{2}, \  \frac{1}{H^{2}}, \  -1\right)$$

## Family Q

Formal scalar chart: $H=h$.

Torus generator order: `n`, `Theta`, `rho`, `Y`.

$$C=\left[\begin{matrix}1 & \frac{1}{H^{2}} & -1 & - \frac{1}{H^{2}}\\H^{2} & 1 & 1 & \frac{1}{H^{2}}\\-1 & 1 & 1 & 1\\- H^{2} & H^{2} & 1 & 1\end{matrix}\right]$$

Coset `0`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `1`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & -1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( - i H, \  H, \  1, \  1\right)$$

## Family R

Formal scalar chart: $H=h$.

Torus generator order: `r`, `X`, `Z1`, `W`.

$$C=\left[\begin{matrix}1 & 1 & 1 & 1\\1 & 1 & \frac{1}{H^{2}} & \frac{1}{H^{4}}\\1 & H^{2} & 1 & 1\\1 & H^{4} & 1 & 1\end{matrix}\right]$$

Coset `00`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `01`:

$$M=\left[\begin{matrix}-1 & 2 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  - H^{2}, \  1, \  1\right)$$

Coset `10`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( -1, \  1, \  \frac{1}{H}, \  - \frac{1}{H^{2}}\right)$$

Coset `11`:

$$M=\left[\begin{matrix}-1 & 2 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  - H^{2}, \  \frac{1}{H}, \  - \frac{1}{H^{2}}\right)$$

## Family V

Formal scalar chart: $H=h$.

Torus generator order: `n`, `Theta`, `rho`, `Y`.

$$C=\left[\begin{matrix}1 & \frac{1}{H^{2}} & 1 & \frac{1}{H^{2}}\\H^{2} & 1 & -1 & - \frac{1}{H^{2}}\\1 & -1 & 1 & 1\\H^{2} & - H^{2} & 1 & 1\end{matrix}\right]$$

Coset `0`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( 1, \  1, \  1, \  1\right)$$

Coset `1`:

$$M=\left[\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & -1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]$$

$$(c_i)=\left( H, \  - i H, \  1, \  1\right)$$

## Family D

Here H=h. Generator order: x1, theta, y1, y2. The birational relations
are certified through Omega and Xi; the model needs no crossed-product factor.

$$C=\left[\begin{matrix}1 & -1 & - \frac{1}{H p} & \frac{1}{H p}\\-1.0 & 1 & \frac{1}{p} & p\\- H p & p & 1 & \frac{1}{p}\\H p & \frac{1}{p} & p & 1\end{matrix}\right]$$

## Family G

Here H=h. Generator order: x1, theta, y1, y2. The birational relations
are certified through Omega and Xi; the model needs no crossed-product factor.

$$C=\left[\begin{matrix}1 & 1 & \frac{1}{H p} & \frac{1}{H p}\\1.0 & 1 & \frac{1}{p} & p\\H p & p & 1 & \frac{1}{p}\\H p & \frac{1}{p} & p & 1\end{matrix}\right]$$

