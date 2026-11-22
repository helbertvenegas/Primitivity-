"""Zhang--Zhang Section 4, with the twist parameter retained.

Source: J. Algebra 322 (2009), 373--409, printed pp. 395, 397--400.
The published Sigma matrices are the defining data. M is also transcribed
independently, rather than generated from Sigma, to check the transcription.
"""
from dataclasses import dataclass
from functools import lru_cache
import sympy as sp

p, q, f, g, h, s, t = sp.symbols('p q f g h s t')
I = sp.I
E2 = sp.eye(2)


@dataclass(frozen=True)
class Family:
    name: str
    P: tuple
    Q: tuple
    sigma: sp.ImmutableMatrix
    printed_m: sp.ImmutableMatrix
    expected_detsigma: sp.ImmutableMatrix
    relations: tuple
    nonzero: tuple
    section: str
    page: int
    pi_condition: str

    @lru_cache(maxsize=65536)
    def scalar(self, expression):
        """Exact rational arithmetic modulo the family's parameter equation."""
        expression = sp.sympify(expression)
        if expression.has(sp.Float):
            raise TypeError('Use exact Rational or algebraic coefficients, not floats')
        expression = sp.cancel(expression)
        for polynomial, variable in self.relations:
            numerator, denominator = sp.fraction(expression)
            numerator = sp.rem(numerator, polynomial, variable)
            denominator = sp.rem(denominator, polynomial, variable)
            if denominator == 0:
                raise ZeroDivisionError('Denominator vanishes in parameter field')
            if variable in denominator.free_symbols:
                numerator *= sp.invert(denominator, polynomial, variable)
                numerator = sp.rem(numerator, polynomial, variable)
                denominator = sp.Integer(1)
            expression = sp.cancel(numerator / denominator)
        return sp.expand(expression)

    def matrix(self, matrix):
        return sp.ImmutableMatrix(matrix).applyfunc(self.scalar)

    def specialize(self, substitutions):
        """Reject provably inadmissible exact specializations; retain symbolic ones."""
        substitutions = {sp.sympify(k): sp.sympify(v) for k, v in substitutions.items()}
        if any(value.has(sp.Float) for value in substitutions.values()):
            raise TypeError('Use exact Rational or algebraic parameters, not floats')
        for equation, variable in self.relations:
            value = sp.simplify(equation.subs(substitutions, simultaneous=True))
            if variable in substitutions and value != 0 and value.equals(0) is not True:
                raise ValueError(f'{self.name}: parameter equation {equation}=0 fails')
        for restriction in self.nonzero:
            if sp.simplify(restriction.subs(substitutions, simultaneous=True)) == 0:
                raise ValueError(f'{self.name}: requires {restriction} != 0')
        rel = tuple((sp.expand(eq.subs(substitutions, simultaneous=True)), var)
                    for eq, var in self.relations
                    if var not in substitutions)
        sub = lambda a: sp.sympify(a).subs(substitutions, simultaneous=True)
        mat = lambda a: sp.ImmutableMatrix(a.subs(substitutions, simultaneous=True))
        return Family(self.name, tuple(map(sub, self.P)), tuple(map(sub, self.Q)),
                      mat(self.sigma), mat(self.printed_m), mat(self.expected_detsigma),
                      rel, tuple(map(sub, self.nonzero)), self.section, self.page,
                      self.pi_condition)

    def reshuffled_m(self):
        order = [0, 2, 1, 3]
        return self.sigma.extract(order, order)

    def detsigma(self):
        """Matrix acting on (x1,x2)^T; composition reverses matrix order.

        If sigma_ij(X)=A_ij X, then (sigma_kl o sigma_ij)(X)
        has matrix A_ij A_kl, since each entry is a scalar.
        This is NOT the ordinary determinant of the 4 by 4 Sigma matrix.
        """
        a11 = self.sigma[:2, :2]
        a12 = self.sigma[:2, 2:]
        a21 = self.sigma[2:, :2]
        a22 = self.sigma[2:, 2:]
        p12, p11 = self.P
        return self.matrix(-p11*a11*a12 + a11*a22 - p12*a21*a12)


def _family(name, P, Q, sigma, m, det, section, page, pi, nonzero=(), relation=None):
    twist = 1 if name == 'N' else h
    restrictions = (() if name == 'N' else (h,)) + tuple(nonzero)
    return Family(name, tuple(map(sp.sympify, P)), tuple(map(sp.sympify, Q)),
                  sp.ImmutableMatrix(sigma)*twist, sp.ImmutableMatrix(m)*twist,
                  sp.ImmutableMatrix(det), (() if relation is None else (relation,)),
                  restrictions, section, page, pi)


def registry():
    """All 26 families, in the original coordinates (in particular for G)."""
    result = {}
    def add(*args, **kw):
        fam = _family(*args, **kw)
        result[fam.name] = fam
    add('A', (1,1), (1,0),
        [[1,0,0,0],[0,1,1,0],[0,0,1,0],[0,-2,-1,1]],
        [[1,0,0,0],[0,1,0,0],[0,1,1,0],[0,-1,-2,1]],
        h**2*E2, '4.1.4',395,'never')
    add('B', (p,0), (p,0),
        [[0,0,0,1],[0,0,1,0],[0,-1,0,0],[1,0,0,0]],
        [[0,0,0,1],[0,0,-1,0],[0,1,0,0],[1,0,0,0]],
        p*h**2*sp.diag(1,-1), '4.2.2',397,'h is a root of unity',
        relation=(p**2+1,p))
    add('C', (p,0), (p,0),
        [[-1,p**2,1,-p],[-p,1,1,-p],[-p,-2*p**2,p,-p],[-p,p**2,1,-1]],
        [[-1,1,p**2,-p],[-p,p,-2*p**2,-p],[-p,1,1,-p],[-p,1,p**2,-1]],
        -3*h**2*sp.diag(p,1),'4.2.2',397,'3*h^2 is a root of unity',
        relation=(p**2+p+1,p))
    add('D', (p,0), (-1,0),
        [[-p,0,0,0],[0,-p**2,1,0],[0,0,p,0],[1,0,0,1]],
        [[-p,0,0,0],[0,p,0,0],[0,1,-p**2,0],[1,0,0,1]],
        -p**2*h**2*E2,'4.2.3',397,'p,h roots of unity and p != +/-1',nonzero=(p,))
    add('E', (p,0), (-1,0),
        [[0,0,1,1],[0,0,1,-1],[-1,1,0,0],[1,1,0,0]],
        [[0,1,0,1],[-1,0,1,0],[0,1,0,-1],[1,0,1,0]],
        2*p*h**2*sp.Matrix([[0,1],[-1,0]]),'4.2.3',397,'2*h^2 is a root of unity',
        relation=(p**2+1,p))
    add('F', (p,0), (-1,0),
        [[-1,-p,1,-1],[-p,1,1,1],[-p,p,p,1],[-p,-p,1,-p]],
        [[-1,1,-p,-1],[-p,p,p,1],[-p,1,1,1],[-p,1,-p,-p]],
        -2*p*h**2*E2,'4.2.3',397,'2*h^2 is a root of unity',relation=(p**2+1,p))
    add('G', (p,0), (1,0),
        [[p,0,0,0],[p,p**2,1,0],[0,0,p,0],[f,0,-1,1]],
        [[p,0,0,0],[0,p,0,0],[p,1,p**2,0],[f,-1,0,1]],
        p**2*h**2*E2,'4.2.4',398,'p,h roots of unity',nonzero=(p,p-1,p+1,f))
    add('H', (-1,0), (1,1),
        [[0,0,1,0],[0,0,f,1],[1,0,0,0],[f,1,0,0]],
        [[0,1,0,0],[1,0,0,0],[0,f,0,1],[f,0,1,0]],
        h**2*sp.Matrix([[1,0],[2*f,1]]),'4.3.1',398,'never')
    add('I', (-1,0), (q,0),
        [[-q,-q,1,-q],[1,1,1,-q],[1,q,q,-q],[-1,-q,1,-1]],
        [[-q,1,-q,-q],[1,q,q,-q],[1,1,1,-q],[-1,1,-q,-1]],
        2*h**2*sp.diag(1,-1),'4.3.2',398,'2*h^2 is a root of unity',
        relation=(q**2+1,q))
    # The printed M for J omits its overall h; tensor reshuffling restores it.
    add('J', (-1,0), (q,0),
        [[0,1,0,1],[-1,0,1,0],[0,1,0,-1],[1,0,1,0]],
        [[0,0,1,1],[0,0,1,-1],[-1,1,0,0],[1,1,0,0]],
        2*h**2*E2,'4.3.2',398,'2*h^2 is a root of unity',relation=(q**2+1,q))
    add('K', (-1,0), (q,0),
        [[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,f,0,0]],
        [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,f,0]],
        h**2*sp.diag(1,f),'4.3.2',398,'q,f,h roots of unity',nonzero=(q,f))
    add('L', (-1,0), (q,0),
        [[0,0,f,0],[0,0,0,1],[f,0,0,0],[0,1,0,0]],
        [[0,f,0,0],[f,0,0,0],[0,0,0,1],[0,0,1,0]],
        h**2*sp.diag(f**2,1),'4.3.2',398,'q,f,h roots of unity',nonzero=(q,f))
    add('M', (-1,0), (-1,0),
        [[0,1,1,0],[f,0,0,-1],[1,0,0,-1],[0,-1,-f,0]],
        [[0,1,1,0],[1,0,0,-1],[f,0,0,-1],[0,-f,-1,0]],
        (1-f)*h**2*E2,'4.3.3',398,'f != 0 and h*(s+/-1) roots of unity',nonzero=(f-1,))
    add('N', (-1,0), (-1,0),
        [[0,-g,0,f],[g,0,f,0],[0,f,0,-g],[f,0,g,0]],
        [[0,0,-g,f],[0,0,f,-g],[g,f,0,0],[f,g,0,0]],
        (f**2-g**2)*E2,'4.3.3',398,'f+g and f-g roots of unity',nonzero=(f**2-g**2,))
    add('O', (-1,0), (-1,0),
        [[1,0,0,f],[0,-1,1,0],[0,f,-1,0],[1,0,0,1]],
        [[1,0,0,f],[0,-1,f,0],[0,1,-1,0],[1,0,0,1]],
        (f-1)*h**2*E2,'4.3.3',398,'f != 0 and h*(1+/-s) roots of unity',nonzero=(f-1,))
    add('P', (-1,0), (-1,0),
        [[0,0,1,f],[0,0,1,1],[1,-f,0,0],[-1,1,0,0]],
        [[0,1,0,f],[1,0,-f,0],[0,1,0,1],[-1,0,1,0]],
        (1-f)*h**2*E2,'4.3.3',399,'f != 0 and h*(1+/-s) roots of unity',nonzero=(f-1,))
    add('Q', (-1,0), (-1,0),
        [[0,0,1,0],[1,1,1,0],[-1,0,0,0],[1,0,-1,1]],
        [[0,1,0,0],[-1,0,0,0],[1,1,1,0],[1,-1,0,1]],
        h**2*sp.diag(-1,1),'4.3.3',399,'h is a root of unity')
    add('R', (-1,0), (-1,0),
        [[1,1,1,0],[0,0,1,0],[0,1,0,0],[0,-1,-1,1]],
        [[1,1,1,0],[0,0,1,0],[0,1,0,0],[0,-1,-1,1]],
        h**2*sp.Matrix([[0,1],[-1,0]]),'4.3.3',399,'h is a root of unity')
    add('S', (-1,0), (-1,0),
        [[-1,1,1,1],[1,-1,1,1],[1,1,-1,1],[1,1,1,-1]],
        [[-1,1,1,1],[1,-1,1,1],[1,1,-1,1],[1,1,1,-1]],
        4*h**2*E2,'4.3.3',399,'2*h is a root of unity')
    add('T', (-1,0), (-1,0),
        [[-1,1,1,1],[1,-1,1,1],[1,1,1,-1],[1,1,-1,1]],
        [[-1,1,1,1],[1,1,1,-1],[1,1,-1,1],[1,-1,1,1]],
        4*h**2*sp.Matrix([[0,1],[1,0]]),'4.3.3',399,'2*h is a root of unity')
    add('U', (-1,0), (-1,0),
        [[-1,1,1,1],[1,1,1,-1],[1,1,-1,1],[1,-1,1,1]],
        [[-1,1,1,1],[1,-1,1,1],[1,1,1,-1],[1,1,-1,1]],
        4*h**2*E2,'4.3.3',399,'2*h is a root of unity')
    add('V', (-1,0), (1,0),
        [[0,1,1,0],[0,1,0,0],[-1,1,0,0],[0,0,0,1]],
        [[0,1,1,0],[-1,0,1,0],[0,0,1,0],[0,0,0,1]],
        h**2*sp.Matrix([[-1,1],[0,1]]),'4.3.4',399,'h is a root of unity')
    add('W', (-1,0), (1,0),
        [[0,f,1,0],[1,0,0,-1],[1,0,0,f],[0,-1,1,0]],
        [[0,1,f,0],[1,0,0,f],[1,0,0,-1],[0,1,-1,0]],
        (f+1)*h**2*E2,'4.3.4',399,'f != 0 and h*(s+/-i) roots of unity',nonzero=(f+1,))
    add('X', (-1,0), (1,0),
        [[0,0,1,0],[0,0,1,1],[1,0,0,0],[1,1,0,0]],
        [[0,1,0,0],[1,0,0,0],[0,1,0,1],[1,0,1,0]],
        h**2*sp.Matrix([[1,0],[2,1]]),'4.3.4',399,'never')
    add('Y', (-1,0), (1,0),
        [[1,0,0,0],[f,-1,1,0],[0,0,1,0],[1,0,f,-1]],
        [[1,0,0,0],[0,1,0,0],[f,1,-1,0],[1,f,0,-1]],
        h**2*E2,'4.3.4',399,'never')
    add('Z', (1,0), (-1,0),
        [[1,0,0,1],[0,1,1,0],[0,f,-1,0],[f,0,0,-1]],
        [[1,0,0,1],[0,-1,f,0],[0,1,1,0],[f,0,0,-1]],
        -(f+1)*h**2*E2,'4.4.4',400,'h*(s+/-i) roots of unity',nonzero=(f,f+1))
    return result


FAMILIES = registry()
