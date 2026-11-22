"""The 14 basis changes and lambda scalars in the manuscript's Table tab:mono."""
from dataclasses import dataclass, replace
import sympy as sp
from .families import FAMILIES, f, g, h, s


@dataclass(frozen=True)
class MonomialDatum:
    name: str
    x_basis: sp.ImmutableMatrix
    y_basis: sp.ImmutableMatrix
    kind: str
    lambdas: tuple
    square_root: bool = False

    @property
    def family(self):
        if self.square_root:
            family = FAMILIES[self.name].specialize({f:s**2})
            return replace(family,nonzero=family.nonzero+(s,))
        return FAMILIES[self.name]

    def transformed_sigma(self):
        change = sp.kronecker_product(self.y_basis,self.x_basis)
        return self.family.matrix(change*self.family.sigma*change.inv())

    def target_sigma(self):
        matrix = sp.zeros(4)
        switch_x = self.kind in ('I','III')
        switch_y = self.kind in ('I','II')
        for a in range(2):
            for b in range(2):
                aa,bb = (1-a if switch_y else a),(1-b if switch_x else b)
                matrix[2*a+b,2*aa+bb] = self.lambdas[2*a+b]
        return sp.ImmutableMatrix(matrix)

    def invariants(self):
        a,b,c,d = self.lambdas
        values = {'0':(a,b,c,d),'I':(a*d,b*c),
                  'II':(a*c,b*d,a/b,c/d),'III':(a*b,c*d,a/c,b/d)}[self.kind]
        return tuple(self.family.scalar(x) for x in values)


def monomial_registry():
    eye = sp.eye(2)
    plus = sp.Matrix([[1,1],[1,-1]])
    complex_plus = sp.Matrix([[1,sp.I],[1,-sp.I]])
    result = {}
    def add(name,X,Y,kind,coefficients,square_root=False):
        twist = 1 if name == 'N' else h
        result[name] = MonomialDatum(name,sp.ImmutableMatrix(X),sp.ImmutableMatrix(Y),
                                    kind,tuple(twist*c for c in coefficients),square_root)
    add('B',eye,eye,'I',(1,1,-1,1))
    add('E',complex_plus,eye,'I',(1+sp.I,1-sp.I,-1+sp.I,-1-sp.I))
    add('J',eye,complex_plus,'I',(1+sp.I,-1+sp.I,1-sp.I,-1-sp.I))
    add('K',eye,[[1,1/s],[1,-1/s]],'0',(1,s,1,-s),True)
    add('L',eye,eye,'II',(f,1,f,1))
    add('M',[[1,sp.I/s],[1,-sp.I/s]],complex_plus,'I',
        (sp.I*(s+1),-sp.I*(s-1),sp.I*(s-1),-sp.I*(s+1)),True)
    add('N',eye,plus,'III',(f-g,f+g,-f-g,g-f))
    add('O',[[1,sp.I*s],[1,-sp.I*s]],complex_plus,'I',(1-s,1+s,1+s,1-s),True)
    add('P',[[1,s],[1,-s]],eye,'II',(1+s,1-s,1-s,1+s),True)
    add('S',plus,plus,'0',(2,-2,-2,-2))
    add('T',plus,complex_plus,'II',(2*sp.I,-2,-2*sp.I,-2))
    add('U',complex_plus,plus,'III',(2*sp.I,-2*sp.I,-2,-2))
    add('W',[[1,sp.I*s],[1,-sp.I*s]],plus,'III',
        (sp.I*(s-sp.I),-sp.I*(s+sp.I),sp.I*(s+sp.I),-sp.I*(s-sp.I)),True)
    add('Z',plus,[[1,sp.I/s],[1,-sp.I/s]],'II',
        (sp.I*(s-sp.I),-sp.I*(s+sp.I),-sp.I*(s+sp.I),sp.I*(s-sp.I)),True)
    return result


MONOMIAL = monomial_registry()
