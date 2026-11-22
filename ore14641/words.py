"""Exact Laurent word calculus for the plane bases used in the manuscript.

Each plane is either a coordinate quantum plane or a balanced plane
(v1**2=v2**2). A word is normalized first by the independently checked mixed
relations and then by the plane relations. Inverses are handled algebraically.
The normal form is r**a e**b rho**c z**d, with all exponents integers.
"""
from dataclasses import dataclass
from functools import lru_cache
import sympy as sp
from .families import FAMILIES, p, q, h, s, f, g
from .monomial import MONOMIAL

alpha, beta, H, S, F = sp.symbols('alpha beta H S F')


class WordModel:
    def __init__(self, scalar, x_kind, y_kind, x_q, y_q, mixing, name=''):
        self.scalar = scalar
        self.x_kind,self.y_kind = x_kind,y_kind
        self.x_q,self.y_q = sp.sympify(x_q),sp.sympify(y_q)
        # (input y,input x) -> (output x,output y,coefficient)
        self.mixing = mixing
        self.name = name
        self.normalize = lru_cache(maxsize=32768)(self._normalize)
        self.e1,self.e2,self.z1,self.z2 = [Word(self,((plane,index,1),))
                                          for plane,index in [(0,0),(0,1),(1,0),(1,1)]]

    def crossing(self,a,b,sy,sx):
        if sy==1 and sx==1:
            j,i,c=self.mixing[a,b]
            return j,i,c
        matches=[]
        for (i,j),(jj,ii,c) in self.mixing.items():
            if sy==1 and sx==-1 and j==b and ii==a:
                matches.append((jj,i,1/c))
            elif sy==-1 and sx==1 and i==a and jj==b:
                matches.append((j,ii,1/c))
            elif sy==-1 and sx==-1 and jj==b and ii==a:
                matches.append((j,i,c))
        if len(matches)!=1:
            raise ValueError('The mixed inverse relation is not uniquely determined')
        return matches[0]

    def plane_normal(self, letters,kind,parameter):
        # r^a e^b with r=e1^-1 e2. Append individual signed generators.
        a,b,c=0,0,sp.Integer(1)
        for _,index,sign in letters:
            if index==0:
                aa,bb,cc=0,sign,sp.Integer(1)
            elif kind=='balanced':
                aa,bb,cc=-1,sign,sp.Integer(1)
            elif sign==1:
                aa,bb,cc=1,1,1/parameter
            else:
                aa,bb,cc=-1,-1,sp.Integer(1)
            if kind=='balanced':
                a,b,c=a+(-1 if b%2 else 1)*aa,b+bb,c*cc
            else:
                c*=cc*parameter**(-b*aa)
                a,b=a+aa,b+bb
        return int(a),int(b),self.scalar(c)

    def _normalize(self, letters):
        letters=list(letters)
        coefficient=sp.Integer(1)
        while True:
            for j in range(len(letters)-1):
                if letters[j][0]==1 and letters[j+1][0]==0:
                    _,a,sy=letters[j];_,b,sx=letters[j+1]
                    bb,aa,c=self.crossing(a,b,sy,sx)
                    letters[j:j+2]=[(0,bb,sx),(1,aa,sy)]
                    coefficient*=c
                    break
            else:
                break
        xs=[a for a in letters if a[0]==0]
        ys=[a for a in letters if a[0]==1]
        a,b,cx=self.plane_normal(xs,self.x_kind,self.x_q)
        c,d,cy=self.plane_normal(ys,self.y_kind,self.y_q)
        return (a,b,c,d),self.scalar(coefficient*cx*cy)

    def identity(self):
        return Word(self,())

    def canonical_word(self,exponents):
        a,b,c,d=map(int,exponents)
        r=self.e1**-1*self.e2
        rho=self.z1**-1*self.z2
        return r**a*self.e1**b*rho**c*self.z1**d


@dataclass(frozen=True)
class Word:
    model: WordModel
    letters: tuple
    coefficient: object = sp.Integer(1)

    def __mul__(self,other):
        if not isinstance(other,Word):
            return Word(self.model,self.letters,self.coefficient*sp.sympify(other))
        if other.model is not self.model:
            raise ValueError('Words belong to different models')
        return Word(self.model,self.letters+other.letters,self.coefficient*other.coefficient)

    __rmul__=__mul__

    def __pow__(self,n):
        if not isinstance(n,int):
            raise TypeError('Use integer Laurent exponents')
        if n<0:
            inverse=Word(self.model,tuple((a,b,-c) for a,b,c in reversed(self.letters)),
                         1/self.coefficient)
            return inverse**(-n)
        if n==0:
            return self.model.identity()
        return Word(self.model,self.letters*n,self.coefficient**n)

    def normal(self):
        exponents,coefficient=self.model.normalize(self.letters)
        return exponents,self.model.scalar(coefficient*self.coefficient)

    def equals(self,other):
        a,c=self.normal();b,d=other.normal()
        return a==b and self.model.scalar(c-d)==0

    def conjugate(self,other):
        return self*other*self**-1

    def commutation(self,other):
        a,c=(self*other).normal();b,d=(other*self).normal()
        if a!=b:
            raise ValueError('These monomials do not commute up to a scalar')
        return self.model.scalar(c/d)

    def certificate(self):
        a,c=self.normal()
        return {'exponents_r_e_rho_z':list(a),'coefficient':str(c)}


def mixing_from_datum(datum,lambdas=None):
    lambdas=datum.lambdas if lambdas is None else lambdas
    sx=datum.kind in ('I','III');sy=datum.kind in ('I','II')
    return {(a,b):(1-b if sx else b,1-a if sy else a,lambdas[2*a+b])
            for a in range(2) for b in range(2)}


def monomial_model(name,formal=False):
    datum=MONOMIAL[name]
    family=datum.family
    coordinate_x={'B':p,'J':q,'K':q,'L':q,'N':-1,'W':1}
    coordinate_y={'B':p,'E':p,'L':-1,'P':-1,'Z':1}
    substitutions={}
    lambdas=datum.lambdas
    if formal:
        if name=='B':
            lambdas=(H,H,-H,H);substitutions={H:h}
        elif name=='E':
            lambdas=(alpha,-sp.I*alpha,sp.I*alpha,-alpha)
            substitutions={alpha:h*(1+sp.I)}
        elif name=='J':
            lambdas=(alpha,sp.I*alpha,-sp.I*alpha,-alpha)
            substitutions={alpha:h*(1+sp.I)}
        elif name=='K':
            lambdas=(H,H*S,H,-H*S);substitutions={H:h,S:s}
        elif name=='L':
            lambdas=(H*F,H,H*F,H);substitutions={H:h,F:f}
        elif name=='M':
            lambdas=(sp.I*alpha,-sp.I*beta,sp.I*beta,-sp.I*alpha)
            substitutions={alpha:h*(s+1),beta:h*(s-1)}
        elif name=='N':
            lambdas=(alpha,beta,-beta,-alpha)
            substitutions={alpha:f-g,beta:f+g}
        elif name=='O':
            lambdas=(beta,alpha,alpha,beta)
            substitutions={alpha:h*(1+s),beta:h*(1-s)}
        elif name=='P':
            lambdas=(alpha,beta,beta,alpha)
            substitutions={alpha:h*(1+s),beta:h*(1-s)}
        elif name=='S':
            lambdas=(H,-H,-H,-H);substitutions={H:2*h}
        elif name=='T':
            lambdas=(sp.I*H,-H,-sp.I*H,-H);substitutions={H:2*h}
        elif name=='U':
            lambdas=(sp.I*H,-sp.I*H,-H,-H);substitutions={H:2*h}
        elif name=='W':
            lambdas=(sp.I*alpha,-sp.I*beta,sp.I*beta,-sp.I*alpha)
            substitutions={alpha:h*(s-sp.I),beta:h*(s+sp.I)}
        elif name=='Z':
            lambdas=(sp.I*alpha,-sp.I*beta,-sp.I*beta,sp.I*alpha)
            substitutions={alpha:h*(s-sp.I),beta:h*(s+sp.I)}
    model=WordModel(family.scalar,'coordinate' if name in coordinate_x else 'balanced',
                    'coordinate' if name in coordinate_y else 'balanced',
                    coordinate_x.get(name,1),coordinate_y.get(name,1),
                    mixing_from_datum(datum,lambdas),name)
    model.substitutions=substitutions
    model.family=family
    return model


def normalizer_model(name):
    """Models justified by the separately checked Theta/Rnorm/Fnorm/Inorm identities."""
    scalar=FAMILIES[name].scalar
    if name in ('Q','V'):
        lambdas=(-sp.I*H,H,sp.I*H,H) if name=='Q' else (H,-sp.I*H,H,sp.I*H)
        mixing={(a,b):(b,a,lambdas[2*a+b]) for a in range(2) for b in range(2)}
        model=WordModel(scalar,'coordinate','balanced',H**2,1,mixing,name)
        model.phase_variables=(H,)
        model.parameter_map={H:h}
        return model
    if name=='R':
        mixing={(0,0):(0,0,H),(0,1):(1,0,H),
                (1,0):(1,1,H),(1,1):(0,1,-H)}
        model=WordModel(scalar,'coordinate','coordinate',-1,1,mixing,name)
        model.phase_variables=(H,)
        model.parameter_map={H:h}
        return model
    # a is the first eigenvalue of conjugation by Z1; b its swap analogue.
    a=H
    b=-H if name=='F' else H
    mixing={(0,0):(0,0,a),(0,1):(1,0,-a),(1,0):(1,1,b),(1,1):(0,1,b)}
    model=WordModel(scalar,'coordinate','balanced',-1,1,mixing,name)
    model.phase_variables=(H,)
    model.parameter_map={H:-h*(1-p) if name=='F' else 1/(h*(1-q))}
    return model


class TorusWords:
    """Extract commutation and action matrices from exact Laurent words."""
    def __init__(self,model,generators,names):
        self.model=model
        self.generators=tuple(generators)
        self.names=tuple(names)
        columns=[sp.Matrix(w.normal()[0]) for w in generators]
        self.basis=sp.ImmutableMatrix.hstack(*columns)
        if self.basis.det()==0:
            raise ValueError('Torus generators must have independent exponent vectors')

    def coordinates(self,word):
        exponents,coefficient=word.normal()
        vector=self.basis.inv()*sp.Matrix(exponents)
        if not all(v.is_Integer for v in vector):
            raise ValueError('The word is outside the specified torus lattice')
        powers=tuple(map(int,vector))
        base=self.model.identity()
        for gen,power in zip(self.generators,powers):
            base=base*gen**power
        base_exponents,base_coefficient=base.normal()
        if base_exponents!=exponents:
            raise ValueError('The exponent chart is not additive on this subgroup')
        return powers,self.model.scalar(coefficient/base_coefficient)

    def commutations(self):
        result=sp.ones(4)
        for i in range(4):
            for j in range(i+1,4):
                result[i,j]=self.generators[i].commutation(self.generators[j])
                result[j,i]=self.model.scalar(1/result[i,j])
        return sp.ImmutableMatrix(result)

    def action(self,conjugator):
        images=[self.coordinates(conjugator.conjugate(v)) for v in self.generators]
        return {'matrix':sp.ImmutableMatrix.hstack(*(sp.Matrix(v) for v,c in images)),
                'scalars':tuple(c for v,c in images)}


def standard_torus(model):
    e,z=model.e1,model.z1
    r=e**-1*model.e2;rho=z**-1*model.z2
    return TorusWords(model,(r,rho,e**2,z**2),('r','rho','X','Y'))


def normalizer_torus(model):
    e,z=model.e1,model.z1
    r=e**-1*model.e2;rho=z**-1*model.z2
    if model.name in ('Q','V'):
        return TorusWords(model,(e,model.e2,rho,z**2),('n','Theta','rho','Y')),(z,)
    if model.name=='R':
        return TorusWords(model,(r,e**2,z,model.z2**2),('r','X','Z1','W')),(e,model.z2)
    return TorusWords(model,(r,e**2,z**2,rho**2),('r','X','Y','P')),(e,z,rho)
