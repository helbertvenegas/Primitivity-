"""Computational ingredients of the primitivity arguments (not ideal-theory formalization)."""
from dataclasses import replace
from itertools import product
from math import gcd
import sympy as sp
from .families import FAMILIES,p,q,h
from .words import Word,WordModel,monomial_model,standard_torus,alpha,beta


class WordSum:
    """Finite sums in the exact Laurent word normal forms."""
    def __init__(self,model,terms=None):
        self.model=model
        self.terms={a:c for a,b in (terms or {}).items() if (c:=model.scalar(b))!=0}

    @classmethod
    def from_word(cls,word):
        a,c=word.normal()
        return cls(word.model,{a:c})

    def as_word(self,exponents,coefficient):
        word=self.model.canonical_word(exponents)
        _,factor=word.normal()
        return word*(coefficient/factor)

    def __add__(self,other):
        result=self.terms.copy()
        for a,c in other.terms.items():result[a]=result.get(a,0)+c
        return WordSum(self.model,result)

    def __neg__(self):return WordSum(self.model,{a:-c for a,c in self.terms.items()})
    def __sub__(self,other):return self+-other

    def __mul__(self,other):
        result={}
        for a,c in self.terms.items():
            for b,d in other.terms.items():
                exponent,coefficient=(self.as_word(a,c)*self.as_word(b,d)).normal()
                result[exponent]=result.get(exponent,0)+coefficient
        return WordSum(self.model,result)

    def conjugate(self,unit):
        result=WordSum(self.model)
        for a,c in self.terms.items():
            result=result+WordSum.from_word(unit.conjugate(self.as_word(a,c)))
        return result

    def nonconstant(self):
        return any(any(a) for a in self.terms)

    def certificate(self):
        return [{'exponents':list(a),'coefficient':str(c)} for a,c in sorted(self.terms.items())]


def orbit_polynomial(model,torus,central_word):
    """Product (T-g(c)) over the four crossed-product cosets; all exact."""
    e,z=model.e1,model.z1
    units=(model.identity(),e,z,e*z)
    orbit=[w.conjugate(central_word) for w in units]
    coefficients=[WordSum.from_word(model.identity())]
    for word in orbit:
        value=WordSum.from_word(word)
        next_coefficients=[WordSum(model) for _ in range(len(coefficients)+1)]
        for j,c in enumerate(coefficients):
            next_coefficients[j]=next_coefficients[j]-c*value
            next_coefficients[j+1]=next_coefficients[j+1]+c
        coefficients=next_coefficients
    witnesses=[]
    all_central=True
    for coefficient in coefficients[:-1]:
        invariant=all(not (coefficient.conjugate(w)-coefficient).terms
                      for w in torus.generators+(e,z))
        all_central &= invariant
        if coefficient.nonconstant():witnesses.append(coefficient.certificate())
    return {'central':all_central,'has_nonconstant_coefficient':bool(witnesses),
            'central_monomial':central_word.certificate(),
            'orbit':[w.certificate() for w in orbit],
            'nonconstant_coefficients':witnesses}


def fixed_family_orbit_certificate(name):
    model=monomial_model(name,True);T=standard_torus(model)
    for gen,power in product(T.generators,(1,2,4)):
        central=gen**power
        if all(central.commutation(v)==1 for v in T.generators):
            return orbit_polynomial(model,T,central)
    raise ArithmeticError('No fixed central monomial found in the stated model')


def conditional_orbit_certificate(name,which,root_order):
    """A non-PI torsion branch with the other scalar still transcendental.

    Works for any specified root order. Reduction modulo Phi_n(U) is symbolic;
    the other free phase parameter beta is not specialized to a root of unity.
    """
    if name not in ('M','N','O','P','W','Z') or which not in ('u','v'):
        raise ValueError('Use one of the six primitive monomial families and u/v')
    if not isinstance(root_order,int) or root_order<1:raise ValueError('Positive root order required')
    if which=='u' and name!='N' and root_order<=2:
        raise ValueError('u=+/-1 is outside the admissible nondegenerate chart')
    U=sp.Symbol('U')
    old=monomial_model(name,True)
    substitution={alpha:(-U*beta if name=='Z' else U*beta)} if which=='u' else {alpha:U/beta}
    family=replace(FAMILIES[name],relations=((sp.cyclotomic_poly(root_order,U),U),))
    mixing={key:(a,b,c.subs(substitution,simultaneous=True))
            for key,(a,b,c) in old.mixing.items()}
    model=WordModel(family.scalar,old.x_kind,old.y_kind,old.x_q,old.y_q,mixing,name)
    T=standard_torus(model);r,rho,X,Y=T.generators
    m=root_order//gcd(root_order,2)
    central=r**m if which=='u' else (X*r if name in ('N','W') else X)**m
    if not all(central.commutation(v)==1 for v in T.generators):
        raise ArithmeticError('The conditional central monomial failed')
    result=orbit_polynomial(model,T,central)
    result.update(family=name,root_parameter=which,root_order=root_order,
                  remaining_parameter='beta is symbolic and unconstrained')
    return result


def dg_mixed_center_certificate(name,torsion,root_order):
    if name not in ('D','G') or torsion not in ('p','q') or root_order<1:
        raise ValueError('Use family D/G, torsion p/q, and a positive root order')
    if torsion=='p' and root_order<=2:
        raise ValueError('The D/G birational chart excludes p=+/-1')
    n=root_order;P,Q=sp.symbols('P Q')
    root=P if torsion=='p' else Q
    family=replace(FAMILIES['S'],relations=((sp.cyclotomic_poly(n,root),root),))
    if torsion=='q':
        a=n if name=='G' else 2*n
        vectors=[(a,0,0,0),(0,n,n,n)]
    else:vectors=[(0,2*n,0,0),(0,0,2*n,-2*n)]
    valid=True
    for a,b,c,d in vectors:
        signs=((-1)**(b+c),(-1)**a,(-1)**a,1) if name=='D' else (1,1,1,1)
        values=(signs[0]*Q**(-c-d),signs[1]*P**(d-c),
                signs[2]*Q**a*P**(b-d),Q**a*P**(c-b))
        valid &= all(family.scalar(v-1)==0 for v in values)
    return {'torsion_parameter':torsion,'symbolic_other_parameter':'Q' if torsion=='p' else 'P',
            'vectors':vectors,'root_order':n,'central':valid,
            'rank':sp.Matrix.hstack(*(sp.Matrix(v) for v in vectors)).rank()}


def primitivity_checks(check):
    for name in ('B','E','J','K','L','S','T','U'):
        result=fixed_family_orbit_certificate(name)
        check(f'{name}.central_orbit','nonprimitive orbit invariants',name,
              result['central'] and result['has_nonconstant_coefficient'],result,
              labels=('thm:primmono',))
    for name in ('M','N','O','P','W','Z'):
        model=monomial_model(name,True);T=standard_torus(model)
        u=-alpha/beta if name=='Z' else alpha/beta
        units=(model.e1,model.z1,model.e1*model.z1)
        witnesses=[];valid=True
        for unit in units:
            action=T.action(unit)
            if action['matrix']!=sp.eye(4):
                witnesses.append({'type':'nonidentity exponent matrix',
                                  'matrix':action['matrix'].tolist()})
                continue
            found=None
            for index in (0,1):
                scalar=action['scalars'][index]
                for odd in (-1,1):
                    sign=model.scalar(scalar/u**odd)
                    if sign not in (1,-1):continue
                    row=[T.commutations()[index,j] for j in range(4)]
                    even=[]
                    for value in row:
                        exponents=[k for k in (-4,-2,0,2,4) if model.scalar(value-u**k)==0]
                        if not exponents:break
                        even.append(exponents[0])
                    if len(even)==4:
                        found={'type':'odd/even exponent obstruction','generator':T.names[index],
                               'outer_scalar':str(scalar),'sign':int(sign),'odd_exponent':odd,
                               'inner_commutator_exponents':even,
                               'consequence_if_inner':'u^(2j-odd)=+/-1, forcing u torsion'}
                        break
                if found:break
            valid &= found is not None
            witnesses.append(found)
        check(f'{name}.outer_actions','primitive localization obstructions',name,valid,
              witnesses,labels=('thm:primmono',))
        # The general orbit construction is recorded above. These extra runs
        # certify the mixed PI/non-PI branches with one parameter still free.
        for which,n in [('u',3),('u',4),('v',1),('v',3)]:
            result=conditional_orbit_certificate(name,which,n)
            check(f'{name}.conditional_orbit.{which}.{n}','conditional orbit regressions',name,
                  result['central'] and result['has_nonconstant_coefficient'],result,
                  labels=('thm:primmono',))

    # The center equations of D/G: elimination in the non-torsion case leaves
    # the following full-rank linear system on (a,b,c,d).
    matrix=sp.Matrix([[0,0,-1,1],[0,0,1,1],[0,2,-1,-1],[1,0,0,0]])
    check('DG.center_nontorsion_elimination','DG center lattice certificates','D/G',
          matrix.rank()==4,{'matrix':matrix.tolist(),'determinant':int(matrix.det()),
                           'hypotheses':'neither p nor q=hp is torsion'},labels=('thm:DG',))
    # For the mixed branches, exhibit two independent central exponent vectors
    # and verify the scalar equations in a cyclotomic coefficient field.
    for name in ('D','G'):
        for torsion in ('p','q'):
            certificate=dg_mixed_center_certificate(name,torsion,5)
            check(f'{name}.mixed_center.{torsion}','DG center lattice certificates',name,
                  certificate['central'] and certificate['rank']==2,certificate,labels=('thm:DG',))
