"""Algebraic certificates for general recurrences, module bases and weights."""
from dataclasses import replace
from functools import lru_cache
import sympy as sp
from .families import FAMILIES,p,q,h,t
from .pbw import Algebra
from .identities import ore_algebra


def algebraic_checks(check):
    def scalar(identifier,expression,label,family='general',reducer=sp.cancel):
        residual=reducer(expression)
        check(identifier,'general algebraic certificates',family,residual==0,
              {'residual':str(residual)},labels=(label,))
    # The twist exponent is a 2-cocycle for all bidegrees, not just test degrees.
    ka,ja,kb,jb,kc,jc=sp.symbols('ka ja kb jb kc jc',integer=True)
    scalar('twist.cocycle',ja*kb+(ja+jb)*kc-jb*kc-ja*(kb+kc),'thm:twist')
    # Scalar identities that connect the independently checked normal/torus
    # relations to the PI criterion, rather than just evaluating its table.
    fam=FAMILIES['C']
    ca,cb=3*(1-p)*h**3,3*(2*p+1)*h**3
    for identifier,expression in [('alpha_p2beta',ca-p**2*cb),
                                   ('alpha_square',ca**2+p*(3*h**2)**3),
                                   ('beta_square',cb**2+(3*h**2)**3)]:
        scalar('C.criterion.'+identifier,expression,'thm:C','C',fam.scalar)
    scalar('F.criterion.c_square',(h*(1-p))**2+2*p*h**2,'thm:F','F',FAMILIES['F'].scalar)
    scalar('I.criterion.c_square',(h*(1-q))**2+2*q*h**2,'thm:I','I',FAMILIES['I'].scalar)
    scalar('EJ.criterion.alpha_square',(h*(1+sp.I))**2-2*sp.I*h**2,'tab:main','E/J')
    # The only sign correction in the L matrix changes no alternating Smith
    # invariants: it changes the sign of one entry, and a34 is zero.
    N,theta,phi,eta=sp.symbols('N theta phi eta')
    from .torus import alternating
    wrong=alternating((2*eta+phi,0,2*theta,N/2,2*eta+2*phi,0))
    correct=alternating((-2*eta-phi,0,2*theta,N/2,2*eta+2*phi,0))
    pf=lambda m:m[0,1]*m[2,3]-m[0,2]*m[1,3]+m[0,3]*m[1,2]
    scalar('L.corrected_pfaffian',pf(correct)-theta*N,'thm:L','L')
    scalar('L.sign_correction_same_pfaffian',pf(wrong)-pf(correct),'thm:L','L')
    k,kp=sp.symbols('k kp')
    for name,upper in [('D',(N/2,-N/2-k-kp,-k-kp,-k,k,-k)),
                       ('G',(0,-k-kp,-k-kp,-k,k,-k))]:
        scalar(name+'.general_pfaffian',pf(alternating(upper))-2*k*(k+kp),
               'prop:DGdeg',name)

    l11,l12,l21,l22,c1,c2,d1,d2,R=sp.symbols('l11 l12 l21 l22 c1 c2 d1 d2 R')
    values=(l11,l12,l21,l22)
    for kind in ('0','I','II','III'):
        sx=kind in ('I','III');sy=kind in ('I','II')
        cs=(c1,c2);ds=(d1,d2)
        changed=tuple(values[2*a+b]*ds[a]*cs[b]/(ds[1-a if sy else a]*cs[1-b if sx else b])
                      for a in range(2) for b in range(2))
        def invariant(v):
            a,b,c,d=v
            return {'0':(a,b,c,d),'I':(a*d,b*c),
                    'II':(a*c,b*d,a/b,c/d),'III':(a*b,c*d,a/c,b/d)}[kind]
        for j,(left,right) in enumerate(zip(invariant(changed),invariant(values))):
            scalar(f'rescale.{kind}.{j}',left-right,'lem:rescale')
    # Normalize type I by r^2*l11*l12=1 and u=(r*l11)^(-1).
    fam=replace(FAMILIES['S'],relations=((l11*l12*R**2-1,R),))
    rr,uu=R,1/(R*l11)
    changed=(rr*uu*l11,rr/uu*l12,uu/rr*l21,l22/(rr*uu))
    for j,(left,right) in enumerate(zip(changed,(1,1,l12*l21,l11*l22))):
        scalar(f'rescale.typeI_normalized.{j}',left-right,'lem:rescale',reducer=fam.scalar)

    # Plane bases: all coefficients are arbitrary, with invertible 2 by 2 change.
    a1,b1,a2,b2=sp.symbols('a1 b1 a2 b2')
    delta=a1*b2-a2*b1
    for sign,name in [(-1,'D'),(1,'G')]:
        A=Algebra(FAMILIES[name]);x1,x2,_,_=A.generators
        e1,e2=a1*x1+b1*x2,a2*x1+b2*x2
        u,v,aa,ap=e1**2,e2**2,e1*e2,e2*e1
        central=(u,v,aa+ap,aa*ap)
        for j,element in enumerate(central):
            for k,gen in enumerate((x1,x2)):
                residual=element.commutator(gen)
                check(f'plane.{sign}.central.{j}.{k}','plane lemma identities','general',
                      not residual,residual.certificate(),labels=('lem:plane',))
        residual=aa*ap-u*v
        check(f'plane.{sign}.product','plane lemma identities','general',not residual,
              residual.certificate(),labels=('lem:plane',))
        if sign==-1:
            for j,residual in enumerate([
                delta**2*x1**2-(b2**2*u+b1**2*v-b1*b2*(aa+ap)),
                delta**2*x2**2-(a2**2*u+a1**2*v-a1*a2*(aa+ap))]):
                check(f'plane.inverse_squared_basis.{j}','plane lemma identities','general',
                      not residual,residual.certificate(),labels=('lem:plane',))
    A0,B0,u,v=sp.symbols('A0 B0 u v')
    scalar('plane.power_sum_recurrence',u**2*A0+v**2*B0-(u+v)*(u*A0+v*B0)+u*v*(A0+B0),'lem:plane')

    # Rank-four module over k[ab,a^2+b^2], including both left and right forms.
    A=Algebra(FAMILIES['D']);a,b,_,_=A.generators
    w1,w2=a*b,a**2+b**2
    identities=[a**3-a*w2-b*w1,a**3-w2*a+w1*b,
                b**3-b*w2+a*w1,b**3-w2*b-w1*a,a**4-w2*a**2-w1**2,
                w1*w2-w2*w1,w1**2+a**2*b**2]
    for j,residual in enumerate(identities):
        check(f'rank4.identity.{j}','finite-module certificates','F/I',not residual,
              residual.certificate(),labels=('lem:rank4',))
    z=sp.Symbol('z')
    scalar('rank4.Hilbert',(1+2*z+z**2)/(1-z**2)**2-1/(1-z)**2,'lem:rank4')

    # Dickson diagonalization is generic in p,c12,c21,w. The final root-of-unity
    # step is represented by P=p^n and W=w^n, so the certificate has no n cutoff.
    c12,c21,w,PP,WW,AA=sp.symbols('c12 c21 w P W A')
    alpha=-c12*c21*p**2/(p**2-1)**2
    tt=w+alpha/w
    l1=p*w+alpha/(p*w);l2=w/p+alpha*p/w
    M=sp.Matrix([[p*tt,c12],[c21,tt/p]])
    V=sp.Matrix([[c12,c12],[l1-p*tt,l2-p*tt]])
    residual=(M*V-V*sp.diag(l1,l2)).applyfunc(sp.cancel)
    check('Dickson.generic_eigenvectors','general Dickson certificate','D/G',not any(residual),
          {'eigenvector_matrix':str(V),'determinant':str(sp.factor(V.det())),
           'residual':str(residual)},labels=('lem:dickson',))
    # The lemma itself permits c12=0; a triangular eigenbasis covers that chart.
    Mt=sp.Matrix([[p*t,0],[c21,t/p]])
    Vt=sp.Matrix([[(p-1/p)*t,0],[c21,1]])
    rt=(Mt*Vt-Vt*sp.diag(p*t,t/p)).applyfunc(sp.cancel)
    check('Dickson.triangular_eigenvectors','general Dickson certificate','general',not any(rt),
          {'c12':0,'eigenvector_matrix':str(Vt),'determinant':str(sp.factor(Vt.det())),
           'residual':str(rt)},labels=('lem:dickson',))
    scalar('Dickson.trace',l1+l2-(p+1/p)*tt,'lem:dickson')
    scalar('Dickson.determinant',l1*l2-tt**2+c12*c21,'lem:dickson')
    scalar('Dickson.all_orders_branch1',(PP*WW+AA/(PP*WW)-WW-AA/WW).subs(PP,1),'lem:dickson')
    scalar('Dickson.all_orders_branch2',(WW/PP+AA*PP/WW-WW-AA/WW).subs(PP,1),'lem:dickson')
    # The q-integer induction only needs the geometric recurrence; Q=q^k.
    qq,Q=sp.symbols('qq Q')
    scalar('geometric.kplus1',(1-qq*Q)/(1-qq)-(1-Q)/(1-qq)-Q,'lem:potY')

    # Exact weight transition identities used in lem:hw and lem:lower.
    B1,B2=sp.symbols('beta1 beta2')
    bb=-B2/(p**2-1);cc=p**2*B1/(p**2-1);aa=bb*cc
    lam=lambda v:v+aa/v
    scalar('weights.aident.1',aa*(p-1/p)**2+B1*B2,'eq:aident')
    scalar('weights.aident.2',aa*(p**2-1)-p**2*B1*bb,'eq:aident')
    scalar('weights.trace',lam(p*w)+lam(w/p)-(p+1/p)*lam(w),'lem:hw')
    scalar('weights.determinant',lam(p*w)*lam(w/p)-lam(w)**2+B1*B2,'lem:hw')
    gamma=p*bb/w
    scalar('weights.hw_gamma',(lam(w/p)-lam(w)/p)/B1-gamma,'lem:hw')
    scalar('weights.hw_Omega',gamma*lam(w/p)-cc*gamma**2-bb,'lem:hw')
    gamma2=bb/w
    scalar('weights.lower_gamma',(p*lam(w)-lam(p*w))/(p*B1)-gamma2,'lem:lower')
    scalar('weights.lower_Omega',lam(w)*gamma2-bb-cc*gamma2**2,'lem:lower')
    scalar('weights.involution',lam(aa/w)-lam(w),'lem:lower')
    P0,Q0=sp.symbols('P0 Q0')
    scalar('weights.distinct_weight_factor',lam(P0*w)-lam(Q0*w)
           -(P0-Q0)*(w-aa/(P0*Q0*w)),'prop:weights')

    # C: explicit module reduction and Hilbert-series certificate for rank nine.
    A=Algebra(FAMILIES['C']);_,_,y1,y2=A.generators
    ca,cb=y1**2*y2,y1**3-y2**3
    identities=[ca*cb-cb*ca,ca**3-y1**6*y2**3,
                ca*y1-p*y1*ca,cb*y1-y1*cb,cb*y2-y2*cb,
                y1**5-cb*y1**2-ca*y2**2]
    for j,residual in enumerate(identities):
        check(f'C.rank9.identity.{j}','C rank-nine module','C',not residual,
              residual.certificate(),labels=('prop:bounds',))
    degrees=[i+j for i,j in C_BASIS]
    scalar('C.rank9.Hilbert',sum(z**d for d in degrees)/(1-z**3)**2-1/(1-z)**2,'prop:bounds','C')
    for total in range(10):
        for i in range(total+1):
            j=total-i;reduction=c_module_reduce(i,j)
            reconstruction=A.element(0)
            for (aexp,bexp,basis),coefficient in reduction.items():
                ii,jj=C_BASIS[basis]
                reconstruction=reconstruction+coefficient*ca**aexp*cb**bexp*y1**ii*y2**jj
            residual=y1**i*y2**j-reconstruction
            check(f'C.rank9.reduce.{i}.{j}','C module reduction regressions','C',not residual,
                  residual.certificate(),labels=('prop:bounds',))


C_BASIS=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(3,0),(1,2),(4,0))


@lru_cache(maxsize=None)
def c_module_reduce(i,j):
    """Unbounded exact left k[a,b]-module reduction for a=y1^2*y2,b=y1^3-y2^3."""
    if i<0 or j<0:raise ValueError('Nonnegative exponents required')
    if (i,j) in C_BASIS:return {(0,0,C_BASIS.index((i,j))):sp.Integer(1)}
    if i>=2 and j>=1:
        summands=[(1,0,p**(2-i),c_module_reduce(i-2,j-1))]
    elif j>=3:
        summands=[(0,0,1,c_module_reduce(i+3,j-3)),(0,1,-1,c_module_reduce(i,j-3))]
    elif i>=5:
        summands=[(0,1,1,c_module_reduce(i-3,0)),(1,0,p**(5-i),c_module_reduce(i-5,2))]
    else:raise ArithmeticError('Incomplete C module recurrence')
    result={}
    scalar=FAMILIES['C'].scalar
    for da,db,coefficient,terms in summands:
        for (a,b,k),value in terms.items():
            key=(a+da,b+db,k)
            result[key]=result.get(key,0)+coefficient*value
    return {key:scalar(value) for key,value in result.items() if scalar(value)!=0}


def right_weight_matrix(n,beta1,beta2):
    """Right t on y1^j*y2^(n-j), with coefficients on the left in k[t]."""
    if not isinstance(n,int) or n<0:raise ValueError('n must be a nonnegative integer')
    matrix=sp.zeros(n+1)
    qint=lambda k,v:sum(v**i for i in range(k))
    for j in range(n+1):
        l=n-j
        matrix[j,j]=p**(j-l)*t
        if j:matrix[j-1,j]=-p**(1-l)*beta1*qint(j,p**2)
        if l:matrix[j+1,j]=-p**(l-2)*beta2*qint(l,p**-2)
    return sp.ImmutableMatrix(matrix)


def weight_matrix_checks(check,max_degree=5):
    B1,B2=sp.symbols('beta1 beta2')
    A=ore_algebra(B1,B2);tt,y1,y2=A.generators
    for n in range(max_degree+1):
        matrix=right_weight_matrix(n,B1,B2)
        valid=True;residuals=[]
        for j in range(n+1):
            actual=y1**j*y2**(n-j)*tt
            expected=A.element(0)
            for i in range(n+1):
                coefficient=sp.Poly(matrix[i,j],t)
                for (degree,),value in coefficient.terms():
                    expected=expected+value*tt**degree*y1**i*y2**(n-i)
            residual=actual-expected
            valid &= not residual
            residuals.append(residual.certificate())
        check(f'weights.matrix.{n}','weight matrix regressions','D/G',valid,
              {'matrix':str(matrix),'residuals':residuals},labels=('prop:weights',))
    # Characteristic polynomials over k(p,a,w), no specialization of p to numbers.
    a,w,z=sp.symbols('a w z')
    for n in range(max_degree+1):
        matrix=right_weight_matrix(n,1,-a*(p-1/p)**2).subs(t,w+a/w)
        actual=matrix.charpoly(z).as_expr()
        expected=sp.prod(z-(p**nu*w+a/(p**nu*w)) for nu in range(-n,n+1,2))
        residual=sp.cancel(actual-expected)
        check(f'weights.spectrum.{n}','weight spectrum regressions','D/G',residual==0,
              {'weights':list(range(-n,n+1,2)),'residual':str(residual)},labels=('prop:weights',))


def dg_central_elements(name,order):
    """Construct the displayed central generators for any primitive order >=3, h=1."""
    from math import lcm
    from .families import f
    if name not in ('D','G') or not isinstance(order,int) or order<3:
        raise ValueError('Use family D/G and primitive order >=3')
    family=replace(FAMILIES[name].specialize({h:1}),relations=((sp.cyclotomic_poly(order,p),p),))
    A=Algebra(family);x1,x2,y1,y2=A.generators
    if name=='G':x2=x2+x1*(1/(p-1))
    alpha=(1 if name=='D' else -f)/(p**2-1)**2
    previous,current=sp.Integer(2),t
    for j in range(2,order+1):previous,current=current,sp.expand(t*current-alpha*previous)
    ell=lcm(2,order) if name=='D' else order
    polynomial=sp.Poly(sp.expand(current**(ell//order)),t)
    phi=A.element(0)
    for (k,),coefficient in polynomial.terms():
        sign=(-1)**(k*(k-1)//2) if name=='D' else 1
        phi=phi+family.scalar(coefficient)*sign*x1**(ell-k)*x2**k
    return A,(x1**ell,y1**ell,y2**order,phi)


def dg_central_checks(check):
    for name in ('D','G'):
        for order in (3,4,5,6):
            A,elements=dg_central_elements(name,order)
            for j,element in enumerate(elements):
                for k,gen in enumerate(A.generators):
                    residual=element.commutator(gen)
                    check(f'{name}.central_generators.{order}.{j}.{k}',
                          'DG central generator regressions',name,not residual,
                          residual.certificate(),labels=('thm:DGPI','lem:potY','lem:dickson'))
