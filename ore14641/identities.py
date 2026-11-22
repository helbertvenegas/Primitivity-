"""Symbolic polynomial and rational-coordinate identities, with manuscript labels."""
from dataclasses import replace
import sympy as sp
from .families import FAMILIES,p,q,f,g,h,s,t
from .pbw import Algebra
from .monomial import MONOMIAL


def ratio_of_elements(left,right):
    if not left or not right or set(left.terms)!=set(right.terms):return None
    word=next(iter(left.terms))
    scalar=left.algebra.scalar(left.terms[word]/right.terms[word])
    return scalar if not left-scalar*right else None


def rational_coordinate_matrix(family,normal_index):
    """Derive Y*t=M(t)Y directly from Sigma, t=n^-1*x'."""
    n=normal_index;other=1-n
    sigma=family.sigma
    S=sp.Matrix([[sigma[2*a+n,2*b+n] for b in range(2)] for a in range(2)])
    cross=sp.Matrix([[sigma[2*a+n,2*b+other] for b in range(2)] for a in range(2)])
    if any(cross):raise ValueError('The selected generator is not normal in this chart')
    U=sp.Matrix([[sigma[2*a+other,2*b+n] for b in range(2)] for a in range(2)])
    V=sp.Matrix([[sigma[2*a+other,2*b+other] for b in range(2)] for a in range(2)])
    return family.matrix(S.inv()*(U+t*V)),family.matrix(S)


def ore_algebra(beta1,beta2,scalar=sp.cancel):
    return Algebra.from_rules({
        (1,0):{(0,1):p,(2,):-p*beta1},
        (2,0):{(0,2):1/p,(1,):-beta2/p},
        (2,1):{(1,2):p}},scalar,('t','y1','y2'))


def extended_polynomial_checks(check):
    def polynomial(identifier,family,residual,label,category='additional PBW identities'):
        check(identifier,category,family,not residual,residual.certificate(),labels=(label,))
    A=Algebra(FAMILIES['L']);x1,x2,y1,y2=A.generators
    uu=x1**2;ww=x1*x2
    for j,y in enumerate((y1,y2)):
        polynomial(f'L.torus.u.{j}','L',y*uu-(h*f)**2*uu*y,'eq:Ltorus')
        polynomial(f'L.torus.w.{j}','L',y*ww-h**2*f*ww*y,'eq:Ltorus')
    polynomial('L.torus.wu','L',ww*uu-q**2*uu*ww,'eq:Ltorus')
    polynomial('L.torus.correct_first_sign','L',ww*y1-(1/(h**2*f))*y1*ww,'thm:L')
    wrong=ww*y1-h**2*f*y1*ww
    check('L.manuscript_first_sign_rejected','manuscript discrepancy controls','L',bool(wrong),
          {'expected':'nonzero residual for the manuscript sign','residual':wrong.certificate(),
           'correction':'a12=-(2*eta+phi)'},labels=('thm:L',))
    for name in ('Q','V'):
        A=Algebra(FAMILIES[name]);x1,x2,y1,y2=A.generators
        theta=(x2*(y2**2-y1**2)-x1*(y1**2-2*y1*y2+y2**2) if name=='Q'
               else x1*(y1**2-y2**2)-x2*(y1**2+y1*y2))
        n=x1 if name=='Q' else x2
        tests=([y1*theta-h*theta*y1,y2*theta-h*theta*y2] if name=='Q'
               else [y1*theta-h*theta*y2,y2*theta+h*theta*y1])
        tests.append(theta*n-h**2*n*theta)
        for j,r in enumerate(tests):polynomial(f'{name}.Theta.{j}',name,r,'lem:Theta')
        for sign in (1,-1):
            z=y1+sign*sp.I*y2
            xn=-sign*sp.I*h if name=='Q' else h
            th=h if name=='Q' else -sign*sp.I*h
            polynomial(f'{name}.model.n.{sign}',name,z*n-xn*n*z,'eq:Qmodel' if name=='Q' else 'eq:Vmodel')
            polynomial(f'{name}.model.theta.{sign}',name,z*theta-th*theta*z,'eq:Qmodel' if name=='Q' else 'eq:Vmodel')

    A=Algebra(FAMILIES['R']);x1,x2,y1,y2=A.generators
    zp,zm=y1+y2,y1-y2;ep,em=x1+x2,x1-x2
    Z1=x1*zp+x2*zm;Z2=x1*y2+x2*y1
    tests=[Z1*x1-h*x1*Z1,Z1*x2-h*x2*Z1,Z2*x1-h*x2*Z2,
           Z2*x2+h*x1*Z2,Z2*Z1-Z1*Z2]
    for j,r in enumerate(tests):polynomial(f'R.normalizers.{j}','R',r,'lem:Rnorm')
    tests=[zp*x1-h*(ep*zp+x2*zm),zp*x2+h*x2*zm,
           zm*x1-h*x1*zp,zm*x2-h*(x1*zp-em*zm)]
    for j,r in enumerate(tests):polynomial(f'R.zbasis.{j}','R',r,'eq:Rz')

    A=Algebra(FAMILIES['F']);x1,x2,y1,y2=A.generators
    gg,gb=x1+p*x2,x1-p*x2;c=h*(1-p)
    Z1=gg*y1-gb*y2;Z2=gg*y1+p*gb*y2
    tests=[Z1*x1+c*x1*Z1,Z1*x2-c*x2*Z1,Z2*x1-c*x2*Z2,
           Z2*x2-c*x1*Z2,Z2**2+p*Z1**2]
    for j,r in enumerate(tests):polynomial(f'F.normalizers.{j}','F',r,'lem:Fnorm')
    tests=[y1*gg-h*(1+p)*gg*y2,y1*gb-h*(-2*gg*y1+(1-p)*gb*y2),
           y2*gb+h*(1+p)*gb*y1,y2*gg-h*((1-p)*gg*y1+2*p*gb*y2),gg**2-gb**2]
    for j,r in enumerate(tests):polynomial(f'F.gbasis.{j}','F',r,'eq:Fg')

    A=Algebra(FAMILIES['I']);x1,x2,y1,y2=A.generators
    w1,w2=y2-q*y1,-y2-q*y1;c=h*(1-q)
    Z1=x1*w1+x2*w2;Z2=x1*w1+q*x2*w2
    tests=[y1*Z1-c*Z1*y1,y2*Z1+c*Z1*y2,y1*Z2-c*Z2*y2,
           y2*Z2-c*Z2*y1,Z2**2-q*Z1**2]
    for j,r in enumerate(tests):polynomial(f'I.normalizers.{j}','I',r,'lem:Inorm')
    tests=[w1*x1+h*(1+q)*x2*w1,w1*x2-h*((1-q)*x1*w1+2*x2*w2),
           w2*x2-h*(1+q)*x1*w2,w2*x1-h*(-2*q*x1*w1+(1-q)*x2*w2),w1**2-w2**2]
    for j,r in enumerate(tests):polynomial(f'I.wbasis.{j}','I',r,'eq:Iw')

    for name in ('D','G'):
        A=Algebra(FAMILIES[name]);x1,x2,y1,y2=A.generators
        if name=='G':x2=x2+x1*(1/(p-1))
        b=(1 if name=='D' else f)/(p**2-1)
        c=(1 if name=='D' else -1)/(p**2-1)
        Xi=x2*y1*y2-x1*(c*y2**2+b/p*y1**2)
        qq=h*p
        weights=[qq**2,-qq**2,-1/qq,1/qq] if name=='D' else [qq**2,qq**2,1/qq,1/qq]
        for j,(gen,weight) in enumerate(zip((x1,x2,y1,y2),weights)):
            polynomial(f'{name}.Xi.normal.{j}',name,Xi*gen-weight*gen*Xi,'rem:DGnormal')

    B1,B2=sp.symbols('beta1 beta2')
    A=ore_algebra(B1,B2);tt,y1,y2=A.generators
    b=-B2/(p**2-1);c=p**2*B1/(p**2-1)
    omega=tt*y1*y2-b/p*y1**2-c*y2**2
    for j,gen in enumerate(A.generators):
        polynomial(f'A.Omega.central.{j}','D/G',omega.commutator(gen),'lem:Omega')
    for word,residual in A.overlaps().items():
        polynomial(f'A.PBW.{word}','D/G',residual,'eq:Arel')

    # Normal elements at which the primitive monomial algebras are localized.
    for name in ('M','N','O','P','W','Z'):
        datum=MONOMIAL[name];A=Algebra(datum.family);x1,x2,y1,y2=A.generators
        es=[datum.x_basis[j,0]*x1+datum.x_basis[j,1]*x2 for j in range(2)]
        zs=[datum.y_basis[j,0]*y1+datum.y_basis[j,1]*y2 for j in range(2)]
        nx=(es[0]**2 if name in ('M','O','P','Z') else es[0]*es[1] if name=='W' else x1*x2)
        ny=(y1*y2 if name=='P' else zs[0]*zs[1] if name=='Z' else zs[0]**2)
        for element_name,element in [('Nx',nx),('Ny',ny)]:
            for j,gen in enumerate(es+zs):
                factor=ratio_of_elements(element*gen,gen*element)
                check(f'{name}.localizer.{element_name}.{j}','localization normality',name,
                      factor is not None,{'factor':str(factor)},labels=('lem:loc',))

    # Anti-isomorphisms: evaluate every defining relation with reversed products.
    target=Algebra(FAMILIES['E']);x1,x2,y1,y2=target.generators
    source=replace(FAMILIES['J'],Q=(1/p,sp.Integer(0)),relations=((p**2+1,p),))
    _anti_relations(check,'EJ',Algebra(source),target,(y1,y2,x1,x2),'eq:EJ')
    target=Algebra(FAMILIES['Z'].specialize({f:s**2}));x1,x2,y1,y2=target.generators
    source=Algebra(FAMILIES['W'].specialize({f:1/s**2,h:s*h}))
    _anti_relations(check,'WZ',source,target,(s*y1+y2,s*(s*y1-y2),x1,x2),'rem:monodeg')


def _anti_relations(check,name,source,target,images,label):
    for word,terms in source.rules.items():
        lhs=images[word[1]]*images[word[0]]
        rhs=target.element(0)
        for v,c in terms.items():
            term=target.element(c)
            for j in reversed(v):term=term*images[j]
            rhs=rhs+term
        residual=lhs-rhs
        check(f'duality.{name}.{word}','anti-isomorphism identities',name,not residual,
              residual.certificate(),labels=(label,))


def rational_center_checks(check):
    """Certify fractions by a common normalizing automorphism for numerator/denominator."""
    for name in ('B','C','F','I','R','K','L','S','H','X'):
        family=FAMILIES[name];A=Algebra(family);x1,x2,y1,y2=A.generators
        if name=='B':
            num=(x1**4+x2**4)**2;den=x1**4*x2**4
            images=(x1,x2,h**8*y1,h**8*y2)
            for j,y in enumerate((y1,y2)):
                for k,(left,right) in enumerate(((x1,x2),(x2,x1))):
                    residual=y*left**4-h**4*right**4*y
                    check(f'B.fourth_power_swap.{j}.{k}','central rational certificates','B',
                          not residual,residual.certificate(),labels=('eq:Bswap',))
        elif name in ('F','R'):
            num=(x1*x2)**2;den=(x1**2+x2**2)**2
            factor=(-2*p*h**2)**2 if name=='F' else h**4
            images=(x1,x2,(1/factor)*y1,(1/factor)*y2)
            # Use right normalization num*y = factor^(-1)*y*num.
        elif name=='I':
            num=(y1*y2)**2;den=(y1**2+y2**2)**2
            factor=(-2*q*h**2)**2
            images=(factor*x1,factor*x2,y1,y2)
        elif name=='C':
            a=y1**2*y2;b=y1**3-y2**3
            num=a**3;den=b**3
            factor=(3*(2*p+1)*h**3)**3
            images=(factor*x1,factor*x2,y1,y2)
        elif name=='K':
            num=f**2*y1**4+y2**4;den=y1**2*y2**2
            images=(h**4*x1,h**4*f**2*x2,y1,y2)
        elif name=='L':
            num=y1**4+y2**4;den=y1**2*y2**2
            images=(h**4*f**4*x1,h**4*x2,y1,y2)
        elif name=='S':
            num=2*(x1**4-6*x1**2*x2**2+x2**4);den=(x1**2+x2**2)**2
            images=(x1,x2,(1/(16*h**4))*y1,(1/(16*h**4))*y2)
        else:
            plus=y1**2+y2**2;minus=y1**2-y2**2
            num=minus**2;den=plus**2
            ff=f if name=='H' else sp.Integer(1)
            for j,(gen,image) in enumerate(zip((x1,x2),(x1,2*ff*x1+x2))):
                for sign,element in [(1,plus),(-1,minus)]:
                    residual=element*gen-sign*h**2*image*element
                    check(f'{name}.S_pm.{sign}.{j}','central rational certificates',name,
                          not residual,residual.certificate(),labels=('prop:HXprim',))
            images=(h**4*x1,h**4*(4*ff*x1+x2),y1,y2)
        label=('prop:Bnotprim' if name=='B' else 'prop:HXprim' if name in ('H','X')
               else 'thm:primmono' if name in ('K','L','S') else 'prop:CFIQRV')
        for which,element in [('numerator',num),('denominator',den)]:
            for j,(gen,image) in enumerate(zip(A.generators,images)):
                # B uses left normalization g*N=phi(g)*N; all others N*g=phi(g)*N.
                residual=(gen*element-element*image if name=='B' else element*gen-image*element)
                check(f'{name}.rational_center.{which}.{j}','central rational certificates',name,
                      not residual,residual.certificate(),labels=(label,))
        residual=num.commutator(den)
        check(f'{name}.rational_center.commuting','central rational certificates',name,
              not residual,residual.certificate(),labels=(label,))


def affine_nonpi_checks(check):
    table=[('X',{},0,sp.eye(2),1),('Y',{},0,sp.Matrix([[f,1],[1,f]]),-1),
           ('D',{p:1},0,sp.Matrix([[0,-1],[1,0]]),1),
           ('D',{p:-1},0,sp.Matrix([[0,1],[-1,0]]),-1),
           ('M',{f:0},1,sp.Matrix([[0,1],[-1,0]]),-1),
           ('P',{f:0},0,sp.diag(-1,1),1),
           ('W',{f:0},0,sp.Matrix([[0,1],[1,0]]),-1)]
    z=sp.Symbol('z')
    for index,(name,subs,normal,C,a) in enumerate(table):
        family=FAMILIES[name].specialize(subs)
        matrix,S=rational_coordinate_matrix(family,normal)
        residual=family.matrix(matrix-a*t*sp.eye(2)-C)
        check(f'{name}.affine.{index}','non-PI affine matrices',name,not any(residual),
              {'derived_matrix':str(matrix),'residual':str(residual)},labels=('prop:noPI',))
        characteristic=sp.factor(C.charpoly(z).as_expr())
        valid=(C.trace()!=0 or C.det()!=0) if a==1 else sp.factor(C.trace()**2-4*C.det())!=0
        check(f'{name}.affine.eigenvalues.{index}','non-PI affine matrices',name,valid,
              {'charpoly':str(characteristic),'a':a},labels=('lem:shiftmat','prop:noPI'))
        square=family.matrix(S*S-h**2*sp.eye(2))
        check(f'{name}.central_power_matrix.{index}','central power certificates',name,
              not any(square),{'S_squared':str(family.matrix(S*S)),
                               'general_exponent':'m=lcm(2,ord(h))'},labels=('prop:HXprim',))
    for name in ('D','G'):
        matrix,_=rational_coordinate_matrix(FAMILIES[name],0)
        if name=='G':
            shift=1/(p-1)
            matrix=matrix.subs(t,t-shift)+shift*sp.eye(2)
        expected=sp.Matrix([[p*t,(-1 if name=='D' else 1)/p],
                            [(1 if name=='D' else f)/p,t/p]])
        residual=FAMILIES[name].matrix(matrix-expected)
        check(f'{name}.rational_matrix_from_Sigma','normal-generator matrices',name,
              not any(residual),str(residual),labels=('eq:MDG',))
    for name,normal,expected in [('Q',0,sp.Matrix([[-1,1-t],[1+t,1]])),
                                 ('V',1,sp.Matrix([[1,t],[1-t,0]]))]:
        matrix,_=rational_coordinate_matrix(FAMILIES[name],normal)
        residual=FAMILIES[name].matrix(matrix-expected)
        check(f'{name}.rational_matrix_from_Sigma','normal-generator matrices',name,
              not any(residual),str(residual),labels=('prop:QV',))
    for name in ('A','Y'):
        _,S=rational_coordinate_matrix(FAMILIES[name],0)
        residual=FAMILIES[name].matrix(S-h*sp.eye(2))
        check(f'{name}.central_power_all_orders','central power certificates',name,
              not any(residual) and FAMILIES[name].Q==(1,0),
              {'normalizing_matrix':str(S),'general_exponent':'ord(h)'},labels=('prop:HXprim',))
