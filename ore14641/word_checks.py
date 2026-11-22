"""Word-level certificates for torus relations, conjugations and crossed products."""
from itertools import product
import sympy as sp
from .families import h,p,q,s
from .monomial import MONOMIAL
from .pbw import Algebra
from .words import monomial_model,normalizer_model,standard_torus,TorusWords,alpha,beta,H
from .degrees import symbolic_model


def word_checks(check):
    for name,datum in MONOMIAL.items():
        model=monomial_model(name)
        formal=monomial_model(name,True)
        e,z=model.e1,model.z1
        es=(e,model.e2);zs=(z,model.z2)
        r=e**-1*model.e2;rho=z**-1*model.z2;X=e**2;Y=z**2
        sx=datum.kind in ('I','III');sy=datum.kind in ('I','II')
        sig=lambda b:1-b if sx else b
        tau=lambda a:1-a if sy else a
        lam=lambda a,b:datum.lambdas[2*a+b]
        def verify(identifier,left,right,label):
            check(name+'.'+identifier,'Laurent identities',name,left.equals(right),
                  {'left':left.certificate(),'right':right.certificate()},labels=(label,))
        for a,b in product(range(2),repeat=2):
            verify(f'A1.z_inverse.{a}.{b}',zs[a]*es[b]**-1,
                   (1/lam(tau(a),b))*es[sig(b)]**-1*zs[tau(a)],'lem:A1')
            verify(f'A1.e_inverse.{a}.{b}',es[b]*zs[a]**-1,
                   lam(tau(a),b)*zs[tau(a)]**-1*es[sig(b)],'lem:A1')
            verify(f'A3.square_e.{a}.{b}',zs[a].conjugate(es[b]**2),
                   lam(a,b)*lam(tau(a),b)*es[sig(b)]**2,'eq:A3')
            verify(f'A3.square_z.{a}.{b}',es[b].conjugate(zs[a]**2),
                   (1/(lam(tau(a),b)*lam(tau(a),sig(b))))*zs[tau(a)]**2,'eq:A3')
        for a,b,c in product(range(2),repeat=3):
            verify(f'A2.ratio_e.{a}.{b}.{c}',zs[a].conjugate(es[b]**-1*es[c]),
                   (lam(tau(a),c)/lam(tau(a),b))*es[sig(b)]**-1*es[sig(c)],'eq:A2')
            verify(f'A2.ratio_z.{a}.{b}.{c}',es[b].conjugate(zs[a]**-1*zs[c]),
                   (lam(tau(a),b)/lam(tau(c),b))*zs[tau(a)]**-1*zs[tau(c)],'eq:A2')
        verify('A3.kappa',rho.conjugate(r),
               ((lam(tau(1),1)/lam(tau(1),0))/(lam(tau(0),1)/lam(tau(0),0)))*r,'lem:A3')
        verify('A3.z_r',z.conjugate(r),(lam(tau(0),1)/lam(tau(0),0))*r**(-1 if sx else 1),'lem:A3')
        verify('A3.e_rho',e.conjugate(rho),(lam(tau(0),0)/lam(tau(1),0))*rho**(-1 if sy else 1),'lem:A3')
        # Check the actual transformed plane relation, so the word calculus does
        # not merely assume that a proposed basis is balanced or coordinate.
        A=Algebra(datum.family);x1,x2,y1,y2=A.generators
        xe=[datum.x_basis[j,0]*x1+datum.x_basis[j,1]*x2 for j in range(2)]
        yz=[datum.y_basis[j,0]*y1+datum.y_basis[j,1]*y2 for j in range(2)]
        for side,kind,param,values in [('x',model.x_kind,model.x_q,xe),('y',model.y_kind,model.y_q,yz)]:
            residual=(values[0]**2-values[1]**2 if kind=='balanced'
                      else values[1]*values[0]-param*values[0]*values[1])
            check(f'{name}.plane_chart.{side}','plane charts',name,not residual,
                  residual.certificate(),labels=('lem:A2','tab:mono'))
        actualT,formalT=standard_torus(model),standard_torus(formal)
        am,fm=actualT.commutations(),formalT.commutations()
        residual=datum.family.matrix(am-fm.subs(formal.substitutions,simultaneous=True))
        check(f'{name}.phase_model_matches_source','formal parameter charts',name,not any(residual),
              str(residual),labels=('thm:monodeg',))
        for j,(aw,fw) in enumerate([(model.e1,formal.e1),(model.z1,formal.z1)]):
            aa,fa=actualT.action(aw),formalT.action(fw)
            residual=datum.family.matrix(sp.Matrix(aa['scalars'])-
                       sp.Matrix(fa['scalars']).subs(formal.substitutions,simultaneous=True))
            check(f'{name}.phase_action_matches_source.{j}','formal parameter charts',name,
                  aa['matrix']==fa['matrix'] and not any(residual),str(residual),labels=('thm:monodeg',))

    # Symbolic tensor splitting, including the coordinate-basis corrections.
    for name in ('B','E','J','M','N','O','P','S','T','U','W','Z'):
        model=monomial_model(name,True);base=standard_torus(model)
        r,rho,X,Y=base.generators
        xp=X*rho if name in ('B','J') else X*r if name in ('N','U','W') else X
        yp=Y*r if name in ('B','E') else Y*rho if name in ('P','T','Z') else Y
        split=TorusWords(model,(r,rho,xp,yp),('r','rho','Xprime','Yprime'))
        matrix=split.commutations()
        expected=sp.ones(4)
        first=(-1 if name in ('B','E','J','S','T','U')
               else (alpha/beta)**2 if name in ('M','P','Z') else (beta/alpha)**2)
        second=(H**4 if name in ('B','S','T','U') else -alpha**4 if name in ('E','J')
                else (alpha*beta)**2)
        expected[0,1]=first;expected[1,0]=1/first
        expected[2,3]=1/second;expected[3,2]=second
        residual=sp.Matrix(matrix-expected).applyfunc(model.scalar)
        change=sp.Matrix.hstack(*(sp.Matrix(base.coordinates(w)[0]) for w in split.generators))
        check(f'{name}.tensor_split','symbolic torus splitting',name,
              not any(residual) and abs(change.det())==1,
              {'commutation_matrix':str(matrix),'basis_change':change.tolist(),
               'pair_parameters':[str(first),str(second)],'residual':str(residual)},
              labels=('thm:monodeg','eq:twoparams'))

    for name in ('B','E','J','K','L','M','N','O','P','S','T','U','W','Z','F','I','Q','R','V'):
        model,T,units,labels=symbolic_model(name)
        identity=model.identity()
        inverse_checks=all((w*w**-1).equals(identity) and (w**-1*w).equals(identity)
                           for w in (model.e1,model.e2,model.z1,model.z2))
        check(f'{name}.Laurent_units','Laurent inverse controls',name,inverse_checks,
              labels=('prop:crossed' if name in MONOMIAL else 'lem:cross',))
        valid=True;entries=[]
        for i,j in product(range(len(units)),repeat=2):
            matches=[]
            for k,unit in enumerate(units):
                try:
                    exponent,scalar=T.coordinates(units[i]*units[j]*unit**-1)
                except ValueError:continue
                matches.append((k,exponent,str(scalar)))
            valid &= len(matches)==1
            entries.append({'left':labels[i],'right':labels[j],'result':matches})
        check(f'{name}.crossed_product_closure','crossed-product closure',name,valid,
              {'index_of_torus_lattice':abs(int(T.basis.det())),
               'coset_count':len(units),'multiplication':entries},
              labels=('prop:crossed' if name in MONOMIAL else 'lem:cross',))
        check(f'{name}.coset_index','crossed-product closure',name,
              abs(int(T.basis.det()))==len(units),labels=('prop:crossed' if name in MONOMIAL else 'lem:cross',))
        # Every entry is displayed in the generated matrix catalogue.
        actions=[T.action(w) for w in units]
        unimodular=all(abs(a['matrix'].det())==1 for a in actions)
        check(f'{name}.action_lattice_automorphisms','torus action matrices',name,unimodular,
              {'torus_generators':list(T.names),'commutations':str(T.commutations()),
               'actions':[{'coset':label,'matrix':a['matrix'].tolist(),
                           'scalars':list(map(str,a['scalars']))} for label,a in zip(labels,actions)]},
              labels=('thm:monodeg' if name in MONOMIAL else 'thm:FIQRV',))
