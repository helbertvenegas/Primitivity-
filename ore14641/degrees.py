"""PI degree certificates for all 21 families with a closed formula.

No root-order cutoff: inputs are exact integer phases in a common cyclic
group. Torus matrices and crossed-product actions are derived from Laurent
words. The center lattice is computed by Hermite normal form and its index
independently by Smith normal form. Theoretical division-ring arguments are
identified in the coverage document, not claimed as machine formalization.
"""
from functools import lru_cache
from itertools import product
from math import gcd,prod,lcm
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form
from .families import p,q
from .words import (monomial_model,standard_torus,normalizer_model,normalizer_torus,
                    alpha,beta,H,S,F)
from .torus import alternating,torus_certificate,dg_closed_formula


def phase(expression,order,assignments):
    expression=sp.factor(sp.cancel(expression))
    if expression==1:return 0
    if expression==-1:
        if order%2:raise ValueError('The phase group must contain -1')
        return order//2
    if expression==sp.I:
        if order%4:raise ValueError('The phase group must contain i')
        return order//4
    if expression in assignments:
        value=sp.sympify(assignments[expression])
        if not value.is_Integer:raise ValueError('Phase exponents must be exact integers')
        return int(value)%order
    if expression.is_Mul:
        return sum(phase(a,order,assignments) for a in expression.args)%order
    if expression.is_Pow and expression.exp.is_Integer:
        return int(expression.exp)*phase(expression.base,order,assignments)%order
    raise ValueError(f'Not a Laurent monomial in the supplied phase variables: {expression}')


def center_lattice(matrix,order):
    """Gamma={a: A*a=0 mod N}, using the dual of Z^d+(A^T/N)Z^d."""
    matrix=sp.Matrix(matrix);size=matrix.rows
    hnf=hermite_normal_form(sp.Matrix.hstack(order*sp.eye(size),matrix.T))
    basis=order*hnf.inv().T
    if not all(a.is_Integer for a in basis):
        raise ArithmeticError('The dual center lattice must be integral')
    if any(int(a)%order for a in matrix*basis):
        raise ArithmeticError('Computed basis does not belong to the modular kernel')
    expected=torus_certificate(matrix,order)['image_size']
    if abs(int(basis.det()))!=expected:
        raise ArithmeticError('Hermite center index disagrees with Smith index')
    return sp.ImmutableMatrix(basis)


def cocycle(matrix,a,b):
    return -sum(int(matrix[i,j])*a[j]*b[i]
                for i in range(len(a)) for j in range(i+1,len(a)))


def monomial_multiply(matrix,order,left,right):
    k,a=left;l,b=right
    return ((k+l+cocycle(matrix,a,b))%order,tuple(x+y for x,y in zip(a,b)))


def monomial_power(matrix,order,monomial,n):
    k,a=monomial
    return ((n*k+n*(n-1)//2*cocycle(matrix,a,a))%order,tuple(n*x for x in a))


def apply_action(matrix,order,action,exponent):
    result=(0,(0,)*matrix.rows)
    for j,n in enumerate(exponent):
        image=(action['phases'][j],tuple(map(int,action['matrix'].col(j))))
        result=monomial_multiply(matrix,order,result,monomial_power(matrix,order,image,int(n)))
    return result


def action_signature(matrix,order,action,center):
    return tuple(apply_action(matrix,order,action,tuple(map(int,center.col(j))))
                 for j in range(center.cols))


@lru_cache(maxsize=None)
def symbolic_model(name):
    if name in ('F','I','Q','R','V'):
        model=normalizer_model(name)
        torus,conjugators=normalizer_torus(model)
    else:
        model=monomial_model(name,formal=True)
        torus=standard_torus(model)
        conjugators=(model.e1,model.z1)
    units=[]
    labels=[]
    for powers in product((0,1),repeat=len(conjugators)):
        word=model.identity()
        for w,k in zip(conjugators,powers):word=word*w**k
        units.append(word)
        labels.append(''.join(map(str,powers)))
    return model,torus,tuple(units),tuple(labels)


def closed_degree(name,order,assignments):
    def value(expr):return phase(expr,order,assignments)
    def root_order(expr):return order//gcd(order,value(expr))
    if name=='B':return 8*root_order(H**4)
    if name in ('E','J'):return 8*root_order(-alpha**4)
    if name in ('S','T','U'):return 8*root_order(H**4)
    if name in ('F','I'):return 8*root_order(H**4)
    if name in ('Q','R','V'):return 4*root_order(H**2)
    if name in ('M','N','O','P','W','Z'):
        u=(-alpha/beta if name=='Z' else alpha/beta)
        v=alpha*beta
        factor=4
        if name in ('O','W','Z') and root_order(u)%2 and root_order(v)%2:
            factor=2
        return factor*root_order(u**2)*root_order(v**2)
    if name in ('K','L'):
        ff=S**2 if name=='K' else F
        denominator=gcd(order,value(-1),value(q),value(ff),value(H))
        n=order//denominator
        if name=='K':return 2*n
        n0=order//gcd(order,value(-1),value(q),value(ff))
        return 2*n if (n//n0)%2 else n
    raise KeyError(name)


def degree_certificate(name,order,assignments,include_action_tables=True):
    """Compute from torus + actions, then compare with the manuscript formula.

    The common order may be any multiple of the actual scalar orders. Including
    i in that common cyclic group is convenient for every monomial family.
    """
    if not isinstance(order,int) or order<1:raise ValueError('Positive integer order required')
    if not all(sp.sympify(v).is_Integer for v in assignments.values()):
        raise ValueError('Phase exponents must be exact integers')
    assignments={sp.sympify(k):int(v)%order for k,v in assignments.items()}
    if name in ('D','G'):
        k,kp=assignments[p],assignments[H]
        if order//gcd(order,k)<=2:raise ValueError('D/G model requires p != +/-1')
        if name=='D' and order%2:raise ValueError('D requires a common group containing -1')
        entries=(order//2,-order//2-k-kp,-k-kp,-k,k,-k) if name=='D' else (0,-k-kp,-k-kp,-k,k,-k)
        matrix=alternating(entries)
        cert=torus_certificate(matrix,order)
        divisor=gcd(order,order//2,k,kp) if name=='D' else gcd(order,k,kp)
        formula=dg_closed_formula(order//divisor,k//divisor,kp//divisor)
        return {'family':name,**cert,'center_lattice':list(map(list,center_lattice(matrix,order).tolist())),
                'center_action_order':1,'algebra_pi_degree':cert['pi_degree'],
                'closed_formula':formula,'agrees':cert['pi_degree']==formula}
    if name in ('B','E'):
        if order%4 or (2*assignments[p])%order!=order//2:
            raise ValueError('p must have square -1')
    if name=='J':
        if order%4 or (2*assignments[q])%order!=order//2:
            raise ValueError('q must have square -1')
    if name in ('M','O','P','W','Z'):
        difference=(assignments[alpha]-assignments[beta])%order
        if difference==0 or (order%2==0 and difference==order//2):
            raise ValueError('These root phases give a singular/degenerate basis or h=0')
    model,torus,units,labels=symbolic_model(name)
    commutation=torus.commutations()
    matrix=alternating([phase(commutation[i,j],order,assignments)
                        for i in range(4) for j in range(i+1,4)])
    certificate=torus_certificate(matrix,order)
    center=center_lattice(matrix,order)
    signatures=[];actions=[]
    for unit,label in zip(units,labels):
        data=torus.action(unit)
        action={'matrix':data['matrix'],
                'phases':tuple(phase(c,order,assignments) for c in data['scalars'])}
        if any(int(v)%order for v in action['matrix'].T*matrix*action['matrix']-matrix):
            raise ArithmeticError('Action does not preserve the torus pairing')
        signatures.append(action_signature(matrix,order,action,center))
        actions.append({'coset':label,'exponent_matrix':data['matrix'].tolist(),
                        'phases':list(action['phases'])})
    image_size=len(set(signatures))
    if len(units)==4:
        if image_size not in (2,4):raise ArithmeticError('Unexpected center action image')
    elif image_size!=len(units):
        raise ArithmeticError('The normalizer model requires distinct center actions')
    computed=image_size*certificate['pi_degree']
    formula=closed_degree(name,order,assignments)
    result={'family':name,**certificate,'center_lattice':center.tolist(),
            'torus_generators':list(torus.names),'coset_count':len(units),
            'center_action_order':image_size,'algebra_pi_degree':computed,
            'closed_formula':formula,'agrees':computed==formula}
    if include_action_tables:result['actions']=actions
    return result


def c_bound_certificate(order,p_exponent,beta_exponent):
    """C remains open in the manuscript: return its certified possible degrees."""
    if order%3 or order//gcd(order,p_exponent)!=3:
        raise ValueError('p must be a primitive cube root')
    k,l=p_exponent,beta_exponent
    matrix=alternating((-k,-l-2*k,-l,-l-2*k,-l,0))
    certificate=torus_certificate(matrix,order)
    n=order//gcd(order,k,l)
    if certificate['pi_degree']!=n:raise ArithmeticError('C subtorus degree disagrees')
    ratios=sorted({ell*b for a in (1,3,9) for b in (1,3,9) for ell in (1,3)
                   if a*b*ell*ell==9})
    return {'family':'C','subtorus':certificate,'module_rank':9,
            'possible_degree_ratios':ratios,'possible_pi_degrees':[n*x for x in ratios],
            'exact_degree_is_open':True}


def phases_of_exact_roots(values):
    """Discrete logs in one exact cyclotomic field, without numerical angles."""
    from .pi import exact_root_order
    orders={key:exact_root_order(value) for key,value in values.items()}
    if any(value in (None,0) for value in orders.values()):
        raise ValueError('All phase variables must be assigned exact roots of unity')
    order=lcm(4,*orders.values())
    zeta=sp.exp(2*sp.pi*sp.I/order)
    field=sp.QQ.algebraic_field(zeta)
    root=field.from_sympy(zeta)
    lookup={};power=field.one
    for k in range(order):
        lookup[tuple(power.to_list())]=k
        power*=root
    assignments={key:lookup[tuple(field.from_sympy(sp.simplify(value)).to_list())]
                 for key,value in values.items()}
    return order,assignments


def degree_from_parameters(name,parameters):
    """Exact p,q,f,g,h interface. C returns its three possibilities, not a guess."""
    from .families import f,g,s,h
    from .pi import pi_decision
    parameters={sp.sympify(k):sp.sympify(v) for k,v in parameters.items()}
    defaults={}
    if name in ('B','E','F') and p not in parameters:defaults[p]=sp.I
    if name in ('I','J') and q not in parameters:defaults[q]=sp.I
    if name=='C' and p not in parameters:defaults[p]=(-1+sp.sqrt(3)*sp.I)/2
    parameters={**defaults,**parameters}
    decision=pi_decision(name,parameters)
    if decision['status'] is not True:
        raise ValueError('Parameters are non-PI or not sufficiently specified')
    value=lambda x:sp.simplify(sp.sympify(x).subs(parameters,simultaneous=True))
    root=sp.sqrt(value(f))
    if name in ('D','G'):values={p:value(p),H:value(h)}
    elif name=='B':values={p:value(p),H:value(h)}
    elif name in ('E','J'):
        values={alpha:value(h)*(1+sp.I),(p if name=='E' else q):value(p if name=='E' else q)}
    elif name=='K':values={q:value(q),S:root,H:value(h)}
    elif name=='L':values={q:value(q),F:value(f),H:value(h)}
    elif name=='M':values={alpha:value(h)*(root+1),beta:value(h)*(root-1)}
    elif name=='N':values={alpha:value(f-g),beta:value(f+g)}
    elif name in ('O','P'):values={alpha:value(h)*(1+root),beta:value(h)*(1-root)}
    elif name in ('S','T','U'):values={H:2*value(h)}
    elif name in ('W','Z'):values={alpha:value(h)*(root-sp.I),beta:value(h)*(root+sp.I)}
    elif name=='F':values={H:-value(h*(1-p))}
    elif name=='I':values={H:1/value(h*(1-q))}
    elif name in ('Q','R','V'):values={H:value(h)}
    elif name=='C':values={p:value(p),beta:value(3*(2*p+1)*h**3)}
    else:raise KeyError(name)
    order,assignments=phases_of_exact_roots(values)
    result=(c_bound_certificate(order,assignments[p],assignments[beta]) if name=='C'
            else degree_certificate(name,order,assignments))
    result['exact_parameters']={str(k):str(v) for k,v in parameters.items()}
    result['default_parameter_branches']={str(k):str(v) for k,v in defaults.items()}
    result['phase_assignments']={str(k):int(v) for k,v in assignments.items()}
    return result


def all_degree_checks(check):
    names=('B','E','J','K','L','M','N','O','P','S','T','U','W','Z','F','I','Q','R','V','D','G')
    label=lambda name:('prop:DGdeg' if name in ('D','G') else 'thm:L' if name=='L'
                       else 'thm:FIQRV' if name in ('F','I','Q','R','V') else 'thm:monodeg')
    for name in names:
        for order in (8,12,20,3988):
            assignments={p:order//4,q:order//4,H:1,S:1,F:1,alpha:1,beta:2}
            if name in ('D','G'):assignments[p]=1 if order!=3988 else 4
            certificate=degree_certificate(name,order,assignments)
            check(f'{name}.complete_degree_model.{order}','all-family degree certificates',name,
                  certificate['agrees'],certificate,labels=(label(name),'tab:degrees'))
    for name in ('O','W','Z'):
        for a,b in ((2,6),(5,7),(4,8)):
            certificate=degree_certificate(name,12,{alpha:a,beta:b})
            check(f'{name}.parity_branch.{a}.{b}','center-action parity branches',name,
                  certificate['agrees'],certificate,labels=('thm:monodeg',))
    for name in ('B','E','J'):
        for branch in (1,3):
            assignment={p:3*branch,q:3*branch,H:0,alpha:1}
            certificate=degree_certificate(name,12,assignment)
            check(f'{name}.degree_branch.{branch}','fourth-root degree branches',name,
                  certificate['agrees'],certificate,labels=(label(name),))
    # L's two parity cases and K's same-parameter distinction.
    for name in ('K','L'):
        for assignment in ({q:0,H:0,F:0,S:0},{q:0,H:1,F:0,S:0},{q:1,H:0,F:0,S:0}):
            certificate=degree_certificate(name,4,assignment)
            check(f'{name}.degree_special.{assignment}','K/L parity certificates',name,
                  certificate['agrees'],certificate,labels=('thm:L','tab:degrees'))
    for order,k,l in ((12,4,1),(12,8,3),(60,20,7)):
        certificate=c_bound_certificate(order,k,l)
        check(f'C.bound_model.{order}.{k}.{l}','C bound certificates','C',
              certificate['possible_degree_ratios']==[1,3,9],certificate,labels=('prop:bounds',))
