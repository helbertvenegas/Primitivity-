"""Exact instantiation of the manuscript's PI criterion (not a new proof)."""
import sympy as sp
from .families import FAMILIES, p, q, f, g, h


def exact_root_order(value):
    """Return n for a root of unity, 0 for certified non-torsion, None if unknown.

    For an exact algebraic number, torsion is detected by its minimal polynomial
    being cyclotomic. No numerical tolerance or arbitrary maximum order is used.
    If phi(n)=d, then n <= 2*d*d, so the subsequent order recovery is exhaustive.
    """
    value = sp.sympify(value)
    if value.has(sp.Float):
        raise TypeError('Floating point cannot certify root-of-unity status')
    if value.free_symbols:
        return None
    if value == 0:
        return 0
    if value.is_algebraic is False:
        return 0
    z = sp.Symbol('z')
    try:
        poly = sp.Poly(sp.minpoly(value,z),z,domain=sp.QQ).monic()
    except (sp.polys.polyerrors.NotAlgebraic,NotImplementedError,ValueError):
        return None
    if not poly.is_cyclotomic:
        return 0
    degree = poly.degree()
    for n in range(1,2*degree*degree+1):
        if sp.totient(n) == degree and poly == sp.Poly(sp.cyclotomic_poly(n,z),z,domain=sp.QQ):
            return n
    raise ArithmeticError('Cyclotomic polynomial has no order in the proved bound')


def pi_decision(name, substitutions):
    family = FAMILIES[name].specialize(substitutions)
    substitutions = {sp.sympify(k):sp.sympify(v) for k,v in substitutions.items()}
    value = lambda x: sp.simplify(sp.sympify(x).subs(substitutions,simultaneous=True))
    result = {'family':name,'method':'Instantiation of manuscript Table tab:main',
              'status':None,'required_scalars':[]}
    if name in ('A','H','X','Y'):
        result['status'] = False
        return result
    if name == 'D' and value(p) in (1,-1):
        result['status'] = False
        return result
    if name in ('M','O','P','W') and value(f) == 0:
        result['status'] = False
        return result
    if name in ('B','Q','R','V'):
        scalars = [h]
    elif name == 'C':
        scalars = [3*h**2]
    elif name in ('E','F','I','J'):
        scalars = [2*h**2]
    elif name in ('D','G'):
        scalars = [p,h]
    elif name in ('K','L'):
        scalars = [q,f,h]
    elif name == 'N':
        scalars = [f+g,f-g]
    elif name in ('S','T','U'):
        scalars = [2*h]
    else:
        root = sp.sqrt(value(f))
        if name == 'M':
            scalars = [value(h)*(root+1),value(h)*(root-1)]
        elif name in ('O','P'):
            scalars = [value(h)*(1+root),value(h)*(1-root)]
        else:
            scalars = [value(h)*(root+sp.I),value(h)*(root-sp.I)]
    orders = []
    for scalar in scalars:
        scalar = value(scalar)
        order = exact_root_order(scalar)
        orders.append(order)
        result['required_scalars'].append({'scalar':str(scalar),'root_order':order})
    result['status'] = (False if 0 in orders else None if None in orders else True)
    return result
