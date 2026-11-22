"""Exact Smith normal form certificates for quantum tori at roots of unity.

Convention: v_i v_j = zeta**A[i,j] v_j v_i, with zeta primitive of order N.
This computes the PI degree of THAT TORUS, not automatically its crossed product.
"""
from itertools import combinations, product
from math import gcd, isqrt, lcm, prod
from functools import reduce
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def alternating(upper, size=4):
    upper = tuple(upper)
    if len(upper) != size*(size-1)//2:
        raise ValueError('Incorrect number of upper triangular entries')
    matrix = sp.zeros(size)
    for (i,j),value in zip(combinations(range(size),2),upper):
        matrix[i,j],matrix[j,i] = value,-value
    return sp.ImmutableMatrix(matrix)


def torus_certificate(matrix, order):
    if not isinstance(order,int) or order < 1:
        raise ValueError('The root order must be a positive integer')
    matrix = sp.Matrix(matrix)
    if matrix.rows != matrix.cols or matrix+matrix.T != sp.zeros(matrix.rows):
        raise ValueError('The commutation matrix must be alternating')
    if not all(v.is_Integer for v in matrix):
        raise ValueError('Use integer exponents, never approximate logarithms')
    smith = smith_normal_form(matrix,domain=sp.ZZ)
    diagonal = [abs(int(smith[i,i])) for i in range(matrix.rows)]
    image_size = prod(order//gcd(order,d) for d in diagonal)
    degree = isqrt(image_size)
    if degree*degree != image_size:
        raise ArithmeticError('An alternating pairing must have square image size')
    return {'root_order':order,'commutation_matrix':[list(map(int,matrix.row(i)))
                                                  for i in range(matrix.rows)],
            'smith_diagonal':diagonal,'image_size':image_size,
            'radical_size_mod_N':order**matrix.rows//image_size,'pi_degree':degree}


def torus_from_rational_phases(upper):
    values = tuple(sp.sympify(a) for a in upper)
    if any(a.has(sp.Float) for a in values):
        raise TypeError('Use exact rational phases, not decimal approximations')
    phases = tuple(sp.Rational(a) for a in values)
    order = lcm(*(int(a.q) for a in phases))
    return torus_certificate(alternating([int(a*order) for a in phases]),order)


def brute_radical_size(matrix, order):
    """Independent finite enumeration for small regression examples only."""
    matrix = sp.Matrix(matrix)
    rows = [tuple(map(int,matrix.row(i))) for i in range(matrix.rows)]
    return sum(all(sum(a*b for a,b in zip(row,v))%order == 0 for row in rows)
               for v in product(range(order),repeat=matrix.cols))


def dg_matrix(name, order, p_exponent, h_exponent):
    k,kp = p_exponent,h_exponent
    if name not in ('D','G'):
        raise ValueError('Only D and G share this torus model')
    if order < 3 or order//gcd(order,k) <= 2:
        raise ValueError('This birational model requires p != +/-1')
    if name == 'D':
        if order%2 or gcd(order,order//2,k,kp) != 1:
            raise ValueError('For D, zeta must generate <-1,p,h>')
        return alternating((order//2,-order//2-k-kp,-k-kp,-k,k,-k))
    if gcd(order,k,kp) != 1:
        raise ValueError('For G, zeta must generate <p,h>')
    return alternating((0,-k-kp,-k-kp,-k,k,-k))


def dg_closed_formula(order,p_exponent,h_exponent):
    return order**2//gcd(order,2*p_exponent*(p_exponent+h_exponent))


def l_matrix(order,q_exponent,f_exponent,h_exponent):
    if order < 2 or order%2 or gcd(order,order//2,q_exponent,f_exponent,h_exponent) != 1:
        raise ValueError('For L, zeta must generate <-1,q,f,h> and N is even')
    th,ph,eta = q_exponent,f_exponent,h_exponent
    # w*y1=(h^2*f)^(-1)*y1*w; the manuscript's first upper entry has a sign typo.
    return alternating((-2*eta-ph,0,2*th,order//2,2*eta+2*ph,0))


def l_closed_formula(order,q_exponent,f_exponent,h_exponent):
    # N0 = order(<-1,q,f>); the crossed product contributes a factor of 2.
    n0 = order//gcd(order,order//2,q_exponent,f_exponent)
    return 2*order if (order//n0)%2 else order


def dickson_matrix_residual(order, family_name):
    """Check D_n(M(t),alpha)=D_n(t,alpha)I in Q(p,f,t)/(Phi_n(p)).

    Each call is a symbolic certificate for primitive roots of this specified
    order. A finite list of calls does not prove the lemma for every order.
    """
    from dataclasses import replace
    from .families import FAMILIES, p, f, t
    if order < 3 or family_name not in ('D','G'):
        raise ValueError('Use n >= 3 and family D or G')
    fam = replace(FAMILIES[family_name],relations=((sp.cyclotomic_poly(order,p),p),))
    c12,c21 = (-1/p,1/p) if family_name == 'D' else (1/p,f/p)
    matrix = fam.matrix([[p*t,c12],[c21,t/p]])
    alpha = fam.scalar(-c12*c21*p**2/(p**2-1)**2)
    prev2,prev1 = 2*sp.eye(2),matrix
    scalar2,scalar1 = sp.Integer(2),t
    for _ in range(2,order+1):
        nxt = fam.matrix(matrix*prev1-alpha*prev2)
        scalar_next = fam.scalar(t*scalar1-alpha*scalar2)
        prev2,prev1 = prev1,nxt
        scalar2,scalar1 = scalar1,scalar_next
    return fam.matrix(prev1-scalar1*sp.eye(2))
