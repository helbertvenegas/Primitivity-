"""Run the exact checks and write a human-readable and machine-readable audit."""
from collections import Counter
from dataclasses import replace
from math import gcd
from pathlib import Path
import argparse
import json
import platform
import sympy as sp
from .families import FAMILIES, p, q, f, g, h, s, t
from .pbw import Algebra
from .monomial import MONOMIAL
from .torus import (torus_certificate, alternating, brute_radical_size,
                    dg_matrix, dg_closed_formula, l_matrix, l_closed_formula,
                    dickson_matrix_residual)
from .pi import pi_decision, exact_root_order
from .identities import extended_polynomial_checks,rational_center_checks,affine_nonpi_checks
from .word_checks import word_checks
from .algebraic import algebraic_checks,weight_matrix_checks,dg_central_checks
from .degrees import all_degree_checks
from .primitivity import primitivity_checks
from .uniform import uniform_checks


def run_all(include_degree_samples=True, include_dickson_samples=True):
    checks = []
    def check(identifier,category,family,passed,detail=None,labels=()):
        if not labels:
            mapping={'source matrices':('sec:prelim',),'PBW overlaps':('sec:prelim',),
                     'det sigma table':('tab:det',),'monomial basis changes':('tab:mono',),
                     'central elements at h=1':('prop:QV',),'normalization of G':('eq:MDG',),
                     'rational coordinate matrices':('prop:QV',),
                     'degenerate presentations':('prop:noPI',),'PI table examples':('tab:main',),
                     'Dickson samples n=3..8':('lem:dickson',),
                     'DG degree samples N<=10':('prop:DGdeg',),'L degree samples N<=8':('thm:L',)}
            labels=mapping.get(category,())
            if category=='normal elements':
                labels=({'F':('lem:Fnormal',),'I':('eq:Inormal',),'C':('eq:Cnormal',)}[family])
        checks.append({'id':identifier,'category':category,'family':family,
                       'passed':bool(passed),'detail':detail,'manuscript_labels':list(labels)})
    for name,family in FAMILIES.items():
        algebra = Algebra(family)
        x1,x2,y1,y2 = algebra.generators
        residual = family.matrix(family.reshuffled_m()-family.printed_m)
        check(f'{name}.M','source matrices',name,not any(residual),str(residual))
        det = family.scalar(family.sigma.det())
        check(f'{name}.detSigma','source matrices',name,det != 0,
              {'determinant':str(sp.factor(det)),
               'nonzero_parameters':[str(v) for v in family.nonzero],
               'meaning':'nonzero in the generic admissible parameter field'})
        residual = family.matrix(family.detsigma()-family.expected_detsigma)
        check(f'{name}.detsigma','det sigma table',name,not any(residual),str(residual))
        det_matrix = family.detsigma()
        u = det_matrix[0,0]*x1+det_matrix[0,1]*x2
        v = det_matrix[1,0]*x1+det_matrix[1,1]*x2
        q12,q11 = family.Q
        residual = v*u-q12*u*v-q11*u**2
        check(f'{name}.det_preserves_plane','det sigma table',name,not residual,residual.certificate())
        for word,residual in algebra.overlaps().items():
            check(f'{name}.overlap.{word}','PBW overlaps',name,not residual,residual.certificate())

    for name,datum in MONOMIAL.items():
        family = datum.family
        transformed = datum.transformed_sigma()
        target = datum.target_sigma()
        for row in range(4):
            residual = family.matrix(transformed[row,:]-target[row,:])
            check(f'{name}.monomial.{row+1}','monomial basis changes',name,not any(residual),str(residual))
        check(f'{name}.basis_invertible','monomial basis changes',name,
              datum.x_basis.det() != 0 and datum.y_basis.det() != 0,
              {'det_x':str(datum.x_basis.det()),'det_y':str(datum.y_basis.det()),
               'invariants':[str(v) for v in datum.invariants()],
               'additional_restriction':'s != 0' if datum.square_root else None})

    # Identities used in the PI proofs, at a symbolic twist h.
    for name in ('F','I','C'):
        algebra = Algebra(FAMILIES[name])
        x1,x2,y1,y2 = algebra.generators
        if name == 'F':
            omega = [x1*x2,x1**2+x2**2]
            residuals = [y*w-(-2*p*h**2)*w*y for y in (y1,y2) for w in omega]
        elif name == 'I':
            omega = [y1*y2,y1**2+y2**2]
            residuals = [w*x-(-2*q*h**2)*x*w for x in (x1,x2) for w in omega]
        else:
            a,b = y1**2*y2,y1**3-y2**3
            residuals = [a*x-3*(1-p)*h**3*x*a for x in (x1,x2)]
            residuals += [b*x-3*(2*p+1)*h**3*x*b for x in (x1,x2)]
        for j,residual in enumerate(residuals,1):
            check(f'{name}.normal.{j}','normal elements',name,not residual,residual.certificate())

    for name in ('Q','V'):
        algebra = Algebra(FAMILIES[name].specialize({h:1}))
        x1,x2,y1,y2 = algebra.generators
        elements = [y1**2+y2**2,y1**2*y2**2]
        elements += ([x1**2+x2**2,x1**2*x2**2+sp.Rational(1,2)*x2**4]
                     if name=='Q' else [x2,x1**2*(x1-x2)**2])
        for j,element in enumerate(elements,1):
            for k,gen in enumerate(algebra.generators,1):
                residual = element.commutator(gen)
                check(f'{name}.central.{j}.{k}','central elements at h=1',name,
                      not residual,residual.certificate())

    algebra = Algebra(FAMILIES['G'])
    x1,x2,y1,y2 = algebra.generators
    e2 = x2+(1/(p-1))*x1
    residuals = [y1*x1-h*p*x1*y1,
                 y1*e2-h*(p**2*e2*y1+x1*y2),
                 y2*x1-h*p*x1*y2,
                 y2*e2-h*(f*x1*y1+e2*y2)]
    for j,residual in enumerate(residuals,1):
        check(f'G.coordinate_shift.{j}','normalization of G','G',not residual,residual.certificate())

    matrices = {'Q':sp.Matrix([[-1,1-t],[1+t,1]]),
                'V':sp.Matrix([[1,t],[1-t,0]])}
    for name,matrix in matrices.items():
        image = (matrix**2*(2*sp.eye(2)-matrix**2) if name=='Q'
                 else (matrix**2-matrix)**2)
        polynomial = t**2*(2-t**2) if name=='Q' else (t**2-t)**2
        residual = sp.simplify(image-polynomial*sp.eye(2))
        check(f'{name}.invariant_matrix','rational coordinate matrices',name,not any(residual),str(residual))

    # Degenerate algebras remain admissible; monomial bases must not be used there.
    for name,value in [('D',1),('D',-1),('M',0),('O',0),('P',0),('W',0)]:
        family = FAMILIES[name].specialize({p:value} if name=='D' else {f:value})
        for word,residual in Algebra(family).overlaps().items():
            check(f'{name}.degenerate.{value}.{word}','degenerate presentations',name,
                  not residual,residual.certificate())
    o0,dminus = FAMILIES['O'].specialize({f:0}),FAMILIES['D'].specialize({p:-1})
    check('O0.equals.Dminus','degenerate presentations','O/D',
          o0.sigma==dminus.sigma and o0.P==dminus.P and o0.Q==dminus.Q)

    # Negative controls: exact, nonzero residuals must be detected.
    bad = FAMILIES['Z'].sigma.as_mutable()
    bad[0,0],bad[0,2],bad[3,0],bad[3,2] = 0,h,0,f*h
    wrong = Algebra(replace(FAMILIES['Z'],sigma=bad.as_immutable()))
    for word,residual in wrong.overlaps().items():
        check(f'Z.printed_relations_rejected.{word}','negative controls','Z',bool(residual),
              {'expected':'nonzero','residual':residual.certificate()})
    for name,subs in [('G',{p:1}),('Z',{f:0}),('K',{f:0}),('B',{p:2}),('C',{p:1}),('E',{h:0})]:
        try:
            FAMILIES[name].specialize(subs)
            rejected = False
        except ValueError:
            rejected = True
        check(f'{name}.invalid.{subs}','negative controls',name,rejected)

    # Exact root-of-unity decisions test the table evaluator; they are not proofs
    # of the theorem assigning the table to all double Ore algebras.
    for name,scale in [('C',1/sp.sqrt(3)),('E',1/sp.sqrt(2)),('F',1/sp.sqrt(2)),
                       ('I',1/sp.sqrt(2)),('J',1/sp.sqrt(2)),('S',sp.Rational(1,2)),
                       ('T',sp.Rational(1,2)),('U',sp.Rational(1,2))]:
        unscaled,scaled = pi_decision(name,{h:1}),pi_decision(name,{h:scale})
        check(f'{name}.twist_example','PI table examples',name,
              unscaled['status'] is False and scaled['status'] is True,
              {'h=1':unscaled,'PI_specialization':scaled,'h':str(scale)})
    for name,subs in [('O',{f:9,h:1/(2*sp.sqrt(2))}),('N',{f:2,g:sp.sqrt(3)})]:
        family = FAMILIES[name].specialize(subs)
        decision = pi_decision(name,subs)
        check(f'{name}.detfinite_not_PI','PI table examples',name,
              family.detsigma()==sp.eye(2) and decision['status'] is False,decision)
    for value,expected in [(1,1),(-1,2),(sp.I,4),(-sp.I,4),
                           ((1+sp.I)/sp.sqrt(2),8),(sp.sqrt(2),0),(sp.Rational(3,2),0)]:
        actual = exact_root_order(value)
        check(f'root_order.{value}','root-of-unity arithmetic','scalar',actual==expected,
              {'value':str(value),'order':actual})

    extended_polynomial_checks(check)
    rational_center_checks(check)
    affine_nonpi_checks(check)
    word_checks(check)
    algebraic_checks(check)
    primitivity_checks(check)
    uniform_checks(check)
    # These APIs take arbitrary orders; the listed calls are regression inputs.
    all_degree_checks(check)
    if include_degree_samples:
        weight_matrix_checks(check)
        dg_central_checks(check)

    if include_dickson_samples:
        for order in range(3,9):
            for name in ('D','G'):
                residual = dickson_matrix_residual(order,name)
                check(f'{name}.Dickson.order{order}','Dickson samples n=3..8',name,
                      not any(residual),str(residual))

    if include_degree_samples:
        for name in ('D','G'):
            for order in range(3,11):
                for k in range(order):
                    for kp in range(order):
                        try:
                            matrix = dg_matrix(name,order,k,kp)
                        except ValueError:
                            continue
                        certificate = torus_certificate(matrix,order)
                        predicted = dg_closed_formula(order,k,kp)
                        check(f'{name}.SNF.{order}.{k}.{kp}','DG degree samples N<=10',name,
                              certificate['pi_degree']==predicted,
                              {'N':order,'k':k,'h_exponent':kp,
                               'smith_diagonal':certificate['smith_diagonal'],
                               'computed':certificate['pi_degree'],'formula':predicted})
        for order in range(2,9,2):
            for th in range(order):
                for ph in range(order):
                    for eta in range(order):
                        try:
                            matrix = l_matrix(order,th,ph,eta)
                        except ValueError:
                            continue
                        certificate = torus_certificate(matrix,order)
                        computed = 2*certificate['pi_degree']
                        predicted = l_closed_formula(order,th,ph,eta)
                        check(f'L.SNF.{order}.{th}.{ph}.{eta}','L degree samples N<=8','L',
                              computed==predicted,
                              {'N':order,'theta':th,'phi':ph,'eta':eta,
                               'computed_algebra_degree':computed,'formula':predicted,
                               'torus_certificate':certificate})
        for order in range(2,6):
            matrix = alternating((1,2,0,1,1,2))
            certificate = torus_certificate(matrix,order)
            actual = brute_radical_size(matrix,order)
            check(f'SNF.independent_enumeration.{order}','independent Smith controls','torus',
                  actual==certificate['radical_size_mod_N'],certificate)

    totals = Counter(check['category'] for check in checks)
    failed = [check['id'] for check in checks if not check['passed']]
    from .coverage import coverage_report
    coverage=coverage_report(checks)
    return {'project':'Double Ore 14641 PI verification','version':'0.3.0',
            'environment':{'python':platform.python_version(),'sympy':sp.__version__},
            'summary':{'total':len(checks),'passed':len(checks)-len(failed),'failed':len(failed),
                       'categories':dict(totals),'families':26,'monomial_families':14},
            'failed_checks':failed,
            'scope':{'symbolic':'26 presentations; 14 monomial data; all torus/action models; normalizers; central rational and orbit certificates; localization elements; module and recurrence certificates',
                     'arbitrary_order_apis':'21 closed PI degree computations, C bounds, center lattices, conditional central orbits, central generators and weight matrices',
                     'sampled':'Finite-order regressions are identified by category; they are not the sole basis for symbolic identities',
                     'not_machine_formalized':['abstract division-ring and ideal-theoretic arguments',
                                               'a proof-assistant formalization of the manuscript'],
                     'open_in_manuscript':['exact PI degree of C','the remaining primitivity question']},
            'coverage':coverage,
            'checks':checks}


def write_report(report,directory):
    directory = Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    (directory/'verification.json').write_text(json.dumps(report,indent=2,ensure_ascii=False,default=lambda a:int(a) if getattr(a,'is_Integer',False) else str(a))+'\n')
    (directory/'coverage.json').write_text(json.dumps(report['coverage'],indent=2,ensure_ascii=False)+'\n')
    lines = ['# Exact verification report','',
             f"Checks passed: **{report['summary']['passed']}/{report['summary']['total']}**.",
             f"Failed checks: **{report['summary']['failed']}**.",'',
             '| Category | Checks |','|---|---:|']
    lines += [f'| {name} | {count} |' for name,count in report['summary']['categories'].items()]
    lines += ['', 'The parameter equations are reduced symbolically; no floating point is used.',
              'For the 26 presentations, all four cubic ambiguities are resolved in the admissible parameter field.',
              'The 26 Sigma matrices and independently transcribed M matrices agree after tensor reshuffling.',
              'The overall factor h missing in the printed M for J is restored explicitly.',
              'All 26 det(sigma) matrices agree with the manuscript table.',
              'All 14 monomial basis changes are checked at symbolic parameters.',
              '', '## Scope', '',
              'All 21 closed degree models are implemented with exact torus matrices, center lattices and crossed-product action images, without a root-order cutoff.',
              'The derived degree is compared independently with the manuscript formula. General recurrence and module identities have separate algebraic certificates.',
              'Finite regression inputs remain explicitly labeled. Abstract division-ring, ideal-theoretic and independence proofs are documented mathematical dependencies, not proof-assistant formalizations.',
              'The PI table evaluator still instantiates Table tab:main; it is not presented as an independent proof.',
              'The exact degree of C and remaining primitivity cases stay open, as in the manuscript.',
              '', '## Manuscript correction found by the completed checks','',
              'In the L torus ordered as w,y1,y1^-1*y2,u, the first upper exponent must be -(2*eta+phi).',
              'The manuscript prints +(2*eta+phi). The relation w*y1=(h^2*f)^(-1)*y1*w verifies the correction.',
              'Its Pfaffian and alternating Smith invariants are unchanged, so the stated degree formula remains unchanged.',
              'See docs/MANUSCRIPT_CORRECTIONS.md for the exact replacement.',
              '', '## Negative control for the printed Z relations', '',
              'Replacing x1*y1 by x1*y2 in the two indicated positions produces nonzero cubic overlap residuals.',
              'The four explicit residuals are retained in verification.json.',
              '',f"Python {report['environment']['python']}; SymPy {report['environment']['sympy']}."]
    if report['failed_checks']:
        lines += ['', '## Failed checks','']+report['failed_checks']
    (directory/'verification.md').write_text('\n'.join(lines)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='reports')
    parser.add_argument('--quick',action='store_true',help='omit finite order degree and Dickson samples')
    args = parser.parse_args()
    report = run_all(not args.quick,not args.quick)
    write_report(report,args.output)
    print(json.dumps(report['summary'],indent=2))
    if report['failed_checks']:
        raise SystemExit('FAILED: '+', '.join(report['failed_checks']))


if __name__ == '__main__':
    main()
