"""Export inspectable inputs and witnesses for the independent checker."""
from pathlib import Path
import argparse
import json
import sympy as sp
from sympy.matrices.normalforms import smith_normal_decomp
from .families import FAMILIES
from .monomial import MONOMIAL
from .uniform import uniform_bundle
from .degrees import center_lattice


def export_bundle(report):
    sources=[]
    for name,family in FAMILIES.items():
        sources.append({'family':name,'P':list(map(str,family.P)),'Q':list(map(str,family.Q)),
           'sigma':[[str(v) for v in row] for row in family.sigma.tolist()],
           'M':[[str(v) for v in row] for row in family.printed_m.tolist()],
           'expected_detsigma':[[str(v) for v in row] for row in family.expected_detsigma.tolist()],
           'relations':[[str(poly),str(var)] for poly,var in family.relations],
           'admissibility':list(map(str,family.nonzero)),'source':{'section':family.section,'page':family.page}})
    cases=[]
    for check in report['checks']:
        c=check.get('detail')
        if check['category'] not in ('all-family degree certificates','center-action parity branches',
                                     'fourth-root degree branches','K/L parity certificates'):continue
        A=sp.Matrix(c['commutation_matrix']);N=c['root_order']
        D,U,V=smith_normal_decomp(A,domain=sp.ZZ)
        Gamma=sp.Matrix(c['center_lattice']);factors=[N//sp.gcd(N,D[i,i]) for i in range(4)]
        kernel=V*sp.diag(*factors);W=kernel.inv()*Gamma
        if not all(v.is_Integer for v in W) or abs(W.det())!=1:raise ArithmeticError('Kernel equivalence must be unimodular')
        integer=lambda mat:[[int(v) for v in row] for row in mat.tolist()]
        cases.append({'id':check['id'],'family':c['family'],'N':N,'A':integer(A),
             'D':integer(D),'U':integer(U),'V':integer(V),'Gamma':integer(Gamma),
             'kernel_change':integer(W),'actions':[
                 {'coset':a['coset'],'exponent_matrix':integer(sp.Matrix(a['exponent_matrix'])),
                  'phases':list(map(int,a['phases']))} for a in c.get('actions',[])],
             'torus_degree':c['pi_degree'],'center_action_order':c['center_action_order'],
             'algebra_degree':c['algebra_pi_degree'],'formula':c['closed_formula']})
    basis_changes=[]
    for name,datum in MONOMIAL.items():
        basis_changes.append({'family':name,'square_root':bool(datum.square_root),
            'E':[[str(v) for v in row] for row in datum.x_basis.tolist()],
            'Z':[[str(v) for v in row] for row in datum.y_basis.tolist()],
            'target':[[str(v) for v in row] for row in datum.target_sigma().tolist()]})
    return {'schema':'ore14641.certificate_bundle.v1','sources':sources,'basis_changes':basis_changes,
            'uniform':uniform_bundle(),'integer_cases':cases,
            'scope':'Inputs and witnesses, never only pass flags or simplified zero residuals.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--reports',default='reports');args=parser.parse_args()
    directory=Path(args.reports)
    report=json.loads((directory/'verification.json').read_text())
    bundle=export_bundle(report)
    (directory/'certificates.json').write_text(json.dumps(bundle,indent=2,ensure_ascii=False)+'\n')
    print('Exported',len(bundle['integer_cases']),'integer certificates and',len(bundle['uniform']['identities']),'generic identities.')


if __name__=='__main__':main()
