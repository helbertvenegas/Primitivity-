"""Claim-level coverage, with hypotheses, inputs, certificate IDs and dependencies."""
from .families import FAMILIES

# Labels sharing an algebraic claim point to the same exact data contract.
GROUPS = [
(('sec:prelim',),'All admissible source parameters; coefficient field characteristic zero.',
 'P,Q,Sigma,M and all critical cubic ambiguities agree with the fixed presentation.',
 ['source matrices','PBW overlaps'],['termination of the degree-lex rewrite order']),
(('tab:det',),'All admissible parameters; determinant composition convention in CONVENTIONS.md.',
 'The determinant automorphism has the displayed two by two matrix and preserves the base-plane relation.',
 ['det sigma table'],['definition of the double-extension determinant']),
(('tab:mono','lem:rescale','lem:plane','lem:A1','eq:A2','eq:A3','lem:A2','lem:A3'),
 'Fourteen monomial families; square-root charts require s^2=f and s!=0.',
 'Invertible basis changes, plane charts, signed inverse rules and rescaling invariants hold exactly.',
 ['monomial basis changes','plane charts','Laurent identities','plane lemma identities','general algebraic certificates'],
 ['PBW independence; no use of a square-root chart at f=0']),
(('prop:noPI',),'A,H,X,Y or the explicitly admissible degenerate specializations.',
 'Affine conjugation matrices and degenerate defining identities used in the non-PI obstruction.',
 ['degenerate presentations','non-PI affine matrices','additional PBW identities'],
 ['characteristic-zero affine conjugation obstruction and Jordan-plane theorem']),
(('eq:MDG','rem:DGnormal','lem:Omega','eq:Arel','eq:aident'),
 'D with p^2!=1, or G with p!=0,+1,-1 and f!=0.',
 'Shifted relations, Omega centrality, Xi normality, and the birational coefficient identities.',
 ['normalization of G','additional PBW identities','normal-generator matrices','general algebraic certificates'],
 ['denominators are nonzero on this parameter chart; PBW independence']),
(('lem:potY','lem:dickson','thm:DGPI'),
 'p^2!=1; for the PI application p has exact finite order and h is torsion.',
 'General eigenvector and geometric recurrences, Dickson root-power identity, and central-generator construction.',
 ['general Dickson certificate','general algebraic certificates','DG central generator regressions','Dickson samples n=3..8'],
 ['Dickson defining identity and induction; finite central-module criterion']),
(('lem:Fnormal','eq:Inormal','eq:Cnormal','thm:F','thm:I','thm:C'),
 'F:p^2=-1; I:q^2=-1; C:p^2+p+1=0; h!=0.',
 'Displayed normal elements have the required mixed scalar relations, including the exact factors of h.',
 ['normal elements','additional PBW identities','general algebraic certificates'],
 ['finite module/quantum plane PI criteria; PBW algebraic independence']),
(('prop:QV','lem:Theta','eq:Qmodel','eq:Vmodel','lem:Rnorm','eq:Rz','lem:Fnorm','eq:Fg','lem:Inorm','eq:Iw'),
 'The relevant Q,V,R,F,I source presentation and all its admissible parameters.',
 'Polynomial normalizers, their squares and mixed relations justify the Laurent model.',
 ['additional PBW identities','central elements at h=1','rational coordinate matrices','anti-isomorphism identities'],
 ['normalizer independence and generation in the division ring']),
(('prop:crossed','lem:cross','thm:monodeg','thm:FIQRV','eq:Ltorus','thm:L','prop:DGdeg','tab:degrees'),
 'A PI specialization; phases are exact integers in a common cyclic group; generated subgroup normalized where needed.',
 'Torus relations, coset closure, center lattices, action images and the 21 closed degree formulas.',
 ['crossed-product closure','torus action matrices','symbolic torus splitting','all-family degree certificates',
  'unbounded algebraic proof objects','unbounded integer ideal certificates'],
 ['center of a quantum torus; alternating Smith theorem; crossed-product degree theorem; docs/UNBOUNDED_PROOFS.md']),
(('prop:bounds','lem:rank4'),
 'C with p of exact order three, or the anti-commuting plane of F/I.',
 'Finite-module reduction identities, Hilbert-series rank, and exact degree-index constraints.',
 ['finite-module certificates','C rank-nine module','C module reduction regressions','C bound certificates'],
 ['PBW dimension comparison; finite division-ring dimension argument; C exact degree remains open']),
(('thm:primmono','prop:Bnotprim','prop:CFIQRV','prop:HXprim','thm:DG','eq:twoparams','lem:loc'),
 'Relevant admissible families; stated torsion/non-torsion hypotheses are retained for each branch.',
 'Central rational/orbit witnesses, localization normality and non-inner-action obstructions.',
 ['central rational certificates','central power certificates','nonprimitive orbit invariants',
  'primitive localization obstructions','conditional orbit regressions','DG center lattice certificates',
  'localization normality','symbolic torus splitting'],
 ['simplicity, localization of prime ideals, rationality and Nullstellensatz arguments in the manuscript']),
(('lem:hw','lem:lower','prop:weights'),
 'p not a root of unity; b*c!=0; w transcendental; Omega regular after flat extension.',
 'Two spectral projectors, nonvanishing branch obstruction and affine interval induction for every n>=0.',
 ['unbounded algebraic proof objects','general algebraic certificates','weight matrix regressions','weight spectrum regressions'],
 ['PBW dimension n+1; domain and flatness; ideal-theoretic consequence separately stated']),
(('thm:twist',),'Fixed mixed matrix and twist h; root-of-unity multiplier zeta.',
 'Bidegree two-cocycle and exact scalar identities supporting the finite-order twist argument.',
 ['general algebraic certificates','additional PBW identities'],
 ['finite Veronese/finite central extension argument in the manuscript']),
(('tab:main',),'All source admissibility conditions; exact algebraic inputs for instantiated examples.',
 'Exact scalar reductions and PI-table instantiations; proofs use the explicitly linked structural ingredients.',
 ['PI table examples','general algebraic certificates'],
 ['the monomial, normalizer, non-PI and C/F/I/D/G arguments; table evaluation alone is not an independent PI proof'])]

CONTRACTS={label:{'hypotheses':hyp,'checked_claim':claim,'certificate_categories':cats,
                  'mathematical_dependencies':deps}
           for labels,hyp,claim,cats,deps in GROUPS for label in labels}


def claim_contract(label,checks):
    selected=[c for c in checks if label in c.get('manuscript_labels',())]
    value=dict(CONTRACTS.get(label,{}))
    if not value:return None
    value['certificate_ids']=[c['id'] for c in selected]
    value['families_actually_checked']=sorted({c['family'] for c in selected})
    value['certificate_locations']=['reports/verification.json#checks','reports/certificates.json','reports/independent_verification.json']
    return value


def validate_contracts(checks,required):
    missing=[label for label in required if label not in CONTRACTS]
    if missing:raise ArithmeticError('Missing claim contracts: '+', '.join(missing))
    ids={c['id'] for c in checks}
    # Merely assigning one check a section label cannot certify source coverage.
    expected={f'{name}.{suffix}' for name in FAMILIES for suffix in ('M','detSigma','detsigma','det_preserves_plane')}
    missing=expected-ids
    if missing:raise ArithmeticError('Missing source obligations: '+', '.join(sorted(missing)))
    for name in FAMILIES:
        overlaps=[c for c in checks if c['category']=='PBW overlaps' and c['family']==name]
        if len(overlaps)!=4:raise ArithmeticError('Four distinct PBW overlaps required for '+name)
    for name in ('B','E','J','K','L','M','N','O','P','S','T','U','W','Z','F','I','Q','R','V','D','G'):
        if not any(c['category']=='all-family degree certificates' and c['family']==name for c in checks):
            raise ArithmeticError('Missing degree model: '+name)
    return True
