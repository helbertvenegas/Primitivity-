"""Traceability to every label in the supplied manuscript snapshot.

The word 'covered' denotes verified computational ingredients, not a formal
proof of the abstract theorem. Open questions remain explicitly open.
"""
from collections import Counter
from .contracts import claim_contract, validate_contracts

REQUIRED_COMPUTATIONS = (
    'sec:prelim','tab:det','tab:mono','tab:degrees','prop:noPI','eq:MDG',
    'lem:Fnormal','eq:Inormal','eq:Cnormal','prop:QV','lem:rescale','lem:plane',
    'lem:A1','eq:A2','eq:A3','lem:A2','lem:A3','prop:crossed','eq:twoparams',
    'thm:monodeg','lem:Theta','eq:Qmodel','eq:Vmodel','lem:Rnorm','eq:Rz',
    'lem:Fnorm','eq:Fg','lem:Inorm','eq:Iw','lem:cross','thm:FIQRV','eq:Ltorus',
    'thm:L','prop:DGdeg','prop:bounds','rem:DGnormal','lem:Omega','eq:Arel',
    'lem:hw','lem:lower','eq:aident','lem:loc','thm:primmono','prop:Bnotprim',
    'prop:CFIQRV','prop:HXprim','thm:DG','thm:twist','thm:C','thm:F','thm:I',
    'tab:main','lem:potY','lem:dickson','lem:rank4','thm:DGPI','prop:weights'
)

DEPENDENCIES = {
    'eq:A1':['lem:A1'],
    'eq:Theta':['lem:Theta'],
    'thm:A':['tab:main','tab:mono','thm:DGPI','thm:F','thm:I','thm:C','prop:QV','prop:noPI'],
    'thm:B':['thm:twist'],
    'thm:mono':['lem:rescale','lem:plane','lem:A1','lem:A3','tab:mono'],
    'cor:mono':['tab:mono'],
    'thm:FIdeg':['lem:Fnorm','lem:Inorm','thm:FIQRV'],
    'thm:QVdeg':['lem:Theta','eq:Qmodel','eq:Vmodel','thm:FIQRV'],
    'thm:Rdeg':['lem:Rnorm','eq:Rz','thm:FIQRV'],
    'cor:det':['tab:det','tab:main'],
    'eq:Cindex':['prop:bounds'],
    'prop:weights':['lem:hw','lem:lower','lem:Omega'],
    'lem:DGmodel':['rem:DGnormal','lem:Omega','eq:MDG'],
    'eq:Dtorus':['rem:DGnormal','lem:Omega','prop:DGdeg'],
    'eq:Gtorus':['rem:DGnormal','lem:Omega','prop:DGdeg'],
    'thm:DGprim':['rem:DGnormal','lem:Omega','prop:weights','lem:lower','thm:DG'],
    'cor:tfinite':['lem:hw'],
    'eq:Bx':['rem:DGnormal','lem:Omega'],
    'eq:Fdec':['lem:rank4','lem:Fnormal'],
    'eq:Fnormal':['lem:Fnormal'],
    'eq:mono':['tab:mono'],
    'eq:Mt':['eq:MDG','prop:QV'],
    'eq:twistscalar':['thm:twist'],
    'sec:normal':['eq:MDG','prop:QV','thm:L'],
    'sec:CFI':['lem:Fnormal','eq:Inormal','eq:Cnormal'],
    'sec:monodeg':['thm:monodeg','prop:crossed'],
    'sec:norm':['thm:FIQRV'],
    'sec:degrees':['prop:bounds','tab:degrees'],
    'sec:det':['tab:det'],
    'sec:noPI':['prop:noPI'],
    'sec:mono':['tab:mono','lem:rescale','lem:plane'],
    'sec:primitive':['thm:primmono','thm:DGprim','prop:CFIQRV'],
    'subsec:HXtorsion':['prop:HXprim'],
    'def:mono':['tab:mono'],
    'rem:monodeg':['thm:monodeg'],
    'rem:bounds':['thm:FIQRV','lem:rank4'],
    'rem:detsigma':['tab:main','tab:det'],
}

THEORETICAL = {
    'lem:fracind':'Clearing Ore denominators and the PBW tensor-product independence argument.',
    'lem:qplane':'Embedding of a quantum plane and its standard PI criterion.',
    'lem:finite':'Posner and a faithful left-regular representation on a finite-dimensional division-ring space.',
    'lem:shift':'Minimal-polynomial argument for affine conjugation on a finite central division algebra.',
    'lem:shiftmat':'Eigenvector reduction to the affine-conjugation lemma; its case matrices are checked.',
    'prop:jordan':'Jordan-plane inclusion for A/H and its characteristic-zero non-PI property.',
    'lem:words':'Dimension/independence from the original PBW basis and invertible normalizers.',
    'lem:centerQT':'The center of the division ring of a quantum torus; implemented lattices use this theorem.',
    'prop:notprimitive':'A primitive affine PI domain over an algebraically closed field cannot have GK dimension four.',
    'prop:rational':'Dixmier--Moeglin rationality implication using the stated Nullstellensatz hypotheses.',
    'rem:LM':'Bibliographic comparison; no separate computational assertion.',
}


def coverage_report(checks):
    validate_contracts(checks,REQUIRED_COMPUTATIONS)
    counts=Counter(label for check in checks for label in check.get('manuscript_labels',()))
    missing=[label for label in REQUIRED_COMPUTATIONS if counts[label]==0]
    if missing:raise ArithmeticError('Missing computational blocks: '+', '.join(missing))
    entries=[]
    for label in SOURCE_LABELS:
        kind=RESULT_KINDS.get(label,'equation/table/section')
        deps=DEPENDENCIES.get(label,[])
        support={name:counts[name] for name in deps if counts[name]}
        if kind=='question':status='open in manuscript; intentionally not assigned a computed answer'
        elif counts[label]:status='computational ingredients checked'
        elif label in THEORETICAL:status='mathematical/theoretical dependency; no machine-proof claim'
        elif support:status='uses checked ingredients; remaining reasoning is mathematical'
        else:status='contextual/bibliographic label; inspect source rather than infer a numerical claim'
        entries.append({'label':label,'kind':kind,'status':status,'direct_checks':counts[label],
                        'supporting_checks':support,'theoretical_dependency':THEORETICAL.get(label),'claim_contract':claim_contract(label,checks)})
    return {'manuscript_label_count':len(SOURCE_LABELS),
            'required_computational_blocks':len(REQUIRED_COMPUTATIONS),
            'missing_required_blocks':missing,'all_21_closed_degree_models_implemented':True,
            'formal_proof_assistant_verification':False,'entries':entries}

SOURCE_LABELS = ['thm:A', 'tab:main', 'thm:B', 'sec:prelim', 'lem:fracind', 'lem:qplane', 'lem:finite', 'lem:shift', 'lem:shiftmat', 'eq:twistscalar', 'thm:twist', 'sec:noPI', 'prop:jordan', 'prop:noPI', 'rem:LM', 'sec:mono', 'def:mono', 'eq:mono', 'lem:rescale', 'lem:plane', 'thm:mono', 'cor:mono', 'tab:mono', 'eq:EJ', 'sec:normal', 'eq:Mt', 'eq:MDG', 'lem:potY', 'lem:dickson', 'thm:DGPI', 'prop:QV', 'thm:L', 'eq:Ltorus', 'sec:CFI', 'lem:rank4', 'lem:Fnormal', 'eq:Fnormal', 'thm:F', 'eq:Fdec', 'thm:I', 'eq:Inormal', 'thm:C', 'eq:Cnormal', 'sec:monodeg', 'prop:crossed', 'lem:A1', 'eq:A1', 'eq:A2', 'eq:A3', 'lem:A2', 'lem:A3', 'thm:monodeg', 'rem:monodeg', 'sec:norm', 'lem:cross', 'lem:words', 'eq:Theta', 'lem:Theta', 'eq:Qmodel', 'eq:Vmodel', 'thm:QVdeg', 'eq:Rz', 'lem:Rnorm', 'thm:Rdeg', 'eq:Fg', 'lem:Fnorm', 'eq:Iw', 'lem:Inorm', 'thm:FIdeg', 'thm:FIQRV', 'rem:bounds', 'sec:degrees', 'prop:bounds', 'eq:Cindex', 'tab:degrees', 'q:degrees', 'sec:det', 'tab:det', 'cor:det', 'rem:detsigma', 'sec:primitive', 'prop:notprimitive', 'prop:rational', 'prop:Bnotprim', 'eq:Bswap', 'prop:CFIQRV', 'eq:Dtorus', 'eq:Gtorus', 'lem:centerQT', 'thm:DG', 'eq:Arel', 'eq:aident', 'eq:Bx', 'lem:Omega', 'lem:DGmodel', 'prop:DGdeg', 'lem:hw', 'cor:tfinite', 'lem:lower', 'prop:weights', 'thm:DGprim', 'rem:DGnormal', 'eq:twoparams', 'lem:loc', 'thm:primmono', 'subsec:HXtorsion', 'prop:HXprim', 'q:remainingprim']

RESULT_KINDS = {'thm:A': 'theorem', 'thm:B': 'theorem', 'lem:fracind': 'lemma', 'lem:qplane': 'lemma', 'lem:finite': 'lemma', 'lem:shift': 'lemma', 'lem:shiftmat': 'lemma', 'thm:twist': 'theorem', 'prop:jordan': 'proposition', 'prop:noPI': 'proposition', 'rem:LM': 'remark', 'lem:rescale': 'lemma', 'lem:plane': 'lemma', 'thm:mono': 'theorem', 'cor:mono': 'corollary', 'lem:potY': 'lemma', 'lem:dickson': 'lemma', 'thm:DGPI': 'theorem', 'prop:QV': 'proposition', 'thm:L': 'theorem', 'lem:rank4': 'lemma', 'lem:Fnormal': 'lemma', 'thm:F': 'theorem', 'thm:I': 'theorem', 'thm:C': 'theorem', 'prop:crossed': 'proposition', 'lem:A1': 'lemma', 'lem:A2': 'lemma', 'lem:A3': 'lemma', 'thm:monodeg': 'theorem', 'rem:monodeg': 'remark', 'lem:cross': 'lemma', 'lem:words': 'lemma', 'lem:Theta': 'lemma', 'thm:QVdeg': 'theorem', 'lem:Rnorm': 'lemma', 'thm:Rdeg': 'theorem', 'lem:Fnorm': 'lemma', 'lem:Inorm': 'lemma', 'thm:FIdeg': 'theorem', 'thm:FIQRV': 'theorem', 'rem:bounds': 'remark', 'prop:bounds': 'proposition', 'q:degrees': 'question', 'cor:det': 'corollary', 'rem:detsigma': 'remark', 'prop:notprimitive': 'proposition', 'prop:rational': 'proposition', 'prop:Bnotprim': 'proposition', 'prop:CFIQRV': 'proposition', 'lem:centerQT': 'lemma', 'thm:DG': 'theorem', 'lem:Omega': 'lemma', 'lem:DGmodel': 'lemma', 'prop:DGdeg': 'proposition', 'lem:hw': 'lemma', 'cor:tfinite': 'corollary', 'lem:lower': 'lemma', 'prop:weights': 'proposition', 'thm:DGprim': 'theorem', 'rem:DGnormal': 'remark', 'lem:loc': 'lemma', 'thm:primmono': 'theorem', 'prop:HXprim': 'proposition', 'q:remainingprim': 'question'}
