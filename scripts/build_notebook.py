"""Build a portable notebook from the readable package sources (standard library)."""
import ast
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def portable_source(source):
    """Remove relative imports and CLI entrypoint, keeping executable definitions."""
    tree = ast.parse(source)
    remove = set()
    for node in ast.walk(tree):
        if isinstance(node,ast.ImportFrom) and node.level:
            remove.update(range(node.lineno,node.end_lineno+1))
    for node in tree.body:
        if (isinstance(node,ast.If) and isinstance(node.test,ast.Compare)
                and isinstance(node.test.left,ast.Name) and node.test.left.id=='__name__'):
            remove.update(range(node.lineno,node.end_lineno+1))
    return ''.join(line for n,line in enumerate(source.splitlines(keepends=True),1) if n not in remove)


def build():
    cells = []
    def add(kind,source,metadata=None):
        source = source.strip()+'\n'
        cell = {'cell_type':kind,'id':hashlib.sha256((str(len(cells))+source).encode()).hexdigest()[:12],
                'metadata':metadata or {},'source':source.splitlines(keepends=True)}
        if kind == 'code':
            cell.update(execution_count=None,outputs=[])
        cells.append(cell)
    def md(source): add('markdown',source)
    def code(source,metadata=None): add('code',source,metadata)

    md(r'''# Exact verification for double Ore extensions of type (14641)

**Companion to the manuscript by James Gómez and Helbert Venegas.**

Source of the 26 presentations: J. J. Zhang and J. Zhang,
*Double extension regular algebras of type (14641)*, J. Algebra 322 (2009),
373--409, Section 4. [Article](https://doi.org/10.1016/j.jalgebra.2009.03.041) ·
[arXiv](https://arxiv.org/abs/0712.2550).

This notebook is self-contained. Upload it to Google Colab and run all cells;
the package files and the source PDF are not needed at runtime. CPU is sufficient.
All mathematical coefficients are exact. The library cells below reproduce the
readable Python sources supplied in the repository and can be collapsed.

Version **0.3.0** adds an independent checker, generic induction and parity
certificates, claim-level hypotheses, and execution in a real Jupyter kernel.
It includes all 21 closed degree models, their torus/action
matrices, normalizers, rational central witnesses, localization identities,
module bases and general recurrence certificates. A coverage ledger links the
checks to all 108 labels in the supplied manuscript. Finite regressions are
identified separately. The exact degree of C and the remaining primitivity
questions remain open, exactly as in the manuscript.

The mathematical division-ring and ideal-theory arguments remain explicit
dependencies, not proof-assistant formalizations. The PI table helper applies
the stated criterion; the degree calculator instead derives a value from the
torus and its center actions and compares it with the closed formula.''')
    code('''import importlib.util
import importlib.metadata
import subprocess
import sys
if importlib.util.find_spec('sympy') is None:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'sympy==1.14.0'])
import sympy as sp
from IPython.display import display, Math, Markdown
print('Python:', sys.version.split()[0], '| SymPy:', sp.__version__)
if sp.__version__ != '1.14.0':
    print('Reference version: SymPy 1.14.0. Install that version for the recorded environment.')''')
    md(r'''## 1. Defining matrices, parameter restrictions, and conventions

We retain the full twist $h$. Family N has no independent $h$.

$$
\begin{pmatrix}y_1x_1\\y_1x_2\\y_2x_1\\y_2x_2\end{pmatrix}
=\Sigma\begin{pmatrix}x_1y_1\\x_2y_1\\x_1y_2\\x_2y_2\end{pmatrix}.
$$

M is obtained by interchanging positions 2 and 3 on both sides. Its entries
were independently transcribed for cross-checking. The missing overall $h$ in
the printed M for J is explicitly restored.

If $\sigma_{ij}(X)=A_{ij}X$, the matrix of $\det\sigma$ is
$$-p_{11}A_{11}A_{12}+A_{11}A_{22}-p_{12}A_{21}A_{12}.$$
Composition reverses the displayed coefficient matrix order. This object is
different from the ordinary scalar determinant of the $4\times4$ matrix $\Sigma$.''')
    modules = [('families.py','Defining data'),('pbw.py','Noncommutative PBW reduction'),
               ('monomial.py','Fourteen monomial basis changes'),
               ('torus.py','Smith normal form and Dickson samples'),
               ('pi.py','Exact instantiation of the PI table'),
               ('words.py','Exact Laurent words and source-derived action matrices'),
               ('identities.py','Normalizer, affine and rational center identities'),
               ('degrees.py','All 21 degree models and exact center lattices'),
               ('word_checks.py','Crossed-product and tensor-splitting certificates'),
               ('algebraic.py','General recurrences, finite modules and weight matrices'),
               ('primitivity.py','Central orbit and outer-action certificates'),
               ('uniform.py','Unbounded spectral, integer ideal and parity certificates'),
               ('contracts.py','Claim-level hypotheses and certificate contracts'),
               ('coverage.py','Manuscript coverage ledger and proof dependencies'),
               ('verify.py','Auditable verification suite'),
               ('certificates.py','Independent certificate bundle exporter')]
    for filename,label in modules:
        source = (ROOT/'ore14641'/filename).read_text()
        code('# '+label+'\n'+portable_source(source),
             {'source_file':'ore14641/'+filename,
              'source_sha256':hashlib.sha256(source.encode()).hexdigest(),
              'jupyter':{'source_hidden':True}})
    md(r'''## 2. Inspect a family

Change `family_name` to any letter A--Z. The displayed matrices are in the
original Section 4 coordinates. In particular, G has not silently been replaced
by its shifted presentation. The x-plane and y-plane relations are determined
by Q and P. The four displayed mixed rules follow directly from the matrix.''')
    code('''family_name = 'Z'
family = FAMILIES[family_name]
print('Family', family_name, '| ZZ09 subcase', family.section, '| printed page', family.page)
print('P =', family.P, '| Q =', family.Q)
print('Nonzero restrictions:', ', '.join(map(str,family.nonzero)))
print('Parameter equations:', [(str(a), str(v)) for a,v in family.relations])
display(Math(r'\\Sigma = '+sp.latex(family.sigma)))
display(Math(r'M = '+sp.latex(family.printed_m)))
display(Math(r'\\det(\\sigma)|_{A_1} = '+sp.latex(family.detsigma())))
print('Ordinary det(Sigma) =', sp.factor(family.scalar(family.sigma.det())))
algebra = Algebra(family)
for lhs,rhs in algebra.rules.items():
    if lhs[0] >= 2 and lhs[1] < 2:
        print(' '.join(('x1','x2','y1','y2')[i] for i in lhs),'=',Element(algebra,rhs))''')
    md(r'''## 3. Run all exact checks

For every family there are four critical cubic ambiguities. All reductions
decrease degree-lex order, so their resolution checks confluence of the six
quadratic rules. Parameter equations are reduced symbolically, covering both
allowed fourth-root or cube-root branches. No random floating point test is used.

The generic Dickson eigenvector and root-power identities are certified
symbolically for the general argument. Additional explicit order checks are
regressions. The degree APIs accept arbitrary common root orders, and derive
center-action images from exact modular lattices rather than a bounded search.
Finite lists of input orders are still labeled as regression checks.''')
    code('''report = run_all()
display(Markdown('| Category | Checks |\\n|---|---:|\\n'+'\\n'.join(
    f'| {name} | {count} |' for name,count in report['summary']['categories'].items())))
print('PASS:',report['summary']['passed'],'/',report['summary']['total'])
print('FAILED:',report['failed_checks'])
if report['failed_checks']:
    raise RuntimeError('Verification failed: '+str(report['failed_checks']))''')
    md(r'''## 3b. Independent checker, induction and parity certificates

This checker uses only exact integers and Python's `Fraction` arithmetic.
The proof objects contain both sides of each identity, signed source words,
and unimodular Smith/kernel witnesses, rather than only zero residuals.
The uniform argument and its mathematical assumptions are documented in
`UNBOUNDED_PROOFS.md` in the companion package.''')
    checker_source=(ROOT/'independent/check.py').read_text()
    code('from pathlib import Path\n'+
         'Path("independent_check.py").write_text('+repr(checker_source)+')\n'+
         '''import json, subprocess, sys
bundle=export_bundle(report)
Path('independent_certificates.json').write_text(json.dumps(bundle,indent=2))
result=subprocess.run([sys.executable,'-I','independent_check.py',
                'independent_certificates.json','independent_report.json'],
                check=True,capture_output=True,text=True)
print(result.stdout)
independent_report=json.loads(Path('independent_report.json').read_text())
assert independent_report['failed']==0
assert 'ore14641' not in sys.modules
print('Primary and independent verification completed without repository imports.')''',
         {'independent_source_file':'independent/check.py',
          'independent_source_sha256':hashlib.sha256(checker_source.encode()).hexdigest()})
    md(r'''## 4. Inspect a monomial basis change

Rows of E and Z give the new x- and y-bases. The code uses
$(Z\otimes E)\Sigma(Z\otimes E)^{-1}$ and compares all four rows with the
claimed monomial relations. For square-root bases, $f=s^2$ and $s\ne0$.
Do not use these bases at the degenerate $f=0$ endpoints.''')
    code('''monomial_family = 'O'
datum = MONOMIAL[monomial_family]
print('Family',monomial_family,'| type',datum.kind)
display(Math('E = '+sp.latex(datum.x_basis)))
display(Math('Z = '+sp.latex(datum.y_basis)))
display(Math(r"\\Sigma' = "+sp.latex(datum.transformed_sigma())))
print('Full lambdas:',datum.lambdas)
print('Invariants:',datum.invariants())
assert not any(datum.family.matrix(datum.transformed_sigma()-datum.target_sigma()))''')
    md(r'''## 5. A normal-element identity in F

With $\omega_1=x_1x_2$, $\omega_2=x_1^2+x_2^2$, the manuscript uses
$y_j\omega_k=(-2ph^2)\omega_k y_j$. The products below are computed in the
noncommutative algebra, rather than with commuting SymPy generator symbols.''')
    code('''A = Algebra(FAMILIES['F'])
x1,x2,y1,y2 = A.generators
omega = (x1*x2,x1**2+x2**2)
for j,y in enumerate((y1,y2),1):
    for k,w in enumerate(omega,1):
        residual = y*w-(-2*p*h**2)*w*y
        print(f'y{j} omega{k} - (-2 p h^2) omega{k} y{j} =',residual)
        assert not residual''')
    md(r'''## 6. Smith normal form and PI degrees

For an integer alternating exponent matrix A and a primitive Nth root,
$|\operatorname{im}(A\bmod N)|=\prod_i N/\gcd(N,d_i)$, where the $d_i$ are
the integer Smith entries. The torus PI degree is the square root of this index.

The first two examples use the manuscript's birational torus identifications
for D and G. The L algebra degree needs a **factor 2** from its crossed product;
the torus degree alone is not the answer for L.''')
    code('''for family_name,N,k,kp in [('G',5,1,2),('D',8,1,3)]:
    matrix = dg_matrix(family_name,N,k,kp)
    certificate = torus_certificate(matrix,N)
    display(Math('A_{'+family_name+'} = '+sp.latex(matrix)))
    print(family_name,certificate)
    assert certificate['pi_degree']==dg_closed_formula(N,k,kp)
N,theta,phi,eta = 4,0,0,1
certificate = torus_certificate(l_matrix(N,theta,phi,eta),N)
print('L torus PI degree:',certificate['pi_degree'])
print('L algebra PI degree:',2*certificate['pi_degree'])
assert 2*certificate['pi_degree']==l_closed_formula(N,theta,phi,eta)''')
    md(r'''## 7. Why the twist cannot be normalized to 1

The table evaluator below applies the manuscript's PI criterion to exact
parameters. Root-of-unity status is determined from cyclotomic minimal
polynomials; this is not a truncated search for an integer power close to 1.''')
    code('''examples = [('C',1/sp.sqrt(3)),('E',1/sp.sqrt(2)),('F',1/sp.sqrt(2)),
            ('I',1/sp.sqrt(2)),('J',1/sp.sqrt(2)),('S',sp.Rational(1,2)),
            ('T',sp.Rational(1,2)),('U',sp.Rational(1,2))]
rows = []
for family_name,scale in examples:
    rows.append(f"| {family_name} | {pi_decision(family_name,{h:1})['status']} | {scale} | {pi_decision(family_name,{h:scale})['status']} |")
display(Markdown('| Family | PI at h=1 | Another h | PI there |\\n|---|---|---|---|\\n'+'\\n'.join(rows)))''')
    md(r'''## 8. Finite-order determinant does not suffice for PI

These exact specializations from the manuscript have $\det\sigma=I$ but
fail its PI criterion. They illustrate why a matrix order check is not itself
a replacement for the general PI argument.''')
    code('''for family_name,params in [('O',{f:9,h:1/(2*sp.sqrt(2))}),('N',{f:2,g:sp.sqrt(3)})]:
    specialized = FAMILIES[family_name].specialize(params)
    print(family_name,params,'det(sigma) =',specialized.detsigma())
    print(pi_decision(family_name,params))
    assert specialized.detsigma()==sp.eye(2)
    assert pi_decision(family_name,params)['status'] is False''')
    md(r'''## 9. Source correction in Z: a negative control

The expanded relations on printed p. 400 replace $x_1y_1$ by $x_1y_2$
in two positions, while the Sigma matrix has $x_1y_1$. The source-based
presentation passes the four cubic checks. The two replacements make the
following overlap residuals nonzero. The wrong presentation is only a negative
control and is never substituted into the family registry.''')
    code('''for item in report['checks']:
    if item['id'].startswith('Z.printed_relations_rejected.'):
        print(item['id'].split('.')[-1])
        print(item['detail']['residual'])
        assert item['passed']''')
    md(r'''## 10. All 21 closed degree models

Each row below is computed from the source-derived torus and its actions on
the center lattice. The comparison with the manuscript's formula is independent
of that computation. The sample order 12 is just a compact display; the API
also runs at arbitrarily specified common orders (including 3988 in the suite).

The exact-parameter interface accepts algebraic $p,q,f,g,h$, with no numerical
angles or tolerances. It records a default algebraic branch if p or q is omitted.''')
    code('''rows = []
for item in report['checks']:
    if item['category']=='all-family degree certificates' and item['id'].endswith('.12'):
        c = item['detail']
        rows.append(f"| {c['family']} | {c['pi_degree']} | {c['center_action_order']} | {c['algebra_pi_degree']} | {c['closed_formula']} |")
display(Markdown('| Family | Torus degree | Center action order | Algebra degree | Formula |\\n|---|---:|---:|---:|---:|\\n'+'\\n'.join(rows)))
exact = degree_from_parameters('Z',{f:3,h:-sp.Rational(1,2)})
print('Z(3,-1/2): PI degree =',exact['algebra_pi_degree'])
display(Math('A = '+sp.latex(sp.Matrix(exact['commutation_matrix']))))
display(Math(r'\\Gamma = '+sp.latex(sp.Matrix(exact['center_lattice']))))
print('Center-action image order:',exact['center_action_order'])''')
    md(r'''## 11. C: a rank-nine module and an open exact degree

The code gives an unbounded reduction of $y_1^iy_2^j$ in a nine-element left
$\Bbbk[a,b]$ basis, where $a=y_1^2y_2$, $b=y_1^3-y_2^3$. Its recurrence
identities and Hilbert-series equality are checked separately.
The subtorus degree and the exhaustive integer constraints yield the three
possibilities in the manuscript; the code does not silently choose one.''')
    code('''print('C module basis exponents:',C_BASIS)
print('Reduction of y1^8 y2^5:',c_module_reduce(8,5))
c_certificate = degree_from_parameters('C',{h:1/sp.sqrt(3)})
print('Possible PI degrees:',c_certificate['possible_pi_degrees'])
assert c_certificate['exact_degree_is_open']''')
    md(r'''## 12. Primitivity: algebraic witnesses and the remaining proof dependencies

Normality of the localization elements is checked in the original PBW algebra.
Central orbit polynomials are checked by exact Laurent reduction. In the next
example u has order 3 while the other parameter remains symbolic, so this is
not merely a PI example with both parameters specialized to roots of unity.

For the primitive branch, the code certifies the nonidentity exponent actions
and the odd/even obstruction to an inner monomial action. The simplicity,
local-closedness and Nullstellensatz arguments remain mathematical proofs.''')
    code('''witness = conditional_orbit_certificate('O','u',3)
print('Orbit coefficients central:',witness['central'])
print('Nonconstant central coefficient exists:',witness['has_nonconstant_coefficient'])
print(witness['remaining_parameter'])
for item in report['checks']:
    if item['id']=='O.outer_actions':
        print('Outer-action certificates:',item['detail'])''')
    md(r'''## 13. Coverage ledger and the manuscript's L sign correction

Every required computational block must have a labeled check, or the verifier
fails. The ledger also records the theoretical dependencies and open questions.

**Correction:** in the L torus ordered as $(w,y_1,y_1^{-1}y_2,u)$, the first
upper exponent is $-(2\eta+\varphi)$, since
$w y_1=(h^2f)^{-1}y_1w$. The supplied manuscript prints the opposite sign.
The Pfaffian and alternating Smith invariants are unchanged, so the closed
degree formula is unaffected. The program retains an explicit negative control
for the old sign. See the accompanying manuscript-correction note.''')
    code('''print('Manuscript labels:',report['coverage']['manuscript_label_count'])
print('Required computation blocks:',report['coverage']['required_computational_blocks'])
print('Missing blocks:',report['coverage']['missing_required_blocks'])
assert not report['coverage']['missing_required_blocks']
inspect_label = 'lem:Fnorm'  # change to any computational label
for item in report['checks']:
    if inspect_label in item['manuscript_labels']:
        print(item['id'],item['passed'],item['detail'])''')
    md(r'''## 14. Export a verification record

The JSON file contains each check identifier, its status, and the zero residual
or exact certificate. The original source files are not needed to run this
notebook and are not redistributed with it.

For future extensions, add further identities or degree models to the readable
Python source, then regenerate this notebook. Mathematical/theoretical
dependencies and open questions are retained in the coverage record.''')
    code('''from pathlib import Path
write_report(report,Path('notebook_reports'))
Path('notebook_verification.json').write_text(Path('notebook_reports/verification.json').read_text())
print('Saved notebook_verification.json')
print('All',report['summary']['total'],'checks passed.')''')

    notebook = {'nbformat':4,'nbformat_minor':5,
                'metadata':{'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},
                            'language_info':{'name':'python','version':'3.12'},
                            'colab':{'provenance':[]}},'cells':cells}
    dest = ROOT/'notebooks'/'Double_Ore_14641_PI_verification.ipynb'
    dest.write_text(json.dumps(notebook,indent=1,ensure_ascii=False)+'\n')

    from ore14641.families import FAMILIES
    import sympy as sp
    lines = ['# Original matrices for all 26 families','',
             'Source: Zhang--Zhang (2009), Section 4. See SOURCE_NOTES.md for source corrections.',
             'All twists are retained. N has no independent h. The registry is the single defining data source.','']
    for name,family in FAMILIES.items():
        lines += [f'## Family {name}', '', f'Subcase {family.section}, printed p. {family.page}.', '',
                  '$$P='+sp.latex(sp.Tuple(*family.P))+r',\qquad Q='+sp.latex(sp.Tuple(*family.Q))+'$$','',
                  '$$\\Sigma='+sp.latex(family.sigma)+'$$','',
                  '$$M='+sp.latex(family.printed_m)+'$$','',
                  r'$$\det\sigma|_{A_1}='+sp.latex(family.detsigma())+'$$','',
                  'Nonzero restrictions: '+', '.join(str(v) for v in family.nonzero)+'.',
                  'Parameter equations: '+(', '.join(str(v)+'=0' for v,_ in family.relations) or 'none')+'.','',
                  'PI condition (manuscript): '+family.pi_condition+'.','']
    (ROOT/'docs'/'MATRICES.md').write_text('\n'.join(lines)+'\n')
    from ore14641.degrees import symbolic_model
    from ore14641.torus import alternating
    from ore14641.words import H,p
    lines=['# Torus commutation and action matrices','',
           'All columns are generator images. If an action has matrix M and scalars c_i,',
           'the i-th generator is sent to c_i times the ordered monomial with exponent column M_i.',
           'These matrices are derived by exact Laurent word reduction, not copied from the closed degree formula.','']
    exported={}
    for name in ('B','E','J','K','L','M','N','O','P','S','T','U','W','Z','F','I','Q','R','V'):
        model,T,units,labels=symbolic_model(name)
        mapping=getattr(model,'substitutions',getattr(model,'parameter_map',{}))
        lines += [f'## Family {name}','',
                  'Formal scalar chart: '+', '.join(f'${sp.latex(k)}={sp.latex(v)}$' for k,v in mapping.items())+'.','',
                  'Torus generator order: '+', '.join('`'+v+'`' for v in T.names)+'.','',
                  '$$C='+sp.latex(T.commutations())+'$$','']
        actions=[]
        for unit,label in zip(units,labels):
            action=T.action(unit)
            lines += [f'Coset `{label}`:','', '$$M='+sp.latex(action['matrix'])+'$$','',
                      '$$(c_i)='+sp.latex(sp.Tuple(*action['scalars']))+'$$','']
            actions.append({'coset':label,'matrix':action['matrix'].tolist(),
                            'scalars':list(map(str,action['scalars']))})
        exported[name]={'chart':{str(k):str(v) for k,v in mapping.items()},'generators':list(T.names),
                        'commutation_matrix':[[str(v) for v in T.commutations().row(i)] for i in range(4)],
                        'actions':actions}
    for name in ('D','G'):
        matrix=sp.ones(4)
        values=((-1,-1/(H*p),1/(H*p),1/p,p,1/p) if name=='D'
                else (1,1/(H*p),1/(H*p),1/p,p,1/p))
        for (i,j),v in zip([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],values):
            matrix[i,j]=v;matrix[j,i]=1/v
        lines += [f'## Family {name}','',
                  'Here H=h. Generator order: x1, theta, y1, y2. The birational relations',
                  'are certified through Omega and Xi; the model needs no crossed-product factor.','',
                  '$$C='+sp.latex(matrix)+'$$','']
        exported[name]={'chart':{'H':'h'},'generators':['x1','theta','y1','y2'],
                        'commutation_matrix':[[str(v) for v in matrix.row(i)] for i in range(4)],'actions':[]}
    (ROOT/'docs'/'TORUS_ACTION_MATRICES.md').write_text('\n'.join(lines)+'\n')
    (ROOT/'reports'/'torus_models.json').write_text(json.dumps(exported,indent=2,default=str)+'\n')
    print('Built',dest.name,'with',len(cells),'cells')


if __name__ == '__main__':
    build()
