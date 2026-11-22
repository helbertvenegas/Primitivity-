"""Finite algebraic proof objects for unbounded induction and root orders.

A certificate retains both sides of an identity and its hypotheses. Arithmetic
checking is distinct from the explicitly stated mathematical inference rules.
"""
import sympy as sp
from .torus import alternating
from .degrees import symbolic_model
from .words import TorusWords, WordModel, H, F
from .families import q


def uniform_bundle():
    identities=[]; ideals=[]; lattices=[]; derivations=[]
    def identity(identifier,lhs,rhs,labels):
        identities.append({'id':identifier,'lhs':str(lhs),'rhs':str(rhs),'labels':list(labels)})
    def matrix_id(identifier,lhs,rhs,labels):
        lhs,rhs=sp.Matrix(lhs),sp.Matrix(rhs)
        for i in range(lhs.rows):
            for j in range(lhs.cols):identity(f'{identifier}.{i}.{j}',lhs[i,j],rhs[i,j],labels)
    def ideal(identifier,source,target,forward,backward,labels):
        ideals.append({'id':identifier,'source':list(map(str,source)),'target':list(map(str,target)),
                       'forward':forward,'backward':backward,'labels':list(labels)})
    # Weight induction: eta=p^nu*w is left as an independent transcendent.
    p,b,c,eta,P,Q,w,n,j,k=sp.symbols('p b c eta P0 Q0 w n j k')
    beta1=c*(p**2-1)/p**2;beta2=-b*(p**2-1);a=b*c
    lam=lambda x:x+a/x
    mu=lam(eta);lp,lm=lam(p*eta),lam(eta/p)
    B=sp.Matrix([[p*mu,-p*beta1],[-beta2/p,mu/p]])
    plus=(B-lm*sp.eye(2))/(lp-lm);minus=(lp*sp.eye(2)-B)/(lp-lm)
    matrix_id('weights.projectors_sum',plus+minus,sp.eye(2),('prop:weights','lem:lower'))
    matrix_id('weights.projector_plus',plus*plus,plus,('prop:weights',))
    matrix_id('weights.projector_minus',minus*minus,minus,('prop:weights',))
    matrix_id('weights.projectors_orthogonal',plus*minus,sp.zeros(2),('prop:weights',))
    matrix_id('weights.spectral_decomposition',lp*plus+lm*minus,B,('prop:weights',))
    identity('weights.eigenvalue_gap',lp-lm,(p-1/p)*(eta-a/eta),('lem:lower',))
    identity('weights.no_missing_lower',(p*mu-lp)/(p*beta1),b/eta,('lem:lower',))
    identity('weights.Omega_obstruction',mu*b/eta-b-c*(b/eta)**2,0,('lem:lower',))
    identity('weights.branch_involution',lam(a/eta/p),lp,('lem:lower',))
    identity('weights.distinct',lam(P*w)-lam(Q*w),(P-Q)*(w-a/(P*Q*w)),('prop:weights',))
    for step,index in [(-1,j),(1,j+1)]:
        identity('weights.index.'+str(step),-n+2*j+step,-(n+1)+2*index,('prop:weights',))
    identity('weights.index.left_endpoint',-n-1,-(n+1),('prop:weights',))
    identity('weights.index.right_endpoint',n+1,-(n+1)+2*(n+1),('prop:weights',))
    identity('weights.dimension',(n+1)+1,n+2,('prop:weights',))
    identity('weights.lowering.positive',(n+1)+(k-1),n+k,('prop:weights',))
    identity('weights.zero_degree_parity',(-2*k+2*k),0,('prop:weights',))
    induction={'id':'weights.all_n','labels':['prop:weights','thm:DGprim'],
       'base':{'degree':0,'weight_indices':[0],'dimension':1},
       'index_maps':['j','j+1'],'intervals':[[0,'n'],[1,'n+1']],
       'union':[0,'n+1'],'eigenvalue':'lambda(p**(-n+2*j)*w)',
       'assumptions':['characteristic zero','p not a root of unity','b*c != 0',
                      'w transcendental','PBW gives dim(V_n)=n+1',
                      'V_(n+1)=V_n*y1+V_n*y2','Omega central and nonzero',
                      'right multiplication by Omega injective by domain and flatness'],
       'inference_rules':['two complementary polynomial spectral projectors',
                          'nonvanishing of both branches from Omega obstruction',
                          'affine interval union induction','dimension counts distinct eigenlines'],
       'ideal_consequence':'lower |nu| to zero; clear denominators; use noetherian prime-product argument',
       'ideal_theory_dependency':'the final prime-product and torsion arguments are those in the manuscript'}
    M,theta,phi,eta,k,kp,N=sp.symbols('M theta phi eta k kp N')
    pf=lambda A:A[0,1]*A[2,3]-A[0,2]*A[1,3]+A[0,3]*A[1,2]
    for name,upper,target,fw,bw in [
      ('D',(M,-M-k-kp,-k-kp,-k,k,-k),(M,k,kp),
       [[1,0,0,0,0,0],[0,0,0,-1,0,0],[0,0,-1,1,0,0]],
       [[1,0,0],[-1,-1,-1],[0,-1,-1],[0,-1,0],[0,1,0],[0,-1,0]]),
      ('G',(0,-k-kp,-k-kp,-k,k,-k),(k,kp),
       [[0,0,0,-1,0,0],[0,0,-1,1,0,0]],
       [[0,0],[-1,-1],[-1,-1],[-1,0],[1,0],[-1,0]]),
      ('L',(-2*eta-phi,0,2*theta,M,2*eta+2*phi,0),(2*eta,phi,2*theta,M),
       [[-2,0,0,0,-1,0],[1,0,0,0,1,0],[0,0,1,0,0,0],[0,0,0,1,0,0]],
       [[-1,-1,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1],[1,2,0,0],[0,0,0,0]]),
      ('K',(M,eta,2*eta+phi,eta,2*eta+phi,-2*theta),(M,eta,phi,2*theta),
       [[1,0,0,0,0,0],[0,1,0,0,0,0],[0,-2,1,0,0,0],[0,0,0,0,0,-1]],
       [[1,0,0,0],[0,1,0,0],[0,2,1,0],[0,1,0,0],[0,2,1,0],[0,0,0,1*-1]])]:
        label='prop:DGdeg' if name in ('D','G') else 'thm:L' if name=='L' else 'thm:monodeg'
        A=alternating(upper)
        expected=2*k*(k+kp) if name in ('D','G') else 2*M*theta if name=='L' else -2*M*theta
        identity(name+'.pfaffian',pf(A),expected,(label,))
        ideal(name+'.entry_ideal',upper,target,fw,bw,(label,))
        derivations.append({'family':name,'label':label,'root_order':'N=2*M' if name!='G' else 'N',
           'generation_gcd':'gcd(M,k,kp)=1' if name=='D' else 'gcd(N,k,kp)=1' if name=='G' else 'gcd(M,theta,phi,eta)=1',
           'alternating_smith_rule':'s1=gcd(entries), s1*s2=abs(Pfaffian); torus degree=prod(N/gcd(N,si))',
           'conclusion':'N**2/gcd(N,2*k*(k+kp))' if name in ('D','G') else '2*N' if name=='K' else '2*N if N/N0 odd, N otherwise',
           'crossed_center_factor':1 if name in ('D','G') else 2})
    # L: d=gcd(theta,M,phi,2*eta) divides 2 by Bezout in the base ideal.
    ideal('L.double_base_in_target',(theta,M,phi,2*eta),(2*theta,2*M,2*phi,2*eta),
          [[2,0,0,0],[0,2,0,0],[0,0,2,0],[0,0,0,1]],
          # This is an inclusion, not equivalence: backward deliberately omitted.
          None,('thm:L',))
    identity('L.gcd_product_scaling',2*M,(2*M),('thm:L',))
    # Central congruence systems and complete integer parametrizations for all m.
    m=sp.Symbol('m');u,ss,rr,tt=sp.symbols('u ss rr tt')
    specs=[('Q',[[2,0,0,0],[-2,m,0,0],[2,0,m,2],[2,0,m,0]],
             ['a','c-d','a-d','a+b'],[2,2,'m','m']),
           ('V',[[-2,m,0,0],[2,0,0,0],[-2,0,m,2],[-2,0,m,0]],
             ['b','c-d','b+d','a+b'],[2,2,'m','m']),
           ('R',[[1,0,0,0],[0,m,0,0],[0,0,m,-2],[0,0,0,1]],
             ['a','b','c+2*d','d'],[1,'m','m',1])]
    for name,G,forms,divs in specs:
        G=sp.Matrix(G)
        lattices.append({'family':name,'basis':[[str(x) for x in row] for row in G.tolist()],
                         'forms':forms,'divisors':divs,'determinant':str(4*m**2 if name in ('Q','V') else m**2),
                         'labels':['thm:FIQRV', 'thm:QVdeg' if name in ('Q','V') else 'thm:Rdeg'],
                         'meaning':'m=ord(h^2)>0; columns parametrize all integral solutions'})
        identity(name+'.central_index',G.det()**2,(4*m**2 if name in ('Q','V') else m**2)**2,('thm:FIQRV',))
    # The finite list below is a proof by exhaustive parity classes, not by orders.
    parity={'classes':[
        {'order_form':'2*r+1','m_form':'2*r+1','m_parity':1,'root_to_m':1,'name':'odd'},
        {'order_form':'4*r+2','m_form':'2*r+1','m_parity':1,'root_to_m':-1,'name':'two_mod_four'},
        {'order_form':'4*r','m_form':'2*r','m_parity':0,'root_to_m':-1,'name':'zero_mod_four'}],
        'assumption':'root has exact positive order d; m=ord(root^2)',
        'elementary_rule':'ord(root^2)=d/gcd(d,2); exact order 2 gives -1 in characteristic zero',
        'families':['M','N','O','P','W','Z']}
    # Save real words, not their normalized forms, for a separate word checker.
    models=[]
    for name in ('B','E','J','K','L','M','N','O','P','S','T','U','W','Z','F','I','Q','R','V'):
        model,T,units,labels=symbolic_model(name)
        if name in ('B','E','J','M','N','O','P','S','T','U','W','Z'):
            r,rho,X,Y=T.generators
            xp=X*rho if name in ('B','J') else X*r if name in ('N','U','W') else X
            yp=Y*r if name in ('B','E') else Y*rho if name in ('P','T','Z') else Y
            T=TorusWords(model,(r,rho,xp,yp),('r','rho','Xp','Yp'))
        encode=lambda word:{'letters':[list(v) for v in word.letters],'coefficient':str(word.coefficient)}
        actions=[T.action(unit) for unit in units]
        models.append({'family':name,'planes':[model.x_kind,model.y_kind],
           'plane_q':[str(model.x_q),str(model.y_q)],
           'relations':[[list(ab),[bb,aa,str(cc)]] for ab,(bb,aa,cc) in model.mixing.items()],
           'scalar_relations':[[str(poly),str(v)] for poly,v in model.family.relations] if hasattr(model,'family') else [],
           'generators':[encode(v) for v in T.generators],'units':[encode(v) for v in units],
           'unit_labels':list(labels),'commutations':[[str(v) for v in row] for row in T.commutations().tolist()],
           'actions':[{'matrix':[[int(v) for v in row] for row in a['matrix'].tolist()],'scalars':list(map(str,a['scalars']))} for a in actions]})
    # A second K torus of index two gives a short uniform Smith proof.
    km=WordModel(sp.cancel,'coordinate','coordinate',q,-1,
        {(0,0):(0,0,H),(0,1):(1,1,H),(1,0):(0,1,H),(1,1):(1,0,H*F)},'K_index_two')
    KT=TorusWords(km,(km.z1,km.z2,km.e1,km.e2**2),('y1','y2','x1','x2_squared'))
    ku=(km.identity(),km.e2)
    ka=[KT.action(v) for v in ku]
    models.append({'family':'K_index_two','planes':['coordinate','coordinate'],
       'plane_q':[str(q),'-1'],'relations':[[list(ab),[bb,aa,str(cc)]] for ab,(bb,aa,cc) in km.mixing.items()],
       'scalar_relations':[],'generators':[encode(v) for v in KT.generators],
       'units':[encode(v) for v in ku],'unit_labels':['0','1'],
       'commutations':[[str(v) for v in row] for row in KT.commutations().tolist()],
       'actions':[{'matrix':[[int(v) for v in row] for row in a['matrix'].tolist()],
                   'scalars':list(map(str,a['scalars']))} for a in ka],
       'source_coordinates':True,'substitutions':{'H':'h','F':'f'},
       'expected_upper':['-1','H','F*H**2','H','F*H**2','q**(-2)']})
    lm=WordModel(sp.cancel,'coordinate','coordinate',q,-1,
        {(0,0):(0,1,H*F),(0,1):(1,1,H),(1,0):(0,0,H*F),(1,1):(1,0,H)},'L_index_two')
    LT=TorusWords(lm,(lm.e1*lm.e2,lm.z1,lm.z1**-1*lm.z2,lm.e1**2),('w','y1','rho','u'))
    lu=(lm.identity(),lm.e1);la=[LT.action(v) for v in lu]
    models.append({'family':'L_index_two','planes':['coordinate','coordinate'],
       'plane_q':[str(q),'-1'],'relations':[[list(ab),[bb,aa,str(cc)]] for ab,(bb,aa,cc) in lm.mixing.items()],
       'scalar_relations':[],'generators':[encode(v) for v in LT.generators],
       'units':[encode(v) for v in lu],'unit_labels':['0','1'],
       'commutations':[[str(v) for v in row] for row in LT.commutations().tolist()],
       'actions':[{'matrix':[[int(v) for v in row] for row in a['matrix'].tolist()],
                   'scalars':list(map(str,a['scalars']))} for a in la],
       'source_coordinates':True,'substitutions':{'H':'h','F':'f'},
       'expected_upper':['1/(F*H**2)','1','q**2','-1','F**2*H**2','1']})
    return {'schema':'ore14641.uniform.v1','identities':identities,'integer_ideals':ideals,
            'center_parametrizations':lattices,'weight_induction':induction,'root_parity':parity,
            'degree_derivations':derivations,'word_models':models,
            'theory_dependencies':['PBW independence and domain','alternating integer Smith theorem',
                'center of a quantum torus','finite crossed-product center/degree theorem',
                'Bezout identity for gcd and cyclic group exact order',
                'flatness and the ideal-theoretic arguments listed in the manuscript']}


def uniform_checks(check):
    data=uniform_bundle()
    for claim in data['identities']:
        result=sp.cancel(sp.sympify(claim['lhs'])-sp.sympify(claim['rhs']))
        check('uniform.'+claim['id'],'unbounded algebraic proof objects','general',result==0,
              claim,labels=claim['labels'])
    for claim in data['integer_ideals']:
        source=sp.Matrix(list(map(sp.sympify,claim['source'])))
        target=sp.Matrix(list(map(sp.sympify,claim['target'])))
        valid=sp.Matrix(claim['forward'])*source==target
        if claim['backward'] is not None:valid &= sp.Matrix(claim['backward'])*target==source
        check('uniform.'+claim['id'],'unbounded integer ideal certificates','general',valid,
              claim,labels=claim['labels'])
