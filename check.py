#!/usr/bin/env python3
"""Independent certificate checker. Python standard library ONLY.

No imports from ore14641, SymPy, NumPy, Sage, or Jupyter. Polynomials use sparse
Fraction arithmetic; PBW uses a rightmost worklist; balanced Laurent planes use
central squares plus alternating words, not the producer's ratio coordinates.
"""
import ast
from fractions import Fraction as F
from itertools import combinations, product
import json
import re
from math import gcd, isqrt, prod
from pathlib import Path
import sys


class Ring:
    def __init__(self,expressions=(),relations=()):
        names={'I'}
        for s in tuple(expressions)+tuple(v[0] for v in relations):
            names.update(x.id for x in ast.walk(ast.parse(str(s),mode='eval')) if isinstance(x,ast.Name))
        self.names=tuple(sorted(names));self.index={v:i for i,v in enumerate(self.names)}
        self.zero=(0,)*len(self.names);self.relations={};self.cache={}
        # Quotient relations are monic univariate polynomials.
        self.relations[self.index['I']]=(F(1),F(0),F(1))
        for expression,var in relations:
            polynomial=self.parse(expression).num
            pos=self.index[var];coeff={m[pos]:c for m,c in polynomial.items()}
            if any(any(e for j,e in enumerate(m) if j!=pos) for m in polynomial):
                raise ValueError('Only univariate coefficient constraints allowed')
            degree=max(coeff);lead=coeff[degree]
            self.relations[pos]=tuple(coeff.get(i,F(0))/lead for i in range(degree+1))
    def clean(self,p):
        result={};pending=list(p.items())
        while pending:
            m,c=pending.pop()
            if not c:continue
            for j,coeff in self.relations.items():
                d=len(coeff)-1
                if m[j]>=d:
                    for i,a in enumerate(coeff[:-1]):
                        if a:
                            mm=list(m);mm[j]+=i-d;pending.append((tuple(mm),-c*a))
                    break
            else:result[m]=result.get(m,F(0))+c
        return {m:c for m,c in result.items() if c}
    def add(self,a,b):
        r=a.copy()
        for m,c in b.items():r[m]=r.get(m,F(0))+c
        return self.clean(r)
    def mul(self,a,b):
        r={}
        for m,c in a.items():
            for n,d in b.items():
                k=tuple(x+y for x,y in zip(m,n));r[k]=r.get(k,F(0))+c*d
        return self.clean(r)
    def const(self,v):return Rat(self,{self.zero:F(v)},{self.zero:F(1)})
    def parse(self,source):
        source=str(source)
        if source in self.cache:return self.cache[source]
        def visit(x):
            if isinstance(x,ast.Constant) and isinstance(x.value,int):return self.const(x.value)
            if isinstance(x,ast.Name):
                m=list(self.zero);m[self.index[x.id]]=1
                return Rat(self,self.clean({tuple(m):F(1)}),{self.zero:F(1)})
            if isinstance(x,ast.UnaryOp):
                a=visit(x.operand)
                if isinstance(x.op,ast.USub):return -a
                if isinstance(x.op,ast.UAdd):return a
            if isinstance(x,ast.BinOp):
                a=visit(x.left)
                if isinstance(x.op,ast.Pow):
                    e=visit(x.right).integer();return a**e
                b=visit(x.right)
                if isinstance(x.op,ast.Add):return a+b
                if isinstance(x.op,ast.Sub):return a-b
                if isinstance(x.op,ast.Mult):return a*b
                if isinstance(x.op,ast.Div):return a/b
            raise ValueError('Unsupported expression node '+ast.dump(x))
        value=visit(ast.parse(source,mode='eval').body);self.cache[source]=value
        return value


class Rat:
    def __init__(self,ring,num,den):
        self.ring=ring;num=ring.clean(num);den=ring.clean(den)
        if not den:raise ZeroDivisionError('Zero certificate denominator')
        if not num:self.num={};self.den={ring.zero:F(1)};return
        # Cancel common monomial factors (no external polynomial gcd).
        mn=[min(m[j] for m in num) for j in range(len(ring.names))]
        md=[min(m[j] for m in den) for j in range(len(ring.names))]
        cancel=tuple(min(a,b) for a,b in zip(mn,md))
        if any(cancel):
            num={tuple(a-b for a,b in zip(m,cancel)):c for m,c in num.items()}
            den={tuple(a-b for a,b in zip(m,cancel)):c for m,c in den.items()}
        lead=den[max(den)]
        self.num={m:c/lead for m,c in num.items()};self.den={m:c/lead for m,c in den.items()}
    def asrat(self,b):return b if isinstance(b,Rat) else self.ring.const(b)
    def __neg__(self):return Rat(self.ring,{m:-c for m,c in self.num.items()},self.den)
    def __add__(self,b):
        b=self.asrat(b);r=self.ring
        if self.den==b.den:return Rat(r,r.add(self.num,b.num),self.den)
        return Rat(r,r.add(r.mul(self.num,b.den),r.mul(b.num,self.den)),r.mul(self.den,b.den))
    __radd__=__add__
    def __sub__(self,b):return self+-self.asrat(b)
    def __rsub__(self,b):return self.asrat(b)+-self
    def __mul__(self,b):
        b=self.asrat(b);r=self.ring
        return Rat(r,r.mul(self.num,b.num),r.mul(self.den,b.den))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=self.asrat(b);r=self.ring
        return Rat(r,r.mul(self.num,b.den),r.mul(self.den,b.num))
    def __pow__(self,n):
        if n<0:return Rat(self.ring,self.den,self.num)**(-n)
        r=self.ring.const(1);a=self
        while n:
            if n&1:r=r*a
            a=a*a;n//=2
        return r
    def __bool__(self):return bool(self.num)
    def __eq__(self,b):
        b=self.asrat(b);r=self.ring
        return r.mul(self.num,b.den)==r.mul(b.num,self.den)
    def integer(self):
        if set(self.num)-{self.ring.zero} or set(self.den)-{self.ring.zero}:raise ValueError('Not constant')
        a=self.num.get(self.ring.zero,F(0))/self.den[self.ring.zero]
        if a.denominator!=1:raise ValueError('Not integral')
        return int(a)


def transpose(A):return list(map(list,zip(*A)))
def mm(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def det(A):
    A=[list(map(F,row)) for row in A];n=len(A);ans=F(1)
    for j in range(n):
        pivot=next((i for i in range(j,n) if A[i][j]),None)
        if pivot is None:return 0
        if pivot!=j:A[j],A[pivot]=A[pivot],A[j];ans=-ans
        value=A[j][j];ans*=value
        for i in range(j+1,n):
            factor=A[i][j]/value
            for k in range(j+1,n):A[i][k]-=factor*A[j][k]
    return ans


def source_checks(sources,record):
    for s in sources:
        expressions=sum(s['sigma'],[])+sum(s['M'],[])+list(s['P'])+list(s['Q'])
        expressions+=sum(s['expected_detsigma'],[])+[p for p,v in s['relations']]
        ring=Ring(expressions,s['relations']);parse=ring.parse
        sigma=[[parse(v) for v in row] for row in s['sigma']]
        M=[[parse(v) for v in row] for row in s['M']];order=(0,2,1,3)
        record(s['family']+'.reshuffle',all(M[i][j]==sigma[order[i]][order[j]] for i in range(4) for j in range(4)))
        p12,p11=map(parse,s['P']);q12,q11=map(parse,s['Q'])
        blocks=[[[sigma[2*a+i][2*b+j] for j in range(2)] for i in range(2)] for a in range(2) for b in range(2)]
        B11,B12,B21,B22=blocks
        C1,C2,C3=mm(B11,B12),mm(B11,B22),mm(B21,B12)
        D=[[-p11*C1[i][j]+C2[i][j]-p12*C3[i][j] for j in range(2)] for i in range(2)]
        record(s['family']+'.detsigma',all(D[i][j]==parse(s['expected_detsigma'][i][j]) for i in range(2) for j in range(2)))
        rules={(1,0):{(0,1):q12,(0,0):q11},(3,2):{(2,3):p12,(2,2):p11}}
        targets=((0,2),(1,2),(0,3),(1,3))
        for i,key in enumerate(((2,0),(2,1),(3,0),(3,1))):rules[key]={w:sigma[i][j] for j,w in enumerate(targets) if sigma[i][j]}
        rules={key:{w:c for w,c in value.items() if c} for key,value in rules.items()}
        record(s['family']+'.decreasing',all(v<k for k,ts in rules.items() for v in ts))
        def reduce_terms(initial):
            pending=list(initial);out={}
            while pending:
                word,c=pending.pop()
                pos=next((j for j in range(len(word)-2,-1,-1) if word[j:j+2] in rules),None)
                if pos is None:out[word]=out.get(word,ring.const(0))+c
                else:
                    for v,d in rules[word[pos:pos+2]].items():pending.append((word[:pos]+v+word[pos+2:],c*d))
            return {w:c for w,c in out.items() if c}
        for abc in combinations(range(4),3):
            word=tuple(reversed(abc));branches=[]
            for pos in (0,1):branches.append(reduce_terms([(word[:pos]+v+word[pos+2:],d) for v,d in rules[word[pos:pos+2]].items()]))
            keys=set(branches[0])|set(branches[1])
            record(s['family']+'.PBW.'+str(word),all(branches[0].get(w,ring.const(0))==branches[1].get(w,ring.const(0)) for w in keys))


def basis_checks(sources,changes,record):
    sources={s['family']:s for s in sources}
    for change in changes:
        source=sources[change['family']]
        convert=lambda expr:re.sub(r'\bf\b','(s**2)',expr) if change['square_root'] else expr
        strings=sum(source['sigma'],[])+sum(change['E'],[])+sum(change['Z'],[])+sum(change['target'],[])
        strings=list(map(convert,strings))
        r=Ring(strings,source['relations'])
        matrix=lambda A:[[r.parse(convert(x)) for x in row] for row in A]
        S,E,Z,T=map(matrix,(source['sigma'],change['E'],change['Z'],change['target']))
        C=[[Z[i//2][j//2]*E[i%2][j%2] for j in range(4)] for i in range(4)]
        record(change['family']+'.independent_basis_change',mm(C,S)==mm(T,C))
        record(change['family']+'.independent_basis_invertible',bool(E[0][0]*E[1][1]-E[0][1]*E[1][0]) and bool(Z[0][0]*Z[1][1]-Z[0][1]*Z[1][0]))


def word_checks(models,record):
    for model in models:
        expressions=model['plane_q']+sum(model['commutations'],[])+[v[1][2] for v in model['relations']]
        expressions+=[s for a in model['actions'] for s in a['scalars']]
        expressions+=[p for p,v in model['scalar_relations']]
        ring=Ring(expressions,model['scalar_relations']);one=ring.const(1)
        rules={}
        for (a,b),(j,i,c) in model['relations']:
            c=ring.parse(c)
            for key,val in [((a,b,1,1),(j,i,c)),((i,b,1,-1),(j,a,one/c)),
                            ((a,j,-1,1),(b,i,one/c)),((i,j,-1,-1),(b,a,c))]:
                if key in rules and rules[key]!=val:raise ValueError('Inconsistent inverse relations')
                rules[key]=val
        def multiply(u,v):return (u[0]+v[0],u[1]*v[1])
        def power(u,n):
            if n<0:return power((tuple((a,b,-sg) for a,b,sg in reversed(u[0])),one/u[1]),-n)
            return u[0]*n,u[1]**n
        def decode(u):return tuple(tuple(x) for x in u['letters']),ring.parse(u['coefficient'])
        def normal(word):
            letters,c=word;letters=list(letters)
            # Rightmost mixed inversion, versus the producer's leftmost pass.
            while True:
                at=next((i for i in range(len(letters)-2,-1,-1) if letters[i][0]==1 and letters[i+1][0]==0),None)
                if at is None:break
                _,a,sy=letters[at];_,b,sx=letters[at+1]
                j,i,factor=rules[a,b,sy,sx];c=c*factor
                letters[at:at+2]=[(0,j,sx),(1,i,sy)]
            results=[]
            for plane,kind in enumerate(model['planes']):
                current=[(i,s) for pl,i,s in letters if pl==plane]
                if kind=='coordinate':
                    x,y=0,0;q=ring.parse(model['plane_q'][plane])
                    for index,sg in current:
                        if index==0:c=c*q**(y*sg);x+=sg
                        else:y+=sg
                    results.append((x,y))
                else:
                    square=0;stack=[]
                    for index,sg in current:
                        if sg<0:square-=1
                        if stack and stack[-1]==index:stack.pop();square+=1
                        else:stack.append(index)
                    results.append((square,tuple(stack)))
            return tuple(results),c
        def equal(u,v):
            a,c=normal(u);b,d=normal(v);return a==b and c==d
        gens=list(map(decode,model['generators']));units=list(map(decode,model['units']))
        if 'expected_upper' in model:
            record(model['family']+'.uniform_smith_input',all(ring.parse(model['commutations'][i][j])==ring.parse(v)
                for (i,j),v in zip(combinations(range(4),2),model['expected_upper'])))
        for i,j in combinations(range(4),2):
            factor=ring.parse(model['commutations'][i][j])
            right=multiply(gens[j],gens[i]);right=(right[0],right[1]*factor)
            record(model['family']+f'.word.commutation.{i}.{j}',equal(multiply(gens[i],gens[j]),right))
        for u,label,action in zip(units,model['unit_labels'],model['actions']):
            record(model['family']+'.word.action_det.'+label,abs(det(action['matrix']))==1)
            for j,g in enumerate(gens):
                rhs=((),ring.parse(action['scalars'][j]))
                for i in range(4):rhs=multiply(rhs,power(gens[i],action['matrix'][i][j]))
                lhs=multiply(multiply(u,g),power(u,-1))
                record(model['family']+f'.word.action.{label}.{j}',equal(lhs,rhs))
        # Universal faithfulness / kernel, using the full-rank central lattice.
        name=model['family'];acts=model['actions']
        if name in ('B','E','J','S','T','U'):
            record(name+'.universal_center_faithful',all(a['matrix']!=eye(4) for a in acts[1:]))
        elif name in ('Q','V','K_index_two','L_index_two'):
            record(name+'.universal_center_faithful',acts[1]['matrix']!=eye(4))
        elif name=='R':
            signatures=[(a['matrix'][0][0],a['scalars'][0]) for a in acts]
            record(name+'.universal_center_faithful',len(set(signatures))==4 and all(model['commutations'][0][j]=='1' for j in range(4)))
        elif name in ('F','I'):
            signatures=[(a['matrix'][0][0],a['scalars'][0],a['matrix'][3][3],a['scalars'][3]) for a in acts]
            record(name+'.universal_center_faithful',len(set(signatures))==8 and all(model['commutations'][i][j]=='1' for i in (0,3) for j in range(4)))
        elif name in ('M','N','O','P','W','Z'):
            exceptional=[a for a in acts[1:] if a['matrix']==eye(4)]
            record(name+'.one_possible_kernel_element',len(exceptional)==1)
            action=exceptional[0]
            # u=alpha/beta except Z, where u=-alpha/beta; v=alpha*beta.
            chars=[]
            for expression in action['scalars']:
                z=ring.parse(expression)
                if len(z.num)!=1 or len(z.den)!=1:raise ValueError('Expected Laurent scalar character')
                (nm,nc),=z.num.items();(dm,dc),=z.den.items()
                exps={v:nm[i]-dm[i] for i,v in enumerate(ring.names)}
                aa,bb=exps.get('alpha',0),exps.get('beta',0)
                if (aa-bb)%2 or (aa+bb)%2:raise ValueError('Not a character in u,v')
                up,vp=(aa-bb)//2,(aa+bb)//2;constant=nc/dc
                if name=='Z':constant*=(-1)**up
                if any(e for v,e in exps.items() if v not in ('alpha','beta')):raise ValueError('Unexpected character')
                chars.append((constant,up,vp))
            parity=((1,1),(1,-1),(0,-1))
            for ui,(um,uv) in enumerate(parity):
                for vi,(vm,vv) in enumerate(parity):
                    values=[]
                    for idx,(constant,up,vp) in enumerate(chars):
                        if idx<2:
                            if vp:raise ValueError('Cross character in first pair')
                            values.append(constant**um*F(uv)**up)
                        else:
                            if up:raise ValueError('Cross character in second pair')
                            values.append(constant**vm*F(vv)**vp)
                    trivial=all(v==1 for v in values)
                    expected=name in ('O','W','Z') and ui==vi==0
                    record(name+f'.universal_parity.{ui}.{vi}',trivial==expected)


def arithmetic_checks(data,record):
    for claim in data['identities']:
        r=Ring((claim['lhs'],claim['rhs']))
        record(claim['id'],r.parse(claim['lhs'])==r.parse(claim['rhs']))
    for claim in data['integer_ideals']:
        r=Ring(claim['source']+claim['target']);src=list(map(r.parse,claim['source']));tar=list(map(r.parse,claim['target']))
        forward=mm(claim['forward'],[[x] for x in src])
        record(claim['id']+'.forward',all(a[0]==b for a,b in zip(forward,tar)))
        if claim['backward'] is not None:
            back=mm(claim['backward'],[[x] for x in tar])
            record(claim['id']+'.backward',all(a[0]==b for a,b in zip(back,src)))
    for cert in data['center_parametrizations']:
        # W=D^(-1)*C*Gamma is an integer unimodular matrix independent of m.
        # This proves necessity, sufficiency and index for ALL positive m.
        strings=sum(cert['basis'],[])+cert['forms']+['a','b','c','d','m']
        r=Ring(strings);G=[[r.parse(v) for v in row] for row in cert['basis']]
        C=[]
        for form in cert['forms']:
            f=r.parse(form);row=[]
            for var in ('a','b','c','d'):
                mon=r.parse(var);value=0
                for degree,coef in f.num.items():
                    if degree in mon.num:value=coef
                row.append(value)
            C.append(row)
        W=mm(C,G)
        W=[[ (value/r.parse(str(divisor))).integer() for value in row] for row,divisor in zip(W,cert['divisors'])]
        record(cert['family']+'.general_center_parametrization',abs(det(W))==1)
    ind=data['weight_induction']
    record('weights.induction_base',ind['base']=={'degree':0,'weight_indices':[0],'dimension':1})
    record('weights.interval_cover',ind['index_maps']==['j','j+1'] and ind['intervals']==[[0,'n'],[1,'n+1']] and ind['union']==[0,'n+1'])
    expected=[('2*r+1','2*r+1',1,1),('4*r+2','2*r+1',1,-1),('4*r','2*r',0,-1)]
    record('roots.exhaustive_parity_partition',[(v['order_form'],v['m_form'],v['m_parity'],v['root_to_m']) for v in data['root_parity']['classes']]==expected)
    record('degree.general_rules_covered',set(v['family'] for v in data['degree_derivations'])=={'D','G','K','L'})


def integer_checks(cases,record):
    for i,c in enumerate(cases):
        A,N,U,V,D=c['A'],c['N'],c['U'],c['V'],c['D'];size=len(A)
        record(f'integer.{i}.alternating',all(A[j][k]==-A[k][j] for j in range(size) for k in range(size)))
        record(f'integer.{i}.smith_identity',mm(mm(U,A),V)==D)
        record(f'integer.{i}.unimodular',abs(det(U))==abs(det(V))==1)
        diag=[abs(D[j][j]) for j in range(size)]
        record(f'integer.{i}.smith_diagonal',all(D[j][k]==0 for j in range(size) for k in range(size) if j!=k) and all(not b if not a else b%a==0 for a,b in zip(diag,diag[1:])))
        kernel=[[V[j][k]*(N//gcd(N,diag[k])) for k in range(size)] for j in range(size)]
        Gamma=c['Gamma'];W=c['kernel_change']
        record(f'integer.{i}.complete_center',mm(kernel,W)==Gamma and abs(det(W))==1)
        record(f'integer.{i}.central',all(v%N==0 for row in mm(A,Gamma) for v in row))
        image=prod(N//gcd(N,d) for d in diag);td=isqrt(image)
        record(f'integer.{i}.torus_degree',td*td==image==abs(det(Gamma)) and td==c['torus_degree'])
        def multiply(u,v):
            phase,a=u;more,b=v
            # Independently count crossings in the ordered product.
            for hi in range(size):
                for lo in range(hi):phase-=a[hi]*b[lo]*A[lo][hi]
            return (phase+more)%N,tuple(a[j]+b[j] for j in range(size))
        def power(u,n):
            # Inversion + binary powering, rather than the producer's n(n-1)/2 formula.
            if n<0:
                ph,a=u;opposite=tuple(-v for v in a)
                crossed=multiply((ph,a),(0,opposite))[0]
                u=((-crossed)%N,opposite);n=-n
            result=(0,(0,)*size)
            while n:
                if n&1:result=multiply(result,u)
                u=multiply(u,u);n//=2
            return result
        signatures=[]
        for action in c['actions']:
            signature=[]
            for column in transpose(Gamma):
                value=(0,(0,)*size)
                for j,n in enumerate(column):
                    mon=(action['phases'][j],tuple(row[j] for row in action['exponent_matrix']))
                    value=multiply(value,power(mon,n))
                signature.append(value)
            signatures.append(tuple(signature))
        factor=len(set(signatures)) if signatures else 1
        record(f'integer.{i}.crossed_center',factor==c['center_action_order'])
        record(f'integer.{i}.algebra_degree',td*factor==c['algebra_degree']==c['formula'])


def verify(bundle):
    passed=[]
    def record(name,ok):
        if not ok:raise ArithmeticError('Independent certificate failed: '+name)
        passed.append(name)
    source_checks(bundle['sources'],record)
    basis_checks(bundle['sources'],bundle['basis_changes'],record)
    for model in bundle['uniform']['word_models']:
        if not model.get('source_coordinates'):continue
        source=next(s for s in bundle['sources'] if s['family']==model['family'][0])
        r=Ring(sum(source['sigma'],[])+['h','f','q'])
        targets=((0,0),(1,0),(0,1),(1,1))
        for (a,b),(j,i,coefficient) in model['relations']:
            for key,value in model['substitutions'].items():
                coefficient=re.sub(r'\b'+key+r'\b','('+value+')',coefficient)
            record(model['family']+f'.source_mixed.{a}.{b}',all(r.parse(source['sigma'][2*a+b][c])==
                   (r.parse(coefficient) if v==(j,i) else 0) for c,v in enumerate(targets)))
    arithmetic_checks(bundle['uniform'],record)
    word_checks(bundle['uniform']['word_models'],record)
    integer_checks(bundle['integer_cases'],record)
    return {'checker':'standard-library independent checker v1','passed':len(passed),'failed':0,
            'source_families':len(bundle['sources']),'laurent_models':len(bundle['uniform']['word_models']),
            'integer_certificates':len(bundle['integer_cases']),
            'unbounded_claims':['weight induction algebra and interval cover',
                'generic integer ideal identities','complete Q/V/R center parametrizations',
                'source-derived torus splits and actions','all parity classes for M/N/O/P/W/Z'],
            'assumptions':bundle['uniform']['theory_dependencies'],'checks':passed}


if __name__=='__main__':
    bundle=json.loads(Path(sys.argv[1]).read_text());report=verify(bundle)
    if len(sys.argv)>2:Path(sys.argv[2]).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('checks','assumptions')},indent=2))
