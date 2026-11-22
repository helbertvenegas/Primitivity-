"""Small exact noncommutative rewriting engine, using tuples of letters.

Letter order: x1 < x2 < y1 < y2. All six rules decrease degree-lex order.
The four cubic overlaps are checked explicitly before using confluence.
No floating point arithmetic, commuting symbol substitution, or degree cutoff.
"""
from functools import lru_cache
from itertools import combinations
import sympy as sp

LETTERS = ('x1', 'x2', 'y1', 'y2')


class Algebra:
    def __init__(self, family):
        self.family = family
        self.scalar = family.scalar
        self.letters = LETTERS
        q12, q11 = family.Q
        p12, p11 = family.P
        self.rules = {
            (1, 0): self.clean({(0, 1): q12, (0, 0): q11}),
            (3, 2): self.clean({(2, 3): p12, (2, 2): p11}),
        }
        lhs = [(2,0),(2,1),(3,0),(3,1)]
        rhs = [(0,2),(1,2),(0,3),(1,3)]
        for row, word in enumerate(lhs):
            self.rules[word] = self.clean({v: family.sigma[row,col]
                                          for col,v in enumerate(rhs)})
        self._initialize()

    def _initialize(self):
        for word, terms in self.rules.items():
            if any((len(v),v) >= (len(word),word) for v in terms):
                raise ValueError('The rewrite rules must decrease degree-lex order')
        # A per-instance cache avoids retaining every past Algebra instance.
        self.reduce_word = lru_cache(maxsize=32768)(self._reduce_word)
        self.generators = tuple(Element(self, {(i,): sp.Integer(1)}) for i in range(len(self.letters)))

    @classmethod
    def from_rules(cls,rules,scalar,letters):
        algebra=object.__new__(cls)
        algebra.scalar=scalar
        algebra.family=None
        algebra.letters=tuple(letters)
        algebra.rules={w:algebra.clean(v) for w,v in rules.items()}
        algebra._initialize()
        return algebra

    def clean(self, terms):
        return {w:c for w,value in terms.items() if (c:=self.scalar(value)) != 0}

    def _reduce_word(self, word):
        for pos in range(len(word)-1):
            pair = word[pos:pos+2]
            if pair in self.rules:
                return tuple(self.first_reduction(word, pos).items())
        return ((word, sp.Integer(1)),)

    def first_reduction(self, word, position):
        """Force one chosen initial step, then reduce deterministically."""
        result = {}
        for replacement, coefficient in self.rules[word[position:position+2]].items():
            new = word[:position]+replacement+word[position+2:]
            for reduced, factor in self.reduce_word(new):
                result[reduced] = result.get(reduced,0)+coefficient*factor
        return self.clean(result)

    def overlaps(self):
        """All critical overlaps; the two starting reductions are independent."""
        residuals = {}
        for ascending in combinations(range(len(self.letters)),3):
            word = tuple(reversed(ascending))
            if word[:2] not in self.rules or word[1:] not in self.rules:
                continue
            left = Element(self,self.first_reduction(word,0))
            right = Element(self,self.first_reduction(word,1))
            residuals[' '.join(self.letters[i] for i in word)] = left-right
        return residuals

    def element(self, value):
        if isinstance(value, Element):
            if value.algebra is not self:
                raise ValueError('Elements belong to different algebras')
            return value
        return Element(self, {(): self.scalar(sp.sympify(value))})


class Element:
    def __init__(self, algebra, terms):
        self.algebra = algebra
        self.terms = algebra.clean(terms)

    def __add__(self, other):
        other = self.algebra.element(other)
        terms = self.terms.copy()
        for word,coefficient in other.terms.items():
            terms[word] = terms.get(word,0)+coefficient
        return Element(self.algebra,terms)

    __radd__ = __add__

    def __neg__(self):
        return Element(self.algebra,{w:-c for w,c in self.terms.items()})

    def __sub__(self,other):
        return self+-self.algebra.element(other)

    def __rsub__(self,other):
        return self.algebra.element(other)+-self

    def __mul__(self,other):
        other = self.algebra.element(other)
        result = {}
        for w,c in self.terms.items():
            for v,d in other.terms.items():
                for word,factor in self.algebra.reduce_word(w+v):
                    result[word] = result.get(word,0)+c*d*factor
        return Element(self.algebra,result)

    def __rmul__(self,other):
        return self.algebra.element(other)*self

    def __pow__(self, exponent):
        if not isinstance(exponent,int) or exponent < 0:
            raise ValueError('Only nonnegative integer powers belong to this engine')
        result = self.algebra.element(1)
        factor = self
        while exponent:
            if exponent % 2:
                result = result*factor
            exponent //= 2
            if exponent:
                factor = factor*factor
        return result

    def commutator(self,other):
        return self*other-other*self

    def __bool__(self):
        return bool(self.terms)

    def certificate(self):
        return {' '.join(self.algebra.letters[i] for i in w) or '1': str(c)
                for w,c in sorted(self.terms.items())}

    def __repr__(self):
        if not self.terms:
            return '0'
        return ' + '.join(f'({c})*'+('*'.join(self.algebra.letters[i] for i in w) or '1')
                          for w,c in sorted(self.terms.items()))
