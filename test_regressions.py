"""Independent and negative checks for the scientific verification machinery."""
from dataclasses import replace
from itertools import product
import json
from pathlib import Path
import sys
import unittest
import sympy as sp
from ore14641 import FAMILIES,p,q,f,g,h,s
from ore14641.pbw import Algebra
from ore14641.monomial import MONOMIAL
from ore14641.pi import pi_decision,exact_root_order
from ore14641.torus import (alternating,torus_certificate,brute_radical_size,
                            torus_from_rational_phases,dg_matrix)

ROOT = Path(__file__).resolve().parents[1]


class RegressionTests(unittest.TestCase):
    def test_registry_complete(self):
        self.assertEqual(list(FAMILIES),list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
        self.assertEqual(len(MONOMIAL),14)

    def test_noncommutative_order_is_not_commutative_sympy(self):
        A=Algebra(FAMILIES['D'])
        x1,x2,_,_=A.generators
        self.assertFalse(x2*x1+x1*x2)
        self.assertTrue(x1*x2)
        self.assertTrue(x1*x2-x2*x1)

    def test_both_parameter_branches(self):
        for value in [sp.I,-sp.I]:
            family=FAMILIES['B'].specialize({p:value,h:1})
            self.assertEqual(family.detsigma(),value*sp.diag(1,-1))
            self.assertTrue(all(not v for v in Algebra(family).overlaps().values()))
        for value in [(-1+sp.sqrt(3)*sp.I)/2,(-1-sp.sqrt(3)*sp.I)/2]:
            family=FAMILIES['C'].specialize({p:value})
            self.assertFalse(any(family.matrix(family.detsigma()-family.expected_detsigma)))

    def test_invalid_parameters_and_floats(self):
        for name,subs in [('B',{p:2}),('G',{p:1}),('Z',{f:-1}),('K',{q:0})]:
            with self.assertRaises(ValueError):FAMILIES[name].specialize(subs)
        with self.assertRaises(ValueError):FAMILIES['B'].specialize({p:s})
        with self.assertRaises(TypeError):FAMILIES['S'].specialize({h:0.5})
        with self.assertRaises(TypeError):exact_root_order(1.0)

    def test_degenerate_member_does_not_use_singular_monomial_basis(self):
        family=FAMILIES['M'].specialize({f:0})
        self.assertTrue(all(not v for v in Algebra(family).overlaps().values()))
        with self.assertRaises(ValueError):MONOMIAL['M'].family.specialize({s:0})
        self.assertIs(pi_decision('M',{f:0,h:1})['status'],False)

    def test_mutated_source_matrix_is_detected(self):
        matrix=FAMILIES['Z'].sigma.as_mutable()
        matrix[0,0],matrix[0,2]=0,h
        family=replace(FAMILIES['Z'],sigma=matrix.as_immutable())
        self.assertTrue(any(Algebra(family).overlaps().values()))

    def test_root_orders_and_unknown_status(self):
        self.assertEqual(exact_root_order((sp.sqrt(3)+sp.I)/2),12)
        self.assertEqual(exact_root_order(sp.real_root(2,3)),0)
        self.assertEqual(exact_root_order(sp.pi),0)
        self.assertIsNone(exact_root_order(h))
        self.assertIsNone(pi_decision('B',{})['status'])
        self.assertIs(pi_decision('G',{p:2})['status'],False)

    def test_smith_independent_brute_force(self):
        # Composite orders and degenerate alternating matrices are intentional.
        for order in [2,3,4,6]:
            for upper in [(1,2,3,4,5,6),(0,0,0,0,0,0),(2,0,0,0,0,2)]:
                matrix=alternating(upper)
                certificate=torus_certificate(matrix,order)
                self.assertEqual(certificate['radical_size_mod_N'],brute_radical_size(matrix,order))
                self.assertEqual(certificate['pi_degree']**2,certificate['image_size'])

    def test_smith_unimodular_basis_and_lift(self):
        matrix=alternating((1,2,3,4,5,6))
        basis=sp.eye(4);basis[0,1]=3;basis[2,3]=-1
        for order in [4,6,9]:
            expected=torus_certificate(matrix,order)['pi_degree']
            self.assertEqual(torus_certificate(basis.T*matrix*basis,order)['pi_degree'],expected)
            self.assertEqual(torus_certificate(matrix+order*alternating((2,0,1,-1,3,0)),order)['pi_degree'],expected)

    def test_rational_phases_and_invalid_torus_inputs(self):
        certificate=torus_from_rational_phases([sp.Rational(1,3),0,0,0,0,sp.Rational(1,4)])
        self.assertEqual(certificate['pi_degree'],12)
        with self.assertRaises(ValueError):torus_certificate(sp.eye(2),3)
        with self.assertRaises(ValueError):dg_matrix('G',6,2,2)
        with self.assertRaises(TypeError):torus_from_rational_phases([0.5,0,0,0,0,0])

    def test_twist_needed(self):
        self.assertIs(pi_decision('S',{h:1})['status'],False)
        self.assertIs(pi_decision('S',{h:sp.Rational(1,2)})['status'],True)

    def test_notebook_sources_synchronized(self):
        sys.path.insert(0,str(ROOT/'scripts'))
        from build_notebook import portable_source
        import hashlib
        notebook=json.loads((ROOT/'notebooks'/'Double_Ore_14641_PI_verification.ipynb').read_text())
        checked=0
        for cell in notebook['cells']:
            if 'source_file' in cell['metadata']:
                source=(ROOT/cell['metadata']['source_file']).read_text()
                self.assertEqual(cell['metadata']['source_sha256'],hashlib.sha256(source.encode()).hexdigest())
                inline=''.join(cell['source'])
                self.assertEqual(inline.split('\n',1)[1].strip(),portable_source(source).strip())
                checked+=1
        self.assertEqual(checked,16)


if __name__ == '__main__':
    unittest.main()
