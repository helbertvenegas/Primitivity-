"""Regression controls for the added Laurent, lattice and degree machinery."""
from dataclasses import replace
import unittest
import sympy as sp
from ore14641 import FAMILIES,p,q,f,g,h
from ore14641.words import (WordModel,monomial_model,normalizer_model,TorusWords,
                            alpha,beta,H,S,F)
from ore14641.degrees import (degree_certificate,degree_from_parameters,center_lattice,
                              monomial_multiply,monomial_power,c_bound_certificate)
from ore14641.torus import alternating,torus_certificate,brute_radical_size,l_matrix
from ore14641.algebraic import c_module_reduce,C_BASIS,dg_central_elements,right_weight_matrix
from ore14641.pbw import Algebra
from ore14641.primitivity import conditional_orbit_certificate


class CompletionTests(unittest.TestCase):
    def test_coordinate_inverse_has_correct_coefficient(self):
        q0=sp.Symbol('q0')
        mixing={(a,b):(b,a,sp.Integer(1)) for a in range(2) for b in range(2)}
        m=WordModel(sp.cancel,'coordinate','coordinate',q0,q0,mixing)
        for word in (m.e1,m.e2,m.z1,m.z2,m.e2*m.z2*m.e1**-1):
            self.assertTrue((word*word**-1).equals(m.identity()))
            self.assertTrue((word**-1*word).equals(m.identity()))
        self.assertEqual((m.e2**-1*m.e1).normal(),((-1,0,0,0),sp.Integer(1)))

    def test_l_first_sign_matches_defining_relations(self):
        m=monomial_model('L')
        w=m.e1*m.e2
        self.assertEqual(w.commutation(m.z1),1/(f*h**2))
        self.assertEqual(l_matrix(12,1,2,3)[0,1],-8)

    def test_center_lattice_with_composite_modulus(self):
        for order in (2,4,6):
            matrix=alternating((2,0,2,0,0,2))
            lattice=center_lattice(matrix,order)
            self.assertFalse(any(int(v)%order for v in matrix*lattice))
            self.assertEqual(abs(int(lattice.det())),order**4//brute_radical_size(matrix,order))

    def test_torus_cocycle_inverse_and_associativity(self):
        matrix=alternating((1,2,-1,3,1,2));N=12
        a=(3,(2,-1,3,0));b=(7,(-1,4,0,2));c=(1,(1,1,-2,3))
        self.assertEqual(monomial_multiply(matrix,N,monomial_multiply(matrix,N,a,b),c),
                         monomial_multiply(matrix,N,a,monomial_multiply(matrix,N,b,c)))
        self.assertEqual(monomial_multiply(matrix,N,a,monomial_power(matrix,N,a,-1)),(0,(0,0,0,0)))

    def test_all_closed_models_work_above_old_cutoffs(self):
        for name in ('B','E','J','K','L','M','N','O','P','S','T','U','W','Z','F','I','Q','R','V','D','G'):
            assignment={p:1 if name in ('D','G') else 15,q:15,H:7,S:3,F:2,alpha:1,beta:5}
            result=degree_certificate(name,60,assignment,False)
            self.assertTrue(result['agrees'],name)

    def test_large_odd_parity_branch_without_enumeration(self):
        result=degree_certificate('O',3988,{alpha:4,beta:8},False)
        self.assertEqual(result['center_action_order'],2)
        self.assertEqual(result['algebra_pi_degree'],2*997**2)

    def test_w_z_equal_parameters_and_duality(self):
        z=degree_from_parameters('Z',{f:3,h:-sp.Rational(1,2)})
        w=degree_from_parameters('W',{f:3,h:-sp.Rational(1,2)})
        dual=degree_from_parameters('W',{f:sp.Rational(1,3),h:-sp.sqrt(3)/2})
        self.assertEqual((z['algebra_pi_degree'],w['algebra_pi_degree'],dual['algebra_pi_degree']),(6,12,6))

    def test_f_i_degree_and_table_bound_are_distinct(self):
        for name in ('F','I'):
            first=degree_from_parameters(name,{h:1/sp.sqrt(2)})
            second=degree_from_parameters(name,{h:(1+sp.I)/2})
            self.assertEqual(first['algebra_pi_degree'],16)
            self.assertEqual(second['algebra_pi_degree'],8)

    def test_c_is_reported_as_open_with_exact_module_reduction(self):
        cert=c_bound_certificate(12,4,3)
        self.assertTrue(cert['exact_degree_is_open'])
        self.assertEqual(cert['possible_degree_ratios'],[1,3,9])
        A=Algebra(FAMILIES['C']);_,_,y1,y2=A.generators
        ca,cb=y1**2*y2,y1**3-y2**3
        for i,j in [(13,0),(0,12),(8,5)]:
            result=A.element(0)
            for (a,b,k),coefficient in c_module_reduce(i,j).items():
                ii,jj=C_BASIS[k]
                result=result+coefficient*ca**a*cb**b*y1**ii*y2**jj
            self.assertFalse(y1**i*y2**j-result)

    def test_dg_central_construction_accepts_arbitrary_order(self):
        for name in ('D','G'):
            A,elements=dg_central_elements(name,7)
            for element in elements:
                self.assertTrue(all(not element.commutator(v) for v in A.generators))

    def test_conditional_orbit_keeps_other_parameter_free(self):
        result=conditional_orbit_certificate('O','u',5)
        self.assertTrue(result['central'])
        self.assertTrue(result['has_nonconstant_coefficient'])
        self.assertEqual(result['remaining_parameter'],'beta is symbolic and unconstrained')
        with self.assertRaises(ValueError):conditional_orbit_certificate('O','u',2)

    def test_invalid_nontorsion_and_inexact_inputs_rejected(self):
        with self.assertRaises(ValueError):degree_from_parameters('S',{h:1})
        with self.assertRaises(ValueError):degree_certificate('S',12,{H:1.5})
        with self.assertRaises(ValueError):degree_certificate('O',12,{alpha:1,beta:1})


if __name__=='__main__':unittest.main()
