"""Independent negative controls: every altered witness must be rejected."""
import ast
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('independent_checker',ROOT/'independent/check.py')
checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.bundle=json.loads((ROOT/'reports/certificates.json').read_text())
    def reject(self,change):
        bundle=deepcopy(self.bundle);change(bundle)
        with self.assertRaises((ArithmeticError,ValueError,ZeroDivisionError)):checker.verify(bundle)
    def test_independent_full_bundle(self):
        report=checker.verify(self.bundle)
        self.assertEqual(report['failed'],0)
        self.assertEqual(report['source_families'],26)
        self.assertEqual(report['laurent_models'],21)
    def test_checker_import_boundary(self):
        tree=ast.parse((ROOT/'independent/check.py').read_text())
        imports={n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)}
        imports|={a.name for n in ast.walk(tree) if isinstance(n,ast.Import) for a in n.names}
        self.assertLessEqual(imports,{'ast','fractions','itertools','json','re','math','pathlib','sys'})
    def test_source_mutation_rejected(self):
        self.reject(lambda b:b['sources'][0]['sigma'][0].__setitem__(0,'2*h'))
    def test_projector_mutation_rejected(self):
        self.reject(lambda b:b['uniform']['identities'][0].__setitem__('rhs','2'))
    def test_smith_witness_mutation_rejected(self):
        self.reject(lambda b:b['integer_cases'][0]['U'][0].__setitem__(0,99))
    def test_kernel_completeness_mutation_rejected(self):
        self.reject(lambda b:b['integer_cases'][0]['kernel_change'][0].__setitem__(0,99))
    def test_laurent_action_scalar_mutation_rejected(self):
        self.reject(lambda b:b['uniform']['word_models'][0]['actions'][1]['scalars'].__setitem__(0,'2'))
    def test_exhaustive_parity_partition_mutation_rejected(self):
        self.reject(lambda b:b['uniform']['root_parity']['classes'].pop())
    def test_general_center_parametrization_mutation_rejected(self):
        self.reject(lambda b:b['uniform']['center_parametrizations'][0]['basis'][0].__setitem__(0,'4'))
    def test_claim_contract_rejects_count_only_coverage(self):
        from ore14641.coverage import coverage_report
        report=json.loads((ROOT/'reports/verification.json').read_text())
        with self.assertRaises(ArithmeticError):
            coverage_report([c for c in report['checks'] if c['id']!='A.M'])
        for entry in report['coverage']['entries']:
            if entry['label']=='prop:weights':
                self.assertTrue(entry['claim_contract']['hypotheses'])
                self.assertTrue(entry['claim_contract']['certificate_ids'])
    def test_quotient_and_fraction_arithmetic(self):
        r=checker.Ring(['p','h'],[('p**2+p+1','p')])
        self.assertEqual(r.parse('p**3'),1)
        self.assertEqual(r.parse('(p+1)**2'),r.parse('p'))
        self.assertNotEqual(r.parse('p'),r.parse('1'))
        self.assertEqual(r.parse('h**(-2)*h**3'),r.parse('h'))


if __name__=='__main__':unittest.main()
