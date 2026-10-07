import tempfile
from pathlib import Path
import unittest
import external_validation as e

cache=tempfile.TemporaryDirectory(prefix='qhat-transfer-tests-')
e.m.init_numerics(Path(cache.name))


class ExternalValidationTests(unittest.TestCase):
    def test_fixed_input_count_and_distinct_families(self):
        specs=e.specifications()
        self.assertEqual(len(specs),8)
        self.assertEqual(len({s['case'] for s in specs}),8)
        self.assertEqual({s['family'] for s in specs},{'HCN','CH2O','H2S','SiH4'})
        self.assertEqual(len(specs)*len(e.m.BASE)*len(e.TIMES)*len(e.STEPS),216)

    def test_filename_electron_metadata(self):
        for spec in e.specifications():
            parsed=e.m.baseline.parse_case_metadata(Path('/tmp')/e.tensor_relative(spec),10)
            self.assertEqual(parsed.active_occupied,4)
            self.assertEqual(parsed.active_vacant,6)
            self.assertEqual(parsed.molecule,spec['family'])

    def test_geometry_scaling(self):
        specs=e.specifications()
        for i in range(0,8,2):
            for (a,x),(b,y) in zip(specs[i]['geometry'],specs[i+1]['geometry']):
                self.assertEqual(a,b)
                e.m.np.testing.assert_allclose(e.m.np.array(x)*1.25,y)

    def test_pairwise_commutation_structure(self):
        coeff={((0,'X'),(1,'X')):.2,((0,'Y'),(1,'Y')):.2}
        snapshot={'raw_pauli_keys':list(coeff),'coefficients':e.m.canonical_terms(coeff),
                  'parent_buckets':[(0,[0,1])],'fallback':[]}
        result=e.structure(snapshot)
        self.assertEqual(result['noncommuting_pairs'],0)
        self.assertAlmostEqual(result['max_Nalpha_Nbeta_N_commutator_l1'][2],0.)
        # This hopping connects opposite spins: total N is conserved, each spin N is not.
        self.assertGreater(result['max_Nalpha_Nbeta_N_commutator_l1'][0],0.)


if __name__=='__main__':unittest.main()
