import unittest
import apns.module_nao.orbgen_kernel.siab_new as amnsn

class TestSIABNew(unittest.TestCase):

    def test_parse_reference(self):

        result = amnsn.parse_reference("apns/module_nao/orbgen_kernel/test/support/SIAB_INPUT")
        self.assertDictEqual(result, 
                             {'environment': '', 
                              'mpi_command': 'mpirun -np 1', 
                              'abacus_command': 'abacus', 
                              'pseudo_dir': '/root/abacus-develop/pseudopotentials/SG15_ONCV_v1.0_upf', 
                              'pseudo_name': 'Si_ONCV_PBE-1.0.upf', 
                              'ecutwfc': 100, 
                              'bessel_nao_rcut': [6, 7], 
                              'smearing_sigma': 0.01, 
                              'optimizer': 'pytorch.SWAT', 
                              'max_steps': [200], 
                              'spillage_coeff': [0.5, 0.5], 
                              'reference_systems': [
                                  {'shape': 'dimer', 'nbands': 8, 'nspin': 1, 'bond_lengths': [1.8, 2.0, 2.3, 2.8, 3.8]}, 
                                  {'shape': 'trimer', 'nbands': 10, 'nspin': 1, 'bond_lengths': [1.9, 2.1, 2.6]}
                               ], 
                               'orbitals': [
                                   {'zeta_notation': 'Z', 'shape': 'dimer', 'nbands_ref': 4, 'orb_ref': 'none'}, 
                                   {'zeta_notation': 'DZP', 'shape': 'dimer', 'nbands_ref': 4, 'orb_ref': 'Z'}, 
                                   {'zeta_notation': 'TZDP', 'shape': 'trimer', 'nbands_ref': 6, 'orb_ref': 'DZP'}
                               ]
                            })
    def test_generate(self):
        result = amnsn.generate(rcuts=[6, 7, 8, 9, 10],
                                ecutwfc=10086,
                                pseudo_name="Fe_ONCV_PBE-1.2.upf",
                                pseudo_dir="./download/pseudopotentials/sg15_oncv_upf_2020-02-06/1.2",
                                fref="./apns/module_nao/orbgen_kernel/data/Fe.SIAB_INPUT")
        self.assertEqual(result["pseudo_name"], "Fe_ONCV_PBE-1.2.upf")
        self.assertEqual(result["pseudo_dir"], "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/1.2")
        self.assertEqual(result["ecutwfc"], 10086)
        self.assertListEqual(result["bessel_nao_rcut"], [6, 7, 8, 9, 10])
              

if __name__ == "__main__":
    unittest.main()