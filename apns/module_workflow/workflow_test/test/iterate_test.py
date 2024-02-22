import unittest
import apns.module_workflow.workflow_test.iterate as amwi

class TestIterate(unittest.TestCase):

    system = ["Er2O3_679"]
    pseudopot_nao_settings = [
            [
                {#                  for element 0, for element 1
                    "pseudopotential": ["sg15_10", "sg15_10"],
                    "numerical_orbital": ["TZDP_6", "TZDP_10"]
                }, # for combination 1
                {
                    "pseudopotential": ["sg15_11", "sg15_10"],
                    "numerical_orbital": ["TZDP_6", "TZDP_10"]
                }, # for combination 2
                {
                    "pseudopotential": ["pd_04", "sg15_10"],
                    "numerical_orbital": ["TZDP_10", "TZDP_10"]
                } # for combination 3
            ] # for system 1
        ]
    calculation_settings = [
            {
                "basis_type": "lcao",
                "dft_functional": "pbe",
                "ecutwfc": 100
            },
            {
                "basis_type": "pw",
                "dft_functional": "pbe",
                "ecutwfc": 200
            },
            {
                "basis_type": "pw",
                "dft_functional": "pbe",
                "ecutwfc": 300
            }
        ]
    extensive = {
        "characteristic_lengths": [0.00],
        "nkpoints_in_line": 0,
    }
    valid_pseudopotentials = {
            "Er": {
                "sg15_10": {
                    "file": "Er_ONCV_PBE-1.0.upf",
                    "kind": "sg15",
                    "version": "10",
                    "appendix": ""
                },
                "sg15_11": {
                    "file": "Er_ONCV_PBE-1.1.upf",
                    "kind": "sg15",
                    "version": "11",
                    "appendix": ""
                },
                "pd_04": {
                    "file": "Er.PD04.PBE.UPF",
                    "kind": "pd",
                    "version": "04",
                    "appendix": ""
                }
            },
            "O": {
                "sg15_10": {
                    "file": "O_ONCV_PBE-1.0.upf",
                    "kind": "sg15",
                    "version": "10",
                    "appendix": ""
                }
            }
        }
    valid_numerical_orbitals = {
            "Er": {
                "TZDP_6@sg15_10": {
                    "file": "Cr_gga_6au_100Ry_3s3p3d2f.orb",
                    "type": "TZDP",
                    "rcut": "6",
                    "appendix": ""
                },
                "TZDP_6@sg15_11": {
                    "file": "Cr_gga_6au_100Ry_3s3p3d2f.orb",
                    "type": "TZDP",
                    "rcut": "6",
                    "appendix": ""
                },
                "TZDP_10@pd_04": {
                    "file": "Cr_gga_10au_100Ry_3s3p3d2f.orb",
                    "type": "TZDP",
                    "rcut": "10",
                    "appendix": ""
                }
            },
            "O": {
                "TZDP_10@sg15_10": {
                    "file": "O_gga_10au_100Ry_3s3p3d2f.orb",
                    "type": "TZDP",
                    "rcut": "10",
                    "appendix": ""
                }
            }
        }
    pspot_archive = {
            "sg15_10": "/home/zhuzhen/Work/psp/sg15/1.0/",
            "sg15_11": "/home/zhuzhen/Work/psp/sg15/1.1/",
            "pd_04": "/home/zhuzhen/Work/psp/pd/04/"
        }
    nao_archive = {
            "TZDP_6@sg15_10": "/home/zhuzhen/Work/nao/sg15/1.0/",
            "TZDP_10@sg15_10": "/home/zhuzhen/Work/nao/sg15/1.0/",
            "TZDP_6@sg15_11": "/home/zhuzhen/Work/nao/sg15/1.1/",
            "TZDP_10@pd_04": "/home/zhuzhen/Work/nao/pd/04/"
        }
    
    def test_iterate_abacus(self):

        result = amwi.iterate(software="abacus",
                              systems=self.system,
                              pseudopot_nao_settings=self.pseudopot_nao_settings,
                              calculation_settings=self.calculation_settings,
                              extensive=self.extensive,
                              valid_pseudopotentials=self.valid_pseudopotentials,
                              valid_numerical_orbitals=self.valid_numerical_orbitals,
                              pspot_archive=self.pspot_archive,
                              nao_archive=self.nao_archive,
                              test_mode=True)
        self.assertGreater(len(result), 0)

    def test_iterate_qespresso(self):

        result = amwi.iterate(software="qespresso",
                              systems=self.system,
                              pseudopot_nao_settings=self.pseudopot_nao_settings,
                              calculation_settings=self.calculation_settings,
                              extensive=self.extensive,
                              valid_pseudopotentials=self.valid_pseudopotentials,
                              valid_numerical_orbitals=self.valid_numerical_orbitals,
                              pspot_archive=self.pspot_archive,
                              nao_archive=self.nao_archive,
                              test_mode=True)
        self.assertGreater(len(result), 0)

if __name__ == "__main__":
    unittest.main()