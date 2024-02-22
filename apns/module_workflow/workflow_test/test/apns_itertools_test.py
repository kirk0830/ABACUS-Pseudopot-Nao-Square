import unittest
import apns.module_workflow.workflow_test.apns_itertools as amwai

class TestApnsItertools(unittest.TestCase):

    system_list = [
        "TiO2_123", "H2PtCl6_456", "ErCl3_789"
    ]
    vpspot = {
        "Ti": {
            "sg15_10": {
                "file": "Ti_ONCV_PBE-1.0.upf",
                "kind": "sg15",
                "version": "10",
                "appendix": ""
            },
            "sg15_11": {
                "file": "Ti_ONCV_PBE-1.1.upf",
                "kind": "sg15",
                "version": "11",
                "appendix": ""
            }
        },
        "O": {
            "sg15_10": {
                "file": "O_ONCV_PBE-1.0.upf",
                "kind": "sg15",
                "version": "10",
                "appendix": ""
            },
            "sg15_11": {
                "file": "O_ONCV_PBE-1.1.upf",
                "kind": "sg15",
                "version": "11",
                "appendix": ""
            }
        },
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
            }
        },
        "Cl": {
            "sg15_10": {
                "file": "Cl_ONCV_PBE-1.0.upf",
                "kind": "sg15",
                "version": "10",
                "appendix": ""
            },
            "sg15_11": {
                "file": "Cl_ONCV_PBE-1.1.upf",
                "kind": "sg15",
                "version": "11",
                "appendix": ""
            }
        },
        "H": {
            "sg15_10": {
                "file": "H_ONCV_PBE-1.0.upf",
                "kind": "sg15",
                "version": "10",
                "appendix": ""
            },
            "sg15_11": {
                "file": "H_ONCV_PBE-1.1.upf",
                "kind": "sg15",
                "version": "11",
                "appendix": ""
            }
        },
        "Pt": {
            "sg15_10": {
                "file": "Pt_ONCV_PBE-1.0.upf",
                "kind": "sg15",
                "version": "10",
                "appendix": ""
            },
            "sg15_11": {
                "file": "Pt_ONCV_PBE-1.1.upf",
                "kind": "sg15",
                "version": "11",
                "appendix": ""
            }
        }
    }
    vnao = {element: {} for element in vpspot.keys()}
    def test_systems(self):
        print(amwai.systems(system_list=self.system_list, 
                            valid_pseudopotentials=self.vpspot,
                            valid_numerical_orbitals=self.vnao))
        
    def test_extensive(self):
        print(amwai.extensive(extensive_settings={
            "characteristic_lengths": [10, 20, 30],
            "nkpoints_in_line": 10,
            "magnetism": "antiferromagnetic"
        }))

if __name__ == "__main__":
    unittest.main()