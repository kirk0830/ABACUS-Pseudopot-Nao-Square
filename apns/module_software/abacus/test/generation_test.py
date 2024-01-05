import unittest
from apns.module_software.abacus import generation as abacus

class TestAbacus(unittest.TestCase):
    """Test the generation of template files for Abacus
    """
    def test_INPUT(self):
        """Test the generation of INPUT file new generation method
        """
        

if __name__ == "__main__":
    print(abacus.INPUT(work_status={
            "calculation": {
                "basis_type": "pw",
                "functionals": ["PBE"],
                "ecutwfc": [100],
                "cal_force": 1,
                "cal_stress": 1,
                "cell_scaling": [0.0]
            },
            "additional_keywords": {

                "kspacing": ["0.5", "0.25", "0.125"]
            }
        }, template=True))