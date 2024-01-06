import unittest
from apns.module_software.abacus import generation as abacus

class TestAbacus(unittest.TestCase):
    """Test the generation of template files for Abacus
    """
    calculation = {
        "basis_type": "pw",
        "functionals": "PBE",
        "ecutwfc": 100,
        "cal_force": 1,
        "cal_stress": 1,
        "cell_scaling": 0.0
    }

    pseudopotentials = {
        "Er": "Er_ONCV_PBE-1.0.upf"
    }

    def test_INPUT(self):
        """Test the generation of template input file for Abacus
        """
        _input = abacus.INPUT(calculation=self.calculation)
        self.assertGreater(len(_input), 0)
        print(_input)

    def test_STRU_ISOLATED(self):
        """Test the generation of template structure file for isolated system
        """
        _stru, _ = abacus._STRU_ISOLATED_(shape="trimer",
                                       pseudopotentials=self.pseudopotentials,
                                       bond_length=3.0)
        self.assertGreater(len(_stru), 0)
        print(_stru)
    

if __name__ == "__main__":
    unittest.main()