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
        "cal_stress": 1
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

    def test_STRU_Molecule(self):
        """Test the generation of template structure file for isolated system
        """
        _stru, _ = abacus.STRU_Molecule(shape="trimer",
                                       pseudopotentials=self.pseudopotentials,
                                       bond_length=3.0)
        self.assertGreater(len(_stru), 0)
        print(_stru)
    
    def test_KLINE(self):
        """Test the generation of template kline file for Abacus
        """
        _kline = abacus._KLINE_(fname="apns_cache/mp-160.cif",
                                nkpts_in_line=5)
        self.assertGreater(len(_kline), 0)
        print(_kline)

    def test_STRU_Pymatgen(self):
        """Test the generation of template structure file for Abacus
        """
        _stru, _ = abacus.STRU_Pymatgen(fname="apns_cache/mp-160.cif",
                                     pseudopotentials={
                                         "B": "B_ONCV_PBE-1.0.upf",
                                     },
                                     numerical_orbitals=None,
                                     cell_scaling=1.0,
                                     starting_magnetization=None)
        self.assertGreater(len(_stru), 0)
        
        # test mp-8.cif, Re2
        starting_magnetization = [-1, 1] # antiferromagnetic
        _stru, _ = abacus.STRU_Pymatgen(fname="apns_cache/mp-8.cif",
                                     pseudopotentials={
                                         "Re": "Re_ONCV_PBE-1.0.upf",
                                     },
                                     numerical_orbitals=None,
                                     cell_scaling=1.0,
                                     starting_magnetization=starting_magnetization)
        self.assertGreater(len(_stru), 0)
        print(_stru)

if __name__ == "__main__":
    unittest.main()