import unittest
import apns.module_software.qespresso.generation as amsqg

class TestGeneration(unittest.TestCase):

    def test_calculation(self):
        kwargs = {"basis_type": "lcao", "dft_functional": "pbe", "ecutwfc": 200}
        result = amsqg._calculation(section="system", ntype=1, natom=10, **kwargs)
        self.assertGreater(len(result), 0)
    def test_ATOMIC_SPECIES(self):
        kwargs = {"mass": {"H": "1.008", "Pb": "157.25"}}
        result = amsqg._ATOMIC_SEPCIES(pseudopotential={"H": "H_ONCV_PBE-1.0.upf", "Pb": "Pb_ONCV_PBE-1.0.upf"}, **kwargs)
        self.assertGreater(len(result), 0)
    def test_CIF(self):
        result = amsqg._CIF("apns_cache/mp-679.cif")
        self.assertGreater(len(result), 0)
    def test_ISOLATED(self):
        print(amsqg._ISOLATED(element="H", shape="dimer", characteristic_length=2.0))
if __name__ == "__main__":
    unittest.main()