import unittest

from apns.module_software.vasp import generation

class TestGeneration(unittest.TestCase):

    fname = "apns_cache/mp-568348.cif"
    def test_poscar(self):
        print(generation.POSCAR(fname=self.fname))

if __name__ == "__main__":
    unittest.main()