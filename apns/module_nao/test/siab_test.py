import unittest
import apns.module_nao.siab as siab

class SiabTest(unittest.TestCase):

    def test_orbitalconfig_tolist(self):
        """single orbital"""
        self.assertEqual(siab.orbitalconfig_tolist("1s"), [1])
        self.assertEqual(siab.orbitalconfig_tolist("2p"), [0, 2])
        self.assertEqual(siab.orbitalconfig_tolist("3d"), [0, 0, 3])
        """multiple orbitals"""
        self.assertEqual(siab.orbitalconfig_tolist("1s2p"), [1, 2])
        self.assertEqual(siab.orbitalconfig_tolist("1s3d"), [1, 0, 3])
        self.assertEqual(siab.orbitalconfig_tolist("2p3d"), [0, 2, 3])
        self.assertEqual(siab.orbitalconfig_tolist("1s2p3d"), [1, 2, 3])
        """multiple orbitals with different number of orbitals"""
        self.assertEqual(siab.orbitalconfig_tolist("1s2p3d4f"), [1, 2, 3, 4])

    def test_zeta_notation_toorbitalconfig(self):
        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZ", [2, 1, 1]), "2s1p1d")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZP", [2, 1, 1]), "2s1p1d1f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZ9P", [2, 1, 1]), "2s1p1d9f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("DZ", [2, 1, 1]), "4s2p2d")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("DZP", [2, 1, 1]), "4s2p2d1f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("TZDP", [2, 1, 1]), "6s3p3d2f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("QZDP", [2, 1, 1]), "8s4p4d2f")

        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZ", [2, 0, 1, 1]), "2s1d1f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZP", [2, 0, 1, 1]), "2s1d1f1g")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZ9P", [2, 0, 1, 1]), "2s1d1f9g")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("DZ", [2, 0, 1, 1]), "4s2d2f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("DZP", [2, 0, 1, 1]), "4s2d2f1g")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("TZDP", [2, 0, 1, 1]), "6s3d3f2g")

if __name__ == "__main__":
    unittest.main()