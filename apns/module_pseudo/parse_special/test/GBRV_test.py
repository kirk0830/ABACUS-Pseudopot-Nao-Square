import unittest
import apns.module_pseudo.parse_special.GBRV_Vanderbilt as GBRV_parser
import json

class TestGBRV(unittest.TestCase):

    def test_valence(self):
        with open("./apns/module_pseudo/special_parser/test/support/Os_pbe_v1.2.uspp.f.json", "r") as f:
            parsed = json.load(f)

        result = GBRV_parser.valence(parsed)
        self.assertGreater(len(result), 0)
        self.assertEqual(len(result[0]), 2)
        self.assertEqual(len(result[1]), 2)
        self.assertEqual(len(result[2]), 1)

    def test_pp_header(self):

        test_str="""0                   Version Number
        As                   Element
        US                  Ultrasoft pseudopotential
        T                  Nonlinear Core Correction
        SLA  PW   PBE  PBE     PBE  Exchange-Correlation functional
        5.00000000000      Z valence
        -39.93046173160      Total energy
        0.00000    0.00000 Suggested cutoff for wfc and rho
        2                  Max angular momentum component
        875                  Number of points in mesh
        2    6             Number of Wavefunctions, Number of Projectors
        Wavefunctions         nl  l   occ
                            4S  0  2.00
                            4P  1  3.00"""

        result = GBRV_parser._PP_HEADER_(test_str)
        self.assertGreater(len(result["attributes"]), 0)

    def test_pp_info(self):

        test_str="""Generated using Vanderbilt code, version   7  3  6                              
        Author: kfg        Generation date:   13    5 2013                              
        Automatically converted from original format                                    
            1        The Pseudo was generated with a Scalar-Relativistic Calculation
        1.50000000000E+00    Local Potential cutoff radius
        nl pn  l   occ               Rcut            Rcut US             E pseu
        4S  4  0  2.00      0.00000000000      1.50000000000     -1.06527890600
        4P  4  1  3.00      0.00000000000      1.50000000000     -0.38182955700"""

        result = GBRV_parser._PP_INFO_(test_str)
        self.assertGreater(len(result["attributes"]), 0)

if __name__ == "__main__":
    unittest.main()