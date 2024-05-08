import unittest
import apns.module_pseudo.parse_special.ONCVPSP_D_R_Hamann as ONCV_test
import json
class TestONCV(unittest.TestCase):

    def test_valence(self):

        with open("./apns/module_pseudo/special_parser/test/support/Ag.json", "r") as f:
            parsed = json.load(f)
        result = ONCV_test.valence(parsed)
        self.assertGreater(len(result), 0)
        self.assertEqual(len(result[0]), 2)
        self.assertEqual(len(result[1]), 1)
        self.assertEqual(len(result[2]), 1)

    def test_(self):
        test_str = """# ATOM AND REFERENCE CONFIGURATION
        # atsym  z   nc   nv     iexc    psfile
        As 33.00    5    3       4      upf
        #
        #   n    l    f        energy (Ha)
        1    0    2.00
        2    0    2.00
        2    1    6.00
        3    0    2.00
        3    1    6.00
        3    2   10.00
        4    0    2.00
        4    1    3.00
        #
        # PSEUDOPOTENTIAL AND OPTIMIZATION
        # lmax
        2
        #
        #   l,   rc,      ep,       ncon, nbas, qcut
        0   1.75000  -0.53265    3    8   5.30000
        1   1.70000  -0.19091    4    8   6.00000
        2   1.80000  -1.49370    4    7   9.00000
        #
        # LOCAL POTENTIAL
        # lloc, lpopt,  rc(5),   dvloc0
        4    5   1.60000      0.00000
        #
        # VANDERBILT-KLEINMAN-BYLANDER PROJECTORs
        # l, nproj, debl
        0    2   4.50000
        1    2   4.80000
        2    2   2.50000
        #
        # MODEL CORE CHARGE
        # icmod, fcfact, rcfact
        3   5.00000   1.40000
        #
        # LOG DERIVATIVE ANALYSIS
        # epsh1, epsh2, depsh
        -12.00   12.00    0.02
        #
        # OUTPUT GRID
        # rlmax, drl
        6.00    0.01
        #
        # TEST CONFIGURATIONS
        # ncnf
        0
        # nvcnf
        #   n    l    f"""

        result = ONCV_test._PP_INPUTFILE_(test_str)
        print(result)

if __name__ == "__main__":
    unittest.main()