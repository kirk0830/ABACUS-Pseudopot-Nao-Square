import unittest
import apns.module_pseudo.parse as ampgp

class TestGeneralParser(unittest.TestCase):

    def test_is_numeric_data(self):

        line1 = "  1.00000000E+00  2.00000000E+00  3.00000000E+00\n  4.00000000E+00  5.00000000E+00  6.00000000E+00\n"
        self.assertTrue(ampgp.is_numeric_data(line1))
        line2 = "Generated using \"atomic\" code by A. Dal Corso  v.5.0.99 svn rev. 10869\n    Author: TM + Modified by ADC\n    Generation date: 10Apr2014\n    Pseudopotential type: NC\n    Element: Ac\n    Functional: PBE\n\n    Suggested minimum cutoff for wavefunctions:   0. Ry\n    Suggested minimum cutoff for charge density:   0. Ry\n    The Pseudo was generated with a Scalar-Relativistic Calculation\n    L component and cutoff radius for Local Potential:  0   3.0000\n\n    Valence configuration: \n    nl pn  l   occ       Rcut    Rcut US       E pseu\n    6D  3  2  1.00      3.000      3.000    -0.165522\n    7S  1  0  2.00      4.134      4.134    -0.289937\n    Generation configuration:\n    7P  2  1  0.00      4.000      4.000     0.100000\n    6D  3  2  1.00      3.000      3.000    -0.165522\n    5F  4  3  0.00      3.000      3.000     0.100000\n    7S  1  0  2.00      3.000      3.000    -0.289938\n\n    Pseudization used: troullier-martins"
        self.assertFalse(ampgp.is_numeric_data(line2))
        line3 = "Generated using \"atomic\" code by A. Dal Corso  v.5.0.99 svn rev. 10869"
        self.assertFalse(ampgp.is_numeric_data(line3))
        line4 = "    6D  3  2  1.00      3.000      3.000    -0.165522"
        self.assertFalse(ampgp.is_numeric_data(line4))
        line5 = "1 2 3\n 5 4 6\n"
        self.assertTrue(ampgp.is_numeric_data(line5))
        line6 = "1"
        self.assertTrue(ampgp.is_numeric_data(line6))
        line7 = "1 2 3"
        self.assertTrue(ampgp.is_numeric_data(line7))
        line8 = "0."
        self.assertTrue(ampgp.is_numeric_data(line8))
        line9 = "0. 0. 0. 0.\n 0. 0. 0. 0.\n 0. 0. 0. 0.\n 0. 0. 0. 0.\n"
        self.assertTrue(ampgp.is_numeric_data(line9))
        line10 = "1.23456e-10"
        self.assertTrue(ampgp.is_numeric_data(line10))

    def test_decompose_data(self):

        line1 = "  1.00000000E+00  2.00000000E+00  3.00000000E+00\n  4.00000000E+00  5.00000000E+00  6.00000000E+00\n"
        self.assertEqual(ampgp.decompose_data(line1), [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        line2 = "Generated using \"atomic\" code by A. Dal Corso  v.5.0.99 svn rev. 10869"
        self.assertRaises(ValueError, ampgp.decompose_data, line2)
        line3 = "    6D  3  2  1.00      3.000      3.000    -0.165522"
        self.assertRaises(ValueError, ampgp.decompose_data, line3)
        line4 = "1 2 3\n 5 4 6\n"
        self.assertEqual(ampgp.decompose_data(line4), [1.0, 2.0, 3.0, 5.0, 4.0, 6.0])
        line5 = "1"
        self.assertEqual(ampgp.decompose_data(line5), 1.0)
        line6 = "1 2 3"
        self.assertEqual(ampgp.decompose_data(line6), [1.0, 2.0, 3.0])
        line7 = "0."
        self.assertEqual(ampgp.decompose_data(line7), 0.0)
        line8 = "0. 0. 0. 0.\n 0. 0. 0. 0.\n 0. 0. 0. 0.\n 0. 0. 0. 0.\n"
        self.assertEqual(ampgp.decompose_data(line8), [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        line9 = "1.23456e-10"
        self.assertEqual(ampgp.decompose_data(line9), 1.23456e-10)

    def test_parse(self):

        fname_oncv = "./download/pseudopotentials/nc-fr-04_pbe_standard/Ag.upf"
        parsed_oncv = ampgp.parse(fname_oncv)
        self.assertEqual(parsed_oncv["PP_HEADER"]["attrib"]["pseudo_type"], "NC")
        self.assertEqual(parsed_oncv["PP_HEADER"]["attrib"]["number_of_proj"], 10)
        self.assertEqual(parsed_oncv["PP_HEADER"]["attrib"]["l_max"], 2)
        self.assertEqual(parsed_oncv["PP_HEADER"]["attrib"]["l_local"], -1)
        self.assertEqual(parsed_oncv["PP_HEADER"]["attrib"]["mesh_size"], 1538)
        self.assertEqual(parsed_oncv["PP_HEADER"]["attrib"]["number_of_wfc"], 6)

        fname_hgh = "./download/pseudopotentials/hgh/Al.pbe-hgh.UPF"
        parsed_hgh = ampgp.parse(fname_hgh)
        self.assertEqual(parsed_hgh["PP_HEADER"]["attrib"]["pseudo_type"], "NC")

    def test_determine_code_1(self):

        fname_oncv = "./download/pseudopotentials/nc-fr-04_pbe_standard/Ag.upf"
        oncv = ampgp.parse(fname_oncv)
        self.assertEqual(ampgp.determine_code(oncv), "ONCVPSP")
    
    def test_determine_code_2(self):
        fname_hgh = "./download/pseudopotentials/hgh/Al.pbe-hgh.UPF"
        hgh = ampgp.parse(fname_hgh)
        self.assertEqual(ampgp.determine_code(hgh), "GTH")
    
    def test_determine_code_3(self):
        fname_adc = "./download/pseudopotentials/pslnc_031/Ac.pbe-n-nc.UPF"
        adc = ampgp.parse(fname_adc)
        self.assertEqual(ampgp.determine_code(adc), "ADC")

    def test_is_compatible(self):

        fname_oncv = "./download/pseudopotentials/nc-fr-04_pbe_standard/Ag.upf"
        self.assertTrue(ampgp.is_compatible(fname_oncv))
        fname_hgh = "./download/pseudopotentials/hgh/Al.pbe-hgh.UPF"
        self.assertFalse(ampgp.is_compatible(fname_hgh))
        fname_adc = "./download/pseudopotentials/pslnc_031/Ac.pbe-n-nc.UPF"
        self.assertTrue(ampgp.is_compatible(fname_adc))

    def test_valence_configuration(self):
        fname_oncv = "./download/pseudopotentials/nc-fr-04_pbe_standard/Ag.upf"
        result = ampgp.valence_configuration(fname_oncv)
        reference = [['5S', '4S'], ['4P'], ['4D']]
        self.assertListEqual(result, reference)
        fname_oncv = "./download/pseudopotentials/NCPP-PD04-PBE/3+_f--core-icmod1/Er3+_f--core-icmod1.PD04.PBE.UPF"
        result = ampgp.valence_configuration(fname_oncv)
        reference = [['6S', '5S'], ['5P'], ['5D']]
        self.assertListEqual(result, reference)
        fname_oncv = "./download/pseudopotentials/NCPP-PD04-PBE/sp/Er-sp.PD04.PBE.UPF"
        result = ampgp.valence_configuration(fname_oncv)
        reference = [['6S', '5S'], ['5P'], ['5D'], ['4F']]
        self.assertListEqual(result, reference)

if __name__ == "__main__":
    unittest.main(buffer=False)