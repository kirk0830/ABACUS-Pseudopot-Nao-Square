import apns.module_pseudo.parse_special.Atomic_A_Dal_Corso as ADC_parser
import unittest
import json
class Test_ADC_parser(unittest.TestCase):

    def test_valence_configuration(self):
        with open("./apns/module_pseudo/special_parser/test/support/Ac.pbe-n-nc.json", "r") as f:
            parsed = json.load(f)
        result = ADC_parser.valence_configuration(parsed)
        self.assertGreater(len(result), 0)
        self.assertEqual(len(result[0]), 1)
        self.assertEqual(len(result[1]), 0)
        self.assertEqual(len(result[2]), 1)

if __name__ == "__main__":
    unittest.main()