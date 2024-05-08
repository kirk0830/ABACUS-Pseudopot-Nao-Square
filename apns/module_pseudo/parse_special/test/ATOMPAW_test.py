import unittest
import apns.module_pseudo.parse_special.ATOMPAW_wentzcovitch as ATOMPAW_parser
import json

class TestATOMPAW(unittest.TestCase):

    def test_valence(self):
        with open("./apns/module_pseudo/special_parser/test/support/Ce.paw.z_12.atompaw.wentzcovitch.v1.2.json", "r") as f:
            parsed = json.load(f)

        result = ATOMPAW_parser.valence(parsed)
        self.assertGreater(len(result), 0)
        self.assertEqual(len(result[0]), 2)
        self.assertEqual(len(result[1]), 1)
        self.assertEqual(len(result[2]), 1)
        self.assertEqual(len(result[3]), 1)

if __name__ == "__main__":
    unittest.main()