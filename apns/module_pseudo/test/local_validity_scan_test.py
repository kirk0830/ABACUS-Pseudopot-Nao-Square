import unittest
import apns.module_pseudo.local_validity_scan as amplvs

class TestLocalValidityScan(unittest.TestCase):

    def test_svp(self):
        elements = ["C", "H", "O"]
        pseudopotentials = {
            "kinds": {
                "C": ["all"],
                "H": ["all"],
                "O": ["all"]
            },
            "versions": {
                "C": ["all"],
                "H": ["all"],
                "O": ["all"]
            },
            "appendices": {
                "C": ["all"],
                "H": ["all"],
                "O": ["all"]
            }
        }
        valid_pseudopotentials = amplvs._svp_(elements, pseudopotentials)
        self.assertEqual(len(valid_pseudopotentials), 3)

if __name__ == "__main__":
    unittest.main()