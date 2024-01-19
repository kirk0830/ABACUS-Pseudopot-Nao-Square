import unittest
import apns.module_analysis.walk_gather as amawg

class WalkGatherTest(unittest.TestCase):

    def test_wash_pseudopot_name(self):
        self.assertEqual(amawg.wash_pseudopot_name("pd04spd"), "PD04-spd")
        self.assertEqual(amawg.wash_pseudopot_name("pd04f3+--icmod"), "PD04-f3+--icmod")
        self.assertEqual(amawg.wash_pseudopot_name("dojo03fr"), "DOJO03 (Full Relativistic)")

    def test_initialize_test_result(self):
        test_result = amawg.initialize_test_result(["energy", "natom"])
        self.assertDictEqual(test_result, {'energy': 0.0, 'natom': 0})

    def test_sort_by_ecutwfc(self):
        result = {
            "test1": {
                "ecutwfc": [10, 30, 500, 20],
                "energy": [1, 3, 5, 2],
                "pressure": [10, 30, 50, 20]
            },
            "test2": {
                "ecutwfc": [50, 764, 31, 65, 10],
                "energy": [10, 20, 30, 40, 50],
                "pressure": [100, 200, 300, 400, 500]
            }
        }
        result = amawg.sort_by_ecutwfc(result)
        self.assertListEqual(list(result["test1"]["ecutwfc"]), [10, 20, 30, 500])
        self.assertListEqual(list(result["test1"]["energy"]), [1, 2, 3, 5])
        self.assertListEqual(list(result["test1"]["pressure"]), [10, 20, 30, 50])
        self.assertListEqual(list(result["test2"]["ecutwfc"]), [10, 31, 50, 65, 764])
        self.assertListEqual(list(result["test2"]["energy"]), [50, 30, 10, 40, 20])
        self.assertListEqual(list(result["test2"]["pressure"]), [500, 300, 100, 400, 200])

if __name__ == '__main__':
    unittest.main()