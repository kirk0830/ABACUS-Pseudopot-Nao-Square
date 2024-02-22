import unittest
import apns.module_workflow.workflow_test.initialize as amwinit

class TestInitialize(unittest.TestCase):

    def test_pspot_software_availability(self):
        inp = {
            "global": {
                "software": "abacus"
            }
        }
        valid_pspot = {
            "Cl": {
                "hgh": {
                    "file": "Cl.pbe-hgh.UPF",
                    "kind": "hgh",
                    "version": "",
                    "appendix": ""
                }
            },
            "C": {
                "hgh": {
                    "file": "C.pbe-hgh.UPF",
                    "kind": "hgh",
                    "version": "",
                    "appendix": ""
                }
            }
        }
        pspot_arch = {
            "hgh": "./download/pseudopotentials/hgh"
        }
        self.assertFalse(amwinit.pspot_software_availability(inp, valid_pspot, pspot_arch))

if __name__ == "__main__":
    unittest.main()