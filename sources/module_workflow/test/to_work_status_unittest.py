"""Please DO NOT ALWAYS EXECUTE THIS UNITTEST
or API will be blocked permenantly someday.
"""

from module_workflow.from_work_status import to_work_status as tws
import json
import unittest
import os

class TestToWorkStatus(unittest.TestCase):

    input_Er = """{
    "global": {
        "test_mode": "pseudopotential",
        "software": "ABACUS",
        "work_dir": "./",
        "pseudo_dir": "./module_pseudo/download",
        "orbital_dir": "./module_nao/download",
        "save_log": true
    },
    "calculation": {
        "basis_type": "pw",
        "functionals": ["PBE"],
        "ecutwfc": [100],
        "cell_scaling": [0.00]
    },
    "systems": ["Er2O3"],
    "pseudopotentials": {
        "kinds": ["sg15", "pd"],
        "versions": ["10", "04"],
        "appendices": [""]
    },
    "numerical_orbitals": {
        "types": ["DZP"],
        "rcuts": [7, 8, 9, 10],
        "appendices": [""]
    }
}"""
    with open("input.json", "w") as f:
        f.write(input_Er)
    def test_to_work_status(self):
        work_status = tws(
                fname="input.json",
                api_key="wV1HUdmgESPVgSmQj5cc8WvttCO8NTHp",
                num_cif=1
            )
        with open("work_status.json", "w") as f:
            json.dump(work_status, f, indent=4)

        self.assertTrue(os.path.exists("work_status.json"))
        self.assertEqual(work_status["calculation"]["basis_type"], "pw")
        self.assertEqual(work_status["calculation"]["functionals"], ["PBE"])
        self.assertEqual(work_status["calculation"]["ecutwfc"], [100])
        self.assertEqual(work_status["calculation"]["cell_scaling"], [0.00])
        self.assertEqual(work_status["systems"][0].startswith("Er2O3"), True)
        self.assertEqual(work_status["pseudopotentials"]["Er"]["kinds"], ["sg15", "pd"])
        self.assertEqual(work_status["pseudopotentials"]["Er"]["versions"], ["10", "04"])
        self.assertEqual(work_status["pseudopotentials"]["Er"]["appendices"], [""])
        self.assertEqual(work_status["pseudopotentials"]["O"]["kinds"], ["sg15", "pd"])
        self.assertEqual(work_status["pseudopotentials"]["O"]["versions"], ["10", "04"])
        self.assertEqual(work_status["pseudopotentials"]["O"]["appendices"], [""])
        self.assertEqual(work_status["numerical_orbitals"]["DZP"]["types"], ["DZP"])
        self.assertEqual(work_status["numerical_orbitals"]["DZP"]["rcuts"], [7, 8, 9, 10])
        self.assertEqual(work_status["numerical_orbitals"]["DZP"]["appendices"], [""])
        self.assertEqual(work_status["global"]["test_mode"], "pseudopotential")
        self.assertEqual(work_status["global"]["software"], "ABACUS")

        os.remove("work_status.json")
        os.remove("input.json")
        os.remove("Er2O3_*.cif")
        
if __name__ == "__main__":

    unittest.main()
