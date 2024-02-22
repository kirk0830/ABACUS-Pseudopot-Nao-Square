import unittest
import apns.module_workflow.workflow_orbgen.initialize as amwwoi

class TestInitialize(unittest.TestCase):

    def test_check(self):
        """written by Github.copilot"""
        with self.assertRaises(ValueError):
            amwwoi.check({"global": {}})
        with self.assertRaises(ValueError):
            amwwoi.check({"global": {"test_mode": ""}})
        with self.assertRaises(ValueError):
            amwwoi.check({"global": {"pseudo_dir": ""}})
        with self.assertRaises(ValueError):
            amwwoi.check({"global": {"orbital_dir": ""}})
        with self.assertRaises(ValueError):
            amwwoi.check({"orbgen": {}})
        with self.assertRaises(ValueError):
            amwwoi.check({"orbgen": {"generator": ""}})
        with self.assertRaises(ValueError):
            amwwoi.check({"orbgen": {"environment": {}}})
        with self.assertRaises(ValueError):
            amwwoi.check({"orbgen": {"mpi_command": ""}})
        with self.assertRaises(ValueError):
            amwwoi.check({"orbgen": {"abacus_command": ""}})
        with self.assertRaises(ValueError):
            amwwoi.check({"abacus": {}})
        with self.assertRaises(ValueError):
            amwwoi.check({"systems": []})
        with self.assertRaises(ValueError):
            amwwoi.check({"pseudopotentials": {}})
        with self.assertRaises(ValueError):
            amwwoi.check({"pseudopotentials": {"kinds": []}})
        with self.assertRaises(ValueError):
            amwwoi.check({"pseudopotentials": {"versions": []}})
        with self.assertRaises(ValueError):
            amwwoi.check({"pseudopotentials": {"appendices": []}})
        with self.assertRaises(ValueError):
            amwwoi.check({"numerical_orbitals": {}})
        with self.assertRaises(ValueError):
            amwwoi.check({"numerical_orbitals": {"rcuts": []}})
        with self.assertRaises(ValueError):
            amwwoi.check({"numerical_orbitals": {"types": []}})
        with self.assertRaises(ValueError):
            amwwoi.check({"global": {"pseudo_dir": "not_exist"}})
        with self.assertRaises(ValueError):
            amwwoi.check({"orbgen": {"generator": "not_exist"}})

    def test_link_pspotlib(self):
        
        inp = {
            "global": {
                "pseudo_dir": "./download/pseudopotentials"
            },
            "pseudopotentials": {
                "kinds": ["sg15", "pd"],
                "versions": ["10", "04"],
                "appendices": ["all"]
            }
        }
        result = amwwoi.link_pspotlib(inp)
        self.assertListEqual(result, 
                             ['sg15_10', 'sg15_10_fr', 'pd_04_3+_f--core', 
                              'pd_04_3+_f--core-icmod1', 'pd_04_d', 'pd_04', 
                              'pd_04_fsp', 'pd_04_high', 'pd_04_low', 'pd_04_s', 
                              'pd_04_s-high', 'pd_04_sp', 'pd_04_sp-exc', 
                              'pd_04_sp-high', 'pd_04_spd', 'pd_04_spd-high'])
        inp = {
            "global": {
                "pseudo_dir": "./download/pseudopotentials"
            },
            "pseudopotentials": {
                "kinds": ["sg15"],
                "versions": ["all"],
                "appendices": ["all"]
            }
        }
        result = amwwoi.link_pspotlib(inp)
        self.assertListEqual(result, 
                             ['sg15_10', 'sg15_10_fr', 'sg15_11', 'sg15_11_fr', 'sg15_12'])
        

if __name__ == '__main__':
    unittest.main()