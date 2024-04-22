def run(finp: str):
    import time
    import os
    import apns.module_workflow.workflow_test.initialize as amwinit
    import apns.module_workflow.workflow_test.apns_itertools as amwai
    import apns.module_workflow.workflow_test.iterate as amwi
    import apns.module_io.compress as amic
    import apns.module_io.abacustest as amia
    import apns.module_io.citation as amicite
    # after import, read input file
    test_attribs, structures, upfs, orbs = amwinit.initialize(finp)

    software, calculation, extensive = test_attribs["global"]["software"].lower(), test_attribs["calculation"], test_attribs["extensive"]
    iterable_settings = dict(zip(["systems", "pseudopotentials", "numerical_orbitals", "calculation_settings", "extensive_settings"],
                                 [structures, upfs, orbs, calculation, extensive]))
    iterate_setting = dict(zip(["upforb_bundles", "calculation_settings", "extensive_settings"], amwai.setup_iterables(**iterable_settings)))
    iterate_setting.update(dict(zip(["software", "systems", "upfs", "orbs"], [software, structures, upfs, orbs])))
    folders = amwi.iterate(**iterate_setting)

    amic.pack(folders, "apns_{}.zip".format(time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())))
    os.system("rm -rf {}".format(" ".join(folders)))
    amia.image_information(software)
    amicite.citation(software)

import unittest
class TestRun(unittest.TestCase):
    def test_run(self):
        print(f"Running {__file__}.... This is a pure workflow function, therefore no test is needed.")

if __name__ == "__main__":
    unittest.main()
