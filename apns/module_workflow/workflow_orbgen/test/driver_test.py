import unittest
import apns.module_workflow.workflow_orbgen.driver as amwwod
import os
class TestDriver(unittest.TestCase):

    def test_exsitu_siab(self):
        
        fscript = amwwod.exsitu_siab("./SIAB/SIAB_nouvelle.py", [
            "SIAB_INPUT_20210101000000.json", "SIAB_INPUT_20210101000001.json"])
        self.assertTrue(os.path.exists(fscript))
        #os.remove(fscript)

if __name__ == "__main__":
    unittest.main()