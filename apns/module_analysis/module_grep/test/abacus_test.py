import unittest
import apns.module_analysis.module_grep.abacus as amaa
import os
class TestAbacus(unittest.TestCase):

    def test_fname_setting(self):
        self.assertEqual(amaa.fname_setting(job_dir="./",
                                            calculation="relax",
                                            suffix="unittest",
                                            include_path=True), 
                        ("./out.log", "./OUT.unittest/running_relax.log", "./INPUT"))
        self.assertEqual(amaa.fname_setting(job_dir="./",
                                            calculation="relax",
                                            suffix="unittest",
                                            include_path=False), 
                        ("out.log", "running_relax.log", "INPUT"))
        
    def test_is_normal_end(self):
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf_notconverge.log"
        self.assertFalse(amaa.is_normal_end(fname))
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf_converge.log"
        self.assertTrue(amaa.is_normal_end(fname))
    
    def test_grep_energy(self):
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf_notconverge.log"
        self.assertFalse(amaa.grep_energy(fname))
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf_converge.log"
        self.assertGreater(abs(amaa.grep_energy(fname) - (-0.000000000000000)), 1e-15)

    def test_grep_pressure(self):
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf_notconverge.log"
        self.assertFalse(amaa.grep_pressure(fname))
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf_converge.log"
        self.assertGreater(abs(amaa.grep_pressure(fname) - (-0.000000000000000)), 1e-15)

    def test_grep_natom(self):
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf_notconverge.log"
        self.assertEqual(amaa.grep_natom(fname), 6)
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf_converge.log"
        self.assertEqual(amaa.grep_natom(fname), 3)

if __name__ == "__main__":
    unittest.main()