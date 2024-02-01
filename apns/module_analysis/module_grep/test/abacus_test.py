import unittest
import apns.module_analysis.module_grep.abacus as amaa
import os
import re
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
        self.assertTrue(amaa.is_normal_end(fname))
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf.log"
        self.assertTrue(amaa.is_normal_end(fname))
    
    def test_grep_energy(self):
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf_notconverge.log"
        self.assertFalse(amaa.grep_energy(fname))
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf.log"
        self.assertGreater(abs(amaa.grep_energy(fname) - (-0.000000000000000)), 1e-15)

    def test_grep_pressure(self):
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf_notconverge.log"
        self.assertFalse(amaa.grep_pressure(fname))
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf.log"
        self.assertGreater(abs(amaa.grep_pressure(fname) - (-0.000000000000000)), 1e-15)

    def test_grep_natom(self):
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf_notconverge.log"
        self.assertEqual(amaa.grep_natom(fname), 2)
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf.log"
        self.assertEqual(amaa.grep_natom(fname), 1)

    def test_grep_band(self):
        float_pattern = r"([0-9\-]+\.[0-9\-]+[eEdD][\+\-][0-9]+)|([0-9\-]+\.[0-9\-]+)"
        
        kpt_pattern = r"^(\s*)([0-9]+)(\s*)(%s)(\s+)(%s)(\s*)(%s)(\s*)(%s)(\s*)$"%(float_pattern, float_pattern, float_pattern, float_pattern)
        band_pattern = r"^(\s*)([0-9]+)(\s*)(%s)(\s+)(%s)(\s*)$"%(float_pattern, float_pattern)
        band_header_pattern = r"^(.*)(\=\s*)(%s)(\s+)(%s)(\s+)(%s)(.*)$"%(float_pattern, float_pattern, float_pattern)
        self.assertTrue(re.match(kpt_pattern, "       1   0.00000000   0.00000000   0.00000000 0.0058"))
        #print(re.match(kpt_pattern, "       1   0.00000000   0.00000000   0.00000000 0.0058").groups())
        self.assertTrue(re.match(band_pattern, "       1       -18.8540     0.00583090"))
        #print(re.match(band_pattern, "       1       -18.8540     0.00583090").groups())
        self.assertTrue(re.match(band_header_pattern, " 1/20 kpoint (Cartesian) = 0.0000 0.0000 0.0000 (44 pws)"))
        #print(re.match(band_header_pattern, " 1/20 kpoint (Cartesian) = 0.0000 0.0000 0.0000 (44 pws)").groups())
        fname = os.getcwd() + "/apns/module_analysis/module_grep/test/support/running_scf.log"
        result = amaa.grep_band(fname)
        bands = result[0]
        nelec = result[1]
        nbands = result[2]
        efermi = result[3]
        self.assertEqual(len(bands), 280)
        self.assertEqual(nelec, 8)
        self.assertEqual(nbands, 14)
        self.assertAlmostEqual(efermi, 1.6277681258, 10)
        print("bands = ", bands)
        print("nelec = ", nelec)
        print("nbands = ", nbands)
        print("efermi = ", efermi)

if __name__ == "__main__":
    unittest.main()