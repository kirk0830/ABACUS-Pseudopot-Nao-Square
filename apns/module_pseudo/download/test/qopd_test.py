import unittest
import apns.module_pseudo.download.qespresso_official_pptable_download as ampdqopd
import re
class TestQespressoOfficialPptableDownload(unittest.TestCase):

    ftests = [
        "H.pbe-kjpaw_psl.1.0.0.UPF",
        "H.pbe-rrkjus_psl.1.0.0.UPF",
        "H.pz-kjpaw_psl.1.0.0.UPF",
        "H.pz-rrkjus_psl.1.0.0.UPF",
        "H.rel-pbe-kjpaw_psl.1.0.0.UPF",
        "H.rel-pbe-rrkjus_psl.1.0.0.UPF",
        "H.rel-pz-kjpaw_psl.1.0.0.UPF",
        "H.rel-pz-rrkjus_psl.1.0.0.UPF",
        "H.pbe-kjpaw_psl.0.1.UPF",
        "H.pbe-rrkjus_psl.0.1.UPF",
        "H.pbesol-kjpaw_psl.0.1.UPF",
        "H.pbesol-rrkjus_psl.0.1.UPF",
        "H.pz-kjpaw_psl.0.1.UPF",
        "H.pz-rrkjus_psl.0.1.UPF",
        "H.rel-pbe-kjpaw_psl.0.1.UPF",
        "H.rel-pbe-rrkjus_psl.0.1.UPF",
        "H.rel-pbesol-kjpaw_psl.0.1.UPF",
        "H.rel-pbesol-rrkjus_psl.0.1.UPF",
        "H.rel-pz-kjpaw_psl.0.1.UPF",
        "H.rel-pz-rrkjus_psl.0.1.UPF",
        "H.pbe-hgh.UPF"
    ]
    def test_re_pattern1(self):
        pattern = r"([A-Z][a-z]?)(\.)(rel\-)?([\w]+)(\-)(.*)(.UPF)$"
        for ftest in self.ftests:
            print(ftest)
            self.assertIsNotNone(re.match(pattern, ftest))
            print(re.match(pattern, ftest).groups())

if __name__ == "__main__":
    unittest.main()