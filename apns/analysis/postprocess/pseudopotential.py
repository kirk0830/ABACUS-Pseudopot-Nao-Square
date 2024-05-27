import json
import apns.pspot.manage as ampm
def testname_pspotid(element: str, testname: str, 
                     pseudo_dir = "download/pseudopotentials/", 
                     fpseudo_db="pseudo_db.json"):
    """return pspotid from testname
    
    testname is the pspot_id.replace("_", "")
    """
    def translate_kind(testname):
        """translate pspotid to real pseudopotential name"""
        dictionary = {"dojo": "PseudoDojo", "pslnc": "PSlibrary (NC)", 
                      "pslrrkjus": "PSlibrary (RRKJUS)", "pslkjpaw": "PSlibrary (PAW)",
                      "gth": "Goedecker-Teter-Hutter", "hgh": "Hartwigsen-Goedecker-Hutter"}
        return dictionary[testname] if testname in dictionary.keys() else testname.upper()
    def translate_version(kind, version):
        if kind.upper() == "GTH":
            return " "+version
        elif kind.upper() == "PD":
            return version
        else:
            #version = ".".join([str(version[i]) for i in range(len(version))])
            return "-v" + version if version != "" else ""
    # TEMPORARY FIX
    # this is because presently only GTH LnPP1 are collected, while there are indeed other versions
    # of GTH pseudopotentials, so we need to fix the name temporarily. Other versions like UZH, LnPP2,
    # ... will be added in the future.
    testname = "gthLnPP1" if testname == "gth" else testname
    with open(pseudo_dir + fpseudo_db, "r") as f:
        pseudo_db = json.load(f)
    for key in pseudo_db[element].keys():
        if key.replace("_", "").replace(".", "") == testname:
            fpseudo_withpath = pseudo_db[element][key]
            attribute = ampm.get_attribute(pseudo_dir=pseudo_dir, fpseudo_withpath=fpseudo_withpath)
            kind, version, appendix = attribute["kind"], attribute["version"], attribute["appendix"]
            kind = translate_kind(kind)
            version = translate_version(kind, version)
            label = kind + version
            label += " (" + appendix + ")" if appendix != "" else ""
            label = label.strip()
            return label, key
    return None, None

import unittest
import json
class TestPseudopotential(unittest.TestCase):
    def test_testname_pspotid(self):
        fpseudo_db = "pseudo_db.json"
        pseudo_dir = "download/pseudopotentials/"
        self.assertEqual(testname_pspotid("O", "pd04", pseudo_dir=pseudo_dir, fpseudo_db=fpseudo_db), 
                         ("PD04", "pd_04"))
        self.assertEqual(testname_pspotid("O", "pd04high", pseudo_dir=pseudo_dir, fpseudo_db=fpseudo_db), 
                         ("PD04 (high)", "pd_04_high"))
        self.assertEqual(testname_pspotid("O", "pslkjpaw031", pseudo_dir=pseudo_dir, fpseudo_db=fpseudo_db), 
                         ("PSlibrary (PAW)-v0.3.1", "pslkjpaw_0.3.1"))
        self.assertEqual(testname_pspotid("O", "pslrrkjus031", pseudo_dir=pseudo_dir, fpseudo_db=fpseudo_db), 
                         ("PSlibrary (RRKJUS)-v0.3.1", "pslrrkjus_0.3.1"))
        self.assertEqual(testname_pspotid("O", "pslnc031", pseudo_dir=pseudo_dir, fpseudo_db=fpseudo_db), 
                         ("PSlibrary (NC)-v0.3.1", "pslnc_0.3.1"))
        self.assertEqual(testname_pspotid("O", "gbrv15", pseudo_dir=pseudo_dir, fpseudo_db=fpseudo_db), 
                         ("GBRV-v1.5", "gbrv_1.5"))
        self.assertEqual(testname_pspotid("O", "dojo05sr", pseudo_dir=pseudo_dir, fpseudo_db=fpseudo_db), 
                         ("PseudoDojo-v0.5 (sr)", "dojo_0.5_sr"))
        self.assertEqual(testname_pspotid("O", "dojo04fr", pseudo_dir=pseudo_dir, fpseudo_db=fpseudo_db), 
                         ("PseudoDojo-v0.4 (fr)", "dojo_0.4_fr"))
        self.assertEqual(testname_pspotid("O", "dojo03sr", pseudo_dir=pseudo_dir, fpseudo_db=fpseudo_db), 
                         ("PseudoDojo-v0.3 (sr)", "dojo_0.3_sr"))
        self.assertEqual(testname_pspotid("O", "pd03", pseudo_dir=pseudo_dir, fpseudo_db=fpseudo_db), 
                         ("PD03", "pd_03"))
        self.assertEqual(testname_pspotid("O", "sg151fr", pseudo_dir=pseudo_dir, fpseudo_db=fpseudo_db), 
                         (None, None)) # there is no such pseudopotential called sg151fr
        
if __name__ == "__main__":
    unittest.main()