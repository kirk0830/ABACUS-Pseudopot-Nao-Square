import json
import apns.module_pseudo.manage as ampm
def testname_pspotid(element: str, testname: str):
    """return pspotid from testname
    
    testname is the pspot_id.replace("_", "")
    """
    def translate_kind(testname):
        """translate pspotid to real pseudopotential name"""
        dictionary = {"dojo": "PseudoDojo", "pslnc": "PSlibrary (NC)", 
                      "pslrrkjus": "PSlibrary (RRKJUS)", "pslncpaw": "PSlibrary (PAW)",
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
    with open("download/pseudopotentials/pseudo_db.json", "r") as f:
        pseudo_db = json.load(f)
    for key in pseudo_db[element].keys():
        if key.replace("_", "").replace(".", "") == testname:
            fpseudo_withpath = pseudo_db[element][key]
            attribute = ampm.get_attribute("download/pseudopotentials", fpseudo_withpath=fpseudo_withpath)
            kind, version, appendix = attribute["kind"], attribute["version"], attribute["appendix"]
            kind = translate_kind(kind)
            version = translate_version(kind, version)
            label = kind + version
            label += " (" + appendix + ")" if appendix != "" else ""
            label = label.strip()
            return label, key
    return None, None