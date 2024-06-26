import os
import json
import re

import apns.outdated.identifier as amwi
def load(pseudo_dir: str):
    print("""
rules.json will be loaded on at beginning of test workflow. Once
new pseudopotentials are added, please update rules.json accordingly.
Contents in rules are organized as:
{
    "rules": [
        {
            "kind": "kind",
            "version": "version",
            "appendix": "appendix",
            "re.folder": "regular expression for folder",
            "re.file": "regular expression for file"
        },
        ...
    ]
}
kind: kind of pseudopotential, e.g. sg15, dojo, gbrv, hgh, pd, psl, ...
version: version of pseudopotential, e.g. 1.0, 1.1, 1.2, 1.3, ...
appendix: appendix of pseudopotential, e.g. fr, sr, spd, spd-high, ...
re.folder: regular expression for folder, e.g. "sg15_oncv_upf_2020-02-06"
           will be matched with "sg15_oncv_upf_2020-02-06"
re.file: regular expression for file, e.g. "Si_ONCV_PBE-1.0.upf"
         will be matched with "^([A-Z][a-z]?)(_ONCV_PBE\\-1\\.0\\.upf)$"
pseudo_db.json will be refreshed after loading rules.json.
""")
    if not os.path.exists(pseudo_dir):
        raise FileNotFoundError(f"Directory {pseudo_dir} does not exist")
    if not os.path.exists(os.path.join(pseudo_dir, "rules.json")):
        raise FileNotFoundError(f"File rules.json does not exist in {pseudo_dir}")
    
    with open(os.path.join(pseudo_dir, "rules.json")) as f:
        rules = json.load(f)

    pseudo_db = {}
    for root, dirs, files in os.walk(pseudo_dir):
        for file in files:
            if file.endswith(".upf") or file.endswith(".UPF"):
                folder = root.replace("/", "/").split("/")[-1]
                for i in range(len(rules["rules"])):
                    match_file = re.match(rules["rules"][i]["re.file"], file)
                    match_folder = re.match(rules["rules"][i]["re.folder"], folder)
                    if match_file and match_folder:
                        element = match_file.group(1).capitalize()
                        id = amwi.pseudopotential(kind=rules["rules"][i]["kind"], 
                                                  version=rules["rules"][i]["version"], 
                                                  appendix=rules["rules"][i]["appendix"])
                        pseudo_db.setdefault(element, {})[id] = os.path.join(root, file)
                        break
                if match_file is None or match_folder is None:
                    raise ValueError(f"No rule found for file: {file} in folder: {folder}")
                
    with open(os.path.join(pseudo_dir, "pseudo_db.json"), "w") as f:
        json.dump(pseudo_db, f, indent=4)

    return pseudo_db

def valid_pseudo(pseudo_dir: str, elements: list, pseudo_setting: dict):
    """from `pseudo_db.json`, get available pseudopotentials for each element.

    Args:
        pseudo_dir (str): directory where pseudopotentials are stored
        elements (list): list of elements
        pseudo_setting (dict): setting for pseudopotentials, including kinds, versions, appendices
    """
    with open(os.path.join(pseudo_dir, "pseudo_db.json")) as f:
        pseudo_db = json.load(f)
    result = {element: {} for element in elements}
    for e in elements:
        if pseudo_setting["kinds"][e] == ["all"]:
            result[e] = pseudo_db[e]
        else:
            for kind in pseudo_setting["kinds"][e]:
                if pseudo_setting["versions"][e] == ["all"]:
                    for id in pseudo_db[e]:
                        if id.startswith(f"{kind}_"):
                            result[e][id] = pseudo_db[e][id]
                else:
                    for version in pseudo_setting["versions"][e]:
                        if pseudo_setting["appendices"][e] == ["all"]:
                            for id in pseudo_db[e]:
                                if id.startswith(f"{kind}_{version}"):
                                    result[e][id] = pseudo_db[e][id]
                        else:
                            for appendix in pseudo_setting["appendices"][e]:
                                id = amwi.pseudopotential(kind=kind, version=version, appendix=appendix)
                                if id in pseudo_db[e]:
                                    result[e][id] = pseudo_db[e][id]
        if len(result[e]) == 0:
            raise ValueError(f"No valid pseudopotential for element {e}.")

    return result

def get_attribute(pseudo_dir: str, **kwargs):
    """
    method1: provide fpseudo, folder
    method2: provide fpseudo_withpath
    
    return: {"kind": kind, "version": version, "appendix": appendix}
    """
    with open(os.path.join(pseudo_dir, "rules.json")) as f:
        rules = json.load(f)
    
    fpseudo = kwargs.get("fpseudo", None)
    fpseudo_withpath = kwargs.get("fpseudo_withpath", None)
    folder = kwargs.get("folder", None)
    if fpseudo is None and fpseudo_withpath is None and folder is None:
        raise ValueError("At least one of fpseudo, fpseudo_withpath, folder should be provided.")
    if fpseudo is not None and folder is not None:
        for i in range(len(rules["rules"])):
            match_file = re.match(rules["rules"][i]["re.file"], fpseudo)
            match_folder = re.match(rules["rules"][i]["re.folder"], folder)
            if match_file and match_folder:
                return {"kind": rules["rules"][i]["kind"], 
                        "version": rules["rules"][i]["version"], 
                        "appendix": rules["rules"][i]["appendix"]}
    elif fpseudo_withpath is not None:
        folder = fpseudo_withpath.replace("/", "/").split("/")[-2]
        fpseudo = fpseudo_withpath.replace("/", "/").split("/")[-1]
        for i in range(len(rules["rules"])):
            match_file = re.match(rules["rules"][i]["re.file"], fpseudo)
            match_folder = re.match(rules["rules"][i]["re.folder"], folder)
            if match_file and match_folder:
                return {"kind": rules["rules"][i]["kind"], 
                        "version": rules["rules"][i]["version"], 
                        "appendix": rules["rules"][i]["appendix"]}
    else:
        raise ValueError("Invalid combination of fpseudo, fpseudo_withpath, folder.")

def export(pseudo_dir: str, out_dir: str = "./", included: list = None, excluded: list = None, pseudo_regex: str = None):
    """export pseudopotentials for included elements"""
    with open(os.path.join(pseudo_dir, "pseudo_db.json")) as f:
        pseudo_db = json.load(f)
    included = included if included else pseudo_db.keys()
    excluded = excluded if excluded else []
    pseudo_regex = pseudo_regex if pseudo_regex else ".*"
    description = {}
    for element in included:
        if element in excluded:
            continue
        os.makedirs(os.path.join(out_dir, element), exist_ok=True)
        description[element] = {}
        # save description for each pseudopotential, key is new name
        # values are, original name, kind, version, appendix
        for abbr in pseudo_db[element].keys():
            if re.match(pseudo_regex, abbr):
                fpseudo = pseudo_db[element][abbr]
                words = abbr.split("_")
                words = [words[0], words[1], "_".join(words[2:])] if len(words) > 2 else words
                words = words + [""] * (3 - len(words))
                kind, version, appendix = words
                new_name = f"{element}_{kind}_{version}_{appendix}.UPF"
                description[element][new_name] = {"original": abbr, "kind": kind, "version": version, "appendix": appendix}
                os.system(f"cp {fpseudo} {os.path.join(out_dir, element, new_name)}")
    with open(os.path.join(out_dir, "description.json"), "w") as f:
        json.dump(description, f, indent=4)
    
import unittest
class TestPseudo(unittest.TestCase):
    def test_load(self):
        pseudo_dir = "./download/pseudopotentials/"
        pseudo_db = load(pseudo_dir)
        self.assertTrue(os.path.exists(os.path.join(pseudo_dir, "pseudo_db.json")))
        self.assertTrue(len(pseudo_db) > 0)
    
    def test_valid_pseudo(self):
        pseudo_dir = "./download/pseudopotentials/"
        elements = ["Si", "Ge"]
        pseudo_setting = {
            "kinds": ["all"],
            "versions": [""],
            "appendices": [""]
        }
        pseudo_setting = {key: {element: pseudo_setting[key] for element in elements} for key in pseudo_setting}
        result = valid_pseudo(pseudo_dir, elements, pseudo_setting)
        self.assertTrue(len(result) > 0)
        self.assertTrue(len(result["Si"]) > 0)
        self.assertTrue(len(result["Ge"]) > 0)
        Si_ref = {
            "gbrv_1.5": "./download/pseudopotentials/GBRV_pbe_UPF_v1.5/si_pbe_v1.uspp.F.UPF",
            "hgh_1.0": "./download/pseudopotentials/hgh/Si.pbe-hgh.UPF",
            "dojo_0.4_fr": "./download/pseudopotentials/nc-fr-04_pbe_standard/Si.upf",
            "dojo_0.4_sr": "./download/pseudopotentials/nc-sr-04_pbe_standard_upf/Si.upf",
            "dojo_0.5_sr": "./download/pseudopotentials/nc-sr-05_pbe_standard_upf/Si.upf",
            "pd_03": "./download/pseudopotentials/NCPP-PD03-PBE/Si.PD03.PBE.UPF",
            "pd_04_sp": "./download/pseudopotentials/NCPP-PD04-PBE/Si-sp.PD04.PBE.UPF",
            "pd_04": "./download/pseudopotentials/NCPP-PD04-PBE/Si.PD04.PBE.UPF",
            "dojo_0.3_sr": "./download/pseudopotentials/pbe_s_sr/Si.upf",
            "pslkjpaw_0.3.1": "./download/pseudopotentials/pslibrary-pbe.0.3.1/PSEUDOPOTENTIALS/Si.pbe-n-kjpaw_psl.0.1.UPF",
            "pslrrkjus_0.3.1": "./download/pseudopotentials/pslibrary-pbe.0.3.1/PSEUDOPOTENTIALS/Si.pbe-n-rrkjus_psl.0.1.UPF",
            "pslnc_0.3.1": "./download/pseudopotentials/pslibrary-pbe.0.3.1/PSEUDOPOTENTIALS_NC/Si.pbe-n-nc.UPF",
            "sg15_1.0": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.0.upf",
            "sg15_1.1": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.1.upf",
            "sg15_1.2": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.2.upf",
            "sg15_1.1_fr": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE_FR-1.1.upf"
        }
        for k, v in Si_ref.items():
            self.assertTrue(k in result["Si"])
            self.assertEqual(result["Si"][k], v)
        Ge_ref = {
            "gbrv_1.5": "./download/pseudopotentials/GBRV_pbe_UPF_v1.5/ge_pbe_v1.4.uspp.F.UPF",
            "hgh_1.0": "./download/pseudopotentials/hgh/Ge.pbe-hgh.UPF",
            "dojo_0.4_fr": "./download/pseudopotentials/nc-fr-04_pbe_standard/Ge.upf",
            "dojo_0.4_sr": "./download/pseudopotentials/nc-sr-04_pbe_standard_upf/Ge.upf",
            "dojo_0.5_sr": "./download/pseudopotentials/nc-sr-05_pbe_standard_upf/Ge.upf",
            "pd_03": "./download/pseudopotentials/NCPP-PD03-PBE/Ge.PD03.PBE.UPF",
            "pd_04_d": "./download/pseudopotentials/NCPP-PD04-PBE/Ge-d.PD04.PBE.UPF",
            "pd_04_low": "./download/pseudopotentials/NCPP-PD04-PBE/Ge-low.PD04.PBE.UPF",
            "pd_04_spd-high": "./download/pseudopotentials/NCPP-PD04-PBE/Ge-spd-high.PD04.PBE.UPF",
            "dojo_0.3_sr": "./download/pseudopotentials/pbe_s_sr/Ge.upf",
            "pslkjpaw_0.3.1": "./download/pseudopotentials/pslibrary-pbe.0.3.1/PSEUDOPOTENTIALS/Ge.pbe-dn-kjpaw_psl.0.3.1.UPF",
            "pslrrkjus_0.3.1": "./download/pseudopotentials/pslibrary-pbe.0.3.1/PSEUDOPOTENTIALS/Ge.pbe-dn-rrkjus_psl.0.3.1.UPF",
            "pslnc_0.3.1": "./download/pseudopotentials/pslibrary-pbe.0.3.1/PSEUDOPOTENTIALS_NC/Ge.pbe-n-nc.UPF",
            "sg15_1.0": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Ge_ONCV_PBE-1.0.upf",
            "sg15_1.2": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Ge_ONCV_PBE-1.2.upf",
            "sg15_1.0_fr": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Ge_ONCV_PBE_FR-1.0.upf"
        }
        for k, v in Ge_ref.items():
            self.assertTrue(k in result["Ge"])
            self.assertEqual(result["Ge"][k], v)
        
        elements = ["Si", "Ge"]
        pseudo_setting = {
            "kinds": ["sg15"],
            "versions": ["all"],
            "appendices": [""]
        }
        pseudo_setting = {key: {element: pseudo_setting[key] for element in elements} for key in pseudo_setting}
        result = valid_pseudo(pseudo_dir, elements, pseudo_setting)
        self.assertTrue(len(result) > 0)
        self.assertTrue(len(result["Si"]) > 0)
        self.assertTrue(len(result["Ge"]) > 0)
        Si_ref = {
            "sg15_1.0": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.0.upf",
            "sg15_1.1": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.1.upf",
            "sg15_1.2": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.2.upf",
            "sg15_1.1_fr": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE_FR-1.1.upf"
        }
        for k, v in Si_ref.items():
            self.assertTrue(k in result["Si"])
            self.assertEqual(result["Si"][k], v)
        Ge_ref = {
            "sg15_1.0": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Ge_ONCV_PBE-1.0.upf",
            "sg15_1.2": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Ge_ONCV_PBE-1.2.upf",
            "sg15_1.0_fr": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Ge_ONCV_PBE_FR-1.0.upf"
        }
        for k, v in Ge_ref.items():
            self.assertTrue(k in result["Ge"])
            self.assertEqual(result["Ge"][k], v)

        elements = ["Si", "Ge"]
        pseudo_setting = {
            "kinds": ["sg15"],
            "versions": ["1.0"],
            "appendices": ["all"]
        }
        pseudo_setting = {key: {element: pseudo_setting[key] for element in elements} for key in pseudo_setting}
        result = valid_pseudo(pseudo_dir, elements, pseudo_setting)
        self.assertTrue(len(result) > 0)
        self.assertTrue(len(result["Si"]) > 0)
        self.assertTrue(len(result["Ge"]) > 0)
        Si_ref = {
            "sg15_1.0": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.0.upf"
        }
        for k, v in Si_ref.items():
            self.assertTrue(k in result["Si"])
            self.assertEqual(result["Si"][k], v)
        Ge_ref = {
            "sg15_1.0": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Ge_ONCV_PBE-1.0.upf",
            "sg15_1.0_fr": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Ge_ONCV_PBE_FR-1.0.upf"
        }
        for k, v in Ge_ref.items():
            self.assertTrue(k in result["Ge"])
            self.assertEqual(result["Ge"][k], v)

        elements = ["Si", "Ge"]
        pseudo_setting = {
            "kinds": ["sg15"],
            "versions": ["1.0"],
            "appendices": ["fr"]
        }
        pseudo_setting = {key: {element: pseudo_setting[key] for element in elements} for key in pseudo_setting}
        # Si will not have valid pseudopotential
        with self.assertRaises(ValueError):
            result = valid_pseudo(pseudo_dir, elements, pseudo_setting)
        
        elements = ["Si", "Ge"]
        pseudo_setting = {
            "kinds": ["sg15", "dojo"],
            "versions": ["all"],
            "appendices": [""]
        }
        pseudo_setting = {key: {element: pseudo_setting[key] for element in elements} for key in pseudo_setting}
        result = valid_pseudo(pseudo_dir, elements, pseudo_setting)
        self.assertTrue(len(result) > 0)
        self.assertTrue(len(result["Si"]) > 0)
        self.assertTrue(len(result["Ge"]) > 0)
        Si_ref = {
            "sg15_1.0": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.0.upf",
            "sg15_1.1": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.1.upf",
            "sg15_1.2": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.2.upf",
            "sg15_1.1_fr": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE_FR-1.1.upf",
            "dojo_0.4_fr": "./download/pseudopotentials/nc-fr-04_pbe_standard/Si.upf",
            "dojo_0.4_sr": "./download/pseudopotentials/nc-sr-04_pbe_standard_upf/Si.upf",
            "dojo_0.5_sr": "./download/pseudopotentials/nc-sr-05_pbe_standard_upf/Si.upf",
            "dojo_0.3_sr": "./download/pseudopotentials/pbe_s_sr/Si.upf"
        }
        for k, v in Si_ref.items():
            self.assertTrue(k in result["Si"])
            self.assertEqual(result["Si"][k], v)
        Ge_ref = {
            "sg15_1.0": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Ge_ONCV_PBE-1.0.upf",
            "sg15_1.2": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Ge_ONCV_PBE-1.2.upf",
            "sg15_1.0_fr": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Ge_ONCV_PBE_FR-1.0.upf",
            "dojo_0.4_fr": "./download/pseudopotentials/nc-fr-04_pbe_standard/Ge.upf",
            "dojo_0.4_sr": "./download/pseudopotentials/nc-sr-04_pbe_standard_upf/Ge.upf",
            "dojo_0.5_sr": "./download/pseudopotentials/nc-sr-05_pbe_standard_upf/Ge.upf",
            "dojo_0.3_sr": "./download/pseudopotentials/pbe_s_sr/Ge.upf"
        }
        for k, v in Ge_ref.items():
            self.assertTrue(k in result["Ge"])
            self.assertEqual(result["Ge"][k], v)

    def test_get_attribute(self):
        pseudo_dir = "./download/pseudopotentials/"
        fpseudo = "Si_ONCV_PBE-1.0.upf"
        folder = "sg15_oncv_upf_2020-02-06"
        result = get_attribute(pseudo_dir, fpseudo=fpseudo, folder=folder)
        self.assertEqual(result["kind"], "sg15")
        self.assertEqual(result["version"], "1.0")
        self.assertEqual(result["appendix"], "")
        
        fpseudo_withpath = "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.0.upf"
        result = get_attribute(pseudo_dir, fpseudo_withpath=fpseudo_withpath)
        self.assertEqual(result["kind"], "sg15")
        self.assertEqual(result["version"], "1.0")
        self.assertEqual(result["appendix"], "")
        
        with self.assertRaises(ValueError):
            result = get_attribute(pseudo_dir, fpseudo=fpseudo)
        
        with self.assertRaises(ValueError):
            result = get_attribute(pseudo_dir, folder=folder)
        
        with self.assertRaises(ValueError):
            result = get_attribute(pseudo_dir)

if __name__ == "__main__":
    unittest.main()