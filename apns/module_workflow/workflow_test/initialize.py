"""initialize is the module should be called everytime will the program starts"""

import apns.module_io.input_translate as amiit
import apns.module_io.abacustest as amia
import apns.module_structure.basic as amsb
# import apns.module_nao.nao_archive as amna
def initialize(finp: str):

    structures, elements = search_structure(finp)
    runtime_settings = amiit.read_apns_inp(finp)
    valid_upfs, valid_orbs = scan_upforb(elements=elements, basis_type=runtime_settings["calculation"]["basis_type"],
                                         pseudo_dir=runtime_settings["global"]["pseudo_dir"],
                                         orbital_dir=runtime_settings["global"]["orbital_dir"],
                                         pspot_settings=runtime_settings["pseudopotentials"],
                                         nao_settings=runtime_settings["numerical_orbitals"])
    abacustest = amia.read_apns_inp(finp)
    return runtime_settings, structures, valid_upfs, valid_orbs, abacustest

import json
import os
import apns.module_structure.manage as amsm
def search_structure(finp: str):

    with open(finp, "r") as f:
        inp = json.load(f)
    # the case systems is defined as a list, means auto-generate structure no matter whether download from Materials
    # project or isolated system.
    if isinstance(inp["systems"], list):
        result = {}
        # isolated is reserved for performing numerical atomic orbital test (development motivation)
        isolated = [s for s in inp["systems"] if s.endswith("_dimer") or s.endswith("_trimer") or s.endswith("_tetramer")]
        for system in isolated:
            element, geometry = system.split("_")
            result.setdefault(element, []).append(("::AUTOGEN::" + geometry.upper(), None))
        # bravis lattice is reserved for performing EOS test, taking advantage of ACWF all-electron calculation result as reference data
        bravis = [s for s in inp["systems"] if s.endswith("_sc") or s.endswith("_bcc") or s.endswith("_fcc") or s.endswith("_diamond")]
        for system in bravis:
            element, geometry = system.split("_")
            result.setdefault(element, []).append(("::AUTOGEN::" + geometry.upper(), None))
        # for practical systems
        practical = [s for s in inp["systems"] if s not in isolated and s not in bravis]
        # number of structures needed to search from Materials Project for each formula
        n_structures = inp["materials_project"]["n_structures"]
        # consider magnetism or not
        consider_magnetism = True if inp["calculation"]["nspin"] == 2 else False
        consider_magnetism = True if consider_magnetism and inp["extensive"]["magnetism"] == "materials_project" else False
        formula_tosearch = [] # always download magnetism information, this would not be harmful
        for system in practical:
            record = amsm.lookup(system, n_structures, with_magmom=consider_magnetism)
            if len(record) < n_structures:
                print("Structure not found in cache, will search from Materials Project: ", system)
                formula_tosearch.append(system)
            else:
                print("Structure found in cache: ", system, ", skip searching from Materials Project.")
                result[system] = record
        if formula_tosearch != []:
            # search from Materials Project
            result.update(amsm.search(api_key=inp["materials_project"]["api_key"],
                                      formula=formula_tosearch,
                                      n_structures=[n_structures] * len(formula_tosearch)))
        return result, amsb.scan_elements(list(result.keys()))
    elif isinstance(inp["systems"], dict):
        # copy this section and modify on it
        result = {formula: [] for formula in inp["systems"].keys()}
        assert id(result) != id(inp["systems"])
        # user-self-defined case, nspin = 2 cannot be automatically supported.
        for formula in inp["systems"].keys():
            # seperate isolated and crystal structures, will give None as magmom
            isolated = [s for s in inp["systems"][formula] if s.endswith("_dimer") or s.endswith("_trimer") or s.endswith("_tetramer")]
            for system in isolated:
                element, geometry = system.split("_")
                result[formula].append(("::AUTOGEN::" + geometry.upper(), None))
            # for practical systems
            crystal = [s for s in inp["systems"][formula] if s not in isolated]
            # check existence of these cif files
            for system in crystal:
                if not os.path.exists(system):
                    raise FileNotFoundError("Cif file not found: ", system)
            result[formula].extend([(system, None) for system in crystal])
        return result, amsb.scan_elements(list(result.keys()))
    else:
        raise ValueError("Unknown format of systems: \n", inp["systems"])

import apns.module_pseudo.manage as ampm
import apns.module_nao.manage as amnm
def scan_upforb(elements: list,
                basis_type: str,
                pseudo_dir: str, 
                orbital_dir: str,
                pspot_settings: dict,
                nao_settings: dict) -> tuple[dict, dict]:
    """scan valid pseudopotential for all elements in input file.
    If not basis_type is "lcao", then return valid_orbs as empty dict."""
    valid_upfs = {element: {} for element in elements}
    valid_upfs = ampm.valid_pseudo(pseudo_dir=pseudo_dir, 
                                   elements=elements, 
                                   pseudo_setting=pspot_settings)
    valid_orbs = {element: {} for element in elements} 
    valid_orbs = amnm.valid_nao(orbital_dir=orbital_dir, 
                                elements=elements, 
                                nao_settings=nao_settings)\
                 if basis_type == "lcao" else valid_orbs
    return valid_upfs, valid_orbs

import unittest
class TestInitializeWokflowTest(unittest.TestCase):
    def test_search_structure(self):
        inp = {
            "systems": {
                "Mn": ["./apns_cache/mp-35.cif", "Mn_dimer"],
                "Pu": ["Pu_trimer", "./apns_cache/mp-582819.cif"]
            }
        }
        with open("test.json", "w") as f:
            json.dump(inp, f)
        result = search_structure("test.json")
        ref = {
            "Mn": [("::AUTOGEN::DIMER", None), ("./apns_cache/mp-35.cif", None)],
            "Pu": [("::AUTOGEN::TRIMER", None), ("./apns_cache/mp-582819.cif", None)]
        }
        self.assertEqual(result, ref)
        os.remove("test.json")

        inp = {
            "systems": {
                "Mn": ["cif_file_in_imagination", "Mn_dimer"],
                "Pu": ["Pu_trimer", "./apns_cache/mp-582819.cif"]
            }
        }
        with open("test.json", "w") as f:
            json.dump(inp, f)
        # catch FileNotFoundError
        with self.assertRaises(FileNotFoundError):
            search_structure("test.json")
        os.remove("test.json")

if __name__ == "__main__":
    print(f"Running test on {__file__}....")
    unittest.main()