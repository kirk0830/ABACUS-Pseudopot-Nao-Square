"""initialize is the module should be called everytime will the program starts"""

import apns.module_io.input_translate as amiit
import apns.module_io.abacustest as amia
import apns.module_structure.basic as amsb
# import apns.module_nao.nao_archive as amna
def initialize(finp: str):
    
    inp = read_apns_inp(finp)
    # step1: get structures and elements, elements are for searching pseudopotentials and naos.
    structures = search_structure(systems=inp["systems"],
                                  n_structures=inp["materials_project"]["n_structures"],
                                  nspin=inp["calculation"]["nspin"],
                                  magnetism=inp["extensive"]["magnetism"],
                                  api_key=inp["materials_project"]["api_key"])
    # step2: get pseudopotentials and naos by tag search
    # return available pseudopotentials and numerical atomic orbitals as dict, the element symbol
    # will be the key of the dict.
    elements = amsb.scan_elements(structures.keys())
    upfs, orbs = search_upforb(elements=elements,
                               pptags=inp["pptags"],
                               naotags=inp["naotags"],
                               pseudo_dir=inp["global"]["pseudo_dir"],
                               orbital_dir=inp["global"]["orbital_dir"])
    abacustest = amia.read_apns_inp(finp)
    return inp, structures, upfs, orbs, abacustest

import json
import os
import apns.module_structure.manage as amsm
def search_structure(systems: dict|list,
                     n_structures: int,
                     nspin: int,
                     magnetism: str,
                     api_key: str):
    """according to the input file, search structures from cache or Materials Project.
    returns a dictionary with formula as key and a list of structures as value.
    along with a list of elements."""
    # the case systems is defined as a list, means auto-generate structure no matter whether download from Materials
    # project or isolated system.
    if isinstance(systems, list):
        result = {}
        # isolated is reserved for performing numerical atomic orbital test (development motivation)
        isolated = [s for s in systems if s.endswith("_dimer") or s.endswith("_trimer") or s.endswith("_tetramer")]
        for system in isolated:
            element, geometry = system.split("_")
            result.setdefault(element, []).append(("::AUTOGEN::" + geometry.upper(), None))
        # bravis lattice is reserved for performing EOS test, taking advantage of ACWF all-electron calculation result as reference data
        bravis = [s for s in systems if s.endswith("_sc") or s.endswith("_bcc") or s.endswith("_fcc") or s.endswith("_diamond")]
        for system in bravis:
            element, geometry = system.split("_")
            result.setdefault(element, []).append(("::AUTOGEN::" + geometry.upper(), None))
        # for practical systems
        practical = [s for s in systems if s not in isolated and s not in bravis]
        # consider magnetism or not
        consider_magnetism = True if nspin == 2 else False
        consider_magnetism = True if consider_magnetism and magnetism == "materials_project" else False
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
            result.update(amsm.search(api_key=api_key, # this is left here to raise error if api_key not specified
                                      formula=formula_tosearch,
                                      n_structures=[n_structures] * len(formula_tosearch)))
        return result
    elif isinstance(systems, dict):
        # copy this section and modify on it
        result = {formula: [] for formula in systems.keys()}
        assert id(result) != id(systems)
        # user-self-defined case, nspin = 2 cannot be automatically supported.
        for formula in systems.keys():
            # seperate isolated and crystal structures, will give None as magmom
            isolated = [s for s in systems[formula] if s.endswith("_dimer") or s.endswith("_trimer") or s.endswith("_tetramer")]
            for system in isolated:
                element, geometry = system.split("_")
                result[formula].append(("::AUTOGEN::" + geometry.upper(), None))
            # for practical systems
            crystal = [s for s in systems[formula] if s not in isolated]
            # check existence of these cif files
            for system in crystal:
                if not os.path.exists(system):
                    raise FileNotFoundError("Cif file not found: ", system)
            result[formula].extend([(system, None) for system in crystal])
        return result
    else:
        raise ValueError("Unknown format of systems: \n", systems)

import apns.module_new.tag_search as amds
def search_upforb(elements: list, pptags: list, naotags: list, pseudo_dir: str, orbital_dir: str):
    """search valid pseudopotential and numerical atomic orbitals for all elements in input file.
    update to 2024-05-01, the logic of searching/indexing pseudopotentials are again iterated,
    now because there are too many pseudopotentials, the search is based on tags, rather than
    three domains that needed user to understand.
    Now with specifying tags (or filters), can user get available upf files. It will also be the
    strategy of numerical atomic orbital."""
    # first to make unique the tags
    pptags = list(set(pptags))
    # because nao strongly bind with pseudopotentials, therefore
    # nao tags should also include pseudopotential tags.
    # But user may explicitly repeat some tags if really wants.
    naotags = list(set(pptags + naotags))
    
    # first search pseudopotentials -> can be futherly encapsulated
    fdb = os.path.join(pseudo_dir, "database.json")
    ppsearcher = amds.TagSearcher(fdb)
    ppsearch = {element: ppsearcher(True, False, element, *pptags) for element in elements}
    # then search numerical atomic orbitals
    fdb = os.path.join(orbital_dir, "database.json")
    naosearcher = amds.TagSearcher(fdb)
    naosearch = {element: naosearcher(True, False, element, *naotags) for element in elements}
    return ppsearch, naosearch

def read_apns_inp(finp: str):
    result = {}
    with open(finp, "r") as f:
        inp = json.load(f)
    section = inp.get("global", {})
    ideal = {
        "test_mode": section.get("test_mode", "pseudopotential"),
        "software": section.get("software", "abacus"),
        "work_dir": section.get("work_dir", os.path.abspath("./")),
        "pseudo_dir": section.get("pseudo_dir", os.path.abspath("./download/pseudopotentials")),
        "orbital_dir": section.get("orbital_dir", os.path.abspath("./download/numerical_orbitals")),
        "save_log": section.get("save_log", True)
    }
    result["global"] = ideal
    section = inp.get("calculation", {})
    ideal = {
        "basis_type": section.get("basis_type", "pw"),
        "ecutwfc": section.get("ecutwfc", 100),
        "nspin": section.get("nspin", 1),
        "symmetry": section.get("symmetry", 0)
    }
    result["calculation"] = ideal
    section = inp.get("extensive", {})
    ideal = {
        "characteristic_lengths": section.get("characteristic_lengths", [0.0]),
        "magnetism": section.get("magnetism", "materials_project"),
        "nkpoints_in_line": section.get("nkpoints_in_line", 0),
    }
    result["extensive"] = ideal
    section = inp.get("systems", [])
    ideal = section
    result["systems"] = ideal
    section = inp.get("materials_project", {})
    ideal = {
        "api_key": section.get("api_key", ""),
        "n_structures": section.get("n_structures", 1),
        "theoretical": section.get("theoretical", False),
        "most_stable": section.get("most_stable", True)
    }
    result["materials_project"] = ideal
    section = inp.get("pptags", [])
    ideal = section
    result["pptags"] = ideal
    section = inp.get("naotags", [])
    ideal = section
    result["naotags"] = ideal

    return result

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
        self.maxDiff = None
        inp = {
            "systems": {
                "Mn": ["./apns_cache/mp-35.cif", "Mn_dimer"],
                "V": ["V_trimer", "./apns_cache/mp-146.cif"]
            }
        }
        with open("test.json", "w") as f:
            json.dump(inp, f)
        result = search_structure("test.json")
        ref = {
            "Mn": [("::AUTOGEN::DIMER", None), ("./apns_cache/mp-35.cif", None)],
            "V": [("::AUTOGEN::TRIMER", None), ("./apns_cache/mp-146.cif", None)]
        }
        self.assertEqual(result, ref)
        os.remove("test.json")

        inp = {
            "systems": ["Fe_bcc", "Fe_diamond", "Ni_fcc"]
        }
        with open("test.json", "w") as f:
            json.dump(inp, f)
        result = search_structure("test.json")
        ref = {
            "Fe": [("::AUTOGEN::BCC", None), ("::AUTOGEN::DIAMOND", None)],
            "Ni": [("::AUTOGEN::FCC", None)]
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