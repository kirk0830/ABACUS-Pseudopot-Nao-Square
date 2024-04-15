"""initialize is the module should be called everytime will the program starts"""

import apns.module_io.input_translate as amiit
import apns.module_structure.basic as amsb
# import apns.module_nao.nao_archive as amna
def initialize(finp: str):
    """initialize the program, return runtime information can be determined before running the workflow
    
    Args:
        finp (str): input file
        
    Raises:
    
    Returns:
        tuple[dict, dict, dict]: input, valid_pseudopotentials, valid_numerical_orbitals
        
    Details:
        input: the translated input json, following aspects are modified compared with user-input:  
               1. systems are changed to system with mpids  
               2. pseudopotentials and numerical_orbitals are expanded from list to dict  
               3. default values are set if not explicitly specified  
        valid_pseudopotentials: valid pseudopotentials for all elements in input file, the first layers keys 
                                are elements, the second layers keys are "identifiers" of pseudopotentials, 
                                the third layers keys are "kind", "version", "appendix" and "file", the first
                                three keys can concatenate as a "pseudopotential identifier".
        valid_numerical_orbitals: valid numerical orbitals for all elements in input file. The first layers
                                  keys are elements, the second are "pseudopotential identifier", the third
                                  are "numerical orbital identifier", the fourth are "type", "rcut", "appendix"
                                  and "file", the first three keys can concatenate as a "numerical orbital
                                  identifier".
    """
    # 1. program file I/O initialization: initialize cache directory, by either creating or checking
    # IS IT FOR PREPARING RUNTIME INFORMATION? yes

    # 2. prepare structure: download structure from remote server or generate by user setting, and return a dict whose keys are system
    #    (1) download structure from Materials Project to cache directory, the cif named as mp-xxx.cif
    #    (2) 
    # IS IT FOR PREPARING RUNTIME INFORMATION? yes, structure, download to cache directory
    # the structures is a dict whose keys are system formula and values are lists of tuples of cif file path (or ::AUTOGEN::[geometry]) and magnetism
    structures = search_structure(finp)
    elements = amsb.scan_elements(list(structures.keys()))
    # 3. translate input file to json and modify it
    # IS IT FOR PREPARING RUNTIME INFORMATION? yes
    runtime_settings = amiit.read_apns_inp(finp)
    # 4. scan valid pseudopotentials and numerical orbitals according to newly generated runtime settings
    # IS IT FOR PREPARING RUNTIME INFORMATION? yes, available pseudopotentials and numerical orbitals can be used in tests
    valid_upfs, valid_orbs = scan_upforb(elements=elements, basis_type=runtime_settings["calculation"]["basis_type"],
                                         pseudo_dir=runtime_settings["global"]["pseudo_dir"],
                                         orbital_dir=runtime_settings["global"]["orbital_dir"],
                                         pspot_settings=runtime_settings["pseudopotentials"],
                                         nao_settings=runtime_settings["numerical_orbitals"])
    # will return a dict, key is element and valus is a dict whose keys are "identifiers" of pseudopotentials
    # and values are the description of the pseudopotential
    # , the same for valid_orbs

    return runtime_settings, structures, valid_upfs, valid_orbs

import json
import os
import apns.module_structure.manage as amsm
def search_structure(finp: str):

    with open(finp, "r") as f:
        inp = json.load(f)
    # the case systems is defined as a list, means auto-generate structure no matter whether download from Materials
    # project or isolated system.
    if isinstance(inp["systems"], list):
        """in this case, the most complicated case would be like
        ```python
        {
            "systems": ["Er", "Er_dimer", ...]
        }
        ```
        This should be converted into a dict like:
        ```python
        {
            "systems": {
                "Er": ["/apns_cache/mp-1234.cif", "::AUTOGEN::DIMER"],
                ...
            }
        }
        ```
        """
        result = {}
        # isolated is reserved for performing numerical atomic orbital test (development motivation)
        isolated = [s for s in inp["systems"] if s.endswith("_dimer") or s.endswith("_trimer") or s.endswith("_tetramer")]
        for system in isolated:
            element, geometry = system.split("_")
            result.setdefault(element, []).append(("::AUTOGEN::" + geometry.upper(), None))
        # for practical systems
        crystal = [s for s in inp["systems"] if s not in isolated]
        # number of structures needed to search from Materials Project for each formula
        n_structures = inp["materials_project"]["n_structures"]
        # consider magnetism or not
        consider_magnetism = True if inp["calculation"]["nspin"] == 2 else False
        consider_magnetism = True if consider_magnetism and inp["extensive"]["magnetism"] == "materials_project" else False
        formula_tosearch = [] # always download magnetism information, this would not be harmful
        for system in crystal:
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
        return result
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
        return result
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
    """scan valid pseudopotential for all elements in input file
    
    Args:
        finp (str): input file
        
    Returns:
        tuple[dict, dict]: valid pseudopotentials and valid numerical orbitals
    """
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

import apns.module_pseudo.parse as ampgp
def pspot_software_availability(inp: dict, valid_pseudopotentials: dict, pseudopot_arch: dict) -> bool:
    """check if the input file is compatible with the software
    
    Args:
        inp (dict): input json
        valid_pseudopotentials (dict): valid pseudopotentials
        
    Raises:
    
    Returns:
        bool: True if compatible, False if not
    """

    """check files one by one"""

    for element in valid_pseudopotentials.keys():
        pseudopot_exceptions = []
        for identifier in valid_pseudopotentials[element].keys():
            fpspot = pseudopot_arch[identifier] + "/" + valid_pseudopotentials[element][identifier]["file"]
            print("Parsing pseudopotential file: ", fpspot)
            if not ampgp.is_compatible(fpspot, software=inp["global"]["software"].lower()):
                print("Compatiblity check failed for pseudopotential file: ", fpspot, " will be skipped.")
                pseudopot_exceptions.append(identifier)
            else:
                print("Compatiblity check OK.")
        for identifier in pseudopot_exceptions:
            del valid_pseudopotentials[element][identifier]
        if valid_pseudopotentials[element] == {}:
            return False
    return valid_pseudopotentials != {}

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
    unittest.main()