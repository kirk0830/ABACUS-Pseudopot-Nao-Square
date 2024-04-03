"""initialize is the module should be called everytime will the program starts"""

import apns.module_io.input_translate as amiit
# import apns.module_nao.nao_archive as amna
def initialize(finp: str, test_mode: bool = False) -> tuple[dict, dict, dict, dict, dict]:
    """initialize the program, return runtime information can be determined before running the workflow
    
    Args:
        finp (str): input file
        
    Raises:
    
    Returns:
        tuple[dict, dict, dict, dict, dict]: input, valid_pseudopotentials, valid_numerical_orbitals, pseudopot_arch, nao_arch
        
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
        pseudopot_arch: pseudopotential archive, a dict whose keys are "identifiers" of pseudopotentials,
                        stored in `pseudo_dir`, values are folders' absolute path.
        nao_arch: numerical orbital archive, a dict whose keys are "numerical orbital identifier", stored in
                  `nao_dir`, values are folders' absolute path.
    """
    # 1. program file I/O initialization: initialize cache directory, by either creating or checking
    # IS IT FOR PREPARING RUNTIME INFORMATION? yes

    # 2. prepare structure: download structure from remote server or generate by user setting, and return a dict whose keys are system
    #    (1) download structure from Materials Project to cache directory, the cif named as mp-xxx.cif
    #    (2) 
    # IS IT FOR PREPARING RUNTIME INFORMATION? yes, structure, download to cache directory
    system_with_mpids = download_structure(finp)
    # 3. translate input file to json and modify it
    # IS IT FOR PREPARING RUNTIME INFORMATION? yes
    runtime_settings = amiit.inp_translate(fname=finp, system_with_mpids=system_with_mpids)
    # 4. scan valid pseudopotentials and numerical orbitals according to newly generated runtime settings
    # IS IT FOR PREPARING RUNTIME INFORMATION? yes, available pseudopotentials and numerical orbitals can be used in tests
    valid_upfs, valid_orbs = scan_upforb(runtime_settings)
    # will return a dict, key is element and valus is a dict whose keys are "identifiers" of pseudopotentials
    # and values are the description of the pseudopotential
    # , the same for valid_orbs

    return runtime_settings, valid_upfs, valid_orbs

import json
import apns.module_structure.materials_project as amsmp
def download_structure(finp: str) -> dict:
    """download structure from remote server and return a dict whose keys are system 
    formula and values are lists of corresponding system_mpids
    
    There is also a special case where system has name like "xxx_dimer", "xxx_trimer"
    or "xxx_tetramer", which means the system is a isolated molecule, in this case,
    the system will not be downloaded from remote server. It will directly assume as 
    a "system" with its "mpid", although the "mpid" now is "dimer", "trimer" or
    "tetramer".

    Args:
        finp (str): input file
    
    Raises:

    Returns:
        dict: a dict whose keys are system formula and values are lists of corresponding system_mpids
    """
    # open 
    with open(finp, "r") as f:
        inp = json.load(f)
    # isolated is reserved for performing numerical atomic orbital test (development motivation)
    isolated = [s for s in inp["systems"] if s.endswith("_dimer") or s.endswith("_trimer") or s.endswith("_tetramer")]
    # routinely used, for practical systems
    crystal = [s for s in inp["systems"] if s not in isolated]
    # do not support specifying two kinds of systems together
    if len(isolated)*len(crystal) != 0:
        raise ValueError("Isolated molecule and crystal cannot be mixed.")
    # support the feature that one chemical formula with many mp-id, say many structures
    system_with_mpids = {}
    if len(crystal) > 0:
        consider_magnetism = True if (inp["calculation"]["nspin"] == 2 and inp["extensive"]["magnetism"] == "materials_project") else False
        system_with_mpids = amsmp.composites(api_key=inp["materials_project"]["api_key"],
                                             formula=crystal,
                                             num_cif=inp["materials_project"]["n_structures"],
                                             theoretical=inp["materials_project"]["theoretical"],
                                             is_stable=inp["materials_project"]["most_stable"],
                                             consider_magnetism=consider_magnetism)
    elif len(isolated) > 0:
        for system in isolated:
            element = system.split("_")[0]
            system_with_mpids.setdefault(element, []).append(system)
    else:
        raise ValueError("No system to download or generate, check your setting.")
    
    return system_with_mpids



import apns.module_structure.basic as amsb
import apns.module_pseudo.manage as ampm
def scan_upforb(finp: str|dict) -> tuple[dict, dict]:
    """scan valid pseudopotential for all elements in input file
    
    Args:
        finp (str): input file
        
    Raises:
        
    Returns:
        tuple[dict, dict]: valid pseudopotentials and valid numerical orbitals
    """
    
    if isinstance(finp, str):
        with open(finp, "r") as f:
            inp = json.load(f)
    elif isinstance(finp, dict):
        inp = finp
    else:
        raise TypeError("finp should be str or dict.")
    
    """get all elements across all systems"""
    elements = amsb.scan_elements(inp["systems"])

    """from elements, get all valid pseudopotentials"""
    valid_upfs = ampm.valid_pseudo(pseudo_dir=inp["global"]["pseudo_dir"], 
                                   elements=elements, 
                                   pseudo_setting=inp["pseudopotentials"])

    """from elements, get all valid numerical orbitals"""
    valid_orbs = {element: {} for element in elements}
    if inp["calculation"]["basis_type"] == "lcao":
        raise NotImplementedError("lcao calculation is not supported yet.")
    
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