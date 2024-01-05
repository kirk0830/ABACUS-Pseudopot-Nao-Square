"""initialize is the module should be called everytime will the program starts"""

import apns.module_io.input_translate as amiit
import apns.module_pseudo.upf_archive as ampua
# import apns.module_nao.nao_archive as amna
def initialize(finp: str) -> tuple[dict, dict, dict, dict, dict]:
    """initialize the program
    
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
    initialize_cache()
    system_with_mpids = download_structure(finp)
    _inp = amiit.inp_translate(fname=finp, system_with_mpids=system_with_mpids)
    valid_pseudopotentials, valid_numerical_orbitals = scan_valid_pseudopot_nao(finp)

    pseudopot_arch = ampua.load(_inp["global"]["pseudo_dir"])
    # nao_arch = amna.load(_inp["global"]["nao_dir"])

    return _inp, valid_pseudopotentials, valid_numerical_orbitals, pseudopot_arch, None

import os
import apns.module_workflow.identifier as amwi
def initialize_cache() -> None:

    print("Current working directory: {}".format(os.getcwd()))
    """change id.TEMPORARY_FOLDER to absolute path"""
    amwi.TEMPORARY_FOLDER = os.path.join(os.getcwd(), amwi.TEMPORARY_FOLDER)
    """create cache directory if not exist"""
    if not os.path.exists(amwi.TEMPORARY_FOLDER):
        os.mkdir(amwi.TEMPORARY_FOLDER)
    else:
        print("Cache directory already exists.")

import json
import apns.module_structure.materials_project as amsmp
def download_structure(finp: str) -> dict:
    """download structure from remote server and return a dict whose keys are system 
    formula and values are lists of corresponding system_mpids"""

    with open(finp, "r") as f:
        inp = json.load(f)

    system_with_mpids = amsmp.composites(api_key=inp["materials_project"]["api_key"],
                                         formula=inp["systems"],
                                         num_cif=inp["materials_project"]["n_structures"],
                                         theoretical=inp["materials_project"]["theoretical"],
                                         is_stable=inp["materials_project"]["most_stable"])

    return system_with_mpids

import apns.module_structure.basic as amsb
import apns.module_pseudo.local_validity_scan as amplvs
import apns.module_nao.local_validity_scan as amnlvs
def scan_valid_pseudopot_nao(finp: str) -> tuple[dict, dict]:
    """scan valid pseudopotential for all elements in input file
    
    Args:
        finp (str): input file
        
    Raises:
        
    Returns:
        tuple[dict, dict]: valid pseudopotentials and valid numerical orbitals
    """

    valid_pseudopotentials = []
    valid_numerical_orbitals = []

    with open(finp, "r") as f:
        inp = json.load(f)

    elements = amsb.scan_elements(inp["systems"])

    valid_pseudopotentials = amplvs._svp_(elements, inp["pseudopotentials"])
    for element in elements:
        if element not in valid_pseudopotentials.keys():
            raise ValueError("No valid pseudopotential for element {}.".format(element))

    if inp["calculation"]["basis_type"] == "lcao":
        raise NotImplementedError("lcao calculation is not supported yet.")
        """TO BE IMPLEMENTED
        valid_numerical_orbitals = amnlvs._svno_(element, valid_pseudopotentials, inp["numerical_orbitals"])
        for element in elements:
            if element not in valid_numerical_orbitals.keys():
                raise ValueError("No valid numerical orbital for element {}.".format(element))
        """
    
    return valid_pseudopotentials, valid_numerical_orbitals
