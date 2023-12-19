"""Design of input for user
see keywords.md for details
"""
import json
import os
from module_structure.basic import scan_elements

def check(inp: dict) -> None:

    if inp["global"]["software"] == "qespresso":
        if inp["calculation"]["basis_type"] == "lcao":
            raise ValueError("Quantum ESPRESSO only supports pw calculation.")
        
def default(inp: dict) -> dict:
    """set default value for unset keywords for global and calculation sections

    Args:
        inp (dict): the read input

    Returns:
        dict: input filled with default values.
    """
    sections = ["global", "calculation"]
    for section in sections:
        for key in DEFAULT_INPUT[section].keys():
            if key not in inp[section].keys():
                inp[section][key] = DEFAULT_INPUT[section][key]
    return inp

def expand(inp: dict, elements: list) -> dict:
    """expand pseudopotentials and numerical_orbitals from list to dict

    Args:
        inp (dict): 
        elements (list): elements in test

    Returns:
        dict: pseudopotential and numerical_orbitals expanded input
    """
    for key in inp["pseudopotentials"]:
        if isinstance(inp["pseudopotentials"][key], list):
            _dict = {}
            for element in elements:
                _dict[element] = inp["pseudopotentials"][key]
            inp["pseudopotentials"][key] = _dict
    for key in inp["numerical_orbitals"]:
        if isinstance(inp["numerical_orbitals"][key], list):
            _dict = {}
            for element in elements:
                _dict[element] = inp["numerical_orbitals"][key]
            inp["numerical_orbitals"][key] = _dict

    return inp

def convert_to_absolute_path(relative_path: str) -> str:
    """convert one path to absolute one, possible to raise NonExistError

    Args:
        relative_path (str): _description_

    Returns:
        str: absolute path
    """
    path_0 = os.getcwd()
    os.chdir(relative_path)
    path = os.getcwd()
    os.chdir(path_0)

    return path

def render(fname: str, **kwargs) -> dict:
    """render input file to a dict
    1. expand pseudopotentials and numerical_orbitals information to element-by-element
    2. set default values if not explicity specified
    3. convert relative path to absolute one

    Args:
        fname (str): input file name

    Returns:
        dict: input dict, see keywords.md for details.
    """
    # read
    with open(fname, "r") as f:
        inp = json.load(f)
    # get all elements as list from inp
    elements = []
    for system in inp["systems"]:
        _elements = scan_elements(system)
        for element in _elements:
            if element not in elements:
                elements.append(element)
    # expand systems
    if "system_mpids" in kwargs:
        for system in kwargs["system_mpids"]:
            if system in inp["systems"]:
                # remove this system from inp["systems"] list
                inp["systems"].remove(system)
            for structure in kwargs["system_mpids"][system]:
                inp["systems"].append(structure)
    # expand pseudopotentials and numerical_orbitals from list to dict
    inp = expand(inp, elements)
    # for unset attributes, use default values
    inp = default(inp)
    # convert ralative path to absolute path
    inp["global"]["work_dir"] = convert_to_absolute_path(inp["global"]["work_dir"])
    inp["global"]["pseudo_dir"] = convert_to_absolute_path(inp["global"]["pseudo_dir"])
    inp["global"]["orbital_dir"] = convert_to_absolute_path(inp["global"]["orbital_dir"])
    # check rationality of input
    check(inp)

    return inp

DEFAULT_INPUT = {
    "global": {
        "test_mode": "pseudopotential",
        "software": "ABACUS",
        "work_dir": "./",
        "pseudo_dir": "../download/pseudopotentials/",
        "orbital_dir": "../download/numerical_orbitals/",
        "save_log": True
    },
    "calculation": {
        "basis_type": "pw",
        "functionals": ["PBE"],
        "ecutwfc": [100],
        "cell_scaling": [0.00]
    },
    "systems": [],
    "pseudopotentials": {
        "kinds": {},
        "versions": {},
        "appendices": {}
    },
    "numerical_orbitals": {
        "types": {},
        "rcuts": {},
        "appendices": {}
    }
}