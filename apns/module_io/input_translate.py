"""Design of input for user
see keywords.md for details
"""
import json
import os
from apns.module_structure.basic import scan_elements

def check(inp: dict) -> None:
    """check reasonality of input parameters

    Args:
        inp (dict): parsed input

    Raises:
        ValueError: unreasonable parameter will raise this error
    """
    """basic"""
    if inp["global"]["software"] == "qespresso":
        if inp["calculation"]["basis_type"] == "lcao":
            raise ValueError("Quantum ESPRESSO only supports pw calculation.")
    
    """on the Feature Request of EOS calculation"""
    # for EOS calculation, force user must specify ecutwfc for each system. But can also give
    # value null (None) to indicate that result from convergence test will be used.

    """on the Feature Request of band structure calculation"""
    if inp["calculation"]["calculation"] == "scf":
        if inp["extensive"]["nkpoints_in_line"] > 0:
            print("calculation: ", inp["calculation"]["calculation"])
            print("nkpoints_in_line: ", inp["extensive"]["nkpoints_in_line"])
            raise ValueError("confused with calculation requested: for nkpoints > 0 specifies a band calculation, not a scf calculation.")
    
    """on the Feature Request of initialization of magnetism"""
    if "nspin" in inp["calculation"].keys():
        if inp["calculation"]["nspin"] == 1:
            if inp["extensive"]["magnetism"] != "nonmagnetic":
                print("Warning: for nspin = 1, specifying magnetism is meaningless.")
        else:
            if inp["extensive"]["magnetism"] == "materials_project":
                print("Will use magnetism from data reported on materials project website.")
            elif inp["extensive"]["magnetism"] not in ["nonmagnetic", "ferromagnetic", "antiferromagnetic"]:
                raise ValueError("for nspin = 2, magnetism must be 'nonmagnetic', 'ferromagnetic' or 'antiferromagnetic', or 'materials_project'.")
    else:
        print("nspin not specified, will use nspin = 1 case by default.")

def default(inp: dict) -> dict:
    """set default value for unset keywords for global and calculation sections

    Args:
        inp (dict): the read input

    Returns:
        dict: input filled with default values.
    """
    sections = ["global", "calculation", "extensive"]
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
    # expand pseudopotentials
    for key in inp["pseudopotentials"]:
        if isinstance(inp["pseudopotentials"][key], list):
            _dict = {}
            for element in elements:
                _dict[element] = inp["pseudopotentials"][key]
            inp["pseudopotentials"][key] = _dict
    # expand numerical_orbitals
    if inp["calculation"]["basis_type"] == "lcao":
        for key in inp["numerical_orbitals"]:
            if isinstance(inp["numerical_orbitals"][key], list):
                _dict = {}
                for element in elements:
                    _dict[element] = inp["numerical_orbitals"][key]
                inp["numerical_orbitals"][key] = _dict

    return inp

def toabspath(relative_path: str) -> str:
    """convert one path to absolute one, possible to raise NonExistError

    Args:
        relative_path (str): _description_

    Returns:
        str: absolute path
    """
    path_0 = os.getcwd().replace("\\", "/")
    os.chdir(relative_path)
    path = os.getcwd().replace("\\", "/")
    os.chdir(path_0)

    return path

def inp_translate(fname: str, **kwargs) -> dict:
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
    # expand systems, sed the system name with system_mpids in read input
    if "system_with_mpids" in kwargs:
        for system in kwargs["system_with_mpids"].keys(): # example: system = Er2O3
            if system in inp["systems"]:
                # remove this system from inp["systems"] list
                inp["systems"].remove(system)
            for structure in kwargs["system_with_mpids"][system]: # example: structure = Er2O3_2460
                if structure not in inp["systems"]:
                    # for `crystal` case, for example Cr crystal, the above remove will remove Cr, 
                    # but Cr_dimer cannot be removed because it is not the key of kwargs["system_with_mpids"]
                    # , it is just the case of `isolated
                    inp["systems"].append(structure)
    else:
        raise ValueError("system_with_mpids not found in kwargs.")
    # expand pseudopotentials and numerical_orbitals from list to dict
    inp = expand(inp, elements)
    # for unset attributes, use default values
    inp = default(inp)
    # convert ralative path to absolute path
    inp["global"]["work_dir"] = toabspath(inp["global"]["work_dir"])
    inp["global"]["pseudo_dir"] = toabspath(inp["global"]["pseudo_dir"])
    inp["global"]["orbital_dir"] = toabspath(inp["global"]["orbital_dir"])
    # check rationality of input
    check(inp)

    return inp # then we call the output dict as work_status, instead of input

DEFAULT_INPUT = {
    "global": {
        "test_mode": "pseudopotential",
        "software": "ABACUS",
        "work_dir": "./",
        "pseudo_dir": "./download/pseudopotentials/",
        "orbital_dir": "./download/numerical_orbitals/",
        "save_log": True
    },
    "calculation": {
        "basis_type": "pw",
        "functionals": ["PBE"],
        "ecutwfc": [100],
        "calculation": "scf",
        "nspin": 1
    },
    "extensive": {
        "characteristic_lengths": [0.00],
        "nkpoints_in_line": 0,
        "magnetism": "nonmagnetic"
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