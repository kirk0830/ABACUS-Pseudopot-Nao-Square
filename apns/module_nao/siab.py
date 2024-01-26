def siab_program_section(hpc_environment_settings: list = None,
                        mpi_command: str = "mpirun -np 1",
                        abacus_command: str = "abacus"):
    
    result = ""
    if hpc_environment_settings is not None:
        result += "EXE_env "
        value = " ".join(hpc_environment_settings)
        result += value
        result += "\n"
    
    result += "EXE_mpi " + mpi_command + "\n"
    result += "EXE_pw " + abacus_command + "\n"
    return result

def siab_electronic_calculation(element: str,
                               ecutwfc: float,
                               rcut: float|list,
                               fpseudo: str,
                               pseudo_dir: str,
                               smearing_sigma: float):
    keys = ["element", "Ecut", "Rcut", "Pseudo_dir", "Pseudo_name", "smearing_sigma"]
    result = ""
    rcut = rcut if isinstance(rcut, float) else " ".join(rcut)

    values = [element, ecutwfc, rcut, pseudo_dir, fpseudo, smearing_sigma]
    for key, value in zip(keys, values):
        result += key + " " + str(value) + "\n"
    return result

def siab_reference_system(reference_systems: list):
    """reference_systems: list of dict, has structure like:
    
    reference_systems = [
        {
            "shape": "dimer",
            "nbands": 8,
            "lmaxmax": 2,
            "nspin": 1,
            "bond_lengths": [1.8, 2.0, 2.3, 2.8, 3.8],
        },
        {
            "shape": "trimer",
            "nbands": 10,
            "lmaxmax": 2,
            "nspin": 1,
            "bond_lengths": [1.9, 2.1, 2.6],
        }
    ]
    """
    def dict2str(d: dict):
        result = ""
        for key, value in d.items():
            if isinstance(value, list):
                value = " ".join([str(item) for item in value])
            result += value + " "
        return result
    for iref, reference_system in enumerate(reference_systems):
        result = "STRU" + str(iref) + " "
        result += dict2str(reference_system)
        result += "\n"
    return result

def siab_siabparams(reference_systems: list,
                    orbital_configurations: list,
                    maxstep: int = 9000):
    """orbital configurations is a list of dict, has structure like:
    [
        {
            "reference_structure": "dimer",
            "nbands_ref": 4,
            "from_to": [None, "SZ"]
        },
        {
            "reference_structure": "dimer",
            "nbands_ref": 4,
            "from_to": ["SZ", "DZP"]
        },
        {
            "reference_structure": "trimer",
            "nbands_ref": 6,
            "from_to": ["DZP", "TZDP"]
        }
    ]
    """
    result = "max_steps " + str(maxstep) + "\n"
    for iorb, orbital_configuration in enumerate(orbital_configurations):
        result += "Level" + str(iorb) + " "
        result += "STRU" + str(reference_systems.index(
            orbital_configuration["reference_structure"])) + " "
        result += str(orbital_configuration["nbands_ref"]) + " "
        inputorb = orbital_configuration["from_to"][0]
        outputorb = orbital_configuration["from_to"][1]

import re
def orbitalconfig_tolist(orbital_config: str):

    symbols = ["s", "p", "d", "f", "g", "h", "i", "k", "l", "m", "n", "o"]
    pattern = r"([0-9]+)([spdfghiklmno]+)"
    result = re.findall(pattern, orbital_config)
    result = [list(item) for item in result]
    for item in result:
        item[0] = int(item[0])
        item[1] = symbols.index(item[1])
    result.sort(key=lambda x: x[1])
    lmax = result[-1][1]
    swap = result
    result = [0 for i in range(lmax+1)]
    for item in swap:
        result[item[1]] = item[0]
    return result

def zeta_notation_toorbitalconfig(zeta_notation: str, minimum_basis: list = None):

    pattern = r"([SDTQPH]Z)([SDTQ5-9]?P)?"
    symbols = ["s", "p", "d", "f", "g", "h", "i", "k", "l", "m", "n", "o"]
    multiplier = {"S": 1, "D": 2, "T": 3, "Q": 4, "5": 5, "6": 6, "7": 7, "8": 8, "9": 9}
    _match = re.match(pattern, zeta_notation)
    if _match is None:
        raise ValueError("zeta_notation is not valid")
    nzeta = multiplier[_match.group(1)[0]]
    basis = [nzeta*i for i in minimum_basis]
    result = ""
    for i in range(len(minimum_basis)):
        if basis[i] != 0:
            result += str(basis[i]) + symbols[i]
    if _match.group(2) is not None:
        if len(_match.group(2)) > 1:
            result += str(multiplier[_match.group(2)[0]]) + symbols[len(minimum_basis)]
        else:
            result += "1" + symbols[len(minimum_basis)]
    return result

def SIAB_INPUT(element: str,
               ecutwfc: float,
               rcut: float|list,
               fpseudo: str,
               pseudo_dir: str,
               smearing_sigma: float,
               reference_systems: list,
               orbital_configurations: list,
               maxstep: int = 9000,
               hpc_environment_settings: list = None,
               mpi_command: str = "mpirun -np 1",
               abacus_command: str = "abacus"):
    """this is a temporarily function for presently refactor-not-complete PTG_dpsi program the 
    input script SIAB_INPUT generation"""
