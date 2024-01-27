"""Generate SIAB_INPUT file for PTG_dpsi program

Main function: SIAB_INPUT
Parameters:
    element (str), ecutwfc (float), nspin (int), rcut (float|list), fpseudo (str), pseudo_dir (str), 
    minimal_basis (list), smearing_sigma (float), reference_systems (list), orbital_configurations (list),
    maxstep (int), hpc_environment_settings (list), mpi_command (str), abacus_command (str)
Interface:
    manually set:
        element, nspin, rcut, smearing_sigma, orbital_configuration, mpi_command, abacus_command
    manually set (partially):
        only shape in reference_system
    module_workflow/initialize.py set:
        ecutwfc, fpseudo, pseudo_dir
    can automatically calculated:
        minimal_basis, reference_systems, orbital_configuration
    seldom used:
        hpc_environment_settings
    already has reasonable default values:
        maxstep
"""

from collections import defaultdict

def merge_dicts(source: list, value_tostr: bool = False, delimiter: str = " "):
    """merge a list of dicts into one dict whose values are lists (certainly)"""

    """exceptions"""
    keys = list(source[0].keys())
    for _dict in source:
        for key in _dict.keys():
            if key not in keys:
                raise ValueError("keys in source dicts are not consistent")
    """functionality"""
    result = defaultdict(list)
    for _dict in source:
        for key, value in _dict.items():
            if value_tostr:
                value = delimiter.join([str(item) for item in value]) if isinstance(value, list) else str(value)
            result[key].append(value)
    return dict(result)

def dict_tostr(source: dict,
               key_sequence: list,
               direction: str,
               comment_sign: str,
               with_key: bool,
               len_placeholder: int) -> str:
    """print a dict to table, can orgainized in vertical or horizontal direction"""
    """initialize"""
    for key, value in source.items():
        if not isinstance(value, list):
            source[key] = [value]
    """exceptions"""
    if len(source) == 0:
        raise ValueError("source is empty")
    if key_sequence is None:
        print("Warning: key_sequence is not provided, will print in random sequence. This may cause error.")
        key_sequence = list(source.keys())
    if direction not in ["vertical", "horizontal"]:
        raise ValueError("direction must be 'vertical' or 'horizontal'")
    len_value = len(source[key_sequence[0]])
    for key in key_sequence:
        if key not in source.keys():
            raise ValueError("key_sequence contains key not in source")
        if len(source[key]) != len_value:
            raise ValueError("length of values in source is not consistent")
    """functionality"""
    result = ""
    if direction == "vertical":
        for key in key_sequence:
            result += f"%-{len_placeholder}s" % key if with_key else ""
            for i in range(len_value):
                result += f"%-{len_placeholder}s" % str(source[key][i])
            result += "\n"
    elif direction == "horizontal":
        titles = key_sequence
        prefix = comment_sign + " " if with_key else ""
        result += prefix
        for title in titles:
            result += f"%{len_placeholder}s" % title
        result += "\n"
        for i in range(len_value):
            result += " "*len(prefix)
            for title in titles:
                result += f"%{len_placeholder}s" % str(source[title][i])
            result += "\n"
    return result

def dicts_tostr(source: list|dict, 
                key_sequence: list = None, 
                direction: str = "vertical",
                comment_sign: str = None,
                with_key: bool = True, 
                len_placeholder: int = 10):
    """print a list of dict to table, can orgainized in vertical or horizontal direction
    
    source: list of dict, has structure like:
    [
        {
            "key1": value1,
            "key2": value2,
            ...
        },
        {
            "key1": value1,
            "key2": value2,
            ...
        },
        ...
    ]

    key_sequence: list of str, the sequence of keys to be printed. NOTE: BECAUSE THE KEY IN
    DICTIONARY IS NOT ORDERED, SO THE ORDER OF KEYS IN THE DICTIONARY IS NOT GUARANTEED.

    direction: str, "vertical" or "horizontal". "vertical" means the table is organized as 
    key value in each row, "horizontal" means the table will has line of titles at only the 
    first row, and the values of each key will be printed in each column.

    comment_sign: str, NECESSARY FOR `direction` == "horizontal", the comment sign will be added
    at the beginning of the line of titles.

    with_key: bool, whether to print the key in the table.

    len_placeholder: int, the length of placeholder for each column.
    """
    """initialize"""
    if isinstance(source, dict):
        source = [source]
    """exceptions"""
    if len(source) == 0:
        raise ValueError("source is empty")
    if key_sequence is None:
        print("Warning: key_sequence is not provided, will print in random sequence. This may cause error.")
        key_sequence = list(source[0].keys())
    for _dict in source:
        if len(_dict) != len(key_sequence):
            raise ValueError("length of key_sequence is not equal to length of dict")
    if direction not in ["vertical", "horizontal"]:
        raise ValueError("direction must be 'vertical' or 'horizontal'")
    if direction == "horizontal" and comment_sign is None:
        raise ValueError("comment_sign is necessary for horizontal direction")
    """functionality"""

    return dict_tostr(source=merge_dicts(source),
                      key_sequence=key_sequence,
                      direction=direction,
                      comment_sign=comment_sign,
                      with_key=with_key,
                      len_placeholder=len_placeholder)
    
def siab_program_section(hpc_environment_settings: list = None,
                         mpi_command: str = "mpirun -np 1",
                         abacus_command: str = "abacus"):
    
    source = {}
    key_sequence = ["EXE_env", "EXE_mpi", "EXE_pw"]
    if hpc_environment_settings is not None:
        source["EXE_env"] = " ".join(hpc_environment_settings)
    else:
        source["#EXE_env"] = ""
        key_sequence[0] = "#EXE_env"
    source.update(dict(zip(["EXE_mpi", "EXE_pw"], [mpi_command, abacus_command])))

    return dict_tostr(source=source,
                      key_sequence=key_sequence,
                      direction="vertical",
                      comment_sign="#",
                      with_key=True,
                      len_placeholder=20)

def siab_electronic_calculation(element: str,
                                ecutwfc: float,
                                rcut: float|list,
                                fpseudo: str,
                                pseudo_dir: str,
                                smearing_sigma: float):
    keys = ["element", "Ecut", "Rcut", "Pseudo_dir", "Pseudo_name", "smearing_sigma"]
    rcut = str(rcut) if not isinstance(rcut, list) else " ".join([str(item) for item in rcut])

    values = [element, ecutwfc, rcut, pseudo_dir, fpseudo, smearing_sigma]
    return dict_tostr(source=dict(zip(keys, values)),
                      key_sequence=keys,
                      direction="vertical",
                      comment_sign="#",
                      with_key=True,
                      len_placeholder=20)

def siab_reference_system(reference_systems: list, 
                          nspin: int = 1, 
                          lmax: int = 2,
                          orbital_configurations: list = None):
    """reference_systems: list of dict, has structure like:
    
    reference_systems = [
        {
            "shape": "dimer",
            "nbands": 8,
            "bond_lengths": [1.8, 2.0, 2.3, 2.8, 3.8],
        },
        {
            "shape": "trimer",
            "nbands": 10,
            "bond_lengths": [1.9, 2.1, 2.6],
        }
    ]
    """
    """statistic lmax for different reference systems from orbital_configurations"""
    shapes = [item["shape"] for item in reference_systems]
    lmaxs = [lmax]*len(shapes)
    for ish, shape in enumerate(shapes):
        for orb_config in orbital_configurations:
            if orb_config["reference_structure"] == shape:
                if orb_config["from_to"][1].endswith("P"):
                    lmaxs[ish] += 1
    key_sequence = ["identifier", "shape", "nbands", "lmax", "nspin", "bond_lengths"]
    for irs, reference_system in enumerate(reference_systems):
        reference_system["identifier"] = "STRU" + str(irs+1)
        reference_system["lmax"] = lmaxs[irs]
        reference_system["nspin"] = nspin
        reference_system["bond_lengths"] = " ".join([str(item) for item in reference_system["bond_lengths"]])
    return dicts_tostr(source=reference_systems,
                       key_sequence=key_sequence,
                       direction="horizontal",
                       comment_sign="#",
                       with_key=True,
                       len_placeholder=-15)

def siab_siabparams(minimal_basis: list,
                    reference_systems: list,
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

    reference_systems = [item["shape"] for item in reference_systems]
    key_sequence = ["orb_id", "stru_id", "nbands_ref", "orb_ref", "orb_config"]
    source = []
    for iorb, orbital_configuration in enumerate(orbital_configurations):
        shape = orbital_configuration["reference_structure"]
        structure_id = "STRU" + str(reference_systems.index(shape)+1)
        from_to = orbital_configuration["from_to"]
        orbital_referred = "none" if from_to[0] is None else "fix"
        config = zeta_notation_toorbitalconfig(from_to[1], minimal_basis)
        source.append({
            "orb_id": "Level" + str(iorb+1),
            "stru_id": structure_id,
            "nbands_ref": orbital_configuration["nbands_ref"],
            "orb_ref": orbital_referred,
            "orb_config": config
        })
    result += dicts_tostr(source=source,
                          key_sequence=key_sequence,
                          direction="horizontal",
                          comment_sign="#",
                          with_key=True,
                          len_placeholder=-15)
    return result

def siab_save(orbital_configurations: list):
    key_sequence = ["save_id", "orb_id", "zeta_notation"]
    source = []
    for iorb, orbital_configuration in enumerate(orbital_configurations):
        from_to = orbital_configuration["from_to"]
        source.append({
            "save_id": "Save" + str(iorb+1),
            "orb_id": "Level" + str(iorb+1),
            "zeta_notation": from_to[1]
        })
    return dicts_tostr(source=source,
                       key_sequence=key_sequence,
                       direction="horizontal",
                       comment_sign="#",
                       with_key=True,
                       len_placeholder=-15)

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

def zeta_notation_toorbitalconfig(zeta_notation: str, minimal_basis: list = None):

    pattern = r"([SDTQPH]Z)([SDTQ5-9]?P)?"
    symbols = ["s", "p", "d", "f", "g", "h", "i", "k", "l", "m", "n", "o"]
    multiplier = {"S": 1, "D": 2, "T": 3, "Q": 4, "5": 5, "6": 6, "7": 7, "8": 8, "9": 9}
    _match = re.match(pattern, zeta_notation)
    if _match is None:
        raise ValueError("zeta_notation is not valid")
    nzeta = multiplier[_match.group(1)[0]]
    basis = [nzeta*i for i in minimal_basis]
    result = ""
    for i in range(len(minimal_basis)):
        if basis[i] != 0:
            result += str(basis[i]) + symbols[i]
    if _match.group(2) is not None:
        if len(_match.group(2)) > 1:
            result += str(multiplier[_match.group(2)[0]]) + symbols[len(minimal_basis)]
        else:
            result += "1" + symbols[len(minimal_basis)]
    return result

import apns.module_pseudo.general_parser as ampgp
def SIAB_INPUT(element: str,                            # element name
               ecutwfc: float,                          # kinetic energy cutoff of planewave basis
               nspin: int,                              # number of spin channels
               rcut: float|list,                        # realspace cutoff for numerical atomic orbitals
               fpseudo: str,                            # pseudopotential file name
               pseudo_dir: str,                         # directory of pseudopotential file
               minimal_basis: list,                     # minimal basis, usually to be the configuration of valence electrons
               smearing_sigma: float,                   # smearing parameter
               reference_systems: list,                 # reference system for learning planewave wavefunction
               orbital_configurations: list,            # zeta and polarization of orbital configurations
               maxstep: int = 9000,                     # maximum number of iteration steps
               hpc_environment_settings: list = None,   # hpc environment settings
               mpi_command: str = "mpirun -np 1",       # mpi command for running ABACUS, conventionally `mpirun -np [ncores]`
               abacus_command: str = "abacus"):         # command to run ABACUS
    """Main program for generating SIAB_INPUT file for PTG_dpsi program

    Args:
        element (str): element name
        ecutwfc (float): kinetic energy cutoff of planewave basis
        nspin (int): number of spin channels, for formal generation task, this is a quantity needed to test
        rcut (float|list): realspace cutoff for numerical atomic orbitals. From ABACUS 3.5.1, all rcut values can
                           be calculated in one shot.
        fpseudo (str): pseudopotential file name
        pseudo_dir (str): directory of pseudopotential file
        minimal_basis (list): minimal basis, usually to be the configuration of valence electrons. If set to None,
                              program will call pseudopotential parser to grep valence electrons information directly
                              from pseudopotential file.
        smearing_sigma (float): smearing parameter, default value is 0.015 (unit Ry)
        reference_systems (list): reference system for learning planewave wavefunction. The structure of this list is
        like: 
        ```json
        [
           {
               "shape": "dimer",
               "nbands": 8,
               "bond_lengths": [1.8, 2.0, 2.3, 2.8, 3.8],
           },
           {
               "shape": "trimer",
               "nbands": 10,
               "bond_lengths": [1.9, 2.1, 2.6],
           }
        ]
        ```

        * the `shape` is always to be "dimer", "trimer", "tetramer". nbands is the number of
           bands to calculate in pw calculation. `bond_lengths` specifies the bond streching,
           representing various chemical environments. This parameter is not recommended to
           set manually, instead, should be set based on the scan of bond lengths, for example,
           provided by APNS project. In this way, `nbands` can also be set automatically based
           on the auto-parse of pseudopotential file. Default is to set nLUMO = nHOMO.

        orbital_configurations (list): zeta and polarization of orbital configurations. The structure of this list is
        like: 
        ```json
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
        ```
        * orbital_configuration specifies how many levels of orbitals wants to generate. For
          each level, one can specify the reference structure, the number of bands (therefore
          number of states set to learn by numerical atomic orbitals, or say number of bands
          to reproduce). During the generation of orbitals, the higher level can always be
          based on the previous level to generate, so the `from_to` parameter is used to
          specify the zeta notation of the previous level. The second value defines the zeta
          notation of present level. The first value can be set to None, which means generate
          orbital from "void".
        
        maxstep (int, optional): maximum number of iteration steps. Defaults to 9000.
        hpc_environment_settings (list, optional): hpc environment settings. Defaults to None.
        mpi_command (str, optional): mpi command for running ABACUS, conventionally `mpirun -np [ncores]`. Defaults to "mpirun -np 1".
        abacus_command (str, optional): command to run ABACUS. Defaults to "abacus".
    Return:
        str: SIAB_INPUT file content
    """
    result = "# Refactor version of SIAB_INPUT of PTG_dpsi for generating numerical atomic orbitals of ABACUS\n"
    result += "# from ABACUS-Planewave DFT calculations - For high throughput auto-generation and test of\n"
    result += "# pseudopotentials and numerical atomic orbitals. This is included in project ABACUS Pseudopot-\n"
    result += "# Nao Square (APNS). Visit related Github repos for more information:\n"
    result += "# ABACUS (deepmodeling) Github repo: https://github.com/deepmodeling/abacus-develop\n"
    result += "# PTG_dpsi (abacusmodeling/ABACUS-orbitals) Github repo: https://github.com/abacusmodeling/ABACUS-orbitals\n"
    result += "# APNS Github repo: https://github.com/kirk0830/ABACUS-Pseudopot-Nao-Square\n"
    result += "# APNS Github Pages: https://kirk0830.github.io/ABACUS-Pseudopot-Nao-Square\n"
    result += "# APNS is mainly developed and maintained by ABACUS-AISI developer team\n\n"
    result += "# PROGRAM CONFIGURATION\n"
    result += siab_program_section(hpc_environment_settings=hpc_environment_settings,
                                   mpi_command=mpi_command,
                                   abacus_command=abacus_command)
    result += "\n"
    result += "# ELECTRONIC STRUCTURE CALCULATION\n"
    result += siab_electronic_calculation(element=element,
                                          ecutwfc=ecutwfc,
                                          rcut=rcut,
                                          fpseudo=fpseudo,
                                          pseudo_dir=pseudo_dir,
                                          smearing_sigma=smearing_sigma)
    result += "\n"
    result += "# REFERENCE SYSTEMS\n"
    if minimal_basis is None:
        minimal_basis = ampgp.valelec_config(pseudo_dir + "/" + fpseudo)
        minimal_basis = [len(item) for item in minimal_basis]
    result += siab_reference_system(reference_systems=reference_systems,
                                    nspin=nspin,
                                    lmax=len(minimal_basis)-1,
                                    orbital_configurations=orbital_configurations)
    result += "\n"
    result += "# SIAB PARAMETERS\n"
    result += siab_siabparams(minimal_basis=minimal_basis,
                              reference_systems=reference_systems,
                              orbital_configurations=orbital_configurations,
                              maxstep=maxstep)
    result += "\n"
    result += "# SAVE\n"
    result += siab_save(orbital_configurations=orbital_configurations)

    return result

