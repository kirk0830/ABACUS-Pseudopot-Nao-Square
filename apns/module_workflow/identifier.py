"""
This file controls all nomenclature convention used in whole project.
"""

# environment variables
TEMPORARY_FOLDER = "apns_cache"

# free functions
def pseudopotential(kind: str, version: str, appendix: str) -> str:
    """Gnerate identifier of pseudopotential

    Args:
        kind (str): pseudo potential kind
        version (str): pseudo potential version
        appendix (str): pseudo potential appendix

    Raises:
        ValueError: kind cannot be empty

    Returns:
        str: identifier of pseudopotential
    """
    if kind == "":
        raise ValueError("kind cannot be empty")
    return "_".join(i for i in [kind, version, appendix] if i != "")

def r_pseudopotential(identifier: str) -> tuple:
    """Parse identifier of pseudopotential

    Args:
        identifier (str): identifier of pseudopotential

    Returns:
        tuple: kind, version, appendix
    """
    words = identifier.split("_")
    if len(words) == 1:
        return words[0], "", ""
    elif len(words) == 2:
        return words[0], words[1], ""
    else:
        return words[0], words[1], "_".join(words[2:])

def numerical_orbital(type: str, rcut: int, appendix: str, version: bool = "old") -> str:
    """Generate identifier of numerical orbital

    Args:
        type (str): numerical orbital type, DZP, TZDP, etc.
        rcut (int): numerical orbital rcut
        appendix (str): numerical orbital appendix

    Raises:
        ValueError: type cannot be empty

    Returns:
        str: identifier of numerical orbital
    """
    if type == "":
        raise ValueError("type cannot be empty")
    if rcut != 0:
        rcut = str(rcut)
    if version == "old":
        return "_".join(i for i in [type, str(rcut), appendix] if i != "")
    elif version == "new":
        return "_".join(i for i in [type[0], str(rcut), appendix] if i != "")

def r_numerical_orbital(identifier: str) -> tuple:
    """Parse identifier of numerical orbital

    Args:
        identifier (str): identifier of numerical orbital

    Returns:
        tuple: type, rcut, appendix
    """
    words = identifier.split("_")
    if len(words) == 1:
        return words[0], "", ""
    elif len(words) == 2:
        return words[0], words[1], ""
    else:
        return words[0], words[1], "_".join(words[2:])

def pseudopot_nao(pseudopotential: list, numerical_orbital: list = []) -> str:

    result = ""
    for pseudo in pseudopotential:
        result += pseudo.replace("_", "")
    if len(numerical_orbital) > 0:
        result += "_"
        for orbital in numerical_orbital:
            words = orbital.split("_")
            result += words[0][0]+words[1]
            if len(words) > 2:
                result += words[2]
    return result

def folder(functional: str, system: str, specific_test: str) -> str:

    functional = functional.lower()
    folder_name_fragments = ["t", functional, system, specific_test]
    return "_".join(folder_name_fragments)

def r_folder(identifier: str) -> tuple:

    return tuple(i for i in identifier.split("_") if i != "")

def _folder_(system: str, pseudo_nao_identifier: str, calculation_identifier: str) -> str:
    
    """new version of Generate identifier of folder

    Args:
        system (str): system name
        pseudo_nao_identifier (str): identifier of pseudo_nao
        calculation_identifier (str): identifier of calculation

    Returns:
        str: identifier of folder
    """
    folder_name_fragments = [system, pseudo_nao_identifier, calculation_identifier]
    return "_".join(folder_name_fragments)

def qespresso(system: str, template: bool = False):

    """Generate identifier of Quantum ESPRESSO input script

    Args:
        system (str): system name
        template (bool, optional): whether the input script is a template. Defaults to False.

    Returns:
        str: identifier of Quantum ESPRESSO input script
    """
    result = "qespresso"
    if template:
        result += "_" + system
    result += ".in"
    return result

def abacus(system: str, template: bool = False):

    """Generate identifier of ABACUS input scripts INPUT, STRU and KPT

    Args:
        system (str): system name
        template (bool, optional): whether the input script is a template. Defaults to False.

    Returns:
        tuple: identifier of ABACUS input scripts INPUT, STRU and KPT
    """
    if template:
        return "INPUT_"+system, "STRU_"+system, "KPT_"+system
    else:
        return "INPUT", "STRU", "KPT"

def shorten_keywords(keyword: str) -> str:

    """Generate short version of ABACUS input keyword

    Args:
        keyword (str): ABACUS input keyword to shorten

    Returns:
        str: shortened version of keyword
    """
    if keyword in ABBR.keys():
        return ABBR[keyword]
    
    if keyword.count("_") == 0:
        _longest = ""
        for key in ABBR.keys():
            if key in keyword:
                if len(key) > len(_longest):
                    _longest = key
        if _longest == "":
            # means no abbreviation found, in this case, return the first and last letter combination
            if keyword not in IRREDUCIBLE:
                return keyword[0] + keyword[-1] if len(keyword) > 1 else keyword
            else:
                return keyword
        else:
            fragments = keyword.replace(_longest, "_")
            return ABBR[_longest].join([shorten_keywords(fragment) for fragment in fragments.split("_")])
    else:
        return "".join([shorten_keywords(fragment) for fragment in keyword.split("_")])

def calculation(param_suite: dict) -> str:

    result = ""
    for param in param_suite.keys():
        result += shorten_keywords(param) + str(param_suite[param])
    return result

def cif(system_with_mpid: str) -> str:

    """Generate identifier of CIF file

    Args:
        system_with_mpid (str): system name with MPID

    Returns:
        str: identifier of CIF file
    """
    return "mp-" + system_with_mpid.split("_")[-1] + ".cif"

def folder_reduce(folder: str) -> str:
    """Remove some redundant words in folder name to make it shorter,
    but keep it clear enough to identify the calculation."""
    for key in FOLDER_REDUICE.keys():
        folder = folder.replace(key, FOLDER_REDUICE[key])
    return folder

ABBR = {"basis": "bs", "cal": "cl", "ecut": "ec", "force": "fs", "stress": "strs", "cell": "c", "scaling": "scal", "kspacing": "kspc",
        "kpoint": "kpt", "gamma": "gm", "centered": "cen", "pw": "pw", "pbe": "pbe", "pbesol": "pbesol", "lda": "lda", "gga": "gga",
        "noncolin": "nc", "spin": "sp", "relativistic": "fr", "pseudo": "ps", "potential": "pot",
        "band": "bnd", "type": "typ", "function": "fn", "mixing": "mix", "efield": "efield", "block": "blk",
        "method": "msd", "hybrid": "hyb", "threshold": "thr", "kinetic": "kin", "weight": "wt", "height": "ht",
        "smearing": "smr", "bessel": "bsl", "descriptor": "dscrptr", "wannier": "wnr", "grad": "grd",
        "thermostat": "thmst", "factor": "fac", "deepks": "dpks", "relax": "rlx", "nonlocal": "nloc",
        "max": "x", "orbital": "orb", "rcut": "rc", "mesh": "msh", "spin": "sp", "lspinorb": "soc", "bndpar": "bndp", "symmetry": "symm",
        "temperature": "temp", "volume": "vol", "press": "prs", "spline": "spl", "freq": "frq", "wfc": "wf",
        "rho": "ro", "ecutwfc": "ecw", "ecutrho": "ecro", "extrap": "extrp", "proj": "prj", "restart": "rst",
        "seed": "sid", "dipole": "dip", "correction": "corr", "noncolin": "nlcc", "functional": "fnl", "functionals": "xc",
        "dft_functional": "xc", "label": "lb", "model": "mdl", "flag": "flg", "down": "dw",
        "up": "up", "charge": "chg", "damping": "dmp", "alpha": "a", "beta": "b", "gamma_only": "k000",
        "file": "f", "nnkp": "nnkp", "lambda": "lmbd", "tolerance": "tol", "step": "stp", "switch": "", 
        "calculation": ""}
IRREDUCIBLE = ["bs", "cl", "ec", "fs", "strs", "c", "scal", "kspc", "kpt", "gm", "cen", "pw", "pbe", "pbesol", "lda", "gga", "nc", "sp", "fr", "ps", "pot"
               "bnd", "typ", "fn", "mix", "efield", "blk", "msd", "hyb", "thr", "kin", "wt", "ht", "smr", "bsl", "dscrptr", "wnr", "grd", "thmst", "fac", "dpks", "rlx", "nloc",
               "x", "orb", "rc", "msh", "sp", "soc", "bndp", "symm", "temp", "vol", "prs", "spl", "frq", "wf",
               "ro", "ecw", "ecro", "extrp", "prj", "rst", "sid", "dip", "corr", "nlcc", "fnl",
               "xc", "lb", "mdl", "flg", "gate", "dw", "up", "chg", "dmp", "a", "b", "k000",
               "f", "nnkp", "lmbd", "tol", "stp", "exx"]
FOLDER_REDUICE = {"xc": "", "bstyp": "", "cscal": "cell"}