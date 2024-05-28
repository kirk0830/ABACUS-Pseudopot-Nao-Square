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

def pseudopot_nao(pseudopotential: list, numerical_orbital: list = None) -> str:

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

def r_folder(identifier: str) -> tuple:

    return tuple(i for i in identifier.split("_") if i != "")

import uuid
def folder(system: str, 
           pseudo_nao_identifier: str, 
           calculation_identifier: str, 
           extensive_identifier: str) -> str:
    
    """new version of Generate identifier of folder

    Args:
        system (str): system name
        pseudo_nao_identifier (str): identifier of pseudo_nao
        calculation_identifier (str): identifier of calculation

    Returns:
        str: identifier of folder
    """
    part1 = "_".join([system, pseudo_nao_identifier])
    part2 = "_".join([calculation_identifier, extensive_identifier])
    part2uuid = uuid.uuid3(uuid.NAMESPACE_DNS, part2).hex
    print(f"""Folder name encoding with uuid.uuid3 (deterministic)
from: {part2} 
to:   {part2uuid}""")
    return "_".join([part1, part2uuid])

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
    if keyword in KEYWORD_ABBR.keys():
        return KEYWORD_ABBR[keyword]
    
    if keyword.count("_") == 0:
        _longest = ""
        for key in KEYWORD_ABBR.keys():
            if key in keyword:
                if len(key) > len(_longest):
                    _longest = key
        if _longest == "":
            # means no abbreviation found, in this case, return the first and last letter combination
            if keyword not in KEYWORD_IRREDUCIBLE:
                return keyword[0] + keyword[-1] if len(keyword) > 1 else keyword
            else:
                return keyword
        else:
            fragments = keyword.replace(_longest, "_")
            return KEYWORD_ABBR[_longest].join([shorten_keywords(fragment) for fragment in fragments.split("_")])
    else:
        return "".join([shorten_keywords(fragment) for fragment in keyword.split("_")])

def calculation(param_suite: dict) -> str:

    result = ""
    for param in param_suite.keys():
        result += shorten_keywords(param).capitalize() + str(param_suite[param])
    return result

def extensive(param_suite: dict) -> str:
    """to render identifier for extensive settings. This function can only be programmed in a 
    case-by-case manner, because some of the extensive settings are global and some are local.
    to iterate. 
    The params to iterate:
    - characteristic_lengths
    However, if only one value is given, then it is a global setting, so it will not be added in
    the identifier.
    
    The params global:
    - nkpoints_in_line
    - magnetism
    """
    result = ""
    iterate_keys = ["characteristic_lengths"]
    for iterkey in iterate_keys:
        if iterkey in param_suite.keys():
            result += shorten_keywords(iterkey).capitalize() + str(param_suite[iterkey])
    return result

def cif(system_with_mpid: str) -> str:

    """Generate identifier of CIF file

    Args:
        system_with_mpid (str): system name with MPID

    Returns:
        str: identifier of CIF file
    """
    return "mp-" + system_with_mpid.split("_")[-1] + ".cif"

import uuid
def folder_reduce(folder: str) -> str:
    """Remove some redundant words in folder name to make it shorter,
    but keep it clear enough to identify the calculation."""
    for key in FOLDER_ABBR.keys():
        folder = folder.replace(key, FOLDER_ABBR[key])
    if len(folder) > 50:
        domains = folder.split("_")
        first_three_domains = "_".join(domains[:3])
        rest = "_".join(domains[3:])
        # use uuid3 because it is deterministic, use the shortest uuid.NAMESPACE
        rest = uuid.uuid3(uuid.NAMESPACE_DNS, rest).hex
        folder = first_three_domains + "_" + rest
        print(f"""WARNING - APNS - Module Workflow - Identifier - folder_reduce:
The length of folder is too long: {len(folder)}.
Will conserve first two domains and replace all the rest to uuid.
Original: {"_".join(domains)}
After: {folder}
""")
    return folder

KEYWORD_ABBR = {"basis": "bs", "cal": "cl", "ecut": "ec", "force": "fs", "stress": "strs", "cell": "c", "scaling": "scal", "kspacing": "kspc",
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
        "calculation": "", "characteristic_lengths": "chlen"}
KEYWORD_IRREDUCIBLE = ["bs", "cl", "ec", "fs", "strs", "c", "scal", "kspc", "kpt", "gm", "cen", "pw", "pbe", "pbesol", "lda", "gga", "nc", "sp", "fr", "ps", "pot"
                       "bnd", "typ", "fn", "mix", "efield", "blk", "msd", "hyb", "thr", "kin", "wt", "ht", "smr", "bsl", "dscrptr", "wnr", "grd", "thmst", "fac", "dpks", "rlx", "nloc",
                       "x", "orb", "rc", "msh", "sp", "soc", "bndp", "symm", "temp", "vol", "prs", "spl", "frq", "wf",
                       "ro", "ecw", "ecro", "extrp", "prj", "rst", "sid", "dip", "corr", "nlcc", "fnl",
                       "xc", "lb", "mdl", "flg", "gate", "dw", "up", "chg", "dmp", "a", "b", "k000",
                       "f", "nnkp", "lmbd", "tol", "stp", "exx"]
FOLDER_ABBR = {"xc": "", "bstyp": "", "cscal": "cell"}
FOLDER_IRREDUCIBLE = ["sp-exc"] # exc is not allowed to be abbreviated due to PD04.sp-exc functional

import unittest
class IdentifierTest(unittest.TestCase):
    """Test the identifier module
    """
    def test_pseudopotential(self):

        self.assertEqual(id.pseudopotential(kind="GGA", version="PBE", appendix=""), "GGA_PBE")
        self.assertEqual(id.pseudopotential(kind="GGA", version="", appendix=""), "GGA")
        self.assertEqual(id.pseudopotential(kind="GGA", version="PBE", appendix="nc"), "GGA_PBE_nc")
        self.assertEqual(id.pseudopotential(kind="GGA", version="", appendix="nc"), "GGA_nc")

    def test_numerical_orbital(self):

        self.assertEqual(id.numerical_orbital(type="DZP", rcut=0, appendix=""), "DZP_0")
        self.assertEqual(id.numerical_orbital(type="DZP", rcut=0, appendix="nc"), "DZP_0_nc")
        self.assertEqual(id.numerical_orbital(type="DZP", rcut=3, appendix="nc"), "DZP_3_nc")
        self.assertEqual(id.numerical_orbital(type="DZP", rcut=3, appendix=""), "DZP_3")
    """
    def test__pseudopotential(self):

        self.assertEqual(id._pseudopotential(identifier="GGA_PBE"), ("GGA", "PBE", ""))
        self.assertEqual(id._pseudopotential(identifier="GGA"), ("GGA", "", ""))
        self.assertEqual(id._pseudopotential(identifier="GGA_PBE_nc"), ("GGA", "PBE", "nc"))
        
    def test__numerical_orbital(self):
            
        self.assertEqual(id._numerical_orbital(identifier="DZP"), ("DZP", "", ""))
        self.assertEqual(id._numerical_orbital(identifier="DZP_3_nc"), ("DZP", "3", "nc"))
        self.assertEqual(id._numerical_orbital(identifier="DZP_3"), ("DZP", "3", ""))
    """
    def test_folder(self):

        self.assertEqual(id.folder(functional="PBE", system="mp-1", specific_test=""), "t_pbe_mp-1_")
        self.assertEqual(id.folder(functional="PBE", system="mp-1", specific_test="test"), "t_pbe_mp-1_test")
        self.assertEqual(id.folder(functional="PBE", system="mp-1", specific_test="test_1"), "t_pbe_mp-1_test_1")

    def test__folder(self):
        pass

    def test_calculation(self):
        self.assertEqual(id.calculation(
            param_suite={
                "functional": "PBE",
                "system": "mp-1",
                "specific_test": "test"
            }
        ), "FnlPBESmmp-1Sctttest")
    
    def test_extensive(self):
        self.assertEqual(id.extensive(param_suite={
            "characteristic_lengths": 2.5,
            "nkpoints_in_line": 10,
            "magnetism": "ferromagnetic"
        }), "Chlen2.5")

if __name__ == "__main__":
    unittest.main()