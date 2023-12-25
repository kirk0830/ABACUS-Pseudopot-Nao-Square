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
    result = ""
    if kind == "":
        raise ValueError("kind cannot be empty")
    result += kind
    if version != "":
        result += "_" + version
    if appendix != "":
        result += "_" + appendix
    return result

def numerical_orbital(type: str, rcut: int, appendix: str) -> str:
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
    result = ""
    if type == "":
        raise ValueError("type cannot be empty")
    result += type
    if rcut != 0:
        result += "_" + str(rcut)
    if appendix != "":
        result += "_" + appendix
    return result

def folder(functional: str, system: str, specific_test: str) -> str:

    functional = functional.lower()
    folder_name_fragments = ["t", functional, system, specific_test]
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
    result += system + ".in"
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
        return "INPUT", "STRU_"+system, "KPT_"+system
    else:
        return "INPUT", "STRU", "KPT"

if __name__ == "__main__":

    print(folder("pbe", "Er2O3", "sg1510pd04"))