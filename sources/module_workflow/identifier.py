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
    if type == "":
        raise ValueError("type cannot be empty")
    if rcut != 0:
        rcut = str(rcut)

    return "_".join(i for i in [type, rcut, appendix] if i != "")

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

if __name__ == "__main__":

    print(folder("pbe", "Er2O3", "sg1510pd04"))