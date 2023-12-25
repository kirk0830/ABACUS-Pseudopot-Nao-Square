import os
import re

"""This is for archiving pseudopotentials element-wise."""
NORM_CONSERVING_FOLDERS = []
ULTRASOFT_FOLDERS = []
PAW_FOLDERS = []

def _scan_(pseudo_dir: str) -> dict:
    """Returns a dictionary of pseudopotentials, element-wise.

    Args:
        pseudo_dir (str): path to the directory containing pseudopotentials

    Returns:
        dict: a dictionary of pseudopotentials, element-wise. Values are lists of pseudopotentials including path.
    """
    result = {}

    for dir in list(os.walk(pseudo_dir)):
        for pseudopotential in dir[2]:
            match = re.match(r"^([A-Za-z]{1,2})([-._]?.*)(.upf)$", pseudopotential, re.IGNORECASE)
            if match:
                element = match.group(1)[0].upper() + match.group(1)[1:].lower()
                if element not in result:
                    result[element] = []
                _dir = os.path.abspath(dir[0])
                result[element].append(
                    _dir+('\\' if _dir.count('\\') > 0 else _dir+'/')+pseudopotential
                )
    return result

print(_scan_(
    "./download/pseudopotentials/"
))