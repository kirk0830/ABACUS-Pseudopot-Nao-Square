"""how to label nao? I also want to know. Need to know what the pseudopotential file is
from the orb file."""
ORBITAL_DIR = "./download/numerical_orbitals"
FDATABASE = ORBITAL_DIR + "/database.json"

import os
import json
import apns.module_nao.parse as amnp
import apns.module_new.tag_search as amds
def initialize(subfolder: str, pptags: list, refresh: bool = False) -> list[str]:
    """initialize will create a database file if not exists, placing all
    numerical atomic orbitals in the same folder would be beneficial, therefore
    the `pptags` can be set easily.
    `pptag` together with `orbtag`, one can uniquely identify the numerical
    atomic orbital.
    
    subfolder: str, the folder where numerical atomic orbitals sharing the same
    `pptags` are placed.
    pptags: list, of pseudopotential, the tags to be added onto ALL the numerical 
    atomic orbitals. One should make sure for pptags specified, with element, only
    one pseudopotential file can be found. If more than one pseudopotential files
    are found, it means numerical atomic orbitals in subfolder correspond to not
    only one kind of pseudopotential, which is unacceptable. The program will 
    raise an error.
    refresh: bool, if True, the database will be updated, otherwise, the database
    will be read and returned.

    return: list, empty.
    """
    if not os.path.exists(FDATABASE):
        with open(FDATABASE, "w") as f:
            json.dump({}, f)
    if not refresh:
        return []
    with open(FDATABASE) as f:
        database = json.load(f)
    subfolder = os.path.abspath(subfolder)
    for root, dirs, files in os.walk(subfolder):
        for file in files:
            if file.lower().endswith(".orb"):
                key = os.path.abspath(os.path.join(root, file))
                orb = amnp.parse(key)
                tags = ["xc", "rcut", "ecutwfc", "nzeta"]
                database.setdefault(key, {"pptag": list(set(pptags + [orb["element"]]))})\
                    .update({"orbtag": [orb.get(tag, None) for tag in tags]})
    
    with open(FDATABASE, "w") as f:
        json.dump(database, f, indent=4)
    return []