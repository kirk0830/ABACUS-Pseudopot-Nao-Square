import apns.module_nao.orbgen_kernel.siab_new as amngsn

import json
def iterate(elements: list,         # iterate range of elements
            pseudopot_ids: list,    # iterate range of pseudopotential ids
            fecutwfc: str,          # link to ecutwfc convergence test result
            fpseudo_archive: str,   # link to pseudopotential archive, it is always the pseudo_dir + "description.json"
            path_ref: str):         # link to path where reference files are stored
    
    """generate a series of numerical orbital generation input in new version format
    from as-generated SG15 pseudopotential-corresponding numerical orbitals provided
    by ABACUS official"""
    
    # read pseudopotential archive
    with open(fpseudo_archive, "r") as f:
        pspot_db = json.load(f)
    with open(fecutwfc, "r") as f:
        ecutwfc_db = json.load(f)

    for pseudopot_id in pseudopot_ids:
        pseudo_dir = pspot_db[pseudopot_id]
        with open(pseudo_dir + "/description.json", "r") as f:
            fpseudos = json.load(f)
        for element in elements:
            fpseudo = fpseudos[element]
            ecutwfc = 100
            shortenkey = pseudopot_id.replace("_", "")
            if shortenkey in ecutwfc_db[element].keys():
                ecutwfc = ecutwfc_db[element][shortenkey]
            