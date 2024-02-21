import json

import apns.module_workflow.identifier as amwi
def find_ecutwfc(element: str,
                 ecutwfc: float|list = None,
                 pspot_id: str|list = None) -> tuple[list[float], list[str]]:
    
    """preset the ecutwfc and pspot_id for the given element for different cases"""
    fecutwfc = amwi.TEMPORARY_FOLDER + "/ecutwfc_convergence.json"
    if ecutwfc is None and pspot_id is None:
        print("Use all values for one element record in the ecutwfc_convergence.json file")
        with open(fecutwfc, "r") as f:
            ecutwfc_db = json.load(f)[element]
        ecutwfc = list(ecutwfc_db.values())
        pspot_id = list(ecutwfc_db.keys())
    elif ecutwfc is None and not pspot_id is None:
        print("Use ecutwfc record in the ecutwfc_convergence.json file for the given pspot_id")
        with open(fecutwfc, "r") as f:
            ecutwfc_db = json.load(f)[element]
        if isinstance(pspot_id, str):
            ecutwfc = [ecutwfc_db[pspot_id.replace("_", "")]]
            pspot_id = [pspot_id]
        elif isinstance(pspot_id, list):
            ecutwfc = [ecutwfc_db[pspot_id[i].replace("_", "")] for i in range(len(pspot_id))]
    elif not ecutwfc is None and pspot_id is None:
        if isinstance(ecutwfc, float) or isinstance(ecutwfc, int):
            print("Use the given ecutwfc for all pspot_id")
            with open(fecutwfc, "r") as f:
                ecutwfc_db = json.load(f)[element]
            pspot_id = list(ecutwfc_db.keys())
            ecutwfc = [ecutwfc] * len(pspot_id)
        elif isinstance(ecutwfc, list):
            with open(fecutwfc, "r") as f:
                ecutwfc_db = json.load(f)[element]
            pspot_id = list(ecutwfc_db.keys())
            if len(ecutwfc) != len(pspot_id):
                raise ValueError("ecutwfc and pspot_id must have the same length")
    # if they both are not None:
    elif (isinstance(ecutwfc, float) or isinstance(ecutwfc, int)) and isinstance(pspot_id, str):
        print(f"One single task: generate {element} pseudopotential with ecutwfc = {ecutwfc} Ry and pspot_id = {pspot_id}")
        pspot_id = [pspot_id]
        ecutwfc = [ecutwfc]
    elif isinstance(ecutwfc, list) and isinstance(pspot_id, list):
        print("Multiple tasks, one-to-one correspondence must hold in this case")
        if len(ecutwfc) != len(pspot_id):
            raise ValueError("ecutwfc and pspot_id must have the same length")
    elif isinstance(ecutwfc, list) and isinstance(pspot_id, str):
        print("Cannot handle this case")
        raise ValueError("Confusing input: multiple ecutwfc with one pspot_id")
    elif isinstance(ecutwfc, float) and isinstance(pspot_id, list):
        print("Use uniform ecutwfc for all pspot_id")
        ecutwfc = [ecutwfc] * len(pspot_id)
    else:
        raise ValueError("Confusing input")

    return ecutwfc, pspot_id

def find_fpseudo(element: str,
                 pspot_id: str,
                 pseudo_dir: str = "./download/pseudopotentials/") -> tuple[str, str]:
    """find the pseudopotential file for the given element and pspot_id"""
    pseudo_dir = pseudo_dir.replace("\\", "/")
    if not pseudo_dir.endswith("/"):
        pseudo_dir += "/"
    with open(pseudo_dir + "description.json", "r") as f:
        pseudo_db = json.load(f)
        keys = list(pseudo_db.keys())
        for key in keys:
            if key.replace("_", "") == pspot_id:
                pseudo_dir = pseudo_db[key]
                break
    pseudo_dir = pseudo_dir.replace("\\", "/")
    if not pseudo_dir.endswith("/"):
        pseudo_dir += "/"
    with open(pseudo_dir + "description.json", "r") as f:
        pseudo_db = json.load(f)
        print(list(pseudo_db["files"].keys()))
        pseudo_name = pseudo_db["files"][element]
    return pseudo_dir, pseudo_name

import apns.module_nao.orbgen_kernel.siab_new as amnoksn
def siab_generator(element: str,
                   rcuts: list = [6, 7, 8, 9, 10],
                   ecutwfc: float|list = None,
                   pspot_id: str|list = None,
                   pseudo_dir: str = "./download/pseudopotentials/"):
    
    # list   list, with the same length, one-to-one correspondence
    ecutwfc, pspot_id = find_ecutwfc(element, ecutwfc, pspot_id)
    # list[tuple[str, str]], the first is pseudo_dir, the second is pseudo_name
    pspot_pack = [find_fpseudo(element, pspot_id[i], pseudo_dir) for i in range(len(pspot_id))]

    for i in range(len(pspot_id)):
        yield amnoksn.generate(rcuts=rcuts,
                               ecutwfc=ecutwfc[i],
                               pseudo_name=pspot_pack[i][1],
                               pseudo_dir=pspot_pack[i][0],
                               fref="./apns/module_nao/orbgen_kernel/data/"+element+".SIAB_INPUT")