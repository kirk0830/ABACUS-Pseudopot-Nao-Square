import json
import apns.outdated.identifier as amwi
def find_ecutwfc(element: str,
                 ecutwfc: float|list = None,
                 pspot_id: str|list = None):
    
    """preset the ecutwfc and pspot_id for the given element for different cases"""
    print(f"Find the ecutwfc and pspot_id for the given element {element}")
    print(f"ecutwfc = {ecutwfc}, pspot_id = {pspot_id}")
    # import ecutwfc database
    fecutwfc = amwi.TEMPORARY_FOLDER + "/apns_ecutwfc_db.json"
    with open(fecutwfc, "r") as f:
        ecutwfc_db = json.load(f)[element]

    if ecutwfc is None and pspot_id is None:
        print("Use all values for one element record in the apns_ecutwfc_db.json file")
        ecutwfc = list(ecutwfc_db.values())
        pspot_id = list(ecutwfc_db.keys())
    elif ecutwfc is None and not pspot_id is None:
        print("Use ecutwfc record in the apns_ecutwfc_db.json file for the given pspot_id")
        if isinstance(pspot_id, str):
            if pspot_id.replace("_", "") in ecutwfc_db.keys():
                ecutwfc = [ecutwfc_db[pspot_id.replace("_", "")]]
                pspot_id = [pspot_id]
            else:
                return None, None
        elif isinstance(pspot_id, list):
            pspot_id = [pspot.replace("_", "") for pspot in pspot_id]
            pspot_id = [pspot for pspot in pspot_id if pspot in ecutwfc_db.keys()]
            ecutwfc = [ecutwfc_db[_pspot_id] for _pspot_id in pspot_id]
    elif not ecutwfc is None and pspot_id is None:
        if isinstance(ecutwfc, float) or isinstance(ecutwfc, int):
            print("Use the given ecutwfc for all pspot_id")
            pspot_id = list(ecutwfc_db.keys())
            ecutwfc = [ecutwfc] * len(pspot_id)
        elif isinstance(ecutwfc, list):
            pspot_id = list(ecutwfc_db.keys())
            if len(ecutwfc) != len(pspot_id):
                raise ValueError("ecutwfc and pspot_id must have the same length")
    # if they both are not None:
    elif (isinstance(ecutwfc, float) or isinstance(ecutwfc, int)) and isinstance(pspot_id, str):
        print(f"One single task: generate {element} pseudopotential with ecutwfc = {ecutwfc} Ry and pspot_id = {pspot_id}")
        if pspot_id.replace("_", "") in ecutwfc_db.keys():
            ecutwfc = [ecutwfc]
            pspot_id = [pspot_id]
        else:
            return None, None
    elif isinstance(ecutwfc, list) and isinstance(pspot_id, list):
        print("Multiple tasks, one-to-one correspondence must hold in this case")
        if len(ecutwfc) != len(pspot_id):
            raise ValueError("ecutwfc and pspot_id must have the same length")
    elif isinstance(ecutwfc, list) and isinstance(pspot_id, str):
        print("Cannot handle this case")
        raise ValueError("Confusing input: multiple ecutwfc with one pspot_id")
    elif (isinstance(ecutwfc, float) or isinstance(ecutwfc, int)) and isinstance(pspot_id, list):
        print("Use uniform ecutwfc for all pspot_id")
        pspot_id = [pspot.replace("_", "") for pspot in pspot_id]
        pspot_id_fromdb = [pspot for pspot in pspot_id if pspot in ecutwfc_db.keys()]
        if len(pspot_id) != len(pspot_id_fromdb):
            print("WARNING WARNING WARNING!\nWARNING: some pspot_id(s) not found in the apns_ecutwfc_db.json file")
        ecutwfc = [ecutwfc] * len(pspot_id)

    else:
        raise ValueError(f"Confusing input: ecutwfc = {ecutwfc}, pspot_id = {pspot_id}")

    return dict(zip(pspot_id, ecutwfc))

import apns.outdated.orbgen.orbgen_kernel.siab_new as amnoksn
def siab_generator(**kwargs):
    """
    Generate SIAB_INPUT.json file based on the given settings
    compulsory input: element
    optional input: rcuts, ecutwfc, fpseudos, pseudo_dir, other_settings

    Args:
    element (str): element symbol
    rcuts (list): list of cutoff radius for Sphbes
    ecutwfc (float): energy cutoff for wavefunctions
    fpseudos (dict): dict, keys are pspot_id, values are pseudo_name with path
    other_settings (dict): other settings, can overwrite the reference settings
    """
    # compulsory input
    element = kwargs["element"]
    # optional, with default values
    rcuts = kwargs.get("rcuts", [6, 7, 8, 9, 10])
    ecutwfc = kwargs.get("ecutwfc", None)
    fpseudos = kwargs.get("fpseudos", None)
    other_settings = kwargs.get("other_settings", None)

    old_pspotids = [pspot_id.replace(".", "").replace("sr", "") for pspot_id in fpseudos.keys()]
    # old pspotids are those dots "." are removed, like sg15_1.0 to sg15_10
    ecutwfcs = find_ecutwfc(element, ecutwfc, old_pspotids)
    # example of ecutwfcs:
    # {"sg15_10": 60, "sg15_20": 80, "sg15_30": 100, "sg15_40": 120, "sg15_50": 140}
    print(list(ecutwfcs.keys()))
    for pspot_id in fpseudos.keys():
        old_pspotid = pspot_id.replace(".", "").replace("sr", "").replace("_", "")
        ecutwfc = ecutwfcs[old_pspotid]
        pseudo_dir_name = fpseudos[pspot_id].replace("\\", "/")
        pseudo_dir, pseudo_name = pseudo_dir_name.rsplit("/", 1)
        yield amnoksn.generate(rcuts=rcuts,
                               ecutwfc=ecutwfc,
                               pseudo_name=pseudo_name,
                               pseudo_dir=pseudo_dir,
                               other_settings=other_settings,
                               fref="./apns/module_nao/orbgen_kernel/data/"+element+".SIAB_INPUT")