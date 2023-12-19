import module_io.work_status_expand as wse
import module_structure.materials_project as mp
import json

def merge_dictionaries(dictionaries: list):

    """
    merge a list of dictionaries into one dictionary

    Args:
        dictionaries (list): a list of dictionaries

    Returns:
        dict: merged dictionary
    """
    merged = {}
    for dictionary in dictionaries:
        merged.update(dictionary)
    return merged

def to_work_status(fname: str, api_key: str, num_cif: int = 1) -> dict:

    with open(fname, "r") as f:
        inp = json.load(f)    

    systems = inp["systems"]
    _structures = mp.composites(api_key=api_key,
                                formula=systems,
                                num_cif=num_cif)
    return wse.render(fname=fname, system_mpids=_structures)