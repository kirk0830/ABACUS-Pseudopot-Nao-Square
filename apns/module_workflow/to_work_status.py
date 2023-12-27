"""From input script to generate initial work_status dictionary"""

import apns.module_io.work_status_expand as wse
import apns.module_structure.materials_project as mp
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

def to(fname: str, num_cif: int = 1) -> dict:
    """ `input.json` -> `work_status` conversion

    Args:
        fname (str): name of json input file
        api_key (str): Materials Project API key, for downloading structures online
        num_cif (int, optional): Number of CIF files download from Materials Project for one specific system. Defaults to 1.

    Returns:
        dict: work_status

    Details:
        1. About structures automatically download from Materials Project
           use `theoretical=False` tag, to filter out all not experimentally-observed structures. Not use `is_stable` tag because
           it filters out too many commonly seen structures. The `is_stable` does not mean filtering out unstable structures,
           but is filter out all not "destination structures". e.g., Anatase TiO2 will transform to Rutile TiO2 at high 
           temperature, but if `is_stable` is `True`, then Anatase will be filtered out.
    """
    with open(fname, "r") as f:
        inp = json.load(f)    

    systems = inp["systems"]
    _structures = mp.composites(api_key=inp["global"]["api_key"],
                                formula=systems,
                                num_cif=inp["n_structure"])
    del inp["global"]["api_key"] # delete api_key from input file
    # _structures is a dict whose keys are system formula and values are lists of corresponding system_mpids.
    return wse.render(fname=fname, system_with_mpids=_structures)