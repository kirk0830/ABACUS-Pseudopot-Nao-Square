"""
Online interface to Materials Project database

Author: Kirk0830
Date: 2023-11-22

Description:
    This file contains functions that can download data from Materials Project database.
    website: https://materialsproject.org/

Functions:
    elemental_substances: download single element structures from materials project
    composites: download composites structures from materials project
    binary_composites: download binary composites structures from materials project
    _oxides: download oxides structures from materials project
    _nitrides: download nitrides structures from materials project
    _hydrides: download hydrides structures from materials project

Notes:
    1. search tags theoretical and is_stable.
       is_stable is a much more tight threshold, for example, if is_stable is true,
       for TiO2, only rutile crystal will return, anatase is not interpreted as a
       stable one. (indeed, experimentally the anatase can transform to rutile at
       high temperatures)
    2. output information utilization
       please do understand what will be returned by function you use. In most cases
       a dict whose keys are system formula or simply concatenation of elements, its
       value is file names of cif downloaded from Materials Project. This information
       is useful to expand the `work_status` dictionary, and helpful to identify
       between different tests.
"""

from mp_api.client import MPRester
import apns.module_workflow.identifier as id
import os
import json

def check_already_exist(mpid: str) -> bool:
    """
    check whether a structure with given mpid already exist in the local database CACHE directory

    Args:

    :param mpid: Materials Project ID

    Returns:

    :return: True if already exist, False if not
    """
    if not mpid.startswith("mp-"):
        mpid = "mp-" + mpid
    fname = mpid + ".cif"
    return fname in os.listdir(id.TEMPORARY_FOLDER)

def recall_search_memory(formuli: str, n_structure: int) -> list:
    """If there is a file named materials_project_memory.json, 
    and the entry with key `formuli` has more than `n_structure` mpids,
    return the first `n_structure` mpids. Otherwise return a empty list."""
    if not os.path.exists(id.TEMPORARY_FOLDER + "/" + "materials_project_memory.json"):
        return []
    else:
        with open(id.TEMPORARY_FOLDER + "/" + "materials_project_memory.json", "r") as f:
            memory = json.load(f)
        if formuli not in memory.keys():
            return []
        elif len(memory[formuli]) >= n_structure:
            return memory[formuli][:n_structure]
        else:
            return []

def memorize(formuli: str, mpids: list) -> None:
    """If there is not a file named materials_project_memory.json, 
    create one and write the first entry. If there is already a file,
    overwrite the file with new entry."""    
    if not os.path.exists(id.TEMPORARY_FOLDER + "/" + "materials_project_memory.json"):
        with open(id.TEMPORARY_FOLDER + "/" + "materials_project_memory.json", "w") as f:
            json.dump({formuli: mpids}, f, indent=4)
    else:
        with open(id.TEMPORARY_FOLDER + "/" + "materials_project_memory.json", "r") as f:
            memory = json.load(f)
        memory[formuli] = mpids
        with open(id.TEMPORARY_FOLDER + "/" + "materials_project_memory.json", "w") as f:
            json.dump(memory, f, indent=4)

def elemental_substances(api_key: str, element: str, num_cif = 1, theoretical = 0):
    """
    down load single element structures from materials project

    Args:

    :param api_key: Materials Project API key
    :param element: element name
    :param num_cif: number of cif files to download, default is 1

    Returns:

    :return: cif_filenames: list of cif filenames
    """
    with MPRester(api_key) as mpr:
        
        docs = mpr.materials.summary.search(elements=[element], num_elements=1, theoretical=theoretical)
        
    mpi_ids = [doc.material_id for doc in docs]
    if len(mpi_ids) == 0:
        raise ValueError("No structure found for element {}".format(element))
    
    num_cif_downloaded = 0
    cif_filenames = []

    for mpi_id in mpi_ids:
        fname = element + "_" + mpi_id.replace("mp-", "") + ".cif"
        if check_already_exist(mpi_id):
            print("Structure with mpid {} already exist in cache directory, skip downloading".format(mpi_id))
        else:
            print("Downloading structure with mpid {} from Materials Project...".format(mpi_id))
            with MPRester(api_key) as mpr:
                structure = mpr.get_structure_by_material_id(mpi_id)
            structure.to(filename=fname, fmt="cif") # better to add symprec = None
            os.system("mv {} {}".format(fname, id.TEMPORARY_FOLDER))

        num_cif_downloaded += 1
        cif_filenames.append(fname.replace(".cif", ""))
        if num_cif_downloaded == num_cif:
            break

    return {element: cif_filenames}

def composites(api_key: str, 
               formula: list, 
               num_cif: int|list, 
               theoretical = False, 
               is_stable = False):
    """
    down load composites structures from materials project

    Args:

    :param api_key: Materials Project API key
    :param formula: list of elements in the composite
    :param num_cif: number of cif files to download, default is 1

    Returns:

    :return: a dict whose keys are system formula and values are list of system_mpid, 
             uniquely identify a structure
    """
    if isinstance(num_cif, list):
        if len(num_cif) != len(formula):
            raise ValueError("num_cif must have the same length as formula")
    elif isinstance(num_cif, int):
        num_cif = [num_cif] * len(formula)
    else:
        raise TypeError("num_cif must be a list or an integer")
    
    result = {}
    
    print("Establishing connection to Materials Project database...")
    with MPRester(api_key) as mpr:
        
        for ifo, formuli in enumerate(formula):
            # first check whether there is a memory file
            mpi_ids = []

            mpids = recall_search_memory(formuli, num_cif[ifo])
            if len(mpids) == num_cif[ifo]:
                mpi_ids = mpids
                print("Found {} system's Materials Project IDs (mpids) in memory file, skip searching. ".format(num_cif[ifo])
                     +"If not satisfied with result, delete file materials_project_memory.json and restart.")
            else:
                docs = mpr.materials.summary.search(formula=formuli, theoretical=theoretical, is_stable=is_stable)
                mpi_ids = [doc.material_id for doc in docs]
                while len(mpi_ids) == 0 and (is_stable or not theoretical):
                    print("Warning: No structure found for formula %s for filter theoretical = %s and is_stable = %s"%(formuli, theoretical, is_stable))
                    if is_stable:
                        print("Automatically switch to is_stable = False, this operation will include more substable structures.")
                        docs = mpr.materials.summary.search(formula=formuli, theoretical=theoretical, is_stable=False)
                        mpi_ids = [doc.material_id for doc in docs]
                    elif not theoretical:
                        print("Automatically switch to theoretical = True, this operation will include also not experimentally observed structures.")
                        docs = mpr.materials.summary.search(formula=formuli, theoretical=True, is_stable=is_stable)
                        mpi_ids = [doc.material_id for doc in docs]
                    else:
                        raise ValueError("No structure found for formula {}".format(formuli))
                    
                memorize(formuli, mpi_ids)

            num_cif_downloaded = 0
            
            result[formuli] = []

            for mpi_id in mpi_ids:
                fname = "mp-" + mpi_id.replace("mp-", "") + ".cif"
                
                if check_already_exist(mpi_id):
                    print("Structure with mpid {} already exist in cache directory, skip downloading".format(mpi_id))
                else:
                    with MPRester(api_key) as mpr:
                        structure = mpr.get_structure_by_material_id(mpi_id)
                    structure.to(filename=fname, fmt="cif") # better to add symprec = None
                    os.system("mv {} {}".format(fname, id.TEMPORARY_FOLDER))
                # save the cif file as mp-id.cif
                num_cif_downloaded += 1

                result[formuli].append(formuli + "_" + mpi_id.replace("mp-", "")) # record "system_mpid" under system formula
                if num_cif_downloaded == num_cif[ifo]:
                    break
    # then result returns like
    # ["Er2O3": ["Er2O3_2460", "Er2O3_1225560", ...], "TiO2": ["TiO2_XXX", "TiO2_YYY", ...]
    print("Connection closed.")
    return result

def binary_composites(api_key: str, element1: str, element2: str, num_cif = 1, theoretical = 0):
    """
    down load binary composites structures from materials project

    Args:

    :param api_key: Materials Project API key
    :param element1: element 1 in the composite
    :param element2: element 2 in the composite
    :param num_cif: number of cif files to download, default is 1

    Returns:

    :return: cif_filenames: list of cif filenames
    """
    with MPRester(api_key) as mpr:
        
        docs = mpr.materials.summary.search(elements=[element1, element2], num_elements=2, theoretical=theoretical)
        
    mpi_ids = [doc.material_id for doc in docs]

    num_cif_downloaded = 0
    cif_filenames = []

    for mpi_id in mpi_ids:
        fname = element1 + element2 + "_" + mpi_id.replace("mp-", "") + ".cif"
        if check_already_exist(mpi_id):
            print("Structure with mpid {} already exist in cache directory, skip downloading".format(mpi_id))
        else:
            with MPRester(api_key) as mpr:
                structure = mpr.get_structure_by_material_id(mpi_id)
            structure.to(filename=fname, fmt="cif") # better to add symprec = None
            os.system("mv {} {}".format(fname, id.TEMPORARY_FOLDER))
            
        num_cif_downloaded += 1
        cif_filenames.append(fname.replace(".cif", ""))
        if num_cif_downloaded == num_cif:
            break

    return {element1 + element2: cif_filenames}

def _oxides_(api_key: str, element: str, num_cif = 1, theoretical = 0):
    """
    down load oxides structures from materials project

    Args:

    :param api_key: Materials Project API key
    :param element: element in the oxide
    :param num_cif: number of cif files to download, default is 1

    Returns:

    :return: cif_filenames: list of cif filenames
    """
    return composites(api_key, [element, "O"], num_cif = num_cif, theoretical=theoretical)

def _nitrides_(api_key: str, element: str, num_cif = 1, theoretical = 0):
    """
    down load nitrides structures from materials project

    Args:

    :param api_key: Materials Project API key
    :param element: element in the nitride
    :param num_cif: number of cif files to download, default is 1

    Returns:

    :return: cif_filenames: list of cif filenames
    """
    return composites(api_key, [element, "N"], num_cif = num_cif, theoretical=theoretical)

def _hydrides_(api_key: str, element: str, num_cif = 1, theoretical = 0):
    """
    down load hydrides structures from materials project

    Args:

    :param api_key: Materials Project API key
    :param element: element in the hydride
    :param num_cif: number of cif files to download, default is 1

    Returns:

    :return: cif_filenames: list of cif filenames
    """
    return composites(api_key, [element, "H"], num_cif = num_cif, theoretical=theoretical)

if __name__ == "__main__":

    api_key = ""
    elements = ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']
    for element in elements:
        elemental_substances(api_key, element, num_cif = 1)
    