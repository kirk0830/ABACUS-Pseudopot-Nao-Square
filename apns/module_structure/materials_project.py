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
import apns.module_workflow.identifier as amwi
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
    return fname in os.listdir(amwi.TEMPORARY_FOLDER)

def recall_mpids_memory(formuli: str, n_structure: int) -> list:
    """If there is a file named formula_mpid.json, 
    and the entry with key `formuli` has more than `n_structure` mpids,
    return the first `n_structure` mpids. Otherwise return a empty list."""
    if not os.path.exists(amwi.TEMPORARY_FOLDER + "/" + "formula_mpid.json"):
        return []
    else:
        with open(amwi.TEMPORARY_FOLDER + "/" + "formula_mpid.json", "r") as f:
            mpids_memory = json.load(f)
        if formuli not in mpids_memory.keys():
            return []
        elif len(mpids_memory[formuli]) >= n_structure:
            return mpids_memory[formuli][:n_structure]
        else:
            return []

def recall_magmoms_memory(mpids: list) -> list:
    """with mpids, get magmom information. Only when all mpids have magmom information,
    return magmom list. Otherwise return a empty list."""
    if not os.path.exists(amwi.TEMPORARY_FOLDER + "/" + "mpid_magmom.json"):
        return []
    else:
        with open(amwi.TEMPORARY_FOLDER + "/" + "mpid_magmom.json", "r") as f:
            magmoms_memory = json.load(f)
        magmoms = []
        for mpid in mpids:
            if mpid not in magmoms_memory.keys():
                return []
            else:
                magmoms.append(magmoms_memory[mpid])
        return magmoms

def recall_memory(formuli: str, n_structure: int, consider_magnetism: bool = False) -> list:
    """If there is a file named formula_mpid.json, 
    and the entry with key `formuli` has more than `n_structure` mpids,
    return the first `n_structure` mpids. Otherwise return a empty list."""
    mpids = recall_mpids_memory(formuli, n_structure)
    magmoms = []
    if len(mpids) == n_structure:
        if consider_magnetism:
            magmoms = recall_magmoms_memory(mpids)
    return mpids, magmoms

def search(mpr: MPRester, 
           structure_criteria: dict = None,
           with_magmom: bool = False) -> tuple[list, list]:
    """with MPRester handle, return all mpids and corresponding magnetism information in lists
    """
    docs = mpr.materials.summary.search(formula=structure_criteria["formuli"], 
                                        theoretical=structure_criteria["theoretical"], 
                                        is_stable=structure_criteria["is_stable"])
    mpi_ids = [doc.material_id for doc in docs]
    if not with_magmom:
        return mpi_ids, []
    
    docs = mpr.magnetism.search(material_ids=mpi_ids)
    magmoms = [doc.magmoms for doc in docs]
    return mpi_ids, magmoms

def memorize(formuli: str, mpids: list, magmoms: list) -> None:
    """If there is not a file named formula_mpid.json, 
    create one and write the first entry. If there is already a file,
    overwrite the file with new entry."""

    """mpid information"""
    if not os.path.exists(amwi.TEMPORARY_FOLDER + "/" + "formula_mpid.json"):
        with open(amwi.TEMPORARY_FOLDER + "/" + "formula_mpid.json", "w") as f:
            json.dump({formuli: mpids}, f, indent=4)
    else:
        with open(amwi.TEMPORARY_FOLDER + "/" + "formula_mpid.json", "r") as f:
            memory = json.load(f)
        memory[formuli] = mpids
        with open(amwi.TEMPORARY_FOLDER + "/" + "formula_mpid.json", "w") as f:
            json.dump(memory, f, indent=4)
            
    """magnetism information"""
    if len(magmoms) == 0:
        return
    elif len(mpids) != len(magmoms):
        raise ValueError("mpids and magmoms must have the same length")
    if not os.path.exists(amwi.TEMPORARY_FOLDER + "/" + "mpid_magmom.json"):
        with open(amwi.TEMPORARY_FOLDER + "/" + "mpid_magmom.json", "w") as f:
            json.dump(dict(zip(mpids, magmoms)), f, indent=4)
    else:
        with open(amwi.TEMPORARY_FOLDER + "/" + "mpid_magmom.json", "r") as f:
            memory = json.load(f)
        memory.update(dict(zip(mpids, magmoms)))
        with open(amwi.TEMPORARY_FOLDER + "/" + "mpid_magmom.json", "w") as f:
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
            os.system("mv {} {}".format(fname, amwi.TEMPORARY_FOLDER))

        num_cif_downloaded += 1
        cif_filenames.append(fname.replace(".cif", ""))
        if num_cif_downloaded == num_cif:
            break

    return {element: cif_filenames}

import re
import apns.module_structure.CifParser_Pymatgen as amscif
def composites(api_key: str, 
               formula: list, 
               num_cif: int|list, 
               theoretical: bool = False, 
               is_stable: bool = False,
               consider_magnetism: bool = False):
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
    
    if len(formula) == 0:
        raise ValueError("empty request to Materials Project")

    result = {}
    mpid_pattern = r"mp-\d+"
    for fomuli in formula:
        mpid_match = re.match(mpid_pattern, fomuli)
        if mpid_match:
            with open(amwi.TEMPORARY_FOLDER + "/" + fomuli + ".cif", "r") as f:
                chemical_formula = f.readlines()[1].strip().split("_")[-1]
            system_mpid = chemical_formula + "_" + fomuli.replace("mp-", "")
            result.setdefault(chemical_formula, []).append(system_mpid)
            formula.remove(fomuli)
    if len(formula) == 0:
        print(result)
        return result
    
    print("Establishing connection to Materials Project database...")
    with MPRester(api_key) as mpr:
        
        """for each formula, search mpids and magnetism information"""
        for ifo, formuli in enumerate(formula):
            # first check whether there is a memory file
            mpi_ids = []

            mpids, magmoms = recall_memory(formuli, num_cif[ifo], consider_magnetism)
            if len(mpids) == num_cif[ifo] and (len(magmoms) == num_cif[ifo] or not consider_magnetism):
                mpi_ids = mpids
                print("Found {} system's Materials Project IDs (mpids) in memory file, skip searching. ".format(num_cif[ifo])
                     +"If not satisfied with result, delete file formula_mpid.json and restart.")
            else:
                structure_criteria = dict(zip(["formuli", "theoretical", "is_stable"], [formuli, theoretical, is_stable]))
                mpi_ids, magmoms = search(mpr, structure_criteria, with_magmom=consider_magnetism)

                while len(mpi_ids) == 0 and (is_stable or not theoretical):
                    print("Warning: No structure found for formula %s for filter theoretical = %s and is_stable = %s"%(formuli, theoretical, is_stable))
                    if is_stable:
                        print("Automatically switch to is_stable = False, this operation will include more substable structures.")
                        structure_criteria["is_stable"] = False
                        mpi_ids, magmoms = search(mpr, structure_criteria, with_magmom=consider_magnetism)
                        is_stable = False
                    elif not theoretical:
                        print("Automatically switch to theoretical = True, this operation will include also not experimentally observed structures.")
                        structure_criteria["theoretical"] = True
                        mpi_ids, magmoms = search(mpr, structure_criteria, with_magmom=consider_magnetism)
                        theoretical = True
                    else:
                        raise ValueError("No structure found for formula {}".format(formuli))
                    
                memorize(formuli, mpi_ids, magmoms)

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
                    os.system("mv {} {}".format(fname, amwi.TEMPORARY_FOLDER))
                # save the cif file as mp-amwi.cif
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
            os.system("mv {} {}".format(fname, amwi.TEMPORARY_FOLDER))
            
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
    