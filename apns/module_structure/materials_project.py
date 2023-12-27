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

def check_already_exist(mpid: str) -> bool:
    """
    check whether a structure with given mpid already exist in the local database

    Args:

    :param mpid: Materials Project ID

    Returns:

    :return: True if already exist, False if not
    """
    if not mpid.startswith("mp-"):
        mpid = "mp-" + mpid
    fname = mpid + ".cif"
    return fname in os.listdir(id.TEMPORARY_FOLDER)

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

    num_cif_downloaded = 0
    cif_filenames = []

    for mpi_id in mpi_ids:
        fname = element + "_" + mpi_id.replace("mp-", "") + ".cif"
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

    return {element: cif_filenames}

def composites(api_key: str, formula: list, num_cif = 1, theoretical = 0):
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

    result = {}
    with MPRester(api_key) as mpr:
        
        for formuli in formula:
            docs = mpr.materials.summary.search(formula=formuli, theoretical=theoretical)
            mpi_ids = [doc.material_id for doc in docs]

            num_cif_downloaded = 0
            
            result[formuli] = []

            for mpi_id in mpi_ids:
                fname = element + "_" + mpi_id.replace("mp-", "") + ".cif"
                if check_already_exist(mpi_id):
                    print("Structure with mpid {} already exist in cache directory, skip downloading".format(mpi_id))
                else:
                    with MPRester(api_key) as mpr:
                        structure = mpr.get_structure_by_material_id(mpi_id)
                    structure.to(filename=fname, fmt="cif") # better to add symprec = None
                    os.system("mv {} {}".format(fname, id.TEMPORARY_FOLDER))

                num_cif_downloaded += 1
                result[formuli].append(fname.replace(".cif", ""))
                if num_cif_downloaded == num_cif:
                    break
    # then result returns like
    # ["Er2O3": ["Er2O3_2460", "Er2O3_1225560", ...], "TiO2": ["TiO2_XXX", "TiO2_YYY", ...]
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
        fname = element + "_" + mpi_id.replace("mp-", "") + ".cif"
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

    api_key = "wV1HUdmgESPVgSmQj5cc8WvttCO8NTHp"
    elements = ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']
    for element in elements:
        elemental_substances(api_key, element, num_cif = 1)
    