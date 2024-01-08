"""
Interface to external files

Author: Kirk0830
Date: 2023-11-22

Description:
    This file contains functions that can read and write files.

Functions:
    check_folder_availability(folder: str) -> None
    parse_cif(cif_file: str, save_json = False) -> tuple[dict, dict]
"""
import numpy as np
import os
import json

"""basic functions"""
def check_folder_availability(folder: str) -> None:
    """
    check if the folder exists, if not, create one
    """
    if os.path.exists(folder):
        print("Folder check: %s, pass."%folder)
    else:
        print("Folder check: %s, not found, exit."%folder)
        exit(1)

def cellparam_to_latvec(cell_parameters: dict):

    vector_a = cell_parameters["a"] * np.array([1, 0, 0])
    vector_b = cell_parameters["b"] * np.array([np.cos(np.deg2rad(cell_parameters["gamma"])), np.sin(np.deg2rad(cell_parameters["gamma"])), 0])
    vector_c_1 = np.cos(np.deg2rad(cell_parameters["beta"]))
    vector_c_2 = (np.cos(np.deg2rad(cell_parameters["alpha"])) - np.cos(np.deg2rad(cell_parameters["gamma"])) * np.cos(np.deg2rad(cell_parameters["beta"]))) / np.sin(np.deg2rad(cell_parameters["gamma"]))
    vector_c_3 = np.sqrt(1 - vector_c_1**2 - vector_c_2**2)
    vector_c = cell_parameters["c"] * np.array([vector_c_1, vector_c_2, vector_c_3])
    lattice_vectors = [
        vector_a.tolist(),
        vector_b.tolist(),
        vector_c.tolist(),
    ]
    return lattice_vectors

def to_json(data: dict, json_file: str):
    with open(json_file, 'w') as f:
        json.dump(data, f, indent = 4)

"""parse functions"""

"""Crsytal information file parser.

This module is used to parse the crystal information file (CIF) and return the
crystal information in dictionaries.

Methodology:
Method 1:
For Materials project CIF files, the information is stored in a standard way.
Method 2:
CIF files format varies from different software.
This module ONLY reads information starts everytime loop_ appears.

loop_ logic:
in loop_, if there is a line starts with "_", it means the title of the
information, if the next line also starts with "_", it means the information
will behind
"""
def read(fname: str, method: int = 1) -> dict:
    if method == 1:
        cif = read_1(fname)
    elif method == 2:
        cif = read_2(fname)
    else:
        raise NotImplementedError
    return cif

def read_1(fname: str, save: bool = False) -> dict:
    cif = {
        "cell_parameters": {
            "a": 0, "b": 0, "c": 0,
            "alpha": 0, "beta": 0, "gamma": 0,
        },
        "atomic_positions": {},
        "atomic_occupancies": {},
    }
    with open(fname, 'r') as f:
        cif_contents = f.readlines()
    read_atomic_position = False
    for line in cif_contents:
        clean_line = line.strip()
        # import lattice parameters information
        if clean_line.startswith("_cell_length_a"):
            cif["cell_parameters"]["a"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_length_b"):
            cif["cell_parameters"]["b"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_length_c"):
            cif["cell_parameters"]["c"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_angle_alpha"):
            cif["cell_parameters"]["alpha"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_angle_beta"):
            cif["cell_parameters"]["beta"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_angle_gamma"):
            cif["cell_parameters"]["gamma"] = float(clean_line.split()[1])
        # import atomic position information
        elif clean_line.startswith("_atom_site_occupancy"):
            read_atomic_position = True
            cif["atomic_positions"] = {}
        elif read_atomic_position:
            if clean_line.startswith("_"):
                read_atomic_position = False
            else:
                element = clean_line.split()[0]
                # if there is any digit and symbol in element, remove
                element = ''.join([i for i in element if i.isalpha()])
                if element not in cif["atomic_positions"].keys():
                    cif["atomic_positions"][element] = []
                    cif["atomic_occupancies"][element] = []
                cif["atomic_positions"][element].append(
                    [float(clean_line.split()[3]), float(clean_line.split()[4]), float(clean_line.split()[5])]
                )
                cif["atomic_occupancies"][element].append(
                    float(clean_line.split()[6])
                )

    if save:
        jsfname = fname.split(".")[0] + ".json"
        to_json(cif, jsfname)
        print("CIF file original information %s has been converted to JSON file %s."%(fname, jsfname))
    return cif

def _method1_atomic_position(cif_contents: dict) -> dict:
    """Parse the atomic position information.

    Args:
        cif_contents: The contents of the CIF file.

    Returns:
        A dictionary that contains the atomic position information.
    """
    atomic_positions = {}
    for atom in cif_contents["atomic_position"].keys(): # O1, Al1 or something like that
        atom_info = cif_contents["atomic_position"][atom]
        if atom_info["element"] not in atomic_positions.keys():
            atomic_positions[atom_info["element"]] = {
                "n": 1,
                "mag": 0.0,
                "positions": [
                    [atom_info["x"], atom_info["y"], atom_info["z"]],
                ]
            }
        else:
            atomic_positions[atom_info["element"]]["n"] += 1
            atomic_positions[atom_info["element"]]["positions"].append(
                [atom_info["x"], atom_info["y"], atom_info["z"]],
            )

    return atomic_positions

def _method1_lattice(cif_contents: dict) -> dict:
    """Parse the lattice information.

    Args:
        cif_contents: The contents of the CIF file.

    Returns:
        A dictionary that contains the lattice information.
    """
    lattice = {}
    lattice_vectors = cellparam_to_latvec(cif_contents["lattice_parameters"])
    lattice["lattice_vectors"] = lattice_vectors
    lattice["lattice_parameters"] = cif_contents["lattice_parameters"]
    return lattice

def read_2(fname: str):

    with open(fname, 'r') as f:
        cif_contents = f.readlines()
    # concatenate the lines into a string
    concat_string = ""
    for line in cif_contents:
        if line.startswith("#"):
            continue # skip comments
        concat_string += line
    concat_string = concat_string.replace("\n", " ")
    # seperate the string into different loop_ sections
    loop_sections = concat_string.split("loop_")
    # parse each loop_ section
    cifs = []
    for section in loop_sections:
        cifs.append(_method2_parse_loop(section))
    return cifs

def _method2_parse_loop(concat_string: str) -> dict:
    """Parse the information after loop_.

    Args:
        concat_string: The string that contains the information in loop_.

    Returns:
        A dictionary that contains the information in loop_.
    """
    scattered_information = concat_string.split()
    if not scattered_information[0].startswith("_"):
        return # it is not a standard loop_ section, may be the super title
    # other normal case
    # 1. recombine contents in quotes into one item

    within_quotes = False
    content_in_quotes = ""
    recombined_information = []
    for word in scattered_information:
        if word.startswith("\'") or word.startswith("\"") or word.startswith(";"):
            within_quotes = True
            content_in_quotes += word.replace("\'", "").replace("\"", "").replace(";", "") + " "
            continue
        
        if word.endswith("\'") or word.endswith("\"") or word.endswith(";"):
            within_quotes = False
            content_in_quotes += word.replace("\'", "").replace("\"", "").replace(";", "")
            recombined_information.append(content_in_quotes)
            content_in_quotes = ""
            continue
            
        if within_quotes:
            content_in_quotes += " "+word
        else:
            recombined_information.append(word)
    #print(recombined_information)
    # 2. identify cases between key-value and key-key-value-value
    attributes = [] # attributes of lines in order
    for word in recombined_information:
        if word.startswith("_"):
            attributes.append("key")
        else:
            attributes.append("value")
    #print(attributes)
    # 3. seperate the key-value and key-key-value-value. assume key-value as a special case of key-key-value-value
    #    subloop is defined as a group of key-...-value
    subloops = {}
    index = 0
    cache_attribute = ""
    for ia, attribute in enumerate(attributes):
        if cache_attribute == attribute:
            if attribute == "key":
                subloops[index]["keys"].append(recombined_information[ia])
            else:
                subloops[index]["values"].append(recombined_information[ia])
        else:
            if attribute == "key": # case value-value-...-value-key
                index += 1
                subloops[index] = {
                    "keys": [recombined_information[ia]],
                    "values": [],
                }
            else: # case key-key-...-key-value
                subloops[index]["values"].append(recombined_information[ia])

            cache_attribute = attribute
    #print(subloops)

    # 4. re-organize subloops
    cif = {}
    for index in subloops.keys():
        subloop = subloops[index] # get one subloop
        nkeys = len(subloop["keys"])
        key = ""
        #print(subloop)
        for iv, value in enumerate(subloop["values"]):
            if int(iv/nkeys) == 0:
                key = subloop["keys"][iv%nkeys]
                cif[key] = []
            cif[key].append(value)
    #print(cif)
    return cif
