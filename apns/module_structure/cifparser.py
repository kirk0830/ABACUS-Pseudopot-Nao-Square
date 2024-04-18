import numpy as np
import os
import json

class CIFParser:
    def __init__(self, fcif: str):
        self.fcif = fcif


def to_latvec(**kwargs):
    """convert lattice parameters to lattice vectors"""
    # read-in
    a, b, c = kwargs.get("a", 0.0), kwargs.get("b", 0.0), kwargs.get("c", 0.0)
    alpha, beta, gamma = kwargs.get("alpha", 0.0), kwargs.get("beta", 0.0), kwargs.get("gamma", 0.0)
    # convert to radian
    alpha, beta, gamma = np.deg2rad(alpha), np.deg2rad(beta), np.deg2rad(gamma)
    # calculate lattice vectors
    e11, e12, e13 = a, 0, 0
    e21, e22, e23 = b * np.cos(gamma), b * np.sin(gamma), 0
    e31 = c * np.cos(beta)
    e32 = (b * c * np.cos(alpha) - e21 * e31) / e22
    e33 = np.sqrt(c**2 - e31**2 - e32**2)
    return np.reshape([e11, e12, e13, e21, e22, e23, e31, e32, e33], (3, 3)).tolist()

def to_latparam(**kwargs):
    """convert lattice vectors to lattice parameters"""
    lattice_vectors = kwargs.get("lattice_vectors", [])
    a = np.linalg.norm(lattice_vectors[0])
    b = np.linalg.norm(lattice_vectors[1])
    c = np.linalg.norm(lattice_vectors[2])
    alpha = np.arccos(np.dot(lattice_vectors[1], lattice_vectors[2]) / (b * c))
    beta = np.arccos(np.dot(lattice_vectors[0], lattice_vectors[2]) / (a * c))
    gamma = np.arccos(np.dot(lattice_vectors[0], lattice_vectors[1]) / (a * b))
    return {
        "a": a, "b": b, "c": c,
        "alpha": np.rad2deg(alpha), "beta": np.rad2deg(beta), "gamma": np.rad2deg(gamma)
    }

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

def read_1(fcif: str, save: bool = False) -> dict:
    result = dict(zip(
        ["cellparam", "positions", "occupancies"],
        [dict(zip(["a", "b", "c", "alpha", "beta", "gamma"], [0.0]*6)), {}, {}]
    ))
    with open(fcif, 'r') as f:
        contents = f.readlines()
    read_atomic_position = False
    for line in contents:
        clean_line = line.strip()
        # import lattice parameters information
        if clean_line.startswith("_cell_length_a"):
            result["cellparam"]["a"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_length_b"):
            result["cellparam"]["b"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_length_c"):
            result["cellparam"]["c"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_angle_alpha"):
            result["cellparam"]["alpha"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_angle_beta"):
            result["cellparam"]["beta"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_angle_gamma"):
            result["cellparam"]["gamma"] = float(clean_line.split()[1])
        # import atomic position information
        elif clean_line.startswith("_atom_site_occupancy"):
            read_atomic_position = True
            result["positions"] = {}
        elif read_atomic_position:
            if clean_line.startswith("_"):
                read_atomic_position = False
            else:
                element = clean_line.split()[0]
                # if there is any digit and symbol in element, remove
                element = ''.join([i for i in element if i.isalpha()])
                if element not in result["positions"].keys():
                    result["positions"][element] = []
                    result["occupancies"][element] = []
                result["positions"][element].append(
                    [float(clean_line.split()[3]), float(clean_line.split()[4]), float(clean_line.split()[5])]
                )
                result["occupancies"][element].append(
                    float(clean_line.split()[6])
                )

    if save:
        jsfname = fcif.split(".")[0] + ".json"
        to_json(result, jsfname)
        print("CIF file original information %s has been converted to JSON file %s."%(fcif, jsfname))
    return result

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
    lattice_vectors = to_latvec(**cif_contents["lattice_parameters"])
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

def volume_to_celldm(bravis: str, volume: float):
    bravis = bravis.lower()
    assert bravis in ["sc", "bcc", "fcc", "diamond"], "Bravis lattice not supported."
    if bravis == "sc":
        celldm = volume**(1/3)
        return celldm
    elif bravis == "bcc":
        # we can image a cubic cell, in which two atoms in it, the volume would be twice
        # of the primitive cell. The edge length would be then a = (2*v)**(1/3), and the edge
        # length relates with the celldm as np.sqrt(3)/2 * a = celldm
        a = (2*volume)**(1/3)
        celldm = np.sqrt(3)/2 * a
        return celldm
    elif bravis == "fcc":
        # similarly, image a cubic cell, there will be 4 atoms in it, the volume would be 4 times
        # of the primitive cell. The edge length would be then a = (4*v)**(1/3), and the edge
        # length relates with the celldm as np.sqrt(2)/2 * a = celldm
        a = (4*volume)**(1/3)
        celldm = np.sqrt(2)/2 * a
        return celldm
    elif bravis == "diamond":
        # for diamond, it is a little bit special, image again a cubic cell, there will be 8 atoms in it,
        # while in primitive cell there are 2. The volume would be 4 times of the primitive cell.
        # The edge length would be then a = (4*v)**(1/3), and the edge length relates with the celldm as
        # np.sqrt(2)/2 * a = celldm
        a = (4*volume)**(1/3)
        celldm = np.sqrt(2)/2 * a
        return celldm
    else:
        raise NotImplementedError

import unittest
class TestCIF(unittest.TestCase):

    delta = 1e-5
    def test_volume_to_celldm(self):
        self.assertEqual(volume_to_celldm("sc", 1), 1)
        self.assertAlmostEqual(volume_to_celldm("bcc", 20.120701), 2.96771, delta=self.delta)
        self.assertAlmostEqual(volume_to_celldm("fcc", 12.113867), 2.57790, delta=self.delta)
        self.assertAlmostEqual(volume_to_celldm("diamond", 40.888291), 3.86697, delta=self.delta)

    def test_tolatvec_tocellparam(self):
        # test simple cubic
        test = {"a": 1.234567, "b": 1.234567, "c": 1.234567, "alpha": 90, "beta": 90, "gamma": 90}
        latvec = to_latvec(**test)
        cellparam = to_latparam(lattice_vectors=latvec)
        for key in test:
            self.assertAlmostEqual(test[key], cellparam[key], delta=self.delta)
        # body centered cubic
        test = {"a": 2.96771, "b": 2.96771, "c": 2.96771, "alpha": 109.47122, "beta": 109.47122, "gamma": 109.47122}
        latvec = to_latvec(**test)
        cellparam = to_latparam(lattice_vectors=latvec)
        for key in test:
            self.assertAlmostEqual(test[key], cellparam[key], delta=self.delta)
        # face centered cubic
        test = {"a": 2.57790, "b": 2.57790, "c": 2.57790, "alpha": 60, "beta": 60, "gamma": 60}
        latvec = to_latvec(**test)
        cellparam = to_latparam(lattice_vectors=latvec)
        for key in test:
            self.assertAlmostEqual(test[key], cellparam[key], delta=self.delta)
        # arbitrary lattice
        test = {"a": 2.57790, "b": 2.96771, "c": 1.234567, "alpha": 60, "beta": 70, "gamma": 80}
        latvec = to_latvec(**test)
        cellparam = to_latparam(lattice_vectors=latvec)
        for key in test:
            self.assertAlmostEqual(test[key], cellparam[key], delta=self.delta)

if __name__ == "__main__":
    unittest.main()