"""this interface is for parsing Cif file download from Materials Project"""
from pymatgen.io.cif import CifParser

def structure(fname: str, as_dict: bool = False):

    parser = CifParser(fname)
    structure = parser.as_dict()
    key = list(structure.keys())[0]
    xs = structure[key]["_atom_site_fract_x"]
    ys = structure[key]["_atom_site_fract_y"]
    zs = structure[key]["_atom_site_fract_z"]
    elements = structure[key]["_atom_site_type_symbol"]

    if as_dict:
        """visit atom by atom the elements, xs, ys and zs, if element has been already in
        the dict, append the position to the list, otherwise, create a new key-value pair"""
        structure_dict = {}
        for i in range(len(elements)):
            if elements[i] in structure_dict:
                structure_dict[elements[i]].append([float(xs[i]), float(ys[i]), float(zs[i])])
            else:
                structure_dict[elements[i]] = [[float(xs[i]), float(ys[i]), float(zs[i])]]
        return structure_dict
    else:
        """directly return lists of elements, xyz positions"""
        return elements, [[float(xs[i]), float(ys[i]), float(zs[i])] for i in range(len(xs))]

import numpy as np
def lattice(fname: str):
    
    parser = CifParser(fname)
    result = parser.as_dict()
    a = result[list(result.keys())[0]]["_cell_length_a"]
    b = result[list(result.keys())[0]]["_cell_length_b"]
    c = result[list(result.keys())[0]]["_cell_length_c"]
    alpha = result[list(result.keys())[0]]["_cell_angle_alpha"]
    beta = result[list(result.keys())[0]]["_cell_angle_beta"]
    gamma = result[list(result.keys())[0]]["_cell_angle_gamma"]

    a = float(a)
    b = float(b)
    c = float(c)
    alpha = float(alpha)
    beta = float(beta)
    gamma = float(gamma)

    vec_a = a*np.array([1, 0, 0])
    vec_b = b*np.array([np.cos(np.radians(gamma)), np.sin(np.radians(gamma)), 0])
    vec_ci = np.cos(np.radians(beta))
    vec_cj = (np.cos(np.radians(alpha))-np.cos(np.radians(beta))*np.cos(np.radians(gamma)))/np.sin(np.radians(gamma))
    vec_ck = np.sqrt(1-vec_ci**2-vec_cj**2)
    vec_c = c*np.array([vec_ci, vec_cj, vec_ck])
    lattice_vectors = [vec_a.tolist(), vec_b.tolist(), vec_c.tolist()]
    #lattice_vectors = lattice.as_dict()["matrix"]
    return a, b, c, alpha, beta, gamma, lattice_vectors