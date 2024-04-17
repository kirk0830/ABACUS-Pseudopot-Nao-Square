"""this interface is for parsing Cif file download from Materials Project

COD database: https://www.crystallography.net/ has been tested and works as well
"""
from pymatgen.io.cif import CifParser

def cast_to_float(x):
    """there are some data in format accompied with erros, like
    0.0000(1), which means the error of the last digit is 1, or
    0.0000(10), which means the error of the last digit is 10.
    """
    result = ""
    for char in x:
        if char.isdigit() or char == "." or char == "-" or char == "+":
            result += char
        else:
            break
    try:
        return float(result)
    except:
        print("Error: input: ", x, "result: ", result)
        return 0

def structure(fname: str, as_dict: bool = False):
    """get elements list and xyz positions list from cif file, if as_dict is True,
    return a dict with elements as keys and xyz positions as values, otherwise, return
    two lists of elements and xyz positions"""
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
                structure_dict[elements[i]].append([cast_to_float(s) for s in (xs[i], ys[i], zs[i])])
            else:
                structure_dict[elements[i]] = [[cast_to_float(s) for s in (xs[i], ys[i], zs[i])]]
        return structure_dict
    else:
        """directly return lists of elements, xyz positions"""
        return elements, [[cast_to_float(s) for s in (xs[i], ys[i], zs[i])] for i in range(len(xs))]

import numpy as np
def lattice(fname: str):
    """get lattice parameters and lattice vectors from cif file, return a tuple of
    (a, b, c, alpha, beta, gamma, lattice_vectors)"""
    parser = CifParser(fname)
    result = parser.as_dict()
    a = result[list(result.keys())[0]]["_cell_length_a"]
    b = result[list(result.keys())[0]]["_cell_length_b"]
    c = result[list(result.keys())[0]]["_cell_length_c"]
    alpha = result[list(result.keys())[0]]["_cell_angle_alpha"]
    beta = result[list(result.keys())[0]]["_cell_angle_beta"]
    gamma = result[list(result.keys())[0]]["_cell_angle_gamma"]

    a, b, c = cast_to_float(a), cast_to_float(b), cast_to_float(c)
    alpha, beta, gamma = cast_to_float(alpha), cast_to_float(beta), cast_to_float(gamma)

    vec_a = a*np.array([1, 0, 0])
    vec_b = b*np.array([np.cos(np.radians(gamma)), np.sin(np.radians(gamma)), 0])
    vec_ci = np.cos(np.radians(beta))
    vec_cj = (np.cos(np.radians(alpha))-np.cos(np.radians(beta))*np.cos(np.radians(gamma)))/np.sin(np.radians(gamma))
    vec_ck = np.sqrt(1-vec_ci**2-vec_cj**2)
    vec_c = c*np.array([vec_ci, vec_cj, vec_ck])
    lattice_vectors = [vec_a.tolist(), vec_b.tolist(), vec_c.tolist()]
    #lattice_vectors = lattice.as_dict()["matrix"]
    return a, b, c, alpha, beta, gamma, lattice_vectors

if __name__ == "__main__":
    import os
    content = """
#------------------------------------------------------------------------------
#$Date: 2023-03-26 11:09:57 +0300 (Sun, 26 Mar 2023) $
#$Revision: 282068 $
#$URL: file:///home/coder/svn-repositories/cod/cif/9/00/90/9009009.cif $
#------------------------------------------------------------------------------
#
# This file is available in the Crystallography Open Database (COD),
# http://www.crystallography.net/. The original data for this entry
# were provided the American Mineralogist Crystal Structure Database,
# http://rruff.geo.arizona.edu/AMS/amcsd.php
#
# The file may be used within the scientific community so long as
# proper attribution is given to the journal article from which the
# data were obtained.
#
data_9009009
loop_
_publ_author_name
'Wyckoff, R. W. G.'
_publ_section_title
;
 Second edition. Interscience Publishers, New York, New York
 Fluorite structure
;
_journal_name_full               'Crystal Structures'
_journal_page_first              239
_journal_page_last               444
_journal_volume                  1
_journal_year                    1963
_chemical_formula_structural     CmO2
_chemical_formula_sum            'Cm O2'
_space_group_IT_number           225
_symmetry_space_group_name_Hall  '-F 4 2 3'
_symmetry_space_group_name_H-M   'F m -3 m'
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                90
_cell_length_a                   5.372
_cell_length_b                   5.372
_cell_length_c                   5.372
_cell_volume                     155.027
_database_code_amcsd             0011687
_exptl_crystal_density_diffrn    11.954
_cod_original_sg_symbol_H-M      'F m 3 m'
_cod_database_code               9009009
_amcsd_formula_title             CmO2
loop_
_space_group_symop_operation_xyz
x,y,z
x,1/2+y,1/2+z
1/2+x,y,1/2+z
1/2+x,1/2+y,z
z,-x,y
z,1/2-x,1/2+y
1/2+z,-x,1/2+y
1/2+z,1/2-x,y
-y,z,-x
-y,1/2+z,1/2-x
1/2-y,z,1/2-x
1/2-y,1/2+z,-x
x,-y,z
x,1/2-y,1/2+z
1/2+x,-y,1/2+z
1/2+x,1/2-y,z
-z,x,-y
-z,1/2+x,1/2-y
1/2-z,x,1/2-y
1/2-z,1/2+x,-y
y,-z,x
y,1/2-z,1/2+x
1/2+y,-z,1/2+x
1/2+y,1/2-z,x
-x,y,-z
-x,1/2+y,1/2-z
1/2-x,y,1/2-z
1/2-x,1/2+y,-z
x,-z,-y
x,1/2-z,1/2-y
1/2+x,-z,1/2-y
1/2+x,1/2-z,-y
-z,y,x
-z,1/2+y,1/2+x
1/2-z,y,1/2+x
1/2-z,1/2+y,x
y,-x,-z
y,1/2-x,1/2-z
1/2+y,-x,1/2-z
1/2+y,1/2-x,-z
-x,z,y
-x,1/2+z,1/2+y
1/2-x,z,1/2+y
1/2-x,1/2+z,y
z,-y,-x
z,1/2-y,1/2-x
1/2+z,-y,1/2-x
1/2+z,1/2-y,-x
-y,x,z
-y,1/2+x,1/2+z
1/2-y,x,1/2+z
1/2-y,1/2+x,z
x,z,y
x,1/2+z,1/2+y
1/2+x,z,1/2+y
1/2+x,1/2+z,y
-z,-y,-x
-z,1/2-y,1/2-x
1/2-z,-y,1/2-x
1/2-z,1/2-y,-x
y,x,z
y,1/2+x,1/2+z
1/2+y,x,1/2+z
1/2+y,1/2+x,z
-x,-z,-y
-x,1/2-z,1/2-y
1/2-x,-z,1/2-y
1/2-x,1/2-z,-y
z,y,x
z,1/2+y,1/2+x
1/2+z,y,1/2+x
1/2+z,1/2+y,x
-y,-x,-z
-y,1/2-x,1/2-z
1/2-y,-x,1/2-z
1/2-y,1/2-x,-z
z,x,-y
z,1/2+x,1/2-y
1/2+z,x,1/2-y
1/2+z,1/2+x,-y
-y,-z,x
-y,1/2-z,1/2+x
1/2-y,-z,1/2+x
1/2-y,1/2-z,x
x,y,-z
x,1/2+y,1/2-z
1/2+x,y,1/2-z
1/2+x,1/2+y,-z
-z,-x,y
-z,1/2-x,1/2+y
1/2-z,-x,1/2+y
1/2-z,1/2-x,y
y,z,-x
y,1/2+z,1/2-x
1/2+y,z,1/2-x
1/2+y,1/2+z,-x
-x,-y,z
-x,1/2-y,1/2+z
1/2-x,-y,1/2+z
1/2-x,1/2-y,z
-z,x,y
-z,1/2+x,1/2+y
1/2-z,x,1/2+y
1/2-z,1/2+x,y
y,-z,-x
y,1/2-z,1/2-x
1/2+y,-z,1/2-x
1/2+y,1/2-z,-x
-x,y,z
-x,1/2+y,1/2+z
1/2-x,y,1/2+z
1/2-x,1/2+y,z
z,-x,-y
z,1/2-x,1/2-y
1/2+z,-x,1/2-y
1/2+z,1/2-x,-y
-y,z,x
-y,1/2+z,1/2+x
1/2-y,z,1/2+x
1/2-y,1/2+z,x
x,-y,-z
x,1/2-y,1/2-z
1/2+x,-y,1/2-z
1/2+x,1/2-y,-z
-x,z,-y
-x,1/2+z,1/2-y
1/2-x,z,1/2-y
1/2-x,1/2+z,-y
z,-y,x
z,1/2-y,1/2+x
1/2+z,-y,1/2+x
1/2+z,1/2-y,x
-y,x,-z
-y,1/2+x,1/2-z
1/2-y,x,1/2-z
1/2-y,1/2+x,-z
x,-z,y
x,1/2-z,1/2+y
1/2+x,-z,1/2+y
1/2+x,1/2-z,y
-z,y,-x
-z,1/2+y,1/2-x
1/2-z,y,1/2-x
1/2-z,1/2+y,-x
y,-x,z
y,1/2-x,1/2+z
1/2+y,-x,1/2+z
1/2+y,1/2-x,z
-x,-z,y
-x,1/2-z,1/2+y
1/2-x,-z,1/2+y
1/2-x,1/2-z,y
z,y,-x
z,1/2+y,1/2-x
1/2+z,y,1/2-x
1/2+z,1/2+y,-x
-y,-x,z
-y,1/2-x,1/2+z
1/2-y,-x,1/2+z
1/2-y,1/2-x,z
x,z,-y
x,1/2+z,1/2-y
1/2+x,z,1/2-y
1/2+x,1/2+z,-y
-z,-y,x
-z,1/2-y,1/2+x
1/2-z,-y,1/2+x
1/2-z,1/2-y,x
y,x,-z
y,1/2+x,1/2-z
1/2+y,x,1/2-z
1/2+y,1/2+x,-z
-z,-x,-y
-z,1/2-x,1/2-y
1/2-z,-x,1/2-y
1/2-z,1/2-x,-y
y,z,x
y,1/2+z,1/2+x
1/2+y,z,1/2+x
1/2+y,1/2+z,x
-x,-y,-z
-x,1/2-y,1/2-z
1/2-x,-y,1/2-z
1/2-x,1/2-y,-z
z,x,y
z,1/2+x,1/2+y
1/2+z,x,1/2+y
1/2+z,1/2+x,y
-y,-z,-x
-y,1/2-z,1/2-x
1/2-y,-z,1/2-x
1/2-y,1/2-z,-x
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Cm 0.00000 0.00000 0.00000
O 0.25000 0.25000 0.25000
loop_
_cod_related_entry_id
_cod_related_entry_database
_cod_related_entry_code
1 AMCSD 0011687

"""
    with open('test.cif', 'w') as f:
        f.write(content)
    parser = CifParser('test.cif')
    os.remove('test.cif')

    structure_ = parser.get_structures()[0]
    print(structure_)
    primitive_ = structure_.get_primitive_structure()
    print(primitive_)

    primitive_.to(filename='test.cif')