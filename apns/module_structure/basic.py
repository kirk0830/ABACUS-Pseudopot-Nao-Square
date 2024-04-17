import re
def scan_elements(system: str|list) -> list:
    """scan elements from a system string or a list of system strings"""
    if isinstance(system, str):
        system = [system]
    elements = []
    for s in system:
        for element in re.findall(r'[A-Z][a-z]*', s):
            if element not in elements and element in ELEMENTS:
                elements.append(element)
    return elements

def expand_atomic_species(symbols: list,
                          atomic_positions: list,
                          pseudopotentials: dict, 
                          numerical_orbitals: dict = None,
                          starting_magnetization: list|dict = None):
    """
    find the real number of atom types for the case magnetization is given.
    return a dict like
    ```python
    {
        "ATOMTYPE1": {
            "pseudopotential": ...
            "numerical_orbitals": ...
            "starting_magnetization": ...
            "atomic_positions": [...]
        }
    }
    ```
    """
    # symbols and atomic_positions must have the same length
    assert len(symbols) == len(atomic_positions)
    # result would be dict
    result = {}
    # if the magmom is given like atom type by atom type, then assign for each atom type directly
    # but it is required magmom has consistency of keys with pseudopotentials
    if isinstance(starting_magnetization, dict):
        """means magnetization is specified for each atom type"""
        assert len(starting_magnetization) == len(pseudopotentials)
        assert set(starting_magnetization.keys()) == set(pseudopotentials.keys())
        # then update the result dict, by either adding new atom type or appending atomic positions
        for ia in range(len(symbols)):
            symbol = symbols[ia]
            if symbol in result:
                result[symbol]["atomic_positions"].append(atomic_positions[ia])
            else:
                result[symbol] = {
                    "pseudopotential": pseudopotentials[symbol],
                    "starting_magnetization": starting_magnetization[symbol],
                    "atomic_positions": [atomic_positions[ia]]
                }
                if numerical_orbitals is not None:
                    result[symbol]["numerical_orbitals"] = numerical_orbitals.get(symbol, None)
    # else if the magmom is given atom by atom, then assign for each atom directly
    elif isinstance(starting_magnetization, list):
        assert len(starting_magnetization) == len(symbols)
        magmom_spec = [(symbols[i], starting_magnetization[i]) for i in range(len(symbols))]
        # delete repeated tuple but keep the order
        magmom_spec = list(dict.fromkeys(magmom_spec))
        # there are cases that same symbol with different magmom, therefore need to distinguish this by
        # adding an index to the symbol
        for ia in range(len(symbols)):
            index = magmom_spec.index((symbols[ia], starting_magnetization[ia]))
            species = symbols[ia] + str(index + 1)
            if species in result:
                result[species]["atomic_positions"].append(atomic_positions[ia])
            else:
                result[species] = {
                    "pseudopotential": pseudopotentials[symbols[ia]],
                    "starting_magnetization": starting_magnetization[ia],
                    "atomic_positions": [atomic_positions[ia]]
                }
                if numerical_orbitals is not None:
                    result[species]["numerical_orbitals"] = numerical_orbitals.get(symbols[ia], None)
    else:
        raise TypeError("starting_magnetization must be a dict or a list. For dict, it is specified atom-type-wise. For list, it is specified atom-wise.")
    
    return result


    
    natoms = dict(zip(["dimer", "trimer", "tetramer"], [2, 3, 4]))
    if structure in natoms.keys():
        if magnetism == "default":
            return [0.0]*natoms[structure]
        elif magnetism == "ferromagnetic":
            return [1.0]*natoms[structure]
        elif magnetism == "antiferromagnetic":
            if natoms[structure] % 2 != 0:
                print("Warning: antiferromagnetic is not strictly possible for odd number of atoms arranged in present pre-set way,",
                      "spin-frustrated/spin-liquid state yield. Therefore, the magnetization is set to zero.")
                return [0.0]*natoms[structure]
            return [1.0, -1.0]*int(natoms[structure]/2)
        else:
            return None

ELEMENTS = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", 
            "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", 
            "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", 
            "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", 
            "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", 
            "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", 
            "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", 
            "Fl", "Mc", "Lv", "Ts", "Og"]

import unittest
class TestStructureBasic(unittest.TestCase):
    def test_scan_elements(self):
        self.assertEqual(scan_elements("Er2O3"), ["Er", "O"])
        self.assertEqual(scan_elements("Er_dimer"), ["Er"])
        self.assertEqual(scan_elements("Er_fcc"), ["Er"])
        self.assertEqual(scan_elements(["Er2O3", "Er_dimer", "Er_fcc"]), ["Er", "O"])
        self.assertEqual(scan_elements(["BaTiO3", "BaPbI3", "SrTiO3"]), ["Ba", "Ti", "O", "Pb", "I", "Sr"])
        self.assertEqual(scan_elements(["Ba2H4Ti8O12F4", "BaPbI3", "SrTiO3"]), ["Ba", "H", "Ti", "O", "F", "Pb", "I", "Sr"])

    def test_expand_atomic_species(self):
        symbols = ["Er", "Er", "O"]
        atomic_positions = [[0, 0, 0], [0, 0, 1], [0, 0, 2]]
        pseudopotentials = {"Er": "Er.psp8", "O": "O.psp8"}
        numerical_orbitals = {"Er": "Er.nao", "O": "O.nao"}
        starting_magnetization = [1.0, 0.0, 0.0] # two Er have different magnetization
        result = expand_atomic_species(symbols, 
                                       atomic_positions, 
                                       pseudopotentials, 
                                       numerical_orbitals, 
                                       starting_magnetization)
        ref = {'Er1': {'pseudopotential': 'Er.psp8', 
                       'starting_magnetization': 1.0, 
                       'atomic_positions': [[0, 0, 0]], 
                       'numerical_orbitals': 'Er.nao'}, 
                'Er2': {'pseudopotential': 'Er.psp8', 
                        'starting_magnetization': 0.0, 
                        'atomic_positions': [[0, 0, 1]], 
                        'numerical_orbitals': 'Er.nao'}, 
                'O3': {'pseudopotential': 'O.psp8', 
                       'starting_magnetization': 0.0,
                        'atomic_positions': [[0, 0, 2]], 
                        'numerical_orbitals': 'O.nao'}}
        self.assertEqual(result, ref)
        # same magmom for two Er
        starting_magnetization = [1.0, 1.0, 0.0]
        result = expand_atomic_species(symbols, 
                                       atomic_positions, 
                                       pseudopotentials, 
                                       numerical_orbitals, 
                                       starting_magnetization)
        ref = {'Er1': {'pseudopotential': 'Er.psp8',
                      'starting_magnetization': 1.0, 
                      'atomic_positions': [[0, 0, 0], [0, 0, 1]], 
                      'numerical_orbitals': 'Er.nao'}, 
                'O2': {'pseudopotential': 'O.psp8', 
                      'starting_magnetization': 0.0, 
                      'atomic_positions': [[0, 0, 2]], 
                      'numerical_orbitals': 'O.nao'}}
        self.assertEqual(result, ref)
        # magmom as dict
        starting_magnetization = {"Er": 1.0, "O": 0.0}
        result = expand_atomic_species(symbols, 
                                       atomic_positions, 
                                       pseudopotentials, 
                                       numerical_orbitals, 
                                       starting_magnetization)
        ref = {'Er': {'pseudopotential': 'Er.psp8',
                      'starting_magnetization': 1.0, 
                      'atomic_positions': [[0, 0, 0], [0, 0, 1]], 
                      'numerical_orbitals': 'Er.nao'}, 
               'O': {'pseudopotential': 'O.psp8', 
                     'starting_magnetization': 0.0, 
                     'atomic_positions': [[0, 0, 2]], 
                     'numerical_orbitals': 'O.nao'}}
        self.assertEqual(result, ref)

if __name__ == "__main__":
    
    unittest.main()