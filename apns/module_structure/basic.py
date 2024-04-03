def scan_elements(system: str|list) -> list:

    if isinstance(system, str):
        pass
    elif isinstance(system, list):
        system = "_".join(system)
    else:
        raise TypeError("system must be str or list")
    
    # 1. remove all digits
    system = ''.join([i for i in system if not i.isdigit()])
    # 2. remove all special characters
    system = ''.join([i for i in system if i.isalpha()])
    # 3. remove all "dimer", "trimer", "tetramer"
    system = system.replace("dimer", "").replace("trimer", "").replace("tetramer", "")
    # 3. split the system word by uppercase letters
    elements = []
    element = ""
    for index, letter in enumerate(system):
        if letter.isupper():
            if index == 0:
                element += letter
            else:
                elements.append(element)
                element = letter
        else:
            element += letter
    elements.append(element)

    # 4. remove duplicates but keeps order
    elements = list(dict.fromkeys(elements))
    # 5. remove empty strings
    elements = [element for element in elements if element != '']
    return elements

def expand_atomic_species(symbols: list,
                          atomic_positions: list,
                          pseudopotentials: dict, 
                          numerical_orbitals: dict = None,
                          starting_magnetization: list|dict = None, 
                          ):
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
    if len(symbols) != len(atomic_positions):
        raise ValueError("symbols and atomic_positions must have the same length")
    
    result = {}
    if isinstance(starting_magnetization, dict) and len(starting_magnetization) == len(pseudopotentials):
        """means magnetization is specified for each atom type"""
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
                    result[symbol]["numerical_orbitals"] = numerical_orbitals[symbol]
    elif isinstance(starting_magnetization, list) and len(starting_magnetization) == len(symbols):
        """means magnetization is specified for each atom"""

        """list of tuples like (symbol, magnetization)"""
        magnetic_species = []
        for ia in range(len(symbols)):
            if (symbols[ia], starting_magnetization[ia]) not in magnetic_species:
                magnetic_species.append((symbols[ia], starting_magnetization[ia]))
        for ia in range(len(symbols)):
            index = magnetic_species.index((symbols[ia], starting_magnetization[ia]))
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
                    result[species]["numerical_orbitals"] = numerical_orbitals[symbols[ia]]
    else:
        raise TypeError("starting_magnetization must be a dict or a list. For dict, it is specified atom-type-wise. For list, it is specified atom-wise.")
    
    return result

import os
import apns.module_workflow.identifier as amwi
import json
def init_magmom(structure: str, magnetism: str = "default") -> list|None:
    """generate starting magnetization for a given structure
    
    Args:
        structure (str): structure name
        magnetism (str, optional): magnetism type. Defaults to "default", can be "default", "ferromagnetic", "antiferromagnetic", "materials_project".

    Note:
        "materials_project" is only available for CIF files. It will read the magnetization from the file "apns_cache/mpid_magmom.json".
        "default" means no magnetization.
        "ferromagnetic" means all spins are up.
        "antiferromagnetic" means spins are up and down alternatively, but not strictly possible for odd number of atoms.
        
    Returns:
        list|None: starting magnetization
    """
    if structure.endswith(".cif") or structure.endswith(".CIF") or structure.startswith("mp-"):
        if magnetism == "default":
            return None
        if magnetism != "materials_project":
            print("Currently only materials_project is supported for Materials Downloaded CIF files.")
            return None
        mpid = structure.replace(".cif", "").replace(".CIF", "")
        fmagmom = amwi.TEMPORARY_FOLDER + "/mpid_magmom.json"
        if not os.path.exists(fmagmom):
            print("Warning: the file 'mpid_magmom.json' does not exist, magnetization is set to zero. Create an empty file.")
            with open(fmagmom, "w") as f:
                json.dump({}, f)
            return None
        with open(fmagmom, "r") as f:
            magmom = json.load(f)
        if mpid in magmom.keys():
            return magmom[mpid]
        else:
            return None
    
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

if __name__ == "__main__":
    
    print(scan_elements("Er2O3"))