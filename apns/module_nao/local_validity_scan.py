import os
import apns.module_database.database as ai

def scan_elements(system: str) -> list:

    # 1. remove all digits
    system = ''.join([i for i in system if not i.isdigit()])
    # 2. remove all special characters
    system = ''.join([i for i in system if i.isalpha()])
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

    # 4. remove duplicates
    elements = list(set(elements))
    # 5. remove empty strings
    elements = [element for element in elements if element != '']
    return elements

def scan_valid_numerical_orbitals(work_status: dict, valid_pseudopotentials: dict):
    """
    Scan for valid numerical orbitals for valid pseudopotentials.

    Args:
        work_status (dict): descriptor of test
        valid_pseudopotentials (dict): element-pseudopotential(identifier)-(kind, version, appendix) hierachy-like stored valid pseudopotential information, 
                                       yielded by function `scan_valid_pseudopotential(work_status: dict) -> dict`.  

    Return:
        valid_numerical_orbitals (dict): element-pseudopotential(identifier)-numerical_orbital(identifier)-(type, rcut, appendix) hierachy-like stored valid
                                         numerical orbital information.

    Notes:
    Imagine the numerical orbitals are stored in this way:
    in [work_folder]/numerical_orbitals/resources/

    pd_04_spd/
        33_As_DZP/
            As_gga_6au_100Ry_2s2p1d.orb
            As_gga_7au_100Ry_2s2p1d.orb
            ...
        33_As_TZDP/
            As_gga_6au_100Ry_2s2p1d.orb
            ...
    sg15_10/
        33_As_DZP/
            As_gga_6au_100Ry_2s2p1d.orb
            ...
            As_gga_6au_100Ry_2s2p1d_osci_rm.orb
        ...
    ...

    Output is organized as:
    {
        "As": {                 // element
            "pd_04_spd: {       // identifier of pseudopotential
                "D6": {         // identifier of numerical_orbital
                    "type": "DZP",
                    "rcut": 6,
                    "appendix": "",
                    "file": "As_gga_6au_100Ry_2s2p1d.orb",
                },
                ...
                "T7": {
                    "type": "TZDP",
                    "rcut": 7,
                    "appendix": "",
                    "file": "As_gga_7au_100Ry_3s3p2d.orb",
                },
                ...
            }
            "sg15_10": {
                "D6": {
                    "type": "DZP",
                    "rcut": 6,
                    "appendix": "",
                    "file": "As_gga_6au_100Ry_2s2p1d.orb",
                },
                ...
                "D6_osci_rm": {
                    "type": "DZP",
                    "rcut": 6,
                    "appendix": "osci_rm",
                    "file": "As_gga_6au_100Ry_2s2p1d_osci_rm.orb",
                },
                ...
            }
        },
        ...
    }
    ---
    Additional notes of identifying the numerical orbitals to test  
    2023/12/19  
    As_gga_6au_100Ry_2s2p1d.orb  
    -----------^^^---^^^^^^-----  
    1. note that it is not always that the ecutwfc to generate numerical orbital is 100Ry, but whether it is 100Ry or other value
       should be determined before orbital generation, say the convergence test with `basis_type pw` calculation on reference 
       structure like dimer or trimer should have determined this value.
    2. the orbital configuration is determined by pseudopotential used, say for one large-core pseudopotential, if f-electron is 
       also pseudized and there are two s, one p and one d orbitals are kept as valence orbital, then DZP should be:
       4s2p2d1f, otherwise it is, 4s2p2d2f1g.
       Therefore it is not in aspect of orbital configuration to enumerate, instead it is in aspect of pseudopotential.
    """
    elements = []
    valid_numerical_orbitals = {}
    # collect all elements needed for test
    for system in work_status["systems"].keys():
        for element in scan_elements(system):
            if element not in elements:
                elements.append(element)
    # check-in
    os.chdir(work_status["numerical_orbitals"]["nao_paths"]["resources"])
    print("enter folder: %s" % os.getcwd())
    # after getting valid pseudopotentials, we scan the numerical orbitals according to it
    for element in elements:
        for pseudopotential in valid_pseudopotentials[element].keys(): # for all valid pseudopotentials of present element
            if os.path.isdir(pseudopotential):
                os.chdir(pseudopotential)
                folder_header = str(ai.get_element_index(element)) + "_" + element
                for nao_type in work_status["numerical_orbitals"]["types"][element]:
                    folder = folder_header + "_" + nao_type
                    if os.path.isdir(folder):
                        files = os.listdir(folder)
                        for rcut in work_status["numerical_orbitals"]["rcuts"][element]:
                            for appendix in work_status["numerical_orbitals"]["appendices"][element]:
                                # example: As_gga_10au_100Ry_2s2p1d.orb
                                nao_startswith = "%s_gga_%sau_100Ry" % (element, rcut)
                                nao_endswith = "%s.orb" % appendix
                                nao_valid = False
                                for file in files:
                                    if file.startswith(nao_startswith) and file.endswith(nao_endswith):
                                        nao_valid = True
                                        break
                                if not nao_valid:
                                    continue
                                if element not in valid_numerical_orbitals.keys():
                                    valid_numerical_orbitals[element] = {}
                                if pseudopotential not in valid_numerical_orbitals[element].keys():
                                    valid_numerical_orbitals[element][pseudopotential] = {}
                                key_nao = nao_type[0] + str(rcut)
                                if appendix != "":
                                    key_nao += "_" + appendix
                                valid_numerical_orbitals[element][pseudopotential][key_nao] = {
                                    "type": nao_type,
                                    "rcut": rcut,
                                    "appendix": appendix,
                                    "file": file,
                                }
                os.chdir("..")
    return valid_numerical_orbitals