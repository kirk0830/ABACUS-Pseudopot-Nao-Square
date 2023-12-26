import os
import module_database.database as ai
import module_pseudo.archive as arch

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

def scan_valid_pseudopotentials(work_status: dict):
    """
    Record valid pseudopotentials specified with kind, version and appendix for all elements.

    Args: 
        work_status (dict): dictionary descripting the test task

    Return: 
        valid_pseudopotential (dict): see below
        
    Notes:  
    Imagine the pseudpotentials are stored in this way: 
    in [work_folder]/psedupotentials/resources/

    pd_04_spd/
        As.PD04.PBE.UPF
        ...
    sg15_10/
        As_ONCV_PBE-1.0.upf
        ...

    Output is organized as:
    {
        "As": {
            "pd_04_spd": {
                file: "As.PD04.PBE.UPF",
                kind: "pd",
                version: "04",
                appendix: "",
            },
            "sg15_10": {
                file: "As_ONCV_PBE-1.0.upf",
                kind: "sg15",
                version: "10",
                appendix: "",
            }
            ...
        },
    }
    """
    elements = []
    valid_pseudopotentials = {}
    # collect all elements needed for test
    for system in work_status["systems"].keys():
        for element in scan_elements(system):
            if element not in elements:
                elements.append(element)
    # check-in
    os.chdir(work_status["pseudopotentials"]["pseudopot_paths"]["resources"])
    # then for each element, scan available pseudopotentials
    for element in elements:
        for pseudopot_kind in work_status["pseudopotentials"]["kinds"][element]:
            for version in work_status["pseudopotentials"]["versions"][element]:
                for appendix in work_status["pseudopotentials"]["appendices"][element]:
                    pseudopotential = pseudopot_kind
                    if version != "":
                        pseudopotential += "_" + version
                    if appendix != "":
                        pseudopotential += "_" + appendix
                    # if find the folder named as pseudopot_kind + "_" + version + "_" + appendix, then it is valid
                    if os.path.isdir(pseudopotential):
                        # check if the pseudopotential file is in the folder, sometimes the file startswith lowercase
                        files = os.listdir(pseudopotential)
                        pseudopot_valid = False
                        for file in files:
                            if file.startswith(element) or file.startswith(element.lower()):
                                pseudopot_valid = True
                                break
                        if not pseudopot_valid:
                            continue
                        # if really find the pseudopotential, then save it to the list
                        if element not in valid_pseudopotentials.keys():
                            valid_pseudopotentials[element] = {}
                        valid_pseudopotentials[element][pseudopotential] = {
                            "file": file,
                            "kind": pseudopot_kind,
                            "version": version,
                            "appendix": appendix,
                        }
                        
    
    return valid_pseudopotentials

def svp(work_status: dict):

    elements = []
    valid_pseudopotentials = {}
    for system in work_status["systems"]:
        for element in scan_elements(system):
            if element not in elements:
                elements.append(element)
                valid_pseudopotentials[element] = {}

    all_available_pseudopotentials = arch.archive()
    for element in elements:
        for pseudopotential in all_available_pseudopotentials[element]:
            description = arch.description(pseudopotential)
            identifier = "_".join([description[key] for key in description.keys() if description[key] != ""])
            valid_pseudopotentials[element][identifier] = description
            valid_pseudopotentials[element][identifier]["file"] = pseudopotential.split("/")[-1] if pseudopotential.count("/") > 0 else pseudopotential.split("\\")[-1]
    
    return valid_pseudopotentials