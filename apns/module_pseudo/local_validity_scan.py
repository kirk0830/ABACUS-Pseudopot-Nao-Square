import os
import apns.module_database.database as ai
import apns.module_pseudo.archive as arch
import apns.module_workflow.identifier as id
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

def scan_orbs(elements: list, pseudopotentials: dict) -> dict:
    """new version of scan valid pseudopotential"""
    valid_pseudopotentials = { element: {} for element in elements }
    all_available_pseudopotentials = arch.archive()

    for element in elements:
        _kinds = pseudopotentials["kinds"][element] if len(pseudopotentials["kinds"][element]) > 0 else ["all"]
        _versions = pseudopotentials["versions"][element] if len(pseudopotentials["versions"][element]) > 0 else ["all"]
        _appendices = pseudopotentials["appendices"] [element] if len(pseudopotentials["appendices"] [element]) > 0 else ["all"]
        for pseudopotential in all_available_pseudopotentials[element]:
            description = arch.description(pseudopotential)

            _b1 = "all" in _kinds
            _b2 = description["kind"] in _kinds
            _b3 = "all" in _versions
            _b4 = description["version"] in _versions
            _b5 = "all" in _appendices
            _b6 = description["appendix"] in _appendices

            if _b1 or (_b2 and (_b3 or (_b4 and (_b5 or _b6)))):
                identifier = id.pseudopotential(kind=description["kind"], 
                                                version=description["version"], 
                                                appendix=description["appendix"])
                valid_pseudopotentials[element][identifier] = description
                valid_pseudopotentials[element][identifier]["file"] = pseudopotential.split("/")[-1] if pseudopotential.count("/") > 0 else pseudopotential.split("\\")[-1]
        if len(valid_pseudopotentials[element]) == 0:
            raise ValueError("No valid pseudopotential for element %s.\nkinds: %s\nversions: %s\nappendices: %s" % (element, _kinds, _versions, _appendices))
        
    return valid_pseudopotentials

def svp(work_status: dict):
    """scan avaiable pseudopotentials
    if specified as "all" in "kind", then all valid will be returned, otherwise will filter.

    This function will return a dictionary like:
    {
        "Er": {
            "identifier1": {
                "kind": "nc",
                "version": "1.0",
                "appendix": "1.0",
                "file": "Er.nc.UPF"
            },
            ...
        }
    }
    """
    elements = []
    valid_pseudopotentials = {}
    for system in work_status["systems"]:
        for element in scan_elements(system):
            if element not in elements:
                elements.append(element)
                valid_pseudopotentials[element] = {}

    all_available_pseudopotentials = arch.archive()

    for element in elements:
        _kinds = work_status["pseudopotentials"]["kinds"][element]
        _versions = work_status["pseudopotentials"]["versions"][element]
        _appendices = work_status["pseudopotentials"]["appendices"] [element] 
        for pseudopotential in all_available_pseudopotentials[element]:
            description = arch.description(pseudopotential)

            _b1 = "all" in _kinds
            _b2 = description["kind"] in _kinds
            _b3 = "all" in _versions
            _b4 = description["version"] in _versions
            _b5 = "all" in _appendices
            _b6 = description["appendix"] in _appendices

            if _b1 or (_b2 and (_b3 or (_b4 and (_b5 or _b6)))):
                identifier = id.pseudopotential(kind=description["kind"], 
                                                version=description["version"], 
                                                appendix=description["appendix"])
                valid_pseudopotentials[element][identifier] = description
                valid_pseudopotentials[element][identifier]["file"] = pseudopotential.split("/")[-1] if pseudopotential.count("/") > 0 else pseudopotential.split("\\")[-1]
        if len(valid_pseudopotentials[element]) == 0:
            raise ValueError("No valid pseudopotential for element {}.".format(element))
    return valid_pseudopotentials