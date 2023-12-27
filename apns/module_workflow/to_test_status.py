"""Generate a test status file to define the work of this package.
A test is defined as:
[functional]_[pseudo_sets](_[nao_sets])_[system]

test_status
"""
from apns.module_pseudo.local_validity_scan import svp as svp
from apns.module_nao.local_validity_scan import scan_valid_numerical_orbitals as svno
from apns.module_structure.basic import scan_elements as se
import itertools as it
import json
import os
import time

def _batch_info_import_(
        test_status: dict, 
        valid_data: dict, 
        system: str, 
        test_name: str, 
        element: str, 
        pseudo_element: str, 
        nao_type_element = "") -> dict:
    """import information to test_status in batch

    Args:
        test_status (dict): dictionary to import information
        valid_data (dict): can be valid_pseudopotentials or valid_numerical_orbitals
        system (str): system to write information
        test_name (str): test identifier, see module_workflow/README.md for details
        element (str): as it is
        pseudo_element (str): _description_
        nao_type_element (str, optional): _description_. Defaults to "".

    Returns:
        dict: modified test_status
    """

    if nao_type_element == "":
        test_status["tests"][system][test_name]["pseudopotentials"]["files"][element] = valid_data[element][pseudo_element]["file"]
        test_status["tests"][system][test_name]["pseudopotentials"]["info"][element] = {
            "kind": valid_data[element][pseudo_element]["kind"],
            "version": valid_data[element][pseudo_element]["version"],
            "appendix": valid_data[element][pseudo_element]["appendix"]
        }
    else:
        test_status["tests"][system][test_name]["numerical_orbitals"]["files"][element] = valid_data[element][pseudo_element][nao_type_element]["file"]
        test_status["tests"][system][test_name]["numerical_orbitals"]["info"][element] = {
            "type": valid_data[element][pseudo_element][nao_type_element]["type"],
            "rcut": valid_data[element][pseudo_element][nao_type_element]["rcut"],
            "appendix": valid_data[element][pseudo_element][nao_type_element]["appendix"]
        }
    return test_status

def _test_template_(elements, ecutwfc: float = -1) -> dict:
    """create one empty test in test_status

    Args:
        elements (list): collection of elements in this test
        ecutwfc (int, optional): ecutwfc. Defaults to -1.

    Returns:
        dict: initialized test
    """
    result = {
        "elements": elements,
        "pseudopotentials":   {"files": {}, "info": {}},
        "numerical_orbitals": {"files": {}, "info": {}},
        }
    if ecutwfc > 0:
        result["ecutwfc"] = ecutwfc
    return result

def _valid_files_scan_(work_status: dict):
    """scan valid files, (optionally) check correspondences

    Args:
        work_status (dict): _description_

    Returns:
        tuple: (valid_pseudopotentials, valid_numerical_orbitals)
    """
    valid_pseudopotentials = svp(work_status)
    valid_numerical_orbitals = {}

    if work_status["calculation"]["basis_type"] == "lcao":
        valid_numerical_orbitals = svno(work_status, valid_pseudopotentials)
        # available pseudopotentials and numerical orbitals should correspond
        if len(list(valid_numerical_orbitals.keys())) != len(list(valid_pseudopotentials.keys())):
            print("Error: number of elements in valid_numerical_orbitals is not equal to number of elements in valid_pseudopotentials.")
            exit(1)

    return valid_pseudopotentials, valid_numerical_orbitals

def _test_initialization_(work_status: dict, system: str, valid_pseudopotentials: dict, valid_numerical_orbitals: dict = {}):

    test = {}
    # convert test name from tuple (["pd_04_spd", "D6"], ["sg15_10", "D6"], ...) to pd04spdsg1510..._D6D6...., each test is equivalent with one tuple
    test_name = ""
    for item in range(len(test[0])): # item can only be 1 or 2, 1 is for pseudopotential without ecutwfc tests, 
                                        #                          2 is for numerical orbitals or pseudopotential with ecutwfc tests
        for axis_element in test: # axis can be many, for instance InP, it is 2, BaTiO3, it is 3
            test_name += axis_element[item].replace("_", "")
        test_name += "_" # add _ between pseudopotential and numerical orbitals

    elements = se(system)
    if work_status["calculation"]["basis_type"] == "pw":
        # add ecutwfc to test name if there are more than one ecutwfc
        for ecutwfc in work_status["calculation"]["ecutwfc"]:
            if len(work_status["calculation"]["ecutwfc"]) > 1:
                test_name += "ecut" + str(ecutwfc)
            else:
                test_name = test_name[:-1] # remove the last _
            print("prepare test: " + test_name + " for system: " + system)
            # initialize test_status for this test
            test = _test_template_(elements, ecutwfc)
            for ie, element in enumerate(elements): # so we can get the index of element
                pseudo_element = test[ie][0]
                # get pseudopotential information and copy to test_status
                test_status = _batch_info_import_(test_status, 
                                                    valid_pseudopotentials, 
                                                    system, 
                                                    test_name, 
                                                    element, 
                                                    pseudo_element)
            test_name = test_name.split("ecut")[0] # remove ecutwfc from test name

    elif work_status["calculation"]["basis_type"] == "lcao":
        test_name = test_name[:-1] # remove the last _
        print("prepare test: " + test_name + " for system: " + system)
        # initialize test_status for this test
        test = _test_template_(elements)
        # loop over elements
        for ie, element in enumerate(elements): # so we can get the index of element
            pseudo_element = test[ie][0]
            nao_type_element = test[ie][1]
            # get pseudopotential information and copy to test_status
            test_status = _batch_info_import_(test_status, 
                                              valid_pseudopotentials, 
                                              system, 
                                              test_name, 
                                              element, 
                                              pseudo_element)
            # get numerical orbital information and copy to test_status
            test_status = _batch_info_import_(test_status, 
                                              valid_numerical_orbitals, 
                                              system, 
                                              test_name, 
                                              element, 
                                              pseudo_element, 
                                              nao_type_element)

def _test_status_backup_(test_status: dict):
        
    os.chdir(test_status["paths"]["work_dir"])
    test_status_json = "test_status_" + time.strftime("%Y%m%d_%H%M%S", time.localtime()) + ".json"
    with open(test_status_json, "w") as f:
        json.dump(test_status, f, indent=4)
    print("test_status saved to " + test_status_json)

def to(work_status: dict) -> dict:

    test_status = {
        "paths": {
            "work_dir": work_status["global"]["work_dir"],
            "pseudo_dir": work_status["global"]["pseudo_dir"],
            "orbital_dir": work_status["global"]["orbital_dir"]
        },
        "tests": {}
    } # test_status["tests"][system]: {"test1": {}, ...}
    valid_pseudopotentials, valid_numerical_orbitals = _valid_files_scan_(work_status)
    for system in work_status["systems"]:
        # goal is to store all possible tests to test_status system by system
        test_status["tests"][system] = {} # create to store all tests for this system
        # allocate memory for storing valid pseudopotentials list element by element
        valid_pseudopots_to_combine = []
        # get elements in present system
        elements = se(system)

        for element in elements:
            if work_status["calculation"]["basis_type"] == "pw":
                valid_pseudopots_to_combine.append([[valid_pseudopot] for valid_pseudopot in list(valid_pseudopotentials[element].keys())])
                # for every element, valid_pseudopots_to_combine will be like:
                # [
                #     [["pd_04_spd"], ["sg15_10"], ...],
                #     [["pd_04_spd"], ["sg15_10"], ...],
                #     ...
                # ]
            elif work_status["calculation"]["basis_type"] == "lcao":
                valid_pseudopot_nao_type_element = []
                for pseudopotential in valid_numerical_orbitals[element].keys():
                    for nao_type in valid_numerical_orbitals[element][pseudopotential].keys():
                        valid_pseudopot_nao_type_element.append([pseudopotential, nao_type])
                valid_pseudopots_to_combine.append(valid_pseudopot_nao_type_element)
                # thus for every element, valid_pseudopots_to_combine will be like:
                # [
                #     pseudopot and corresponding numerical orbital in each []
                #     [["pd_04_spd", "D6"], ["sg15_10", "D6"], ...], # element 1
                #     [["pd_04_spd", "D6"], ["sg15_10", "D6"], ...], # element 2
                #     ...
                # ]
        # make Cartesian product of all valid pseudopotentials
        # test combinations (per) system
        test_combinations_system = list(it.product(*valid_pseudopots_to_combine))
        # such that each combination will be a list of tuples of list:
        # for lcao, it will be like:
        # [
        #             element1             element2           element3
        #     test1: (["pd_04_spd", "D6"], ["sg15_10", "D6"], ...),
        #     test2: (["pd_04_spd", "D6"], ["sg15_10", "D6_osci_rm"]),
        #     ...
        # ]
        # for pw, it will be like:
        # [
        #             element1       element2     element3
        #     test1: (["pd_04_spd"], ["sg15_10"], ...),
        #     ...
        # ]
        # loop over all combinations
        for test in test_combinations_system:
            # convert test name from tuple (["pd_04_spd", "D6"], ["sg15_10", "D6"], ...) to pd04spdsg1510..._D6D6...., each test is equivalent with one tuple
            test_name = ""
            for item in range(len(test[0])): # item can only be 1 or 2, 1 is for pseudopotential without ecutwfc tests, 
                                             #                          2 is for numerical orbitals or pseudopotential with ecutwfc tests
                for axis_element in test: # axis can be many, for instance InP, it is 2, BaTiO3, it is 3
                    test_name += axis_element[item].replace("_", "")
                test_name += "_" # add _ between pseudopotential and numerical orbitals
            if work_status["calculation"]["basis_type"] == "pw":
                # add ecutwfc to test name if there are more than one ecutwfc
                for ecutwfc in work_status["calculation"]["ecutwfc"]:
                    if len(work_status["calculation"]["ecutwfc"]) > 1:
                        test_name += "ecut" + str(ecutwfc)
                    else:
                        test_name = test_name[:-1] # remove the last _
                    print("prepare test: " + test_name + " for system: " + system)
                    # initialize test_status for this test
                    test_status["tests"][system][test_name] = _test_template_(elements, ecutwfc)
                    for ie, element in enumerate(elements): # so we can get the index of element
                        pseudo_element = test[ie][0]
                        # get pseudopotential information and copy to test_status
                        test_status = _batch_info_import_(test_status, 
                                                          valid_pseudopotentials, 
                                                          system, 
                                                          test_name, 
                                                          element, 
                                                          pseudo_element)
                    test_name = test_name.split("ecut")[0] # remove ecutwfc from test name

            elif work_status["calculation"]["basis_type"] == "lcao":
                test_name = test_name[:-1] # remove the last _
                print("prepare test: " + test_name + " for system: " + system)
                # initialize test_status for this test
                test_status["tests"][system][test_name] = _test_template_(elements)
                # loop over elements
                for ie, element in enumerate(elements): # so we can get the index of element
                    pseudo_element = test[ie][0]
                    nao_type_element = test[ie][1]
                    # get pseudopotential information and copy to test_status
                    test_status = _batch_info_import_(test_status, 
                                                      valid_pseudopotentials, 
                                                      system, 
                                                      test_name, 
                                                      element, 
                                                      pseudo_element)
                    # get numerical orbital information and copy to test_status
                    test_status = _batch_info_import_(test_status, 
                                                      valid_numerical_orbitals, 
                                                      system, 
                                                      test_name, 
                                                      element, 
                                                      pseudo_element, 
                                                      nao_type_element)
    
    if work_status["global"]["save_log"] or "save_log" not in work_status["global"].keys():
        # save test_status to json file, name marked with time stamp
        # note that the test_status is saved in the work folder
        print("save test_status to json file")
        _test_status_backup_(test_status)
    else:
        print("warning: test_status is necessary for analyzing test results, but it is not saved according to input.")
    # change back to root folder
    os.chdir(test_status["paths"]["work_dir"])
    return test_status