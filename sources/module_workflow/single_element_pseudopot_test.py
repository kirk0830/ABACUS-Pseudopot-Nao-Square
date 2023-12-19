"""
Single-element pseudopotential-numerical atomic orbital test workflow

Author: Kirk0830
Date: 2023-11-22

Description:
    This file contains functions that can prepare template input files for single element pseudopotential-numerical atomic orbital test.

Functions:
    prepare_template_input(api_key: str, elements: list, softare = "ABACUS", move = False, **kwargs) -> None
    system_decompose(system: str) -> list
    systems_decompose(systems: list) -> list
    generate_global_section(test_mode: str, analysis_items: list, softare = "ABACUS", save_test_status = False) -> dict
    generate_calculation_section(basis_type: str, functionals: list, ecutwfc: list, stress_deformation_ratios: list = [], **kwargs)
    generate_systems_section(systems: list, **kwargs)
    generate_pseudopotentials_section(work: str, resources: str, elements: list)
    generate_numerical_orbitals_section(work: str, resources: str, elements: list)
    generate_input_json(work_folder = ".", test_mode = "pseudopotential", analysis_items = [], software = "ABACUS", save_test_status = False, basis_type = "pw", functionals = ["PBE"], ecutwfc = [100], systems = [], pseudopot_folder = "./pseudopotentials", pseudopot_resources = "./pseudopotentials/resources", numerical_orbital_folder = "./numerical_orbitals", numerical_orbital_resources = "./numerical_orbitals/resources")
"""
import sources.module_software.qespresso.generation as qe
import sources.module_software.abacus.generation as ab
import sources.module_structure.materials_project as mp
import sources.module_io.batch_set as bs
from main import prepare_test
import os
import json

def prepare_template_input(api_key: str, elements: list, softare = "ABACUS", move = False, **kwargs) -> None:
    """
    Prepare template input files for single element pseudopotential-numerical atomic orbital test.

    Args:
    :param api_key: Materials Project API key
    :param elements: list of elements
    :param softare: software to use, default is ABACUS
    :param move: move template input files to pseudopotential folder, default is False
    :param kwargs: additional arguments for ABACUS

    """
    filenames = []
    for element in elements:
        # download cif file HERE
        cif_filenames = mp.elemental_substances(api_key, element, num_cif = 1)
        if len(cif_filenames) == 0:
            continue
        elif len(cif_filenames) == 1:
            cif_filename = cif_filenames[0]
            if softare == "Quantum ESPRESSO":
                qe.cif_to_quantum_espresso(cif_filename)
                os.rename(cif_filename.replace(".cif", ".in"), f"template_{element}.in")
                filenames.append(f"template_{element}.in")
            elif softare == "ABACUS":
                if "pseudopotentials_in" not in kwargs.keys():
                    raise ValueError("pseudopotentials_in is not provided")
                if "numerical_orbitals_in" not in kwargs.keys():
                    raise ValueError("numerical_orbitals_in is not provided")
                ab.cif_to_STRU(
                    cif_filename, 
                    pseudopotentials_in = kwargs["pseudopotentials_in"], 
                    numerical_orbitals_in = kwargs["numerical_orbitals_in"]
                    )
                filenames.append(f"STRU_{element}.in")
        else:
            index_cif = 0
            for cif_filename in cif_filenames:
                if softare == "Quantum ESPRESSO":
                    qe.cif_to_quantum_espresso(cif_filename)
                    os.rename(cif_filename.replace(".cif", ".in"), f"template_{element}_{index_cif}.in")
                    filenames.append(f"template_{element}_{index_cif}.in")
                elif softare == "ABACUS":
                    if "pseudopotentials_in" not in kwargs.keys():
                        raise ValueError("pseudopotentials_in is not provided")
                    if "numerical_orbitals_in" not in kwargs.keys():
                        raise ValueError("numerical_orbitals_in is not provided")
                    ab.cif_to_STRU(
                        cif_filename, 
                        pseudopotentials_in = kwargs["pseudopotentials_in"], 
                        numerical_orbitals_in = kwargs["numerical_orbitals_in"]
                        )
                    filenames.append(f"STRU_{element}_{index_cif}.in")
                index_cif += 1
    if move:
        for filename in filenames:
            if softare == "Quantum ESPRESSO":
                os.rename(filename, f"./abacus_pseudopotential_square/pseudopotentials/templates/Quantum_ESPRESSO/{filename}")
            elif softare == "ABACUS":
                os.rename(filename, f"./abacus_pseudopotential_square/pseudopotentials/templates/ABACUS/{filename}")

def system_decompose(system: str) -> list:
    """
    Decompose system into elements.

    Args:
    :param system: system name

    Returns:
    :return: elements: list of elements
    """
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

def systems_decompose(systems: list) -> list:
    """
    Decompose systems into elements.

    Args:
    :param systems: list of system names

    Returns:
    :return: elements: list of elements
    """
    elements = []
    for system in systems:
        _elements = system_decompose(system)
        elements += _elements
    elements = list(set(elements))
    return elements

def generate_global_section(
        test_mode: str,
        analysis_items: list,
        softare = "ABACUS",
        save_test_status = False
        ) -> dict:
    """
    Generate global section of input json.

    Args:
    :param test_mode: test mode, default is pseudopotential
    :param analysis_items: list of analysis items
    :param softare: software to use, default is ABACUS
    :param save_test_status: save test status, default is False

    Returns:
    :return: global_section: global section of input json
    """
    return {
        "test_mode": test_mode,
        "analysis_items": analysis_items,
        "software": softare,
        "save_test_status": save_test_status
    }

def generate_calculation_section(
        basis_type: str,
        functionals: list,
        ecutwfc: list,
        stress_deformation_ratios: list = [],
        **kwargs
        ):
    """
    Generate calculation section of input json.

    Args:
    :param basis_type: basis type, default is pw
    :param functionals: list of functionals
    :param ecutwfc: list of ecutwfc
    :param stress_deformation_ratios: list of stress deformation ratios
    :param kwargs: additional arguments

    Returns:
    :return: calculation_section: calculation section of input json
    """
    return_dict = {
        "basis_type": basis_type,
        "functionals": functionals,
        "ecutwfc": ecutwfc,
        "stress_deformation_ratios": stress_deformation_ratios
    }
    for key, value in kwargs.items():
        return_dict[key] = value
    return return_dict

def generate_systems_section(
        systems: list,
        **kwargs
        ):
    """
    Generate systems section of input json.

    Args:
    :param systems: list of systems
    :param kwargs: additional arguments

    Returns:
    :return: systems_section: systems section of input json
    """
    return_dict = {}
    for system in systems:
        return_dict[system] = {
            "experimental_values": {}
        }
        for key, value in kwargs.items():
            if key == system:
                return_dict[system] = value
    return return_dict

def generate_pseudopotentials_section(
        work: str,
        resources: str,
        elements: list
):
    """
    Generate pseudopotentials section of input json.

    Args:
    :param work: pseudopotential folder
    :param resources: pseudopotential resources folder
    :param elements: list of elements

    Returns:
    :return: pseudopotentials_section: pseudopotentials section of input json
    """
    return_dict = {
        "pseudopot_paths": {
            "work": work,
            "resources": resources
        }
    }
    return_dict["kinds"] = {}
    return_dict["versions"] = {}
    return_dict["appendices"] = {}
    for element in elements:
        return_dict["kinds"][element] = []
        return_dict["versions"][element] = []
        return_dict["appendices"][element] = []
    return return_dict

def generate_numerical_orbitals_section(
        work: str,
        resources: str,
        elements: list
):
    """
    Generate numerical orbitals section of input json.

    Args:
    :param work: numerical orbital folder
    :param resources: numerical orbital resources folder
    :param elements: list of elements

    Returns:
    :return: numerical_orbitals_section: numerical orbitals section of input json
    """
    return_dict = {
        "nao_paths": {
            "work": work,
            "resources": resources
        }
    }
    return_dict["types"] = {}
    return_dict["rcuts"] = {}
    return_dict["appendices"] = {}
    for element in elements:
        return_dict["types"][element] = []
        return_dict["rcuts"][element] = []
        return_dict["appendices"][element] = []
    return return_dict

def generate_input_json(
        work_folder = ".",
        test_mode = "pseudopotential",
        analysis_items = [],
        software = "ABACUS",
        save_test_status = False,
        basis_type = "pw",
        functionals = ["PBE"],
        ecutwfc = [100],
        systems = [],
        pseudopot_folder = "./pseudopotentials",
        pseudopot_resources = "./pseudopotentials/resources",
        numerical_orbital_folder = "./numerical_orbitals",
        numerical_orbital_resources = "./numerical_orbitals/resources"
):
    """
    Generate input json.

    Args:
    :param work_folder: work folder, default is "."
    :param test_mode: test mode, default is pseudopotential
    :param analysis_items: list of analysis items
    :param software: software to use, default is ABACUS
    :param save_test_status: save test status, default is False
    :param basis_type: basis type, default is pw
    :param functionals: list of functionals
    :param ecutwfc: list of ecutwfc
    :param systems: list of systems
    :param pseudopot_folder: pseudopotential folder
    :param pseudopot_resources: pseudopotential resources folder
    :param numerical_orbital_folder: numerical orbital folder
    :param numerical_orbital_resources: numerical orbital resources folder

    Returns:
    :return: input_dict: input json
    """
    elements = []
    # decompse system into elements
    for index, system in enumerate(systems):
        _elements = system_decompose(system)
        # add elements to elements list
        elements += _elements
    # remove duplicates
    elements = list(set(elements))
    
    # generate input json
    input_dict = {"work_folder": work_folder}
    input_dict["global"] = generate_global_section(
        test_mode = test_mode,
        analysis_items = analysis_items,
        softare = software,
        save_test_status = save_test_status
    )
    input_dict["calculation"] = generate_calculation_section(
        basis_type = basis_type,
        functionals = functionals,
        ecutwfc = ecutwfc
    )
    input_dict["systems"] = generate_systems_section(
        systems = systems,
    )
    input_dict["pseudopotentials"] = generate_pseudopotentials_section(
        work = pseudopot_folder,
        resources = pseudopot_resources,
        elements = elements
    )
    input_dict["numerical_orbitals"] = generate_numerical_orbitals_section(
        work = numerical_orbital_folder,
        resources = numerical_orbital_resources,
        elements = elements
    )
    return input_dict

if __name__ == "__main__":

    # 0. matrials project api key
    api_key = "wV1HUdmgESPVgSmQj5cc8WvttCO8NTHp"
    # 1. define systems
    #systems = ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']
    systems = ['La'] # does not have pd_04_core, pd_04_icmod
    # 2. prepare template input
    prepare_template_input(api_key, systems, softare = "Quantum ESPRESSO", move = True)
    # 3. generate input json skeleton
    input_dict = generate_input_json(software="Quantum ESPRESSO",systems=systems)
    with open("input.json", "w") as f:
        json.dump(input_dict, f, indent=4)
    # 4. adjust template input according to input json
    pseudopotentials = [
        "pd_04",
        "paw_10",
        "gthpp1",
        "atompaw.wentzcovitch"
        ]
    for pseudopotential in pseudopotentials:
        bs.set_pseudopotential(
            kind = pseudopotential,
            version = "",
            appendix = ""
        )
        prepare_test("input.json")
