"""
Multi-dimensional to one-dimensional batch set functions

Author: Kirk0830
Date: 2023-11-22

Description:
    This file contains functions that can set pseudopotential and numerical orbital information in input.json for all elements uniformly.

Useful cases:
    1. test pseudopotential over systems, rather than elements in each system.
"""
import json

def set_pseudopotential(kind: str, version: str, appendix: str):
    """
    set pseudopotential information in input.json for all elements uniformly

    Args:

    :param kind: pseudopotential kind
    :param version: pseudopotential version
    :param appendix: pseudopotential appendix

    Returns:

    :return: None
    """
    with open("input.json", "r") as f:
        input_json = json.load(f)
    for element in input_json["pseudopotentials"]["kinds"].keys():
        input_json["pseudopotentials"]["kinds"][element] = [kind]
        input_json["pseudopotentials"]["versions"][element] = [version]
        input_json["pseudopotentials"]["appendices"][element] = [appendix]

    with open("input.json", "w") as f:
        json.dump(input_json, f, indent=4)

def set_numerical_orbital(type: str, rcut: int, appendix: str):
    """
    set numerical orbital information in input.json for all elements uniformly

    Args:

    :param type: numerical orbital type
    :param rcut: numerical orbital rcut
    :param appendix: numerical orbital appendix

    Returns:

    :return: None
    """
    with open("input.json", "r") as f:
        input_json = json.load(f)
    for element in input_json["numerical_orbitals"]["types"].keys():
        input_json["numerical_orbitals"]["types"][element] = [type]
        input_json["numerical_orbitals"]["rcuts"][element] = [rcut]
        input_json["numerical_orbitals"]["appendices"][element] = [appendix]

    with open("input.json", "w") as f:
        json.dump(input_json, f, indent=4)
