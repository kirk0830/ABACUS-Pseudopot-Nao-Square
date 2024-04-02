import json
def read(finp: str):
    """read input file"""
    with open(finp, "r") as f:
        return check(json.load(f))

import os
def check(inp: dict):
    """check correctness and completeness of input file"""
    print("Checking input file...")
    # global section
    if "global" not in inp:
        raise ValueError("global section not found in input file")
    if "test_mode" not in inp["global"]:
        raise ValueError("test_mode not found in global section")
    if "pseudo_dir" not in inp["global"]:
        raise ValueError("pseudo_dir not found in global section")
    if not os.path.exists(inp["global"]["pseudo_dir"]):
        raise ValueError("pseudo_dir not found")
    if "orbital_dir" not in inp["global"]:
        raise ValueError("orbital_dir not found in global section")
    if not os.path.exists(inp["global"]["orbital_dir"]):
        print("Warning: orbital_dir not found, creating one. This folder is the target for orbital generation.")
        os.makedirs(inp["global"]["orbital_dir"], exist_ok=True)
    # orbgen section
    if "orbgen" not in inp:
        raise ValueError("orbgen section not found in input file")
    if "generate_mode" not in inp["orbgen"]:
        raise ValueError("generate_mode not found in orbgen section")
    if inp["orbgen"]["generate_mode"] not in ["in-situ", "ex-situ"]:
        raise ValueError("generate_mode should be in-situ or ex-situ")
    if inp["orbgen"]["generate_mode"] == "in-situ":
        print("To run orbgen workflow in-situ, additional check proceed...")
        if "generator" not in inp["orbgen"]:
            raise ValueError("generator not found in orbgen section")
        if not os.path.exists(inp["orbgen"]["generator"]):
            raise ValueError(f"orbgen generator not found: {inp['orbgen']['generator']}")
        print("orbgen generator found.")
    else:
        if "generator" not in inp["orbgen"]:
            print("Warning: generator not found in orbgen section")
        elif not os.path.exists(inp["orbgen"]["generator"]):
            print(f"Warning: orbgen generator not found: {inp['orbgen']['generator']}")
        else:
            print("orbgen generator found.")
    if "environment" not in inp["orbgen"]:
        raise ValueError("environment not found in orbgen section")
    if "mpi_command" not in inp["orbgen"]:
        raise ValueError("mpi_command not found in orbgen environment")
    if "abacus_command" not in inp["orbgen"]:
        raise ValueError("abacus_command not found in orbgen environment")
    if "rcuts" not in inp["orbgen"]:
        raise ValueError("rcuts not found in orbgen section")
    if len(inp["orbgen"]["rcuts"]) == 0:
        raise ValueError("rcuts section is empty")
    if "types" not in inp["orbgen"]:
        raise ValueError("types not found in orbgen section")
    if len(inp["orbgen"]["types"]) == 0:
        raise ValueError("types section is empty")
    # abacus section
    if "abacus" not in inp:
        raise ValueError("abacus section not found in input file")
    # systems section
    if "systems" not in inp:
        raise ValueError("systems section not found in input file")
    if len(inp["systems"]) == 0:
        raise ValueError("systems section is empty")
    # pseudopotentials section
    if "pseudopotentials" not in inp:
        raise ValueError("pseudopotentials section not found in input file")
    if "kinds" not in inp["pseudopotentials"]:
        raise ValueError("kind not found in pseudopotentials section")
    if len(inp["pseudopotentials"]["kinds"]) == 0:
        raise ValueError("kinds section is empty")
    if "versions" not in inp["pseudopotentials"]:
        raise ValueError("versions not found in pseudopotentials section")
    if len(inp["pseudopotentials"]["versions"]) == 0:
        raise ValueError("versions section is empty")
    if "appendices" not in inp["pseudopotentials"]:
        raise ValueError("appendices not found in pseudopotentials section")
    if len(inp["pseudopotentials"]["appendices"]) == 0:
        raise ValueError("appendices section is empty")

    return inp

import apns.module_pseudo.manage as ampm
import apns.module_io.input_translate as amiit
import apns.module_structure.basic as amsb
def initialize(finp: str):
    """setup the orbgen workflow, return the checked inp file contents 
    along with available pspot_ids selected by the input file
    
    Args:
    finp (str): fname of e.g., input_orbgen.json
    
    Returns:
    dict: setting in "orbgen" section
    list: elements
    dict: valid pseudopotentials whose keys are element and values are key-value pair of pspotid
    and fname"""
    inp = read(finp)
    elements = amsb.scan_elements(inp["systems"])
    inp = amiit.expand(inp, elements)
    vupfs = ampm.valid_pseudo(inp["global"]["pseudo_dir"], elements, inp["pseudopotentials"])
    # vupfs will have keys as elements 
    # and values the key-value pairs of pspot_ids and file name with path

    # design note: 
    # the first contains information about: "what is needed for generating orbitals?"
    # the second contains: "what kind of DFT calculation will perform?"
    # the third contains: "for what elements? (although the original "system" key can specify
    # composites, but eventually it would be the question of single element)"
    # the fourth contains: what pseudopotentials
    return inp["orbgen"], inp["abacus"], elements, vupfs

import unittest
class TestInitialize(unittest.TestCase):

    def test_check(self):
        """written by Github.copilot"""
        with self.assertRaises(ValueError):
            check({"global": {}})
        with self.assertRaises(ValueError):
            check({"global": {"test_mode": ""}})
        with self.assertRaises(ValueError):
            check({"global": {"pseudo_dir": ""}})
        with self.assertRaises(ValueError):
            check({"global": {"orbital_dir": ""}})
        with self.assertRaises(ValueError):
            check({"orbgen": {}})
        with self.assertRaises(ValueError):
            check({"orbgen": {"generator": ""}})
        with self.assertRaises(ValueError):
            check({"orbgen": {"environment": {}}})
        with self.assertRaises(ValueError):
            check({"orbgen": {"mpi_command": ""}})
        with self.assertRaises(ValueError):
            check({"orbgen": {"abacus_command": ""}})
        with self.assertRaises(ValueError):
            check({"abacus": {}})
        with self.assertRaises(ValueError):
            check({"systems": []})
        with self.assertRaises(ValueError):
            check({"pseudopotentials": {}})
        with self.assertRaises(ValueError):
            check({"pseudopotentials": {"kinds": []}})
        with self.assertRaises(ValueError):
            check({"pseudopotentials": {"versions": []}})
        with self.assertRaises(ValueError):
            check({"pseudopotentials": {"appendices": []}})
        with self.assertRaises(ValueError):
            check({"numerical_orbitals": {}})
        with self.assertRaises(ValueError):
            check({"numerical_orbitals": {"rcuts": []}})
        with self.assertRaises(ValueError):
            check({"numerical_orbitals": {"types": []}})
        with self.assertRaises(ValueError):
            check({"global": {"pseudo_dir": "not_exist"}})
        with self.assertRaises(ValueError):
            check({"orbgen": {"generator": "not_exist"}})

if __name__ == '__main__':
    unittest.main()