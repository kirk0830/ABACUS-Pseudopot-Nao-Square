import json
def read(finp: str):
    """read input file"""
    with open(finp, "r") as f:
        return json.load(f)

import os
def check(inp: dict):
    """check correctness and completeness of input file"""
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
    if "generator" not in inp["orbgen"]:
        raise ValueError("generator not found in orbgen section")
    if not os.path.exists(inp["orbgen"]["generator"]):
        raise ValueError("orbgen generator not found")
    if "environment" not in inp["orbgen"]:
        raise ValueError("environment not found in orbgen section")
    if "mpi_command" not in inp["orbgen"]["environment"]:
        raise ValueError("mpi_command not found in orbgen environment")
    if "abacus_command" not in inp["orbgen"]["environment"]:
        raise ValueError("abacus_command not found in orbgen environment")
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
    # numerical_orbitals section
    if "numerical_orbitals" not in inp:
        raise ValueError("numerical_orbitals section not found in input file")
    if "rcuts" not in inp["numerical_orbitals"]:
        raise ValueError("rcuts not found in numerical_orbitals section")
    if len(inp["numerical_orbitals"]["rcuts"]) == 0:
        raise ValueError("rcuts section is empty")
    if "types" not in inp["numerical_orbitals"]:
        raise ValueError("types not found in numerical_orbitals section")
    if len(inp["numerical_orbitals"]["types"]) == 0:
        raise ValueError("types section is empty")
    return inp

def link_pspotlib(inp: dict):
    """according to locally available pseudopotential library, link the pspot_ids to the input file
    , and return the available pspot_ids selected by the input file"""
    pspot_ids = []
    with open(inp["global"]["pseudo_dir"] + "/description.json", "r") as f:
        pspotlib = json.load(f)
    
    if inp["pseudopotentials"]["kinds"] == ["all"]:
        pspot_ids = list(pspotlib.keys())
    else:
        for kind in inp["pseudopotentials"]["kinds"]:
            for pspot_id in pspotlib.keys():
                if pspot_id.startswith(kind):
                    if inp["pseudopotentials"]["versions"] == ["all"]:
                        pspot_ids.append(pspot_id)
                    else:
                        for version in inp["pseudopotentials"]["versions"]:
                            if pspot_id.startswith(kind + "_" + version):
                                if inp["pseudopotentials"]["appendices"] == ["all"]:
                                    pspot_ids.append(pspot_id)
                                else:
                                    for appendix in inp["pseudopotentials"]["appendices"]:
                                        if pspot_id.startswith(kind + "_" + version + "_" + appendix):
                                            pspot_ids.append(pspot_id)
    if len(pspot_ids) == 0:
        raise ValueError("no pseudopotential found in range selected by input file. Please check the input file.")
    """NOTE: we did not check the existence of the pseudopotential files here,
    because we plan to support not only one element generating numerical orbitals in one shot,
    therefore the detailed element-specific pseudopotential file existence check will be done 
    in function apns/module_nao/orbgen/siab_generator.py:siab_generator() instead."""
    return pspot_ids

def initialize(finp: str):
    """setup the orbgen workflow, return the inp file contents along with available pspot_ids selected by the input file"""
    inp = read(finp)
    check(inp)
    pspot_ids = link_pspotlib(inp)
    return inp, pspot_ids