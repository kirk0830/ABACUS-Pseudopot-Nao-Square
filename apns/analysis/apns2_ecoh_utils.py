####
# utilities for the calculation of cohesive energy
####

def desc_equal_bulk_vs_atom(desc1: dict, desc2: dict):
    """for cohesive energy calculation task, only admit the difference
    in CellGenerator["identifier"] and CellGenerator["config"]"""
    from apns.analysis.apns2_utils import cal_desc_diff
    diff = cal_desc_diff(desc1, desc2)
    if set(diff.keys()) != {"CellGenerator"}:
        return False
    if set(diff["CellGenerator"].keys()) != {"identifier", "config"}:
        return False
    if set(diff["CellGenerator"]["identifier"]) != {"cif", "molecule"}:
        return False
    return True

def pair(data: list):
    """pair the bulk and isolated atom data from
    collect_jobs function returned value"""
    
    paired = []
    # all (bulk, atom) pair should have identical sections
    # "ParamSet" and "AtomSpecies"
    # the former ensures the same calculation settings
    # the latter ensures the identical pseudopotential and
    # numerical atomic orbitals (if any)
    bulk = [(desc, eks) for desc, eks in data if desc["CellGenerator"]["identifier"] == "cif"]
    atom = [(desc, eks) for desc, eks in data if desc["CellGenerator"]["identifier"] == "molecule"]
    # then pair the bulk and atom data
    for desc_b, e_b in bulk:
        for desc_a, e_a in atom:
            if desc_equal_bulk_vs_atom(desc_b, desc_a):
                paired.append((desc_b, e_b, e_a))
                break
    return paired

def cal_e_cohesive(e_b: float, e_a: float, natom: int):
    """calculate the cohesive energy"""
    return (e_b - e_a * natom) / natom

def plot():
    pass