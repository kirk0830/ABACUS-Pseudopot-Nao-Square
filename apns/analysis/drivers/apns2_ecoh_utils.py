####
# utilities for the calculation of cohesive energy
####

def paramset_equal(paramset1: dict, paramset2: dict, excluded: list) -> bool:
    """check if two paramsets are equal. Due to it is unnecessary to keep the 
    same value for nspin between bulk and atom calculation (in principle, the
    calculation of atom should/must be spin-polarized). Therefore, nspin should
    always be excluded from the comparison."""
    return {k: v for k, v in paramset1.items() if k not in excluded} == \
        {k: v for k, v in paramset2.items() if k not in excluded}

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
            if paramset_equal(desc_b["ParamSet"], desc_a["ParamSet"], ["nspin"]):
                if desc_b["AtomSpecies"] == desc_a["AtomSpecies"]:
                    paired.append((desc_b, e_b, e_a))
    return paired

def cal_e_cohesive(e_b: float, e_a: float, natom: int):
    """calculate the cohesive energy"""
    return (e_b - e_a * natom) / natom

def plot():
    pass