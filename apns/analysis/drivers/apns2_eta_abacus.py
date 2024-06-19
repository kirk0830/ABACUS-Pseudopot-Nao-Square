""""""

def collect_jobs(folder: str):
    """Collect APNS jobs for band structure similarity calculation"""
    import os, json, re
    from apns.analysis.postprocess.read_abacus_out import read_istate, read_kpoints
    print("* * * Collect APNS jobs * * *".center(100))
    result = []
    for root, _, files in os.walk(folder):
        for file in files:
            if re.match(r"(running_)(\w+)(\.log)", file):
                parent = os.path.dirname(root)
                with open(os.path.join(parent, "description.json"), "r") as f:
                    desc = json.load(f)
                # we will have istate.info in this folder
                fistate = os.path.join(root, "istate.info")
                if not fistate:
                    print(f"Warning: {parent} does not have istate.info")
                istate, kpoints = read_istate(fistate)
                fkpoints = os.path.join(root, "kpoints")
                kref = read_kpoints(fkpoints)
                # because the occ information in istate is multiplied by kpoints weight,
                # need divide it.
                kwt = lookup_kwt(kpoints, kref)
                istate = clean_kwt_from_istate(istate, kwt)
                ekb, occ = decompose_istate(istate)
                result.append((desc, ekb, occ, kwt))
                # note! desc is a dict, ekb and occ are nested like [ispin][ik][iband]
                # kwt is a list of kpoints weight
    return result

def lookup_kwt(kpoints, kref):
    """lookup the kpoints weight by comparing between the coordinate extracted
    from istate.info with the kpoints file extracted kref.
    kpoints has the simple format like: 
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.1], ...]
    kref is a little bit nested like:
    [[1, 0.0, 0.0, 0.0, 0.1], [2, 0.0, 0.0, 0.1, 0.1], ...]
    in which the first index is just the index of kpoints, the last is the weight.
    
    will return a list of kpoints weight, in the same order as kpoints read in 
    the istate.info file."""
    kwt = []
    for k in kpoints:
        for i in range(len(kref)):
            if k == kref[i][1:4]:
                kwt.append(kref[i][-1])
                break
    return kwt
    
def clean_kwt_from_istate(istate, kwt):
    """clean the kpoints weight from istate, because in istate.info file the occpuation
    is already multiplied by the weight of kpoints, kwt"""
    for i in range(len(istate)): # loop over spin
        for j in range(len(istate[i])): # loop over kpoints
            istate[i][j][-1] /= kwt[j]
    return istate

def decompose_istate(istate):
    """seperate the istate into energy and occupation,
    indexed by spin and kpoint, respectively.
    istate is the data has structure like:
    [ispin][ik][iband][icol],
    icol = 0: meaningless index
    icol = 1: energy in Ry
    icol = 2: occupation number
    """
    energy, occ = [], []
    for i in range(len(istate)): # loop over spin
        energy.append([])
        occ.append([])
        for j in range(len(istate[i])): # loop over kpoints
            energy[i].append([istate[i][j][k][1] for k in range(len(istate[i][j]))])
            occ[i].append([istate[i][j][k][2] for k in range(len(istate[i][j]))])
    return energy, occ



def pair(collected: list, excluded: dict):
    """pair results from collect_jobs function, store all "the same" data
     together in one list. The "same" means all parameters except defined
    keys in the dict excluded are the same.

    - For PW vs LCAO, it should be:
        excluded = {"ParamSet": ["basis_type"], "AtomSpecies": ["nao"]}
    - For PW ecutwfc convergence test, it should be:
        excluded = {"ParamSet": ["ecutwfc"]}
      or if ultrasoft pseudopotential is used:
        excluded = {"ParamSet": ["ecutwfc", "ecutrho"]}
    """
