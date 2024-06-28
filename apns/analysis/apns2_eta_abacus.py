"""Band structure similarity test for ABACUS"""

def collect(folder: str):
    """Collect APNS jobs for band structure similarity calculation"""
    import os, json, re
    print("* * * Collect APNS jobs * * *".center(100))
    result = []
    for root, _, files in os.walk(folder):
        for file in files:
            if re.match(r"(running_)(\w+)(\.log)", file):
                parent = os.path.dirname(root)
                print(f"Collecting from {parent}")
                try:
                    with open(os.path.join(parent, "description.json"), "r") as f:
                        desc = json.load(f)
                except FileNotFoundError:  
                    print(f"File not found in {parent}")  
                    continue  
                except json.JSONDecodeError:  
                    print(f"Invalid JSON format in {parent}")  
                    continue
                temp = _parse_folder(root)
                if not temp: continue
                ekb, occ, kwt = temp
                result.append((desc, ekb, occ, kwt))
                # note! desc is a dict, ekb and occ are nested like [ispin][ik][iband]
                # kwt is a list of kpoints weight
    return result

def _parse_folder(folder: str):
    """for eta-calculation specific, parse the folder to get the ekb, occ, kwt data.
    
    Args:
        - folder: str, the folder to be parsed, should contain the following files:
            - istate.info
            - kpoints
    
    Returns:
        - ekb: list, the energy of bands, in the form of [ispin][ik][iband]
        - occ: list, the occupation number, in the form of [ispin][ik][iband]
        - kwt: list, the weight of kpoints, in the form of [ik]
    """
    import os, json
    from apns.analysis.postprocess.read_abacus_out import read_istate, read_kpoints
    fistate = os.path.join(folder, "istate.info")

    try:
        temp = read_istate(fistate) # istate to be indiced by nspin, then k
    except FileNotFoundError:
        print(f"istate.info file not found in {folder}")
        return None
    except json.JSONDecodeError:
        print(f"Invalid format in istate.info file in {folder}")
        return None
    
    if temp is None:
        return None
    istate, kpoints = temp
    fkpoints = os.path.join(folder, "kpoints")
    kref = read_kpoints(fkpoints)[0]
    # because the occ information in istate is multiplied by kpoints weight,
    # need divide it.
    kwt = _lookup_kwt(kpoints, kref)
    nspin = len(istate)
    assert nspin in [1, 2], f"nspin should be 1 or 2, but got {nspin}"
    #istate = _clean_kwt_from_istate(istate, kwt)
    ekb, occ = _decompose_istate(istate)
    return ekb, occ, kwt

def _lookup_kwt(kpoints, kref):
    """lookup the kpoints weight by comparing between the coordinate extracted
    from istate.info with the kpoints file extracted kref.
    kpoints has the simple format like: 
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.1], ...]
    kref is a little bit nested like:
    [[1, 0.0, 0.0, 0.0, 0.1], [2, 0.0, 0.0, 0.1, 0.1], ...]
    in which the first index is just the index of kpoints, the last is the weight.
    
    will return a list of kpoints weight, in the same order as kpoints read in 
    the istate.info file."""
    def eq(k1, k2):
        return abs(k1[0] - k2[0]) < 1e-5 and abs(k1[1] - k2[1]) < 1e-5 and abs(k1[2] - k2[2]) < 1e-5
    # convert kref to nested list from list of ndarray
    kref_ = [k[1:4] for k in kref]
    wref_ = [k[4] for k in kref]
    # then compare the kpoints with kref and find the weight
    out = []
    for i, k in enumerate(kpoints):
        if eq(k, kref_[i]): # ideal shortcut
            out.append(wref_[i])
            continue
        else: # otherwise search all in kref_
            j = -1
            for j_, k_ in enumerate(kref_):
                if eq(k, k_):
                    j = j_
                    break
            if j != -1:
                out.append(wref_[j])
                continue
        raise ValueError(f"Cannot find the kpoints weight for {k}")
    assert len(out) == len(kpoints), "kpoints weight should be the same as the number of kpoints"
    return out
    
def _clean_kwt_from_istate(istate, kwt):
    """clean the kpoints weight from istate, because in istate.info file the occpuation
    is already multiplied by the weight of kpoints, kwt"""
    # check length
    nks = [len(ist) for ist in istate]
    assert all([n == nks[0] for n in nks]), "nks should be the same for all spins"
    nks = nks[0]
    assert nks == len(kwt), "kpoints weight should be the same as the number of kpoints"
    for i in range(len(istate)): # loop over spin
        for j in range(len(istate[i])): # loop over kpoints
            w = kwt[j]
            for k in range(len(istate[i][j])): # loop over bands
                _, e, occ = istate[i][j][k]
                if abs(w) < 1e-4 and abs(occ) > 1e-3:
                    raise ValueError(f"Weight is 0, but occupation is not 0: {occ}")
                istate[i][j][k][2] = occ / w if abs(w) > 1e-4 else 0.0
    return istate

def _decompose_istate(istate):
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

def _desc_equal_pw_vs_lcao(desc1: dict, desc2: dict) -> bool:
    """calculate the difference between two description.json contents,
    only following are admitted to be different:
    desc["DFTParamSet"]["basis_type"]
    desc["AtomSpecies"][i]["nao"]
    """
    from apns.analysis.apns2_utils import cal_desc_diff
    diff = cal_desc_diff(desc1, desc2)
    # I doubt whether it is really needed to confine the AtomSpecies
    # must be different in nao for pw and lcao calculation...
    # I would like to skip this check, presently
    # if set(diff.keys()) != {"DFTParamSet", "AtomSpecies"}:
    #     return False
    if set(diff.keys()) > {"DFTParamSet", "AtomSpecies"}:
        return False # but at least other keys should not be different
    if set(diff["DFTParamSet"].get("basis_type", None)) != {"pw", "lcao"}:
        return False
    # for i in range(len(diff["AtomSpecies"])):
    #     if set(diff["AtomSpecies"][i].keys()) != {"nao"}:
    #         return False
    return True

def _desc_equal_between_pw(desc1: dict, desc2: dict) -> bool:
    """calculate the difference between two description.json contents,
    only following are admitted to be different:
    desc["DFTParamSet"]["ecutwfc"] (optionally, ecutrho)
    """
    from apns.analysis.apns2_utils import cal_desc_diff
    diff = cal_desc_diff(desc1, desc2)
    if set(diff.keys()) != {"DFTParamSet"}:
        return False
    if set(diff["DFTParamSet"].keys()) + {"ecutrho"} != {"ecutwfc", "ecutrho"}:
        return False
    return True

def _desc_equal_between_lcao(desc1: dict, desc2: dict) -> bool:
    """calculate the difference between two description.json contents,
    only following are admitted to be different:
    desc["AtomSpecies"][i]["nao"]
    """
    from apns.analysis.apns2_utils import cal_desc_diff
    diff = cal_desc_diff(desc1, desc2)
    if set(diff.keys()) != {"AtomSpecies"}:
        return False
    for i in range(len(diff["AtomSpecies"])):
        if set(diff["AtomSpecies"][i].keys()) != {"nao"}:
            return False
    return True

def pair_pw_lcao(collected: list):
    """pair the pw and lcao data from collect function returned value"""
    paired = []
    pw = [(desc, ekb, occ, kwt) for desc, ekb, occ, kwt in collected if desc["DFTParamSet"].get("basis_type", "pw") == "pw"]
    lcao = [(desc, ekb, occ, kwt) for desc, ekb, occ, kwt in collected if desc["DFTParamSet"].get("basis_type", "pw") == "lcao"]
    print(f"# of pw tests: {len(pw)}")
    print(f"# of lcao tests: {len(lcao)}")
    # then pair the pw and lcao data
    for pw_i in pw:
        for lcao_i in lcao:
            desc_pw, desc_lcao = pw_i[0], lcao_i[0]
            if _desc_equal_pw_vs_lcao(desc_pw, desc_lcao):
                paired.append([pw_i, lcao_i])
    print(f"# of paired tests: {len(paired)}")
    return paired

def pair_between_pw(collected: list):
    """pair the pw data from collect function returned value"""
    paired = []
    pw = [(desc, ekb, occ, kwt) for desc, ekb, occ, kwt in collected if desc["DFTParamSet"].get("basis_type", "pw") == "pw"]
    # then pair the pw data
    paired.append([pw[0]])
    for i in range(1, len(pw)):
        for j in range(len(paired)):
            if _desc_equal_between_pw(pw[i][0], paired[j][0][0]):
                paired[j].append(pw[i])
                break
        else:
            paired.append([pw[i]])
    return paired

def pair_between_lcao(collected: list):
    """pair the lcao data from collect function returned value"""
    paired = []
    lcao = [(desc, ekb, occ, kwt) for desc, ekb, occ, kwt in collected if desc["DFTParamSet"].get("basis_type", "pw") == "lcao"]
    # then pair the lcao data
    paired.append([lcao[0]])
    for i in range(1, len(lcao)):
        for j in range(len(paired)):
            if _desc_equal_between_lcao(lcao[i][0], paired[j][0][0]):
                paired[j].append(lcao[i])
                break
        else:
            paired.append([lcao[i]])
    return paired

def cal_eta_on_pair(pair: list, smear: str, sigma: float, ef_shift: float):
    """calculate the eta value for the pair of jobs. Will return a matrix of eta values"""
    from apns.analysis.apns2_eta_utils import delta_band, cal_nelec
    import numpy as np
    eta, desc = [], []
    for i in range(len(pair)):
        d, ekb, occ, kwt = pair[i]
        nelec = cal_nelec(occ, kwt)
        row = [delta_band(ekb, pair[j][1], nelec, kwt, smear, sigma, ef_shift) \
               for j in range(i + 1, len(pair))]
        row = [0.0] * i + row
        eta.append(row)
        desc.append(d)
    eta = np.array(eta)
    return desc, eta + eta.T

import unittest
class APNS2EtaABACUSTest(unittest.TestCase):
    def test_collect_jobs(self):
        pass

    def test_lookup_kwt(self):
        kpoints = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.1], [0.0, 0.0, 0.2], [0.0, 0.0, 0.3]]
        kref = [[1, 0.0, 0.0, 0.0, 0.1], [2, 0.0, 0.0, 0.1, 0.1], [3, 0.0, 0.0, 0.2, 0.1], [4, 0.0, 0.0, 0.3, 0.1]]
        self.assertEqual(_lookup_kwt(kpoints, kref), [0.1, 0.1, 0.1, 0.1])

    def test_clean_kwt_from_istate(self):
        istate = [
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]],
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]]
        ]
        kwt = [0.1, 0.1, 0.1, 0.1]
        self.assertEqual(_clean_kwt_from_istate(istate, kwt), [
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]],
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]]
        ])

    def test_decompose_istate(self):
        istate = [
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]],
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]]
        ]
        self.assertEqual(_decompose_istate(istate), (
            [[[0.0, 0.0], [0.0, 0.0]], [[0.0, 0.0], [0.0, 0.0]]],
            [[[0.0, 0.1], [0.0, 0.1]], [[0.0, 0.1], [0.0, 0.1]]]
        ))

    def test_desc_equal_pw_vs_lcao(self):
        desc1 = {
            "DFTParamSet": {"basis_type": "pw"},
            "AtomSpecies": [{"nao": 1}, {"nao": 2}]
        }
        desc2 = {
            "DFTParamSet": {"basis_type": "lcao"},
            "AtomSpecies": [{"nao": 1}, {"nao": 2}]
        }
        self.assertTrue(_desc_equal_pw_vs_lcao(desc1, desc2))

    def test_desc_equal_between_pw(self):
        desc1 = {
            "DFTParamSet": {"ecutwfc": 100},
            "AtomSpecies": [{"nao": 1}, {"nao": 2}]
        }
        desc2 = {
            "DFTParamSet": {"ecutwfc": 200},
            "AtomSpecies": [{"nao": 1}, {"nao": 2}]
        }
        self.assertTrue(_desc_equal_between_pw(desc1, desc2))

    def test_desc_equal_between_lcao(self):
        desc1 = {
            "DFTParamSet": {"basis_type": "lcao"},
            "AtomSpecies": [{"nao": 1}, {"nao": 2}]
        }
        desc2 = {
            "DFTParamSet": {"basis_type": "lcao"},
            "AtomSpecies": [{"nao": 1}, {"nao": 3}]
        }
        self.assertTrue(_desc_equal_between_lcao(desc1, desc2))

    def test_pair_pw_lcao(self):
        collected = [
            (
                {"DFTParamSet": {"basis_type": "pw"}, "AtomSpecies": [{"nao": 1}, {"nao": 2}]},
                [[0.0, 0.0], [0.0, 0.0]],
                [[0.0, 0.1], [0.0, 0.1]],
                [0.1, 0.1]
            ),
            (
                {"DFTParamSet": {"basis_type": "lcao"}, "AtomSpecies": [{"nao": 1}, {"nao": 2}]},
                [[0.0, 0.0], [0.0, 0.0]],
                [[0.0, 0.1], [0.0, 0.1]],
                [0.1, 0.1]
            )
        ]
        self.assertEqual(pair_pw_lcao(collected), [[collected]])

    def test_pair_between_pw(self):
        collected = [
            (
                {"DFTParamSet": {"ecutwfc": 100}, "AtomSpecies": [{"nao": 1}, {"nao": 2}]},
                [[0.0, 0.0], [0.0, 0.0]],
                [[0.0, 0.1], [0.0, 0.1]],
                [0.1, 0.1]
            ),
            (
                {"DFTParamSet": {"ecutwfc": 200}, "AtomSpecies": [{"nao": 1}, {"nao": 2}]},
                [[0.0, 0.0], [0.0, 0.0]],
                [[0.0, 0.1], [0.0, 0.1]],
                [0.1, 0.1]
            )
        ]
        self.assertEqual(pair_between_pw(collected), [[collected[0], collected[1]]])

    def test_pair_between_lcao(self):
        collected = [
            (
                {"DFTParamSet": {"basis_type": "lcao"}, "AtomSpecies": [{"nao": 1}, {"nao": 2}]},
                [[0.0, 0.0], [0.0, 0.0]],
                [[0.0, 0.1], [0.0, 0.1]],
                [0.1, 0.1]
            ),
            (
                {"DFTParamSet": {"basis_type": "lcao"}, "AtomSpecies": [{"nao": 1}, {"nao": 3}]},
                [[0.0, 0.0], [0.0, 0.0]],
                [[0.0, 0.1], [0.0, 0.1]],
                [0.1, 0.1]
            )
        ]
        self.assertEqual(pair_between_lcao(collected), [[collected[0], collected[1]]])

    def test_cal_eta_on_pair(self):
        pass

if __name__ == "__main__":
    unittest.main()