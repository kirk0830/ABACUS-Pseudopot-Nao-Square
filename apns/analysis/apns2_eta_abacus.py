"""Band structure similarity test for ABACUS"""

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

def desc_equal_pw_vs_lcao(desc1: dict, desc2: dict) -> bool:
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

def desc_equal_between_pw(desc1: dict, desc2: dict) -> bool:
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

def desc_equal_between_lcao(desc1: dict, desc2: dict) -> bool:
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

def pair_pw_vs_lcao(collected: list):
    """pair the pw and lcao data from collect_jobs function returned value"""
    paired = []
    pw = [(desc, ekb, occ, kwt) for desc, ekb, occ, kwt in collected if desc["DFTParamSet"].get("basis_type", "pw") == "pw"]
    lcao = [(desc, ekb, occ, kwt) for desc, ekb, occ, kwt in collected if desc["DFTParamSet"].get("basis_type", "pw") == "lcao"]
    # then pair the pw and lcao data
    for pw_i in pw:
        for lcao_i in lcao:
            if desc_equal_pw_vs_lcao(pw_i[0], lcao_i[0]):
                paired.append([pw_i, lcao_i])
                break
    return paired

def pair_between_pw(collected: list):
    """pair the pw data from collect_jobs function returned value"""
    paired = []
    pw = [(desc, ekb, occ, kwt) for desc, ekb, occ, kwt in collected if desc["DFTParamSet"].get("basis_type", "pw") == "pw"]
    # then pair the pw data
    paired.append([pw[0]])
    for i in range(1, len(pw)):
        for j in range(len(paired)):
            if desc_equal_between_pw(pw[i][0], paired[j][0][0]):
                paired[j].append(pw[i])
                break
        else:
            paired.append([pw[i]])
    return paired

def pair_between_lcao(collected: list):
    """pair the lcao data from collect_jobs function returned value"""
    paired = []
    lcao = [(desc, ekb, occ, kwt) for desc, ekb, occ, kwt in collected if desc["DFTParamSet"].get("basis_type", "pw") == "lcao"]
    # then pair the lcao data
    paired.append([lcao[0]])
    for i in range(1, len(lcao)):
        for j in range(len(paired)):
            if desc_equal_between_lcao(lcao[i][0], paired[j][0][0]):
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
        self.assertEqual(lookup_kwt(kpoints, kref), [0.1, 0.1, 0.1, 0.1])

    def test_clean_kwt_from_istate(self):
        istate = [
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]],
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]]
        ]
        kwt = [0.1, 0.1, 0.1, 0.1]
        self.assertEqual(clean_kwt_from_istate(istate, kwt), [
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]],
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]]
        ])

    def test_decompose_istate(self):
        istate = [
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]],
            [[[0, 0.0, 0.0, 0.0], [0, 0.0, 0.0, 0.1]], [[0, 0.0, 0.0, 0.2], [0, 0.0, 0.0, 0.3]]]
        ]
        self.assertEqual(decompose_istate(istate), (
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
        self.assertTrue(desc_equal_pw_vs_lcao(desc1, desc2))

    def test_desc_equal_between_pw(self):
        desc1 = {
            "DFTParamSet": {"ecutwfc": 100},
            "AtomSpecies": [{"nao": 1}, {"nao": 2}]
        }
        desc2 = {
            "DFTParamSet": {"ecutwfc": 200},
            "AtomSpecies": [{"nao": 1}, {"nao": 2}]
        }
        self.assertTrue(desc_equal_between_pw(desc1, desc2))

    def test_desc_equal_between_lcao(self):
        desc1 = {
            "DFTParamSet": {"basis_type": "lcao"},
            "AtomSpecies": [{"nao": 1}, {"nao": 2}]
        }
        desc2 = {
            "DFTParamSet": {"basis_type": "lcao"},
            "AtomSpecies": [{"nao": 1}, {"nao": 3}]
        }
        self.assertTrue(desc_equal_between_lcao(desc1, desc2))

    def test_pair_pw_vs_lcao(self):
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
        self.assertEqual(pair_pw_vs_lcao(collected), [[collected]])

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