def load(src: str):
    from apns.analysis.apns2_eta_abacus import collect
    return collect(src)

def make_pair(pw, lcao):
    from apns.analysis.apns2_eta_abacus import pair_pw_lcao
    to_pair = pw + lcao
    print(f"pairing between pw and lcao: {len(to_pair)} tests")

    return pair_pw_lcao(to_pair)

def cal_wrt_pw(pw, lcao, smearing: str = "gaussian", sigma: float = 0.01, efermi_shift: float = 0):
    """calculate eta value with respect to pw results
    
    Args:
    - pw: list, list of pw data
    - lcao: list, list of lcao data
    - smearing: str, smearing method
    - sigma: float, sigma value for smearing
    - efermi_shift: float, efermi shift value
    
    Returns:
    - dict, the data structure that stores the eta values, has following structure:
    ```json
    {
        "ppcases": [[pp1, pp2, ...], [pp3, pp4, ...], ...],
        "pptests": [
            {
                "orbcases": [[orb1, orb2, ...], [orb3, orb4, ...], ...],
                "orbtests": [(eta1, eta_max1), (eta2, eta_max2), ...]
            },
            {
                "orbcases": [[orb5, orb6, ...], [orb7, orb8, ...], ...],
                "orbtests": [(eta3, eta_max3), (eta4, eta_max4), ...]
            },
            ...
        ]
    }
    ```
    """

    nelec_thr = 1e-6

    # Python basic module
    import os
    # APNS module
    from apns.analysis.apns2_eta_utils import delta_band, cal_nelec
    from apns.analysis.apns2_utils import stru_rev_map

    out = {}

    paired = make_pair(pw, lcao)
    # mpid to system reverse mapping
    sysrevmap_ = stru_rev_map("./apns_cache/structures.json", True)

    # calculate each pair
    for pw_, lcao_ in paired:
        desc_pw, ekb_pw, occ_pw, kwt = pw_
        desc_lcao, ekb_lcao, occ_lcao, _ = lcao_
        assert desc_pw["DFTParamSet"]["basis_type"] == "pw", \
            f"basis type should be pw: {desc_pw['DFTParamSet']['basis_type']}"
        assert desc_lcao["DFTParamSet"]["basis_type"] == "lcao", \
            f"basis type should be lcao: {desc_lcao['DFTParamSet']['basis_type']}"
        
        # calculate eta
        nelec_pw = cal_nelec(occ_pw)
        nelec_lcao = cal_nelec(occ_lcao)
        assert abs(nelec_lcao - nelec_pw) <= nelec_thr, \
            f"number of electrons should be the same (at precision {nelec_thr}): {nelec_pw} != {nelec_lcao}"

        eta, eta_max = delta_band(ekb_pw, ekb_lcao, nelec_pw, kwt, smearing, sigma, efermi_shift)

        # save eta
        # basic information
        atom_species, cellgen = desc_lcao["AtomSpecies"], desc_lcao["CellGenerator"]
        system = os.path.basename(cellgen["config"])
        system = sysrevmap_[system]

        # the system layer
        if system not in out:
            out[system] = {}

        # the pp layer
        pps = [as_["pp"] for as_ in atom_species]
        ipps = -1
        if pps not in out[system].setdefault("ppcases", []):
            out[system].setdefault("ppcases", []).append(pps)
        else:
            ipps = out[system]["ppcases"].index(pps)
        if ipps == -1:
            out[system].setdefault("pptests", []).append({"orbcases": [], "orbtests": []})
            ipps = len(out[system]["ppcases"]) - 1

        # the orb layer: one pp layer may have multiple orb
        orbs = [as_["nao"] for as_ in atom_species]
        iorbs = -1
        if orbs not in out[system]["pptests"][ipps]["orbcases"]:
            out[system]["pptests"][ipps]["orbcases"].append(orbs)
        else:
            iorbs = out[system]["pptests"][ipps]["orbcases"].index(orbs)
            
        if iorbs == -1:
            out[system]["pptests"][ipps]["orbtests"].append([(eta, eta_max)])
        else:
            out[system]["pptests"][ipps]["orbtests"][iorbs].append((eta, eta_max))

    return out

def print_postdft(result: dict):
    from apns.analysis.apns2_utils import convert_forb_to_orbid, convert_fpp_to_ppid
    for system, data in result.items():
        print(f"System {system}")
        for i, ppcase in enumerate(data["ppcases"]):
            ppcase = [convert_fpp_to_ppid(pp) for pp in ppcase]
            ppcase = ", ".join(ppcase)
            print(f"  Pseudopotential case {ppcase}")
            for j, orbs in enumerate(data["pptests"][i]["orbcases"]):
                orbs = [convert_forb_to_orbid(orb) for orb in orbs]
                orbs = ", ".join(orbs)
                print(f"    Orbital case {orbs}")
                etas = data["pptests"][i]["orbtests"][j]
                etas, eta_maxs = zip(*etas)
                etas = [f"{eta*1e3:>.4f} meV/atom" for eta in etas]
                eta_maxs = [f"{eta_max*1e3:>.4f} meV/atom" for eta_max in eta_maxs]
                etas = ", ".join(etas)
                eta_maxs = ", ".join(eta_maxs)
                print(f"      Eta: {etas}")
                print(f"      Eta_max: {eta_maxs}")

import unittest
class TestPostdftEta(unittest.TestCase):

    def est_cal_wrt_pw_minimal(self):
        from apns.analysis.apns2_eta_abacus import _parse_folder
        from apns.analysis.apns2_eta_utils import delta_band, cal_nelec
        lcao = _parse_folder("./test_files/lcao/")
        pw = _parse_folder("./test_files/pw/")
        
        ekb_lcao, occ_lcao, kwt = lcao
        ekb_pw, occ_pw, _ = pw
        nelec_lcao = cal_nelec(occ_lcao)
        nelec_pw = cal_nelec(occ_pw)
        self.assertAlmostEqual(nelec_lcao, nelec_pw, delta=1e-5)
        eta, eta_max = delta_band(ekb_pw, ekb_lcao, nelec_pw, kwt, "gaussian", 0.01, 0)
        print(eta, eta_max)

if __name__ == "__main__":
    
    unittest.main(exit=False)

    import os
    # select jobgroup
    group = "sus"

    # import data
    root = "/root/documents/simulation/orbgen/apns-orbgen-project/"
    pw = load(os.path.join(root, f"pw/{group}_eostest_pw"))
    lcao = load(os.path.join(root, f"lcao-v1.0/{group}_eostest_lcao-v1.0"))

    out = cal_wrt_pw(pw, lcao, "gaussian", 0.01, 10)
    print_postdft(out)
    #print(out)