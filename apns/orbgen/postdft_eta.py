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

    nelec_thr = 1e-4

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
        try:
            eta, eta_max = delta_band(ekb_pw, ekb_lcao, nelec_pw, kwt, smearing, sigma, efermi_shift)
        except TypeError as e:
            nbnd_pw = len(ekb_pw[0][0])
            nbnd_lcao = len(ekb_lcao[0][0])
            if nbnd_pw > nbnd_lcao:
                print(f"""{e} in calculating eta for
PW: {desc_pw['CellGenerator']['config']}, with 
LCAO: {desc_lcao['CellGenerator']['config']}
calculation. This is always because of default nbands setting strategy inconsistent between PW and LCAO that
nbands must be smaller than or equal to the number of basis functions. Will truncate PW nbands to match LCAO
nbands.""")
                ekb_pw = [[band[:nbnd_lcao] for band in spin] for spin in ekb_pw]
                eta, eta_max = delta_band(ekb_pw, ekb_lcao, nelec_pw, kwt, smearing, sigma, efermi_shift)
                continue
        eta, eta_max = float(eta), float(eta_max)

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
            out[system]["pptests"][ipps]["orbtests"].append([[eta, eta_max]])
        else:
            out[system]["pptests"][ipps]["orbtests"][iorbs].append([eta, eta_max])

    return out

def print_postdft(result: dict):
    from apns.analysis.apns2_utils import convert_forb_to_orbid, convert_fpp_to_ppid
    for system, data in result.items():
        print(f"System {system}")
        for i, ppcase in enumerate(data["ppcases"]):
            ppcase = [": ".join(convert_fpp_to_ppid(pp)) for pp in ppcase]
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

def _chessboard(fname: str, ppcase: str, pptest: dict, subplots: list = ["eta", "etamax"], scale: float = 1):
    """draw a chessboard plot for the eta and eta_max values for one ppcase.
    will first get all numerical orbital information and extract rcut, orbconf two dimensional data, 
    then eta and eta_max will be represented by color

    Args:
        ppcase: dict, orbtest data of one pseudopotential, in which orbcases and orbtests must be included  
        scale: float, scale of the plot, default is 1
    
    Returns:
        None
    
    """
    import os, re
    import numpy as np
    import matplotlib.pyplot as plt

    assert "orbcases" and "orbtests" in pptest, "pptest should have orbcases and orbtests"
    orbcases = pptest["orbcases"]
    orbtests = pptest["orbtests"]
    assert len(orbcases) == len(orbtests), "orbcases and orbtests should have the same length"
    assert len(orbcases) > 0, "orbcases should not be empty" # then length of orbtests can be indirectly asserted

    orbcases = [[os.path.basename(orb) for orb in orbs] for orbs in orbcases] # will only support the case that len(orbs) == 1
    assert all(len(orbs) == 1 for orbs in orbcases), "only support the case that len(orbs) == 1"
    orbcases = [orb[0] for orb in orbcases]
    orbpat = r"[A-Z][a-z]?_gga_(\d+)au_.*Ry_(\d+\w+)\.orb"
    orbcases = [re.match(orbpat, orb) for orb in orbcases]

    assert all(orb_ is not None for orb_ in orbcases), "orbconf should match the pattern"
    orbcases = [(int(orb_.group(1)), orb_.group(2)) for orb_ in orbcases]

    # use the following to allocate a figure
    nrcuts = len(set([orb[0] for orb in orbcases]))
    norbconfs = len(set([orb[1] for orb in orbcases]))
    # build map from rcut to ircut, orbconf to iorbconf
    rcut_map = {rcut: i for i, rcut in enumerate(sorted(set([orb[0] for orb in orbcases])))}
    orbconf_map = {orbconf: i for i, orbconf in enumerate(sorted(set([orb[1] for orb in orbcases])))}
    
    # build data
    # data eta is indexed by orbtests[iorbtest][itest][0], first shrink to [iorbtest][0] by taking average
    eta = [[vt[0] for vt in orbt] for orbt in orbtests]
    eta = [sum(orbt)/len(orbt)*1e3 for orbt in eta]
    eta_matrix = np.zeros((nrcuts, norbconfs))
    for iorbcase, (rcut, orbconf) in enumerate(orbcases):
        eta_matrix[rcut_map[rcut], orbconf_map[orbconf]] = eta[iorbcase]

    eta_max = [[vt[1] for vt in orbt] for orbt in orbtests]
    eta_max = [sum(orbt)/len(orbt)*1e3 for orbt in eta_max]
    eta_max_matrix = np.zeros((nrcuts, norbconfs))
    for iorbcase, (rcut, orbconf) in enumerate(orbcases):
        eta_max_matrix[rcut_map[rcut], orbconf_map[orbconf]] = eta_max[iorbcase]

    # plot
    fig, ax = plt.subplots(1, len(subplots), figsize=(5*len(subplots)*scale, 6*scale))
    # set suptitle
    ppcase = [os.path.basename(pp) for pp in ppcase]
    ppcase = "_AND_".join(ppcase)
    fig.suptitle(ppcase, fontsize=14*scale)
    
    if "eta" in subplots:
        ax_ = ax if len(subplots) == 1 else ax[0]
        im = ax_.imshow(eta_matrix*1e3, cmap="coolwarm", interpolation="nearest")
        ax_.set_title("$\eta$ (meV)", fontsize=12*scale)
        ax_.set_xticks(np.arange(norbconfs))
        ax_.set_yticks(np.arange(nrcuts))
        ax_.set_xticklabels(sorted(orbconf_map.keys()), fontsize=10*scale)
        ax_.set_yticklabels(sorted(rcut_map.keys()), fontsize=10*scale)
        ax_.set_xlabel("Orbital Configuration", fontsize=12*scale)
        ax_.set_ylabel("rcut (au)", fontsize=12*scale)
        plt.setp(ax_.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        for irow in range(nrcuts):
            for icol in range(norbconfs):
                text = f"{eta_matrix[irow, icol]:.4f}" if eta_matrix[irow, icol] != 0 else "NO DATA"
                ax_.text(icol, irow, text, ha="center", va="center", color="black", fontsize=10*scale)
        # block edge
        for irow in range(nrcuts+1):
            ax_.axhline(irow-0.5, color="black", lw=1)
        for icol in range(norbconfs+1):
            ax_.axvline(icol-0.5, color="black", lw=1)
    
    if "etamax" in subplots:
        ax_ = ax if len(subplots) == 1 else ax[1]
        im = ax_.imshow(eta_max_matrix*1e3, cmap="coolwarm", interpolation="nearest")
        ax_.set_title("$\eta_{max}$ (meV)", fontsize=12*scale)
        ax_.set_xticks(np.arange(norbconfs))
        ax_.set_yticks(np.arange(nrcuts))
        ax_.set_xticklabels(sorted(orbconf_map.keys()), fontsize=10*scale)
        ax_.set_yticklabels(sorted(rcut_map.keys()), fontsize=10*scale)
        ax_.set_xlabel("Orbital Configuration", fontsize=12*scale)
        ax_.set_ylabel("rcut (au)", fontsize=12*scale)
        plt.setp(ax_.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        for irow in range(nrcuts):
            for icol in range(norbconfs):
                text = f"{eta_max_matrix[irow, icol]:.4f}" if eta_max_matrix[irow, icol] != 0 else "NO DATA"
                ax_.text(icol, irow, text, ha="center", va="center", color="black", fontsize=10*scale)
        # block edge
        for irow in range(nrcuts+1):
            ax_.axhline(irow-0.5, color="black", lw=1)
        for icol in range(norbconfs+1):
            ax_.axvline(icol-0.5, color="black", lw=1)

    fig.tight_layout()
    plt.savefig(fname)
    plt.close()

def plot_chessboard(result: dict, subplots: list = ["eta", "etamax"]):
    for system, data in result.items():
        for i, ppcase in enumerate(data["ppcases"]):
            ppcase = [os.path.basename(pp) for pp in ppcase]
            pptest = data["pptests"][i]
            _chessboard(f"{system}_{i}.png", ppcase, pptest, subplots)

import unittest
class TestPostdftEta(unittest.TestCase):

    def test_cal_wrt_pw_minimal(self):
        from apns.analysis.apns2_eta_abacus import _parse_folder
        from apns.analysis.apns2_eta_utils import delta_band, cal_nelec
        lcao = _parse_folder("/root/documents/simulation/abacus/Ag_eos_v2.1/lcao/OUT.ABACUS")
        pw = _parse_folder("/root/documents/simulation/abacus/Ag_eos_v2.1/pw/OUT.ABACUS")
        
        ekb_lcao, occ_lcao, kwt = lcao
        ekb_pw, occ_pw, _ = pw
        nelec_lcao = cal_nelec(occ_lcao)
        nelec_pw = cal_nelec(occ_pw)
        self.assertAlmostEqual(nelec_lcao, nelec_pw, delta=1e-5)
        eta, eta_max = delta_band(ekb_pw, ekb_lcao, nelec_pw, kwt, "gaussian", 0.01, 0)
        print(eta, eta_max)

if __name__ == "__main__":
    
    unittest.main(exit=False)

    import os, json
    # select jobgroup
    group = "short-sp.part2"
    version_lcao = "2.1"

    fjson = f"eta10_lcao-v{version_lcao}.json"
    if not os.path.exists(fjson):
        data = {}
    else:
        with open(fjson, "r") as f:
            data = json.load(f)

    # uncomment the following two lines to plot chessboard, comment out to parse source data
    plot_chessboard(data, ["eta"])
    exit()
    
    # import data
    root = "/root/documents/simulation/orbgen/apns-orbgen-project/eos_test"
    pw = load(os.path.join(root, f"pw/{group}_eostest_pw"))
    lcao = load(os.path.join(root, f"lcao-v{version_lcao}/{group}_eostest_lcao-v{version_lcao}"))

    out = cal_wrt_pw(pw, lcao, "gaussian", 0.01, 10)
    data.update(out)
    with open(fjson, "w") as f:
        json.dump(data, f)
    
    print_postdft(out)
    