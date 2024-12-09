def load(src: str):
    """load collected data from a json file
    Args:
        - src: str, the path to the json file that stores the result from run
    
    Returns:
        - dict, the data structure that stores the result
    """
    out = {}

    import json
    from apns.analysis.apns2_utils import stru_rev_map

    with open(src, "r") as f:
        data = json.load(f)

    sysrevmap_ = stru_rev_map("./apns_cache/structures.json", True)
    for system, v in data.items():
        # dangerous operation: will be data collision if the same system name is used in different structure
        k = sysrevmap_[system]
        assert k not in out, f"system name collision: {k}"
        out[k] = v

    return out

def _transp_vol_ener(pptest: list, sort: bool = True, as_dict: bool = False):
    """merge data points yield from multiple times of tests. The `pptest` has the following
    format:
    ```python
    [{"eks": 1, "volume": 1, "natom": 1}, {"eks": 2, "volume": 2, "natom": 2}]
    ```
    Args:
    - pptest: list, list of data points
    - sort: bool, whether to sort the data points along volume

    Returns:
    - volume: list, list of volumes
    - eks: list, list of energies
    - natom: int, number of atoms
    """
    import numpy as np
    eks = [p["eks"] for p in pptest]
    volume = [p["volume"] for p in pptest]
    natom = [p["natom"] for p in pptest]
    # assert natom to have all the same value
    assert all([n == natom[0] for n in natom]), "natom should be the same for all tests"
    natom = natom[0]
    # sort along volume
    if sort:
        idx = np.argsort(volume)
        volume = np.array(volume)[idx].tolist()
        eks = np.array(eks)[idx].tolist()

    out = volume, eks, natom if not as_dict else {"volume": volume, "eks": eks, "natom": natom}
    return out

def _cal_delta(bm1: dict, bm2: dict, natom: int, vmin: float, vmax: float):
    """calculate the delta value between two Birch-Murnaghan fits"""
    from apns.analysis.apns2_eos_utils import delta_value as delta
    if bm1 is None:
        print("bm1 is None, cannot calculate delta")
        return None
    if bm2 is None:
        print("bm2 is None, cannot calculate delta")
        return None
    return delta(bm1, bm2, vmin, vmax, natom)

def _cal_demin(e1, e2):
    """calculate the difference between two minimum energies"""
    if isinstance(e1, list):
        e1 = min(e1)
    if isinstance(e2, list):
        e2 = min(e2)
    return e1 - e2

def cal_wrt_pw(pw: dict, lcao: dict):
    """calculate the delta and demin values for the LCAO results with respect to the PW results
    Args:
        - pw: dict, the data structure that stores the result
        - lcao: dict, the data structure that stores the result
    
    Returns:
        - dict, the data structure that stores the delta and demin values
    """
    from apns.analysis.apns2_eos_utils import fit_birch_murnaghan as fit_bm
    out = {s: {} for s in pw if s in lcao} # initialize the output
    for system, data_pw in pw.items():
        ntests_pw = len(data_pw["pptests"])
        assert len(data_pw["ppcases"]) == ntests_pw, "number of tests should be the same as number of cases"
        if system not in lcao:
            continue
        data_lcao = lcao[system]
        for i in range(ntests_pw):         # for one pseudopotential case
            pps, _ = data_pw["ppcases"][i] # for one pseudopotential case
            out[system].setdefault("ppcases", []).append(pps)
            out[system].setdefault("pptests", []).append({"orbcases": [], "orbtests": []})
            vol, eks_pw, natom = _transp_vol_ener(data_pw["pptests"][i])
            bm_pw = fit_bm(vol, eks_pw, True)
            idxs = [j for j, pporbs_ in enumerate(data_lcao["ppcases"])\
                    if set(pporbs_[0]) == set(pps)]
            vmin, vmax = min(vol), max(vol) # for delta value calculation
            for j in idxs: # for each orbital case
                _, orbs = data_lcao["ppcases"][j]
                _, eks_lcao, _ = _transp_vol_ener(data_lcao["pptests"][j])
                if len(vol) != len(eks_lcao):
                    print(f"number of data points should be the same: {len(vol)} != {len(eks_lcao)}, from test {system}: {orbs}")
                    delta_, demin = None, None
                else:
                    bm_lcao = fit_bm(vol, eks_lcao, True)
                    delta_ = _cal_delta(bm_pw, bm_lcao, natom, vmin, vmax)
                    demin = _cal_demin(eks_lcao, eks_pw) # the smaller, the better
                out[system]["pptests"][-1]["orbcases"].append(orbs)
                out[system]["pptests"][-1]["orbtests"].append({"delta": delta_, "demin": demin})
    return out

def print_postdft(result: dict):
    """print the post-DFT results. This function is written by Github.copilot"""
    from apns.analysis.apns2_utils import convert_forb_to_orbid, convert_fpp_to_ppid
    for system, data in result.items():
        print(f"System {system}")
        for ppcase, pptest in zip(data["ppcases"], data["pptests"]):
            ppcase = [": ".join(convert_fpp_to_ppid(pp)) for pp in ppcase]
            ppcase = ", ".join(ppcase)
            print(f"  Pseudopotential case {ppcase}")
            for orbcase, orbtest in zip(pptest["orbcases"], pptest["orbtests"]):
                orbs = [convert_forb_to_orbid(orb) for orb in orbcase]
                orbs = ", ".join(orbs)
                print(f"    Orbital case {orbs}")
                delta, demin = orbtest["delta"], orbtest["demin"]
                if delta:
                    print(f"      Delta: {delta*1e3:>.4f} meV/atom")
                else:
                    print("      Delta: WARNING: delta is not successfully calculated, this may be due to the failure of fitting the Birch-Murnaghan equation!")
                if demin:
                    print(f"      Demin: {demin:>.4f} eV")
                else:
                    print("      Demin: WARNING: demin is not successfully calculated, this may be due to the failure of fitting the Birch-Murnaghan equation!")

def barplot_postdft(key: str, val: dict):
    """plot the bar plot for the post-DFT results of one system. The system has the following structure:
    ```python
    key: {
        "ppcases": [[...], [...], ...],
        "pptests": {
            "orbcases": [[...], [...], [...], ...],
            "orbtests": [{"delta": ..., "demin": ...}, {...}, {...}]
        }
    }
    ```
    will plot the bar plot for the given key and val. Different ppcases are plot as different subplot,
    then two barplot in each subplot, one for delta, one for demin. x axis are the orbcases (use 
    os.path.basename to wash).
    """
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    from apns.analysis.apns2_utils import convert_forb_to_orbid, convert_fpp_to_ppid

    fig, axs = plt.subplots(len(val["ppcases"]), 2, figsize=(20, 10*len(val["ppcases"])), squeeze=False)
    # suptitle as the key
    fig.suptitle(key, fontsize=20)
    for i, ppcase in enumerate(val["ppcases"]): # for each pseudopotential case
        # build up title
        ppcase = [": ".join(convert_fpp_to_ppid(pp)) for pp in ppcase]
        ppcase = " + ".join(ppcase)

        # extract data, for checking if there is None in orbtests
        delta = [d["delta"] for d in val["pptests"][i]["orbtests"]]
        demin = [d["demin"] for d in val["pptests"][i]["orbtests"]]
        invalid = [d["delta"] is None or d["demin"] is None for d in val["pptests"][i]["orbtests"]]
        if all(invalid):
            print(f"all invalid data for {ppcase}, skip")
            continue

        # delete the invalid data
        delta = [d * 1e3 for j, d in enumerate(delta) if not invalid[j]]
        demin = [d for j, d in enumerate(demin) if not invalid[j]]

        # update orbcases label
        orbcases = [" + ".join([convert_forb_to_orbid(os.path.basename(orb)) for orb in orbs])\
                    for j, orbs in enumerate(val["pptests"][i]["orbcases"]) if not invalid[j]]
        x = range(len(orbcases))
        # the orbcase will be in the format like "Xau, YRy (Z)", first sort by X and then by Z,
        # also sort the delta and demin correspondingly
        orbcases, delta, demin = zip(*sorted(zip(orbcases, delta, demin), 
                                     key=lambda x: (float(x[0].split("au")[0]), 
                                                   x[0].split("au")[1].split(" ")[-1])))
        ######
        # setup the left subplot
        ######
        ax = axs[i, 0]
        ax.set_title(f"Pseudopotential case: {ppcase}", fontsize=17.5)
        # draw, with gradient color, x vertical, y horizontal
        ax.bar(x, delta, color=plt.cm.rainbow(np.linspace(0, 1, len(delta))), log=True)
        ax.set_xticks(x)
        ax.set_xticklabels(orbcases, rotation=45, ha="right")
        ax.set_ylabel("$\Delta$ (meV/atom)", fontsize=15)
        # set black edge
        for bar in ax.patches:
            bar.set_edgecolor("black")
            bar.set_linewidth(0.5)
        ######
        # setup the right subplot
        ######
        ax = axs[i, 1]
        ax.set_title(f"Pseudopotential case: {ppcase}", fontsize=17.5)
        # draw, with gradient color, log scale
        ax.bar(x, demin, color=plt.cm.rainbow(np.linspace(0, 1, len(demin))), log=True)
        ax.set_xticks(x)
        ax.set_xticklabels(orbcases, rotation=45, ha="right")
        ax.set_ylabel("Basis completeness error (eV)", fontsize=15)
        for bar in ax.patches:
            bar.set_edgecolor("black")
            bar.set_linewidth(0.5)

    # readjust the layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"{key}.png")
    plt.close()

    return f"{key}.png"

if __name__ == "__main__":
    
    import json

    # to save data: call apns/analysis/apns2_eos_abacus.py: collect, to get the data can be loaded by function load()
    # Parse
    # flcao = "/root/documents/simulation/orbgen/apns-orbgen-project/eos_test/lcao-v1.0"
    # out = collect(flcao, "scf")
    # with open("lcao-v1.0.json", "w") as f:
    #     json.dump(out, f)
    # flcao = "/root/documents/simulation/orbgen/apns-orbgen-project/eos_test/lcao-v2.0"
    # out = collect(flcao, "scf")
    # with open("lcao-v2.0.json", "w") as f:
    #     json.dump(out, f)
    # flcao = "/root/documents/simulation/orbgen/apns-orbgen-project/eos_test/lcao-v2.1"
    # out = collect(flcao, "scf")

    # from apns.analysis.apns2_eos_abacus import collect
    # flcao = "/root/documents/simulation/orbgen/ZSM5-TPA/test/"
    # out = collect(flcao, "scf")
    # with open("lcao-v3.0-20241028.json", "w") as f:
    #     json.dump(out, f)
    # exit()

    # Post processng
    pw = load("/root/abacus-develop/apns-orbgen-project/eos_test/out/pw.json")
    lcao = load("lcao-v3.0-20241028.json")
    out = cal_wrt_pw(pw, lcao)
    for title, data in out.items():
        barplot_postdft(title, data)
    print_postdft(out)
    with open("SIABv3.0-EOS-CompletenessError.json", "w") as f:
        json.dump(out, f)