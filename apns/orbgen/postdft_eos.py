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
        print(f"bm1 is None, cannot calculate delta")
        return None
    if bm2 is None:
        print(f"bm2 is None, cannot calculate delta")
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
    out = {s: {} for s in pw.keys() if s in lcao.keys()} # initialize the output
    for system, data_pw in pw.items():
        ntests_pw = len(data_pw["pptests"])
        assert len(data_pw["ppcases"]) == ntests_pw, "number of tests should be the same as number of cases"
        if system not in lcao:
            continue
        data_lcao = lcao[system]
        for i in range(ntests_pw):         # for one pseudopotential case
            pps, _ = data_pw["ppcases"][i] # for one pseudopotential case
            out[system].setdefault("ppcases", []).append(pps)
            vol, eks_pw, natom = _transp_vol_ener(data_pw["pptests"][i])
            bm_pw = fit_bm(vol, eks_pw, True)
            idxs = [j for j, pporbs_ in enumerate(data_lcao["ppcases"])\
                    if set(pporbs_[0]) == set(pps)]
            vmin, vmax = min(vol), max(vol) # for delta value calculation
            for j in idxs:
                _, orbs = data_lcao["ppcases"][j]
                _, eks_lcao, _ = _transp_vol_ener(data_lcao["pptests"][j])
                bm_lcao = fit_bm(vol, eks_lcao, True)
                delta_ = _cal_delta(bm_pw, bm_lcao, natom, vmin, vmax)
                demin = _cal_demin(eks_lcao, eks_pw) # the smaller, the better
                out[system].setdefault("pptests", {}).setdefault("orbcases", []).append(orbs)
                out[system].setdefault("pptests", {}).setdefault("orbtests", []).append({"delta": delta_, "demin": demin})
    return out

def print_postdft(result: dict):
    """print the post-DFT results"""
    from apns.analysis.apns2_utils import convert_forb_to_orbid, convert_fpp_to_ppid
    for system, data in result.items():
        print(f"System {system}")
        for i, ppcase in enumerate(data["ppcases"]):
            ppcase = [convert_fpp_to_ppid(pp) for pp in ppcase]
            ppcase = ", ".join(ppcase)
            print(f"  Pseudopotential case {ppcase}")
            for j, orbs in enumerate(data["pptests"]["orbcases"]):
                orbs = [convert_forb_to_orbid(orb) for orb in orbs]
                orbs = ", ".join(orbs)
                print(f"    Orbital case {orbs}")
                delta, demin = data["pptests"]["orbtests"][j]["delta"], data["pptests"]["orbtests"][j]["demin"]
                if delta:
                    print(f"      Delta: {delta*1e3:>.4f} meV/atom")
                print(f"      Demin: {demin:>.4f} eV")

if __name__ == "__main__":
    pw = load("/root/documents/simulation/orbgen/apns-orbgen-project/out/eos_pw.json")
    lcao = load("/root/documents/simulation/orbgen/apns-orbgen-project/out/eos_lcao-v1.0.json")
    out = cal_wrt_pw(pw, lcao)
    print_postdft(out)