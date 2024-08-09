def load(src: str):
    import os
    from apns.analysis.apns2_ecoh_utils import pair, cal_e_cohesive
    from apns.analysis.apns2_ecoh_abacus import collect

    out = {}
    data = collect(src)
    paired = pair(data)
    for desc, e_b, e_a in paired:
        natom = len(desc["Cell"]["coords"])
        e_coh = cal_e_cohesive(e_b, e_a, natom)
        system = os.path.basename(desc["CellGenerator"]["config"])
        if system not in out:
            out[system] = {}
        
        pps = [os.path.basename(a["pp"]) for a in desc["AtomSpecies"]]
        ipps = -1
        if pps not in out[system].setdefault("ppcases", []):
            out[system].setdefault("ppcases", []).append(pps)
        else:
            ipps = out[system]["ppcases"].index(pps)

        # if it is pw calculation, the orbs is not necessary to be available
        basis_type = desc.get("DFTParamSet").get("basis_type")
        assert basis_type != None, "Severe error: basis_type cannot be found in description.json"

        if ipps == -1:
            initializer = {"orbcases": [], "orbtests": []} if basis_type == "lcao" else e_coh
            out[system].setdefault("pptests", []).append(initializer)
            ipps = len(out[system]["ppcases"]) - 1
        if basis_type == "pw":
            continue

        orbs = [os.path.basename(a["nao"]) for a in desc["AtomSpecies"]]
        assert all(orbs), "Severe error: nao cannot be found in description.json"
        iorbs = -1
        if orbs not in out[system]["pptests"][ipps]["orbcases"]:
            out[system]["pptests"][ipps]["orbcases"].append(orbs)
        else:
            iorbs = out[system]["pptests"][ipps]["orbcases"].index(orbs)
        
        if iorbs == -1:
            initializer = e_coh
            out[system]["pptests"][ipps]["orbtests"].append(initializer)
        else:
            out[system]["pptests"][ipps]["orbtests"][iorbs] = e_coh
    
    return out

def cal_wrt_pw(pw, lcao):
    out = {}
    for system in lcao:
        out.setdefault(system, {"ppcases": [], "pptests": {
            "orbcases": [], "orbtests": []
        }})
        for ipps, ppcase in enumerate(lcao[system]["ppcases"]):
            if ppcase not in out[system]["ppcases"]:
                out[system]["ppcases"].append(ppcase)
            ipps_ = pw[system]["ppcases"].index(ppcase)
            for iorbs, orbcase in enumerate(lcao[system]["pptests"][ipps]["orbcases"]):
                if orbcase not in out[system]["pptests"]["orbcases"]:
                    out[system]["pptests"]["orbcases"].append(orbcase)
                val_lcao = lcao[system]["pptests"][ipps]["orbtests"][iorbs]
                val_pw = pw[system]["pptests"][ipps_]
                out[system]["pptests"]["orbtests"].append(val_lcao - val_pw)
    return out

if __name__ == "__main__":
    """
    ncores = 16, therefore all calculation with orbitals with less than nbase = 16 are
    failed and should be rerun.
    SZ: 1s(1), 1s1p(4), 1s1p1d(9)
    DZP: 2s1p(5), 2s1p1d(10), 2s2p1d(13), 4s2p1d(15)
    TZDP: 3s2p(11)
    """
    
    import json
    # basis = "lcao-v2.1"

    # src = f"/root/documents/simulation/orbgen/apns-orbgen-project/ecoh_test/{basis}/"
    # fout = f"ecoh_{basis}.json"
    # try:
    #     with open(fout, "r") as f:
    #         out = json.load(f)
    # except FileNotFoundError:
    #     out = {}
    # out.update(load(src))
    # with open(fout, "w") as f:
    #     json.dump(out, f, indent=4)

    