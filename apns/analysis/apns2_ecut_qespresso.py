def is_outdir(folder: str):
    import os
    b = os.path.exists(os.path.join(folder, "out.log"))
    b = b and os.path.exists(os.path.join(folder, "description.json"))
    b = b and os.path.exists(os.path.join(folder, "pwscf.in"))
    return b

def collect(folder: str):
    print("* * * Collect QESPRESSO result * * *".center(100))
    import apns.analysis.postprocess.read_qespresso_out as read_qe
    from apns.analysis.apns2_utils import read_apnsjob_desc, convert_fpp_to_ppid
    from apns.pspot.parse import z_valence
    import os
    import json
    result = {}
    for root, _, files in os.walk(folder):
        if is_outdir(root):
            natom = read_qe.read_natom(os.path.join(root, "out.log"))
            eks = read_qe.read_e(os.path.join(root, "out.log"))
            press = read_qe.read_press(os.path.join(root, "out.log"))
            bs = read_qe.read_istate(os.path.join(root, "out.log"))
            if None in [eks, press, bs]:
                print(f"WARNING: Present APNS job is broken: {root}")
                continue
            atom_species, cellgen = read_apnsjob_desc(os.path.join(root, "description.json"))
            system = os.path.basename(cellgen["config"])
            with open(os.path.join(root, "description.json"), "r") as f:
                desc = json.load(f)
            ecutwfc = desc["DFTParamSet"]["system"]["ecutwfc"]
            pps = [a["pp"] for a in atom_species]
            ppids = [convert_fpp_to_ppid(pp) for pp in pps]
            zvals = [z_valence(os.path.join(root, pp)) for pp in pps]
            s = "\n".join(ppids)
            print(f"""In folder {root}
Structure tested: {system}
Number of atoms: {natom}
ecutwfc: {ecutwfc}
Final Kohn-Sham energy: {eks}
Pressure: {press}
Pseudopotentials are used:\n{s}
""")
            data = {"ecutwfc": ecutwfc, "eks": eks, "pressure": press, "istate": bs,
                    "natom": natom, "z_valence": zvals}
            idx = -1 if result.get(system, None) is None \
                or result[system].get("ppcases", None) is None \
                    or result[system]["ppcases"].count(pps) == 0 \
                else result[system]["ppcases"].index(pps)
            if idx == -1:
                result.setdefault(system, {"ppcases": [], "pptests": []}).setdefault("ppcases", []).append(pps)
                result[system]["pptests"].append([data])
            else:
                result[system]["pptests"][idx].append(data)

    return result

if __name__ == "__main__":
    # import apns.analysis.external_frender.htmls as amaeh
    # element = "Yb"
    # html = amaeh.pseudopotentials(element=element, 
    #                               xc_functional="PBE", 
    #                               software="ABACUS",
    #                               fconv=f"{element}.svg",
    #                               fconvlog=f"{element}_logplot.svg")
    # with open(f"{element}.md", "w") as f:
    #     f.write(html)
    # exit()
    from apns.analysis.apns2_ecut_utils import build_sptc_from_nested
    collected = collect("12588486")
    system_and_stpcs = build_sptc_from_nested(collected)
    result = {}
    for s, stpcs in system_and_stpcs.items():
        for stpc in stpcs:
            pp, data = stpc()
            result[(s, pp)] = data
            ecut_conv = stpc.ecuts[stpc.iconv]
            pp = stpc.pp(as_list=True)
            assert len(pp) == 1, "The pseudopotential should be unique for each test case"
            #update_ecutwfc(pp[0], ecut_conv)
    from apns.analysis.apns2_ecut_utils import plot_log, plot_stack
    flogs = plot_log(result)
    fstacks = plot_stack(result)