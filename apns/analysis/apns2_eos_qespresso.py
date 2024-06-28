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
    import apns.pspot.parse as ppparse
    import os, json
    result = {}
    for root, _, files in os.walk(folder):
        if is_outdir(root):
            natom = read_qe.read_natom(os.path.join(root, "out.log"))
            eks = read_qe.read_e(os.path.join(root, "out.log"))
            vol = read_qe.read_volume(os.path.join(root, "out.log"))
            atom_species, cellgen = read_apnsjob_desc(os.path.join(root, "description.json"))
            system = os.path.basename(cellgen["config"])
            pps = [a["pp"] for a in atom_species]
            ppids = [convert_fpp_to_ppid(pp) for pp in pps]
            s = "\n".join(ppids)
            print(f"""In folder {root}
Structure tested: {system}
Number of atoms: {natom}
Final Kohn-Sham energy: {eks}
Volume: {vol} A^3
Pseudopotentials are used:\n{s}
""")
            data = {"eks": eks, "volume": vol, "natom": natom}
            idx = -1 if result.get(system) is None \
                or result[system].get("ppcases") is None \
                    or result[system]["ppcases"].count(pps) == 0 \
                else result[system]["ppcases"].index(pps)
            if idx == -1:
                result.setdefault(system, {"ppcases": [], "pptests": []}).setdefault("ppcases", []).append(pps)
                result[system]["pptests"].append([data])
            else:
                result[system]["pptests"][idx].append(data)
    return result

def prepare(folder: str):
    from apns.analysis.apns2_eos_utils import EOSSingleCase, read_acwf_refdata
    jobs = collect(folder)
    result = {}
    for system, data in jobs.items():
        result.setdefault(system, {}).setdefault("ppcases", system)
        token = EOSSingleCase.tokenize(system)
        ref = read_acwf_refdata(token)
        result[system].setdefault("AEref", ref)
        for i, ppcase in enumerate(data["ppcases"]):
            case = EOSSingleCase(system, ppcase, data["pptests"][i])
            pp, bmfit, delta = case()
            result[system][pp] = bmfit
            result[system][pp].update({"delta": delta, "volume": case.volumes, "energy": case.energies})
    return result

def main(folder: str):
    from apns.analysis.apns2_eos_utils import plot
    to_plot = prepare(folder)
    feos = plot(to_plot)
    return feos

