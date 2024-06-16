def is_outdir(files: list, fromcal: str):
    return f"running_{fromcal}.log" and\
           "STRU_ION_D" in files and\
           "INPUT" in files and\
           "STRU_READIN_ADJUST.cif" in files and\
           "istate.info" in files and\
           "kpoints" in files

def collect_jobs(folder: str):
    """collect APNS EOS test jobs data in structure like:
    ```json
    {
        "system": {
            "ppcases": [
                ["pp1", "pp2", ...],
                //...
            ],
            "pptests": [
                [ // pp1
                    {"eks": 1.0, "volume": 10.0, "natom": 1}, // point 1
                    // other volume-energy points...
                ],
                // the other pps data...
            ]
        },
        // other systems...
    }
    ```
    """
    print("* * * Collect ABACUS result * * *".center(100))
    import apns.analysis.postprocess.read_abacus_out as read_abacus_out
    import os
    from apns.analysis.drivers.apns2_utils import read_apnsjob_desc, convert_fpp_to_ppid
    result = {}
    for root, _, files in os.walk(folder):
        if is_outdir(files, "relax"):
            natom = read_abacus_out.read_natom_fromlog(os.path.join(root, "running_relax.log"))
            eks = read_abacus_out.read_e_fromlog(os.path.join(root, "running_relax.log"))
            vol = read_abacus_out.read_volume_fromstru(os.path.join(root, "STRU_ION_D"), "A")
            parent = os.path.dirname(root)
            atom_species, cellgen = read_apnsjob_desc(os.path.join(parent, "description.json"))
            system = os.path.basename(cellgen["config"])
            pps = [a["pp"] for a in atom_species]
            ppids = [convert_fpp_to_ppid(pp) for pp in pps]
            s = "\n".join(ppids)
            print(f"""In folder {parent}
Structure tested: {system}
Number of atoms: {natom}
Final Kohn-Sham energy: {eks}
Volume: {vol} A^3
Pseudopotentials are used:\n{s}
""")
            data = {"eks": eks, "volume": vol, "natom": natom}
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

def main(folder: str):
    import os
    from apns.analysis.drivers.apns2_eos_utils import EOSSingleCase, read_acwf_refdata
    jobs = collect_jobs(folder)
    for system, data in jobs.items():
        token = EOSSingleCase.tokenize(system)
        ref = read_acwf_refdata(token)
        for i, ppcase in enumerate(data["ppcases"]):
            case = EOSSingleCase(system, ppcase, data["pptests"][i])
            pp, bmfit, delta = case()