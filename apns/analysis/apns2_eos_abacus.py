def is_outdir(files: list, fromcal: str):
    return f"running_{fromcal}.log" and\
           "INPUT" in files and\
           "STRU_READIN_ADJUST.cif" in files and\
           "istate.info" in files and\
           "kpoints" in files

def collect(folder: str, calc_type: str = "relax"):
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
    from apns.analysis.apns2_utils import read_apnsjob_desc, convert_fpp_to_ppid
    result = {}
    for root, _, files in os.walk(folder):
        if is_outdir(files, calc_type):
            natom = read_abacus_out.read_natom_fromlog(os.path.join(root, f"running_{calc_type}.log"))
            eks = read_abacus_out.read_e_fromlog(os.path.join(root, f"running_{calc_type}.log"))
            parent = os.path.dirname(root)
            vol = read_abacus_out.read_volume_fromstru(os.path.join(parent, "STRU"), "A")
            atom_species, cellgen = read_apnsjob_desc(os.path.join(parent, "description.json"))
            system = os.path.basename(cellgen["config"])
            pps = [a["pp"] for a in atom_species]
            orbs = [a["nao"] for a in atom_species]
            ppids = [convert_fpp_to_ppid(pp) for pp in pps]
            ppstr = "\n".join(ppids)
            orbstr = "\n".join(orbs) if all([orb is not None for orb in orbs]) else "none"
            if None in [natom, eks, vol]:
                print(f"Error in {root}")
                continue
            print(f"""In folder {parent}
Structure tested: {system}
Number of atoms: {natom}
Final Kohn-Sham energy: {eks}
Volume: {vol} A^3
Pseudopotentials are used:\n{ppstr}
Numerical atomic orbitals are used:\n{orbstr}
""")        
            pps = [pps, orbs]
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

def prepare(folder: str, calc_type: str = "relax"):
    from apns.analysis.apns2_eos_utils import EOSSingleCase, read_acwf_refdata
    jobs = collect(folder, calc_type)
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
    to_plot = prepare(folder, "scf")
    feos = plot(to_plot)
    return feos

if __name__ == "__main__":
    feos = main("/root/documents/simulation/abacus/12891738")