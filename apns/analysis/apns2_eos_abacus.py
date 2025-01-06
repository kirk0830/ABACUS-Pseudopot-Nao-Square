def is_outdir(files: list, fromcal: str):
    """for determining whether a OUT.ABACUS dir is really valid, there are files must present:
    - running_{fromcal}.log
    - INPUT
    - STRU_READIN_ADJUST.cif or STRU.cif
    - istate.info
    - kpoints
    .
    
    Args:
        files (list): list of filenames in the dir
        fromcal (str): the calculation type, e.g. "scf", "relax", "md"

    Returns:
        bool: whether the dir is valid
    """
    ref = {f"running_{fromcal}.log", "INPUT", "STRU_READIN_ADJUST.cif", "STRU.cif", "istate.info", "kpoints"}
    delta = ref - set(files)
    out = delta == {"STRU_READIN_ADJUST.cif"} or delta == {"STRU.cif"}
    # due to the file name is changed in the latest version of ABACUS
    if not out and delta != ref:
        print(f"Not a valid OUT.ABACUS dir, files diff: {delta}")
    return out

def collect(folder: str, calc_type: str = "relax"):
    """Given folder, extract energy, volume, and number of atoms from ABACUS calculation results.
    The method to determining whether a dir is valid is by checking the presence of files:
    running_{calc_type}.log, INPUT, STRU_READIN_ADJUST.cif or STRU.cif, istate.info, kpoints.
    The data will be stored in a dict, with the key is the system name, and the value is a dict with keys:
    - ppcases: list of pseudopotentials used in the calculation
    - pptests: list of tests corresponding to the pseudopotentials, each test is a dict, either with keys
    `energy`, `volume`, `natom`, or `orbcases` and `orbtests` if the calculation is LCAO. If it is the 
    latter case, then the `orbcases` and `orbtests` will have the same structure as `ppcases` and `pptests`,
    except this time in `orbtests` there are dicts with keys `energy`, `volume`, `natom`.

    Args:
        folder (str): the folder contains the ABACUS calculation results
        calc_type (str, optional): the calculation type. Defaults to "relax".
    
    Returns:
        dict: the dict contains the extracted data
    
    Example:
    1. PW calculation:
    ```json
    {
        "system1": {
            "ppcases": [["pp1", "pp2"], ["pp3", "pp4"]],
            "pptests": [
                [{"energy": 1, "volume": 2, "natom": 3}, {"energy": 4, "volume": 5, "natom": 6}],
                [{"energy": 7, "volume": 8, "natom": 9}, {"energy": 10, "volume": 11, "natom": 12}]
            ]
        },
    }
    ```
    2. LCAO calculation:
    ```json
    {
        "system1": {
            "ppcases": [["pp1", "pp2"], ["pp3", "pp4"]],
            "pptests": [
                [{"orbcases": [["orb1", "orb2"], ["orb3", "orb4"]], "orbtests": [
                    [{"energy": 1, "volume": 2, "natom": 3}, {"energy": 4, "volume": 5, "natom": 6}],
                    [{"energy": 7, "volume": 8, "natom": 9}, {"energy": 10, "volume": 11, "natom": 12}]
                ]},
                {"orbcases": [["orb5", "orb6"], ["orb7", "orb8"]], "orbtests": [
                    [{"energy": 13, "volume": 14, "natom": 15}, {"energy": 16, "volume": 17, "natom": 18}],
                    [{"energy": 19, "volume": 20, "natom": 21}, {"energy": 22, "volume": 23, "natom": 24}]
                ]}
            ]
        },
    }
    ```
    A mixture of PW and LCAO calculation is also supported.
    """
    print("* * * Collect ABACUS result * * *".center(100))
    import apns.analysis.postprocess.read_abacus_out as read_abacus_out
    import os
    from apns.analysis.apns2_utils import read_apnsjob_desc, convert_fpp_to_ppid
    result = {}
    for root, _, files in os.walk(folder):
        if is_outdir(files, calc_type):
            ##############
            # data parse #
            ##############
            natom = read_abacus_out.read_natom_fromlog(os.path.join(root, f"running_{calc_type}.log"))
            eks = read_abacus_out.read_e_fromlog(os.path.join(root, f"running_{calc_type}.log"))
            parent = os.path.dirname(root)
            vol = read_abacus_out.read_volume_fromstru(os.path.join(parent, "STRU"), "A")
            atom_species, cellgen = read_apnsjob_desc(os.path.join(parent, "description.json"))
            system = os.path.basename(cellgen["config"])
            #########################
            # pp and orb info-print #
            #########################
            pps = [a["pp"] for a in atom_species]
            orbs = [a["nao"] for a in atom_species]
            ppids = [": ".join(convert_fpp_to_ppid(pp)) for pp in pps]
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
            ###############
            # data append #
            ###############
            # because we are APNS2 now and want to achieve the universal workflow, here the PW and LCAO calculation
            # needed to be carefully treated.
            data = {"energy": eks, "volume": vol, "natom": natom} # no matter PW or LCAO calculation, the data of interest are the same

            # there will always be pseudopotential, therefore add it to the dict, if needed
            result.setdefault(system, {"ppcases": [], "pptests": []})
            pp_orb = [pps, orbs]
            if pp_orb not in result[system]["ppcases"]:
                result[system]["ppcases"].append(pp_orb)
                result[system]["pptests"].append([])
            icase = result[system]["ppcases"].index(pp_orb)
            result[system]["pptests"][icase].append(data)
    return result

def fit(raw: dict):
    """get the ACWF reference data and fit the EOS for each system in the raw data, calculate
    the delta and store the results in a dict."""
    from apns.analysis.apns2_eos_utils import EquationOfStateSingleTestCase, read_acwf_refdata,\
        ACWF_DATAPATH, ACWF_OXIDE_AEREF, ACWF_UNARIES_AEREF
    import re, os
    unaries = r"[A-Z][a-z]?-X/(SC|BCC|FCC|Diamond)"
    oxides = r"[A-Z][a-z]?-X(O|O2|O3|2O|2O3|2O5)"
    out = []
    for name, data in raw.items():
        # each system corresponds to a collection of "case", but first to tokenize the system
        # to get the refernce data from ACWF
        token = EquationOfStateSingleTestCase.tokenize(name)
        assert re.match(unaries, token) or re.match(oxides, token), f"Error in {name}, token {token} is not valid!"
        acwf = ACWF_UNARIES_AEREF if re.match(unaries, token) else ACWF_OXIDE_AEREF
        acwf = os.path.join(ACWF_DATAPATH, acwf)
        ref = read_acwf_refdata(token, acwf)
        temp = {"name": token, "fcif": name, "data": [{"pp": "AE (averaged over Wien2K and FLEUR)", "orb": None}]}
        if ref is None:
            print(f"Error in {name} (token {token}), no ACWF reference data, skip!")
            continue
        temp["data"][0].update(ref)

        # then do the fit Birch-Murnaghan EOS for each case, it should be noted that the raw
        # data may either be PW or LCAO, so the structure of the data is somewhat complex...
        for case, test in zip(data["ppcases"], data["pptests"]):
            # while it is safe that the "ppcase" will not have different structure...
            # but the "pptest" may have different structure, so we need to check it first (everytime)
            # but it is also possible that there is a mixture of PW and LCAO calculation, in this case
            # there will be a dict with keys "orbcases" and "orbtests" in the "pptest", instead of
            # the dict with keys "energy", "volume", "natom".

            # we first try to collect all {"energy", "volume", "natom"} data, it is a single case of
            # PW calculation, then the rest is LCAO calculation, we will handle with them later.
            pw = [a for a in test if set(a.keys()) == {"energy", "volume", "natom"}] # collect pw cases
            if len(pw) > 0:
                # means there are no data with keys {"energy", "volume", "natom"}, so it is purely the
                # LCAO calculation, we will handle with them later
                pwcase = EquationOfStateSingleTestCase(name, pw, case[0])
                pp, orb, bmfit, delta = pwcase(acwf)
                assert orb is None, "The data is not PW calculation!"
                if bmfit is not None:
                    temp["data"].append({"pp": pp, "orb": orb, "delta": delta, 
                                         "volume": pwcase.volumes, "energy": pwcase.energies, "natoms": pwcase.natom})
                    temp["data"][-1].update(bmfit)
                else:
                    print(f"Error in {name} with {pp}, no Birch-Murnaghan fit!")
            lcao = [a for a in test if set(a.keys()) == {"orbcases", "orbtests"}]
            if len(lcao) > 0:
                for l in lcao:
                    for orbcase, orbtest in zip(l["orbcases"], l["orbtests"]):
                        lcaocase = EquationOfStateSingleTestCase(name, orbtest, case[0], case[1])
                        pp, orb, bmfit, delta = lcaocase(acwf)
                        if bmfit is not None:
                            temp["data"].append({"pp": pp, "orb": orb, "delta": delta, 
                                                 "volume": lcaocase.volumes, "energy": lcaocase.energies, "natoms": pwcase.natom})
                            temp["data"][-1].update(bmfit)
                        else:
                            print(f"Error in {name} with {pp} and {orb}, no Birch-Murnaghan fit!")
        out.append(temp)
    return out

def main(source: str):
    import os
    from apns.analysis.apns2_eos_utils import plot
    raw = collect(source, "relax")
    fitted = fit(raw)
    fdb = os.path.basename(source) + ".json"
    with open(fdb, "w") as f:
        import json
        json.dump(fitted, f, indent=4)
    exit()
    feos = plot(to_plot)
    return feos

if __name__ == "__main__":
    feos = main("/path/to/your/abacus/job/folder")