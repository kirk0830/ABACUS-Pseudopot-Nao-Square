"""this file is for greping and rendering svg plots of ecutwfc convergence test for new version of APNS"""

def build_sptc_from_nested(cases: dict):
    """build instances of SystemPspotTestCase from a nested dict. Example input:
    ```python
    cases = {
        "mp-1234.cif": {
            "ppcases": [["Ag.pbe-spn-rrkjus_psl.0.2.3-tot.UPF", "O.pbe-n-rrkjus_psl.0.2.3-tot.UPF"],
                        ["Ag_ONCV_PBE-1.0_fr.upf", "O_ONCV_PBE-1.0_fr.upf"]],
            "pptests": [
                [{"ecutwfc": 30, "eks": -1.0, "pressure": 0.0, "istate": 0.0, "natom": 1, "z_valence": [11, 6]},
                 {"ecutwfc": 40, "eks": -2.0, "pressure": 0.0, "istate": 0.0, "natom": 1, "z_valence": [11, 6]},
                 {"ecutwfc": 50, "eks": -3.0, "pressure": 0.0, "istate": 0.0, "natom": 1, "z_valence": [11, 6]}]
                [{"ecutwfc": 40, "eks": -1.0, "pressure": 0.0, "istate": 0.0, "natom": 1, "z_valence": [11, 6]},
                 {"ecutwfc": 50, "eks": -2.0, "pressure": 0.0, "istate": 0.0, "natom": 1, "z_valence": [11, 6]},
                 {"ecutwfc": 60, "eks": -3.0, "pressure": 0.0, "istate": 0.0, "natom": 1, "z_valence": [11, 6]}]
                ]
        },
        "mp-5678.cif": {...}
    }
    ```
    This function will return a dict whose keys are different systems, and values are lists of SystemPspotTestCase instances.
    Each instance represents a pseudopotential test case for a specific system."""
    result = {}
    for system, data in cases.items():
        ppcases = data["ppcases"]
        for ipt, pptests in enumerate(data["pptests"]):
            pps = ppcases[ipt]
            sptc = SystemPspotTestCase(system, pps, pptests)
            result.setdefault(system, []).append(sptc)
    return result

class SystemPspotTestCase:
    """the class owning test data varies ecutwfc"""
    system: str
    ecuts: list
    natom: int
    zvals: list
    energies: list
    pressures: list
    istates: list
    pps: list

    iconv: int = None
    def __init__(self, system: str, pps: list, cases: list):
        """init from a nested dict, data can be get like:
        ```python
        scratch = [{"ecutwfc": 30, "eks": -1.0, "pressure": 0.0, "istate": 0.0, "natom": 1, "z_valence": [11, 6]},
                   {"ecutwfc": 40, "eks": -2.0, "pressure": 0.0, "istate": 0.0, "natom": 1, "z_valence": [11, 6]},
                   {"ecutwfc": 50, "eks": -3.0, "pressure": 0.0, "istate": 0.0, "natom": 1, "z_valence": [11, 6]}]
        ```"""
        self.system = system
        self.pps = pps
        self.ecuts = [float(c["ecutwfc"]) for c in cases]
        self.natom = cases[0]["natom"]
        assert all([c["natom"] == self.natom for c in cases]), "The number of atoms should be consistent for all ecutwfc tests"
        self.zvals = cases[0]["z_valence"]
        assert all([c["z_valence"] == self.zvals for c in cases]), "The valence of atoms should be consistent for all ecutwfc tests"
        self.energies = [c["eks"] for c in cases]
        self.pressures = [c["pressure"] for c in cases]
        self.istates = [c["istate"] for c in cases]

    def sort(self):
        """sort the data according to ecutwfc"""
        import numpy as np
        idx = np.argsort(self.ecuts)
        self.ecuts = np.array(self.ecuts)[idx].tolist()
        self.energies = np.array(self.energies)[idx].tolist()
        self.pressures = np.array(self.pressures)[idx].tolist()
        self.istates = np.array(self.istates)[idx].tolist()

    def calc_conv(self, ethr: float = 1e-3, pthr: float = 0.1, bsthr: float = 1e-2):
        """calculate the converged ecutwfc, return the index of the converged ecutwfc"""
        from apns.analysis.postprocess.conv.kernel import default_calculator
        from apns.analysis.postprocess.conv.ecutwfc_istate import calculator
        energies = [e/self.natom for e in self.energies]
        de = default_calculator(energies, energies[-1])
        dp = default_calculator(self.pressures, self.pressures[-1])
        dbs = calculator(self.istates, self.istates[-1])
        
        iconv_e = [i for i, e in enumerate(de) if abs(e) <= ethr][0]
        iconv_p = [i for i, p in enumerate(dp) if abs(p) <= pthr][0]
        iconv_bs = [i for i, bs in enumerate(dbs) if abs(bs) <= bsthr][0]
        iconv = max(iconv_e, iconv_p, iconv_bs)
        self.iconv = iconv
        return {"de": de, "dp": dp, "dbs": dbs, "iconv": iconv, "ecutwfc": self.ecuts, "zvals": self.zvals}
    
    def pp(self, as_list: bool = False):
        """return the pseudopotential string"""
        from apns.analysis.drivers.apns2_utils import convert_fpp_to_ppid
        return self.pps if as_list else "|".join([convert_fpp_to_ppid(pp) for pp in self.pps])
    
    def zval(self, as_list: bool = False):
        """return the valence of atoms"""
        return self.zvals if as_list else "|".join(self.zvals)

    def __call__(self):
        self.sort()
        return self.pp(), self.calc_conv()

def update_ecutwfc(pp: str, ecutwfc: float, cache_dir: str = "./apns_cache/ecutwfc.json"):
    import json
    import os
    if os.path.exists(cache_dir):
        with open(cache_dir, "r") as f:
            cache = json.load(f)
    else:
        cache = {}
    cache[pp] = ecutwfc
    with open(cache_dir, "w") as f:
        json.dump(cache, f)

def collect_apnsjob_data(folder: str):
    print("* * * Collect ABACUS result * * *".center(100))
    import apns.analysis.postprocess.read_abacus_out as read_abacus_out
    import apns.pspot.parse as ppparse
    import os, re
    from apns.analysis.drivers.apns2_utils import read_apnsjob_desc, convert_fpp_to_ppid
    result = {}
    for root, _, files in os.walk(folder):
        for file in files:
            if re.match(r"(running_)(\w+)(\.log)", file): # reach the inner most folder like OUT.ABACUS
                natom = read_abacus_out.read_natom_fromlog(os.path.join(root, file))
                eks = read_abacus_out.read_e_fromlog(os.path.join(root, file))
                pressure = read_abacus_out.read_pressure_fromlog(os.path.join(root, file))
                bs = read_abacus_out.read_istate(os.path.join(root, "istate.info"))
                # continue if there is None among eks, pressure and bs
                parent = os.path.dirname(root)
                if None in [eks, pressure, bs]:
                    print(f"""WARNING: Present APNS job is broken: {parent}""")
                    continue
                atom_species, cellgen = read_apnsjob_desc(os.path.join(parent, "description.json"))
                system = os.path.basename(cellgen["config"])
                abacus_input = read_abacus_out.read_keyvals_frominput(os.path.join(parent, "INPUT"))
                ecutwfc = abacus_input["ecutwfc"]
                pps = [a["pp"] for a in atom_species]
                ppids = [convert_fpp_to_ppid(pp) for pp in pps]
                zvals = [ppparse.z_valence(os.path.join(parent, pp)) for pp in pps]
                s = "\n".join(ppids)
                print(f"""In folder {parent}
Structure tested: {system}
Number of atoms: {natom}
ecutwfc: {ecutwfc}
Final Kohn-Sham energy: {eks}
Pressure: {pressure}
Pseudopotentials are used:\n{s}
""")
                data = {"ecutwfc": ecutwfc, "eks": eks, "pressure": pressure, "istate": bs, "natom": natom, "z_valence": zvals}
                # band structure is not easy to print, therefore omitted
                idx = -1 if result.get(system, None) is None \
                    or result[system].get("ppcases", None) is None \
                        or result[system]["ppcases"].count(pps) == 0 \
                    else result[system]["ppcases"].index(pps)
                if idx == -1:
                    result.setdefault(system, {"ppcases": [], "pptests": []}).setdefault("ppcases", []).append(pps)
                    result[system]["pptests"].append([data])
                else:
                    result[system]["pptests"][idx].append(data)
                #result[(system, "|".join(pps), ecutwfc)] = (natom, zvals, eks, pressure, bs)
    return result

def repair_apnsjob(folder: str):
    """this function is written because one old abacustest version will drop description.json auto-generated by
    apns workflow. Now this has been fixed, so do not call this function in any cases.
    Call example:
    ```python
    repair_apnsjob("12310698")
    ```
    """
    import os
    import re
    import shutil
    for root, _, files in os.walk(folder):
        for file in files:
            if re.match(r"(running_)(\w+)(\.log)", file): # reach the inner most folder like OUT.ABACUS
                parent = os.path.dirname(root)
                f = os.path.basename(os.path.dirname(os.path.dirname(root)))
                shutil.copy2(os.path.join(f"../Yb_ecutwfc_test/{f}/description.json"), os.path.join(parent, "description.json"))

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

    collected = collect_apnsjob_data("12551436")
    system_and_stpcs = build_sptc_from_nested(collected)
    result = {}
    for s, stpcs in system_and_stpcs.items():
        for stpc in stpcs:
            pp, data = stpc()
            result[(s, pp)] = data
            ecut_conv = stpc.ecuts[stpc.iconv]
            pp = stpc.pp(as_list=True)
            assert len(pp) == 1, "The pseudopotential should be unique for each test case"
            update_ecutwfc(pp[0], ecut_conv)
    from apns.analysis.drivers.apns2_ecut_utils import plot_log, plot_stack
    flogs = plot_log(result)
    fstacks = plot_stack(result)