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

#########
# utils #
#########
def read_apnsjob_desc(fdesc: str):
    import json
    with open(fdesc, 'r') as f:
        desc = json.load(f)

    print("Read APNS job description from file: ", fdesc)
    atom_species_symbols = [as_["symbol"] for as_ in desc["AtomSpecies"]]
    s = ", ".join(atom_species_symbols)
    print(f"Atom species: {s}")
    pps = [as_["pp"] for as_ in desc["AtomSpecies"]]
    s = "\n".join(pps)
    print(f"Pseudopotentials bond with:\n{s}")
    # naos = [as_["nao"] for as_ in desc["AtomSpecies"]]
    # s = ", ".join(naos)
    # print(f"Number of atomic orbitals: {s}")
    cellgen = desc["CellGenerator"]
    print(f"""CellGenerator (where the cell generated from)
identifier: {cellgen["identifier"]}
source: {cellgen["config"]}
""")
    return desc["AtomSpecies"], desc["CellGenerator"]

def convert_fpp_to_ppid(fpp: str):
    import re
    import os
    if "hgh" in fpp:
        family, version = "Hartwigsen-Goedecker-Hutter", ""
        appendix = re.match(r"([A-Z][a-z]?\.pbe\-)(.*)(hgh\.UPF)", os.path.basename(fpp)).group(2)
        appendix = "" if appendix is None else appendix
    elif "NCPP-PD04-PBE" in fpp:
        family, version = "PD04", ""
        appendix = re.match(r"([A-Z][a-z]?)([\d\+\-\_\w]*)(\.PD04\.PBE\.UPF)", os.path.basename(fpp)).group(2)
        appendix = "" if appendix is None else appendix[1:] if appendix[0] in ["+", "-"] else appendix
    elif "GBRV_pbe_UPF_v1.5" in fpp:
        family, version, appendix = "GBRV", "1.5", ""
    elif "nc-sr-05_pbe_standard_upf" in fpp:
        family, version, appendix = "PseudoDojo", "0.5", "sr"
    elif "nc-fr-04_pbe_standard" in fpp:
        family, version, appendix = "PseudoDojo", "0.4", "fr"
    elif "pbe_s_sr" in fpp:
        family, version, appendix = "PseudoDojo", "0.3", "sr"
    elif "NCPP-PD03-PBE" in fpp:
        family, version, appendix = "PD03", "", ""
    elif "sg15_oncv_upf_2020-02-06" in fpp:
        family = "SG15"
        match_ = re.match(r"([A-Z][a-z]?(_ONCV_PBE)(_)?(FR)?(\-)(\d\.\d)(\.upf))", os.path.basename(fpp))
        version = match_.group(6)
        appendix = "fr" if match_.group(4) is not None else "sr"
    elif "nc-sr-04_pbe_standard_upf" in fpp:
        family, version, appendix = "PseudoDojo", "0.4", "sr"
    elif "nc-sr-04-3plus_pbe_standard_upf" in fpp:
        family, version, appendix = "PseudoDojo", "0.4", "sr"
    elif "pseudos_ac_she" in fpp:
        family, version = "PseudoDojo", "1.0"
        appendix = "fr" if fpp.endswith("_r.upf") else "sr"
    elif "gth" in fpp:
        family, version = "Goedecker-Teter-Hutter", ""
        appendix = os.path.basename(fpp).split("_")[-1].split(".")[0]
    elif "psl" in fpp:
        match_ = re.match(r"([A-Z][a-z]?)(\.)(rel-)?(pbe|pz)(-\w+)?(-)(rrkjus|kjpaw)(_psl\.)([\.\d]+)(\.UPF)", os.path.basename(fpp))
        family = "PSlibrary"
        version = match_.group(9)
        apps = []
        if match_.group(7): apps.append(match_.group(7).upper())
        if match_.group(3): apps.append("fr")
        if match_.group(5): apps.append(match_.group(5)[1:])
        appendix = ", ".join(apps)
    else:
        raise ValueError(f"Unrecognized pseudopotential file: {fpp}")
    return f"{family} v{version} ({appendix})".replace("v ", "").replace("()", "")

def collect_apnsjob_data(folder: str):
    print("* * * Collect ABACUS result * * *".center(100))
    import apns.analysis.postprocess.read_abacus_out as read_abacus_out
    import apns.pspot.parse as ppparse
    import os
    import re
    result = {}
    for root, _, files in os.walk(folder):
        for file in files:
            if re.match(r"(running_)(\w+)(\.log)", file): # reach the inner most folder like OUT.ABACUS
                natom = read_abacus_out.read_natom_fromlog(os.path.join(root, file))
                eks = read_abacus_out.read_efin_fromlog(os.path.join(root, file))
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

def plot_log(conv_result: dict):
    import matplotlib.pyplot as plt
    import apns.analysis.drivers.driver_EcutwfcConv_20240319 as outdated
    plt.rcParams["font.family"] = "Arial"

    # merge again that indexed like [system][pps]
    merged = {}
    for key, val in conv_result.items():
        system, pps = key
        merged.setdefault(system, {})[pps] = val
    figures = {s: f"{s}_logplog.svg" for s in merged.keys()}
    for s, r in merged.items(): # s stands for system and r stands for result
        # result would be dict indexed by different pps
        pps = list(r.keys())
        xs = [[r[pp]["ecutwfc"] for pp in pps], [r[pp]["ecutwfc"] for pp in pps], [r[pp]["ecutwfc"] for pp in pps]]
        ys = [[r[pp]["de"] for pp in pps], [r[pp]["dp"] for pp in pps], [r[pp]["dbs"] for pp in pps]]
        logplot_style = {"highlight_ys": [1e-3, 0.1, 1e-2], "nrows": 1, 
                         "xtitle": "Planewave kinetic energy cutoff (ecutwfc, in Ry)", 
                         "ytitle": ["Absolute Kohn-Sham energy difference per atom (eV)", 
                                    "Absolute pressure difference (kbar)",
                                    "Band structure difference (eV)"], 
                         "ysymbols": ["$|\Delta E_{KS}|$", "$|\Delta P|$", "$|\eta_{all, 00}|$"],
                         "suptitle": s, 
                         "supcomment": "NOTE: Absence of data points result from SCF convergence failure or walltime limit.",
                         "labels": pps, "fontsize": 19}
        fig, ax = outdated.discrete_logplots(xs, ys, **logplot_style)
        plt.savefig(figures[s])
        plt.close()

    return figures

def plot_stack(conv_result: dict):
    import matplotlib.pyplot as plt
    import apns.analysis.drivers.driver_EcutwfcConv_20240319 as outdated
    plt.rcParams["font.family"] = "Arial"

    # merge again that indexed like [system][pps]
    merged = {}
    for key, val in conv_result.items():
        system, pps = key
        merged.setdefault(system, {})[pps] = val
    figures = {s: f"{s}.svg" for s in merged.keys()}
    for s, r in merged.items(): # s stands for system and r stands for result
        pps = list(r.keys())
        xs = [[r[pp]["ecutwfc"]]*3 for pp in pps]
        ys = [[r[pp]["de"], r[pp]["dp"], r[pp]["dbs"]] for pp in pps]
        lineplot_style = {"highlight_xs": [(pp, r[pp]["ecutwfc"][r[pp]["iconv"]]) for pp in pps], "ncols": 1, 
                          "subtitles": pps, 
                          "z_vals": [int(float(r[pp]["zvals"][0])) for pp in pps], 
                          "grid": True,
                          "xtitle": "Planewave kinetic energy cutoff (ecutwfc, in Ry)", 
                          "ytitle": ["Kohn-Sham energy difference per atom (eV)", 
                                     "Pressure difference (kbar)",
                                     "Band structure difference (eV)"],
                          "suptitle": s, 
                          "supcomment": "NOTE: The red circle indicates the converged ecutwfc wrt. ecutwfc$_{max}$\
 with precision threshold (1.0 meV/atom, 0.1 kbar, 10 meV) respectively.\n \
Absence of data points result from SCF convergence failure or walltime limit.",
                          "fontsize": 13, "alpha": 0.8}
        shift_style = {"shifts": [5, 500, 10], "ld": "pseudopotential",
                       "ysymbols": ["$\Delta E_{KS}$", "$\Delta P$", "$\eta_{all, 00}$"],
                      }
        fig, ax = outdated.shift_lineplots(xs=xs, ys=ys, **shift_style, **lineplot_style)
        plt.savefig(figures[s])
        plt.close()

    return figures

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

    collected = collect_apnsjob_data("12508504")
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
    flogs = plot_log(result)
    fstacks = plot_stack(result)