'''this is a collection of utility functions for ecutwfc convergence test
'''
# built-in modules
import json
import os

# third-party modules
import numpy as np
import matplotlib.pyplot as plt

# local modules
from apns.analysis.apns1_ecut_abacus import discrete_logplots, shift_lineplots
from apns.analysis.postprocess.conv.kernel import cal_diff
from apns.analysis.postprocess.conv.ecutwfc_istate import cal_band_diff
from apns.analysis.apns2_utils import convert_fpp_to_ppid

def build_sptc_from_nested(cases: dict):
    """build instances of KineticEnergyCutoffTestCase from a nested dict. Example input:
    ```python
    cases = {
        "mp-1234.cif": {
            "ppcases": [["Ag.pbe-spn-rrkjus_psl.0.2.3-tot.UPF", 
                         "O.pbe-n-rrkjus_psl.0.2.3-tot.UPF"],
                        ["Ag_ONCV_PBE-1.0_fr.upf", 
                         "O_ONCV_PBE-1.0_fr.upf"]],
            "pptests": [
                [{"ecutwfc": 30, 
                  "eks": -1.0, "pressure": 0.0, "istate": 0.0, "natom": 1, "z_valence": [11, 6]},
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
    This function will return a dict whose keys are different systems, and values are lists of KineticEnergyCutoffTestCase instances.
    Each instance represents a pseudopotential test case for a specific system."""
    result = {}
    for system, data in cases.items():
        ppcases = data["ppcases"]
        for ipt, pptests in enumerate(data["pptests"]):
            pps = ppcases[ipt]
            sptc = KineticEnergyCutoffTestCase(system, pps, pptests)
            result.setdefault(system, []).append(sptc)
    return result

class KineticEnergyCutoffTestCase:
    """Definition: A single case in ecutwfc convergence test is:
    1. one system
    2. many ecutwfc
    3. (therefore) many energies, istates and pressures
    4. one (combination of) pseudopotential(s)"""
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
        # sort and remove duplicates
        idx = np.unique(self.ecuts, return_index=True)[1]
        self.ecuts = np.array(self.ecuts)[idx].tolist()
        self.energies = np.array(self.energies)[idx].tolist()
        self.pressures = np.array(self.pressures)[idx].tolist()
        self.istates = np.array(self.istates)[idx].tolist()

    def calc_conv(self, ethr: float = 1e-3, pthr: float = 0.1, bsthr: float = 1e-2):
        """calculate the converged ecutwfc, return the index of the converged ecutwfc"""

        energies = [e/self.natom for e in self.energies]
        de = cal_diff(energies, energies[-1])
        dp = cal_diff(self.pressures, self.pressures[-1])
        dbs = cal_band_diff(self.istates, self.istates[-1])
        assert len(de) == len(dp) == len(dbs) == len(self.ecuts), \
            "The length of the calculated differences should be the same as the number of ecutwfc tests"
        
        def is_converged(diff, thr):
            return abs(diff) <= thr
        
        i = len(self.ecuts) - 1
        while i > 0:
            if is_converged(de[i], ethr) and \
               is_converged(dp[i], pthr) and \
               is_converged(dbs[i], bsthr):
                i -= 1 # the safest search is to start from the last one
            else:
                break
        if i + 1 == len(self.ecuts):
            print(f"WARNING: The convergence test for {self.system} with pseudopotential {self.pp()} is not converged.")

        self.iconv = i + 1
        return {"de": de, 
                "dp": dp, 
                "dbs": dbs, 
                "iconv": i + 1,
                "ecutwfc": self.ecuts, 
                "zvals": self.zvals}
    
    def pp(self, as_list: bool = False):
        """return the pseudopotential string"""
        
        return self.pps if as_list else "|".join([": ".join(convert_fpp_to_ppid(pp)) for pp in self.pps])
    
    def zval(self, as_list: bool = False):
        """return the valence of atoms"""
        return self.zvals if as_list else "|".join(self.zvals)

    def __call__(self):
        self.sort()
        return self.pp(), self.calc_conv()

def update_ecutwfc(pp: str, ecutwfc: float, cache_dir: str = "./apns_cache/ecutwfc.json"):

    if os.path.exists(cache_dir):
        with open(cache_dir, "r") as f:
            cache = json.load(f)
    else:
        cache = {}
    cache[pp] = ecutwfc
    with open(cache_dir, "w") as f:
        json.dump(cache, f)

def plot_log(conv_result: dict, fmt = 'svg'):

    plt.rcParams["font.family"] = "Arial"

    if fmt not in ['svg', 'pdf', 'png']:
        raise ValueError("fmt should be one of 'svg', 'pdf', 'png'")
    # merge again that indexed like [system][pps]
    merged = {}
    
    for result in conv_result:
        system = result["name"]
        pps = result["pp"]
        merged.setdefault(system, {})[pps] = result
    
    figures = {s: f"{s}-logscale.{fmt}" for s in merged.keys()}
    for s, r in merged.items(): # s stands for system and r stands for result
        # result would be dict indexed by different pps
        pps = list(r.keys())
        pps.sort(key=lambda x: x.split(': ')[-1])
        xs = [[r[pp]["ecutwfc"] for pp in pps]]*3
        ys = [[r[pp]["de"] for pp in pps], 
              [r[pp]["dp"] for pp in pps], 
              [r[pp]["dbs"] for pp in pps]]
        logplot_style = {"highlight_ys": [1e-3, 0.1, 1e-2], "nrows": 1, 
                         "xtitle": "Planewave kinetic energy cutoff (ecutwfc, in Ry)", 
                         "ytitle": ["Absolute Kohn-Sham energy difference per atom (eV)", 
                                    "Absolute pressure difference (kbar)",
                                    "Band structure difference (eV)"], 
                         "ysymbols": ["$|\Delta E_{KS}|$", "$|\Delta P|$", "$|\eta_{all, 00}|$"],
                         "suptitle": s, 
                         "supcomment": "NOTE: Absence of data points result from SCF convergence failure or walltime limit.",
                         "labels": pps, "fontsize": 22.5}
        fig, ax = discrete_logplots(xs, ys, **logplot_style)
        plt.savefig(figures[s])
        plt.close()

    return figures

def plot_stack(conv_result: dict, figfmt = 'svg'):

    plt.rcParams["font.family"] = "Arial"

    if figfmt not in ['svg', 'pdf', 'png']:
        raise ValueError("figfmt should be one of 'svg', 'pdf', 'png'")
    
    # merge again that indexed like [system][pps]
    merged = {}

    for result in conv_result:
        system = result["name"]
        pps = result["pp"]
        merged.setdefault(system, {})[pps] = result
    figures = {s: f"{s}.{figfmt}" for s in merged.keys()}
    figure_style = {"figsize": (20, 10)}
    for s, r in merged.items(): # s stands for system and r stands for result
        pps = list(r.keys())
        pps.sort(key=lambda x: x.split(': ')[-1])
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
        fig, ax = shift_lineplots(xs=xs, ys=ys, **shift_style, **lineplot_style, **figure_style)
        plt.tight_layout()
        plt.savefig(figures[s])
        plt.close()

    return figures
