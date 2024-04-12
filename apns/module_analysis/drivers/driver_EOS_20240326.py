ACWF_DATAPATH = "/root/documents/acwf-verification-scripts-main/acwf_paper_plots/code-data/"
ACWF_REFJSON = "results-unaries-verification-PBE-v1-AE-average.json"

import json
def read_acwf_refdata(token: str, 
                      refdata_path: str = ACWF_DATAPATH + ACWF_REFJSON,
                      domain: str = "BM_fit_data"):
    """Read data from ACWF reference json file
    see Github repo:
    https://github.com/aiidateam/acwf-verification-scripts/tree/main
    
    and their paper:
    Bosoni, E., Beal, L., Bercx, M. et al. 
    How to verify the precision of density-functional-theory implementations via reproducible 
    and universal workflows. Nat Rev Phys 6, 45-58 (2024). 
    https://doi.org/10.1038/s42254-023-00655-3

    This project is referred in Standard Solid State Pseudopotential (SSSP) Library 2.0, which
    is maintained on website Materials Cloud:
    https://www.materialscloud.org/discover/sssp
    """
    with open(refdata_path, "r") as f:
        data = json.load(f)
    return data[domain][token]

import re
import apns.module_analysis.postprocess.eos as amape
def cal_delta_wrtacwf(element: str, bmfit: dict, vmin: float, vmax: float,
                        refdata_path: str = ACWF_DATAPATH + ACWF_REFJSON) -> float:
    """Search the best fit data from ACWF reference json file,
    because the structural data provided by Materials Project
    does not explicitly include the conventional name of crystal
    phase like BCC, FCC, Diamond, ..."""
    delta = 1e10
    syspatn = r"([A-Z][a-z]*)(-X/)([(SC)(FCC)(BCC)(Diamond)])"
    with open(refdata_path, "r") as f:
        data = json.load(f)
    data = data["BM_fit_data"]

    phase = ""
    for key in data.keys():
        _match = re.match(syspatn, key)
        if _match is not None:
            if _match.group(1) == element:
                print("Scan crystal phase from ACWF reference data:", key)
                result = data[key]
                value = amape.delta_value(bm_fit1=bmfit, bm_fit2=result, vmin=vmin, vmax=vmax)
                if value < delta:
                    delta = value
                    phase = key

    return delta, phase

import os
import apns.module_analysis.postprocess.read_abacus_out as amapr
def is_outdir(files: list):
    return "running_cell-relax.log" in files and\
           "STRU_ION_D" in files and\
           "INPUT" in files and\
           "STRU_READIN_ADJUST.cif" in files and\
           "istate.info" in files and\
           "kpoints" in files

def search(path: str, fromcal: str = "cell-relax"):
    
    result = {}
    for root, dirs, files in os.walk(path):
        if is_outdir(files):
            natom = amapr.read_natom_fromlog(f"{root}/running_{fromcal}.log")
            eks = amapr.read_etraj_fromlog(f"{root}/running_{fromcal}.log")[-1]
            v = amapr.read_volume_fromstru(f"{root}/STRU_ION_D", "A")
            _, sys, mpid, pnid, _ = amapr.read_testconfig_fromBohriumpath(root)
            result.setdefault((sys, mpid, pnid), []).append((v, eks, natom))
    
    # sort by volume
    for key in result.keys():
        result[key] = sorted(result[key], key=lambda x: x[0])
        result[key] = list(zip(*result[key]))
        result[key][0] = list(result[key][0]) # volume
        result[key][1] = list(result[key][1]) # energy
        result[key][2] = list(result[key][2]) # natom

    return result

import matplotlib.pyplot as plt
import numpy as np
import apns.module_analysis.postprocess.pseudopotential as amapp
import apns.module_analysis.external_frender.styles as amefs
def plot_eos(sys_mpid_pnid, vs, eks, bm_fit, ref_bm_fit, delta, **kwargs):
    """draw the EOS plot for the system
    
    Args:
        sys_mpid_pnid: tuple, (system, MPID and PNID)
        vs: list, the volume list
        eks: list, the energy list
        bm_fit: dict, the fit result
        ref_bm_fit: dict, the reference fit result
        delta: float, the delta value, REMEMBER TO DIVIDE BY NATOMS
    
    Returns:
        None
    """
    ###################
    # Styles setting  #
    ###################
    fontsize = kwargs.get("fontsize", 15)
    colors = kwargs.get("colors", amefs.styles_factory(property="color", ndim=2))
    marker = kwargs.get("marker", amefs.styles_factory(property="marker", ndim=1))[0]
    markersize = kwargs.get("markersize", 10)
    ###################
    # Data operations #
    ###################
    vs = np.array(vs)
    # shift the energy to 0
    eks = np.array(eks) - min(eks)
    # do inter/extrapolation on vs to let it evenly distributed in 200 points
    # v_eksmin = vs[np.argmin(eks)]
    # deltav_rhs = max(vs) - v_eksmin
    # deltav_lhs = v_eksmin - min(vs)
    # deltav = max(deltav_lhs, deltav_rhs)
    # vs_interp = np.linspace(min(vs) - deltav, max(vs) + deltav, 200)
    vs_interp = np.linspace(min(vs)*0.995, max(vs)*1.01, 200)

    # shift the energy to 0
    ref_interp = amape.birch_murnaghan(vs_interp, ref_bm_fit["E0"], ref_bm_fit["bulk_modulus_ev_ang3"], ref_bm_fit["bulk_deriv"], ref_bm_fit["min_volume"])
    ref_interp = ref_interp - min(ref_interp)
    # shift the energy to 0
    fit_interp = amape.birch_murnaghan(vs_interp, bm_fit["E0"], bm_fit["bulk_modulus_ev_ang3"], bm_fit["bulk_deriv"], bm_fit["min_volume"])
    fit_interp = fit_interp - min(fit_interp)

    pspotid, _ = amapp.testname_pspotid(sys_mpid_pnid[0], sys_mpid_pnid[2])

    fig, ax = plt.subplots(figsize=(10, 8))
    plt.plot(vs_interp, ref_interp, "-", label="AE (averaged over Wien2K and FLEUR)", color=colors[0])
    plt.plot(vs, eks, "o", label=f"DFT: {pspotid}", markersize=markersize, markeredgecolor=colors[1], markerfacecolor="none",
             markeredgewidth=1.5, linestyle="None")
    #plt.plot(vs_interp, fit_interp, "-", label=f"DFT: {pspotid}", color=colors[1])
    # turn on grid
    plt.grid()
    # set titles
    plt.xlabel("Volume ($\AA^3$)", fontsize=fontsize)
    plt.ylabel("Kohn-Sham energy (Shifted, $\Delta E_{KS}$) (eV)", fontsize=fontsize)
    plt.legend(fontsize=fontsize)
    # add annotation the delta value, with unit meV/atom
    plt.annotate(f"$\Delta$ = {delta*1e3:.2f} meV/atom", xy=(0.5, 0.5), xycoords="axes fraction", fontsize=fontsize)
    # set fontsize of xticks and yticks
    plt.xticks(fontsize=fontsize*0.9)
    plt.yticks(fontsize=fontsize*0.9)
    # set font as Arial
    plt.rcParams["font.family"] = "Arial"
    # set super title
    plt.suptitle(f"{sys_mpid_pnid[0]} EOS", fontsize=fontsize*1.2)
    # save the figure
    plt.savefig("eos_" + "_".join(sys_mpid_pnid) + ".png")
    plt.close()

import apns.module_analysis.postprocess.eos as amape
def calculate(sys_mpid_pnid_veks: dict):
    """calculate EOS V0, E0, B0, B0', delta refer to All Electron provided by ACWF,
    plot figure for (system, MPID, PNID) tuple-organized data"""

    result = {}
    for key in sys_mpid_pnid_veks.keys(): # loop over all (system, MPID, PNID) tuple
        print("Processing test:", key)
        # get source data
        ve = sys_mpid_pnid_veks[key][0]
        eks = sys_mpid_pnid_veks[key][1]
        natoms = sys_mpid_pnid_veks[key][2]
        assert len(set(natoms)) == 1, "The number of atoms in the cell is not consistent."
        natoms = natoms[0]
        # fit the data
        bm_fit = amape.birch_murnaghan_eos(ve, eks, as_dict=True)
        # calculate delta
        delta, phase = cal_delta_wrtacwf(key[0], bm_fit, min(ve), max(ve))
        delta = delta/natoms
        # read reference data
        bm_fitref = read_acwf_refdata(phase, domain="BM_fit_data")
        # write-back the data
        sys_mpid_pnid_veks[key] = bm_fit
        sys_mpid_pnid_veks[key]["delta"] = delta
        # save to result
        element, mpid, pnid = key
        if element not in result.keys():
            result[element] = {"mpid": mpid, "AEref": bm_fitref}
        result[element][pnid] = sys_mpid_pnid_veks[key]
        result[element][pnid]["volume"] = ve
        result[element][pnid]["energy"] = eks
        # plot the EOS
        # plot_eos(key, ve, eks, bm_fit, bm_fitref, delta)
        # will save to "eos_[element]_[mpid]_[pnid].png"

    return result

def plot(testresult: dict, ncols: int = 2, **kwargs):
    """plot all in one-shot, from the output of calculate function
    example:
    ```json
    {
        "Si": {
            "mpid": "149",
            "AEref": {"E0": 0, "bulk_deriv": ..., "bulk_modulus_ev_ang3": ..., "min_volume": ..., "residuals": ...},
            "dojo05": {"E0": ..., "bulk_deriv": ..., "bulk_modulus_ev_ang3": ..., "min_volume": ..., "residuals": ..., 
                       "delta": ...,
                       "volume": ..., "energy": ...},
            "pd04": {"E0": ..., "bulk_deriv": ..., "bulk_modulus_ev_ang3": ..., "min_volume": ..., "residuals": ..., 
                     "delta": ...,
                     "volume": ..., "energy": ...},
            ...
        },
        "Al": {...}
    }
    ```
    """
    ###################
    # Styles setting  #
    ###################
    fontsize = kwargs.get("fontsize", 18)
    colors = kwargs.get("colors", amefs.styles_factory(property="color", ndim=2))
    marker = kwargs.get("marker", amefs.styles_factory(property="marker", ndim=1))[0]
    markersize = kwargs.get("markersize", 10)
    
    subplot_height = kwargs.get("subplot_height", 10)
    subplot_width = kwargs.get("subplot_width", 10)

    # each element will create one figure
    for element in testresult.keys():
        print("Plotting element:", element)
        suptitle = f"{element} EOS (mpid: {testresult[element]['mpid']})"
        # calculate number of subplots and their layout
        npspots = len(testresult[element].keys()) - 2
        # however, if the number of subplots is less than ncols, we will set ncols to npspots
        ncols = min(npspots, ncols)
        # calculate nrows
        nrows = npspots // ncols + (npspots % ncols > 0)
        # calculate the figure size
        figsize = (subplot_width*ncols, subplot_height*nrows)
        # create figure
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, squeeze=False)

        pnids = [pnid for pnid in testresult[element].keys() if pnid != "mpid" and pnid != "AEref"]
        for i, pnid in enumerate(pnids):
            row = i // ncols
            col = i % ncols
            print("Plotting:", element, pnid, "at", row, col)
            ax = axes[row, col]
            bm_fit = testresult[element][pnid]
            bm_fitref = testresult[element]["AEref"]
            ve = bm_fit["volume"]
            eks = bm_fit["energy"]
            delta = bm_fit["delta"]
            ###################
            # Data operations #
            ###################
            vs = np.array(ve)
            # shift the energy to 0
            eks = np.array(eks) - min(eks)
            # do inter/extrapolation on vs to let it evenly distributed in 200 points
            vs_interp = np.linspace(min(vs)*0.995, max(vs)*1.01, 200)
            # shift the energy to 0
            ref_interp = amape.birch_murnaghan(vs_interp, bm_fitref["E0"], bm_fitref["bulk_modulus_ev_ang3"], bm_fitref["bulk_deriv"], bm_fitref["min_volume"])
            ref_interp = ref_interp - min(ref_interp)
            # shift the energy to 0
            fit_interp = amape.birch_murnaghan(vs_interp, bm_fit["E0"], bm_fit["bulk_modulus_ev_ang3"], bm_fit["bulk_deriv"], bm_fit["min_volume"])
            fit_interp = fit_interp - min(fit_interp)
            pspotid, _ = amapp.testname_pspotid(element, pnid)
            # plot
            ax.plot(vs_interp, ref_interp, "-", label="AE (averaged over Wien2K and FLEUR)", color=colors[0])
            ax.plot(vs, eks, "o", label=f"DFT: {pspotid}", markersize=markersize, markeredgecolor=colors[1], markerfacecolor="none",
                    markeredgewidth=1.5, linestyle="None")
            ax.grid()
            ax.set_xlabel("Volume ($\AA^3$)", fontsize=fontsize)
            ax.set_ylabel("Kohn-Sham energy (Shifted, $\Delta E_{KS}$) (eV)", fontsize=fontsize)
            ax.legend(fontsize=fontsize)
            ax.annotate(f"$\Delta$ = {delta*1e3:.2f} meV/atom", xy=(0.5, 0.5), xycoords="axes fraction", fontsize=fontsize)
            # set fontsize of xticks and yticks
            ax.tick_params(axis="both", labelsize=fontsize*0.9)
            # set font as Arial
            plt.rcParams["font.family"] = "Arial"
            # set super title
            fig.suptitle(suptitle, fontsize=fontsize*1.2)

        plt.savefig("eos_" + element + ".png")
        plt.close()



import argparse
def entry():
    parser = argparse.ArgumentParser(description='Calculate EOS from abacus output')
    # add -i
    parser.add_argument("-i", "--input", type=str, help="input json file")
    parser.add_argument("-f", "--fromcal", type=str, default="cell-relax", help="from calculation type")
    args = parser.parse_args()
    return args.input, args.fromcal


def run():

    path, fromcal = entry()
    print("Search in path:", path, "from calculation:", fromcal)
    ve_data = search(path, fromcal=fromcal)
    
    data = calculate(ve_data)
    plot(data)

import unittest
class TestEOS(unittest.TestCase):
    def test_is_outdir(self):
        files = ["running_relax.log", "STRU_ION_D", "INPUT", "STRU_SIMPLE.cif", "STRU_READIN_ADJUST.cif"]
        self.assertTrue(is_outdir(files))
        files = ["running_relax.log", "STRU_ION_D", "STRU_SIMPLE.cif", "STRU_READIN_ADJUST.cif"]
        self.assertFalse(is_outdir(files))
        files = ["running_relax.log", "STRU_ION_D", "STRU_SIMPLE.cif", "STRU_READIN_ADJUST.cif", "STRU_READIN_ADJUST.cif"]
        self.assertFalse(is_outdir(files))
        files = ["running_relax.log", "STRU_ION_D", "STRU_SIMPLE.cif", "STRU_READIN_ADJUST.cif", "STRU_READIN_ADJUST.cif", "STRU_READIN_ADJUST.cif"]
        self.assertFalse(is_outdir(files))
    def test_calculate(self):
        data = [[37.9098, -214.406018],
                [38.7163, -214.441355],
                [39.5229, -214.4649915],
                [40.3295, -214.4780729],
                [41.1361, -214.4816059],
                [41.9427, -214.4765565],
                [42.7493, -214.4637563]]
        smpv = {
            ("Si", "149", "dojo05"): ([x[0] for x in data], [x[1] for x in data])
        }
        result = calculate(smpv)
        self.assertAlmostEqual(
            result[("Si", "149", "dojo05")]["E0"], 
            -214.48165631975377, 
            delta=1e-4)
        self.assertAlmostEqual(
            result[("Si", "149", "dojo05")]["bulk_deriv"], 
            4.220117478624572, 
            delta=1e-4)
        self.assertAlmostEqual(
            result[("Si", "149", "dojo05")]["bulk_modulus_ev_ang3"], 
            0.5476332146301505, 
            delta=1e-4)
        # it is not values really fully reliable, just for testing
        # more reliable all-electron data can refer to https://doi.org/10.1038/s42254-023-00655-3
        """
        "Si-X/Diamond": {
            "E0": 0,
            "bulk_deriv": 4.31178461988603,
            "bulk_modulus_ev_ang3": 0.5524442002451444,
            "min_volume": 40.914946909495136,
            "residuals": 0
        },
        """
    
if __name__ == "__main__":
    # path = "../job-abacustest-v0.3.102-4edc4d"
    # result = run(path)
    # print(result)
    #unittest.main()
    run()