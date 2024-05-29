"""this file is for greping and rendering svg plots of ecutwfc convergence test for new version of APNS"""

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
    elif "pslibrary-pbe.0.3.1" in fpp:
        family, version = "PSlibrary", "0.3.1"
        if "kjpaw" in fpp:
            appendix = "PAW"
        elif "rrkjus" in fpp:
            appendix = "RRKJUS"
        elif "nc" in fpp:
            appendix = "NC"
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
    else:
        raise ValueError(f"Unrecognized pseudopotential file: {fpp}")
    return f"{family} v{version} ({appendix})".replace("v ", "")

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
                eks = read_abacus_out.read_etraj_fromlog(os.path.join(root, file))[-1]
                pressure = read_abacus_out.read_pressure_fromlog(os.path.join(root, file))
                bs = read_abacus_out.read_istate(os.path.join(root, "istate.info"))
                # continue if there is None among eks, pressure and bs
                parent = os.path.dirname(root)
                if None in [eks, pressure, bs]:
                    print(f"""WARNING: Present APNS job is broken: {parent}""")
                    continue
                atom_species, cellgen = read_apnsjob_desc(os.path.join(parent, "description.json"))
                system = os.path.basename(cellgen["config"])
                #system = "Yb"
                abacus_input = read_abacus_out.read_keyvals_frominput(os.path.join(parent, "INPUT"))
                ecutwfc = abacus_input["ecutwfc"]
                #pps = [fpp for fpp in os.listdir(parent) if fpp.lower().endswith(".upf")]
                #ppids = pps
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
                # band structure is not easy to print, therefore omitted
                pps = ", ".join(ppids)
                result[(system, pps, ecutwfc)] = (natom, zvals, eks, pressure, bs)
    return result

def merge_collected(collected: dict):
    """the collected result would have the key like:
    (system, (pp1, pp2, pp3), ecutwfc) -> (natom, eks, pressure, bs).
    The task of this function is to merge the collected data along axis of ecutwfc, so that the result would be:
    (system, (pp1, pp2, pp3)) -> {ecutwfc: (natom, eks, pressure, bs), ...}."""
    print("* * * Merge collected job results * * *".center(100))
    result = {}
    for key, value in collected.items():
        system, pps, ecutwfc = key
        result.setdefault((system, pps), {})[ecutwfc] = value
    print(f"""Merge summary\nSystem with pp: ecutwfc -> (natom, eks, pressure, bs)""")
    for k in result.keys():
        print(f"{k}: {len(result[k])} tests")
    return result

def sort_merged(merged: dict):
    """sort the merged result according to ecutwfc, return a sorted one and all properties are stored in lists, except `natom`.
    Note: lists in this step will be converted to np.ndarray data type!"""
    print("* * * Sort the merged data * * *".center(100))
    import numpy as np
    result = {}
    for key in merged.keys():
        ecuts = list(merged[key].keys())
        natoms = [merged[key][ecut][0] for ecut in ecuts]
        natoms_avg = np.mean(natoms)
        assert all([n == natoms_avg for n in natoms]), "The number of atoms should be consistent for all ecutwfc tests"
        natom = natoms_avg
        # zval is the second, should be uniformly the same for all ecutwfc tests
        zvals = [merged[key][ecut][1] for ecut in ecuts]
        assert all([\
            len(zvals[0])==len(zval) and all([zi == zj for zi, zj in zip(zvals[0], zval)]) \
                for zval in zvals])
        zvals = zvals[0]
        # the following part is hard coded temporarily
        eks = [merged[key][ecut][2] for ecut in ecuts]
        press = [merged[key][ecut][3] for ecut in ecuts]
        istate = [merged[key][ecut][4] for ecut in ecuts]
        # sort according to ecutwfc
        ecuts = [float(ecut) for ecut in ecuts]
        idx = np.argsort(ecuts)
        ecuts = np.array(ecuts)[idx].tolist()
        eks = np.array(eks)[idx].tolist()
        press = np.array(press)[idx].tolist()
        istate = np.array(istate)[idx].tolist()
        result[key] = dict(zip(["natom", "zvals", "ecuts", "energies", "pressures", "istates"], \
                               [natom, zvals, ecuts, eks, press, istate]))
    return result

def calc_sorted_conv(sorted_: dict, ethr: float = 1e-3, pthr: float = 0.1, bsthr: float = 1e-2):
    print("* * * Calculate converged ecutwfc * * *".center(100))
    import apns.analysis.postprocess.conv.kernel as amack
    import apns.analysis.postprocess.conv.ecutwfc_istate as istate
    result = {}
    for key, val in sorted_.items():
        system, pps = key
        natom = val["natom"]
        energies = [e/natom for e in val["energies"]]
        de = amack.default_calculator(energies, energies[-1])
        dp = amack.default_calculator(val["pressures"], val["pressures"][-1])
        dbs = istate.calculator(val["istates"], val["istates"][-1])
        
        iconv_e = [i for i, e in enumerate(de) if abs(e) <= ethr][0]
        iconv_p = [i for i, p in enumerate(dp) if abs(p) <= pthr][0]
        iconv_bs = [i for i, bs in enumerate(dbs) if abs(bs) <= bsthr][0]
        iconv = max(iconv_e, iconv_p, iconv_bs)
        print(f"""For system {system}
with pseudopotentials: {pps}
Energy converges to the threshold {ethr*1e3} meV/atom at ecutwfc = {val["ecuts"][iconv_e]} Ry
Pressure converges to the threshold {pthr} kbar at ecutwfc = {val["ecuts"][iconv_p]} Ry
Band structure converges to the threshold {bsthr*1e3} meV at ecutwfc = {val["ecuts"][iconv_bs]} Ry

Overall converged ecutwfc is {val["ecuts"][iconv]} Ry
""")
        result[(system, pps)] = {"de": de, "dp": dp, "dbs": dbs, "iconv": iconv, "ecutwfc": val["ecuts"],
                                 "zvals": val["zvals"]}
        
    print("* * * Calculation end * * *".center(100))
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
    import apns.analysis.external_frender.htmls as amaeh
    element = "Yb"
    html = amaeh.pseudopotentials(element=element, 
                                    xc_functional="PBE", 
                                    software="ABACUS",
                                    fconv=f"{element}.svg",
                                    fconvlog=f"{element}_logplot.svg")
    with open(f"{element}.md", "w") as f:
        f.write(html)
    exit()
    repair_apnsjob("12310698")
    #exit()
    collected = collect_apnsjob_data("12310698")
    merged = merge_collected(collected)
    calculated = calc_sorted_conv(sort_merged(merged))
    plot_log(calculated)
    plot_stack(calculated)