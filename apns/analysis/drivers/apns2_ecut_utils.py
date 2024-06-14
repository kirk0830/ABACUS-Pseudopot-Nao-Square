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