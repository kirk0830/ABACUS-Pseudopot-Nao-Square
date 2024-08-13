"""
Driver Name: EcutwfcConv_20240319

Driver description:
This is driver for APNS pages pseudopotential convergence test
Two kinds of data are post-processed:
1. ecutwfc - energy per atom
2. ecutwfc - pressure of whole cell
"""

import argparse
def entry():
    parser = argparse.ArgumentParser(description="APNS pseudopotential convergence test")
    # add -i
    parser.add_argument("-i", "--input", type=str, help="input json file")
    args = parser.parse_args()
    return args.input

import json
import apns.pspot.parse as ampp
def z_valence(element: str, pspotid: str):
    """pspotid here should be the "key" returned by function testname_pspotid"""
    with open("download/pseudopotentials/pseudo_db.json", "r") as f:
        pseudo_db = json.load(f)
    return ampp.z_valence(pseudo_db[element][pspotid])

def categorize_byelement(labels, data):
    """the data will arrange like:
    ```python
    (
        [('Ag', '8566', 'dojo05', 90.0), ('Ag', '8566', 'pd03', 90.0), 
        ('Ag', '8566', 'sg1512', 60.0), ('Ag', '8566', 'sg1510', 60.0), 
        ...
        ],
        [
            [
                [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 200.0, 300.0], 
                [111.00118765492425, 33.39993995924942, 9.262190618749628, 2.1097572079243037, 0.3720193069993911, 0.041972779749812617, 
                0.002425373099868011, 0.0005534015244847978, 0.00044674614946416114, 0.00019317457463330356, 4.175727463007206e-05, 0.0]
            ], 
            [
                [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 200.0, 300.0], 
                [111.14227980914984, 33.32067861827545, 9.212816671300061, 2.113927723249617, 0.3799892410997927, 0.04445111024961079, 
                0.0027109032753287465, 0.000579334774556628, 0.00048263027474604314, 0.00022845432522444753, 4.262575021130033e-05, 0.0]
            ], 
            ...
        ],
        [('Ag', '8566', 'dojo05', 90.0), ('Ag', '8566', 'pd03', 90.0), 
        ('Ag', '8566', 'sg1512', 60.0), ('Ag', '8566', 'sg1510', 60.0), 
        ...
        ],
        [
            [
                [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 200.0, 300.0], 
                [111.00118765492425, 33.39993995924942, 9.262190618749628, 2.1097572079243037, 0.3720193069993911, 0.041972779749812617, 
                0.002425373099868011, 0.0005534015244847978, 0.00044674614946416114, 0.00019317457463330356, 4.175727463007206e-05, 0.0]
            ], 
            [
                [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 200.0, 300.0], 
                [111.14227980914984, 33.32067861827545, 9.212816671300061, 2.113927723249617, 0.3799892410997927, 0.04445111024961079, 
                0.0027109032753287465, 0.000579334774556628, 0.00048263027474604314, 0.00022845432522444753, 4.262575021130033e-05, 0.0]
            ], 
            ...
        ]
    )
    ```
    each time import two elements, the first is labels, the second is data.

    This function will categorize the data by element, 
    pspotids first indiced by element, second by pspotid,
    the same for decomposed_x and decomposed_y

    Return:
    elements: list of elements
    pspotids: list of list of pseudopotential ids
    decomposed_x: list of list of x data
    decomposed_y: list of list of y data
    """
    assert len(labels) == len(data)
    for label in labels:
        assert len(label) == 4 # element, mpid, pspotid, ecutwfc
        element, mpid, pspotid, ecutwfc = label
        assert isinstance(element, str)
        assert isinstance(mpid, str)
        assert isinstance(pspotid, str)
        assert isinstance(ecutwfc, float)
    for piece in data:
        assert len(piece) == 2 # ecutwfc, data
        assert isinstance(piece[0], list)
        assert isinstance(piece[1], list)
        assert len(piece[0]) == len(piece[1])

    elements = [] # is designed to only collect unique elements
    pspotids = [] # for each element, it would be a list of pspotids
    decomposed_x = [] # corresponds to pspotids, for each element, there would be several tests, and each test has a series of ecutwfc
    decomposed_y = [] # corresponds to pspotids, for each element, there would be several tests, and each test has a series of data
    for itst in range(len(labels)): # for each test, loop. one test is defined as fixed pseudopotential, series of SCF at different ecutwfc
        # there always are duplicates for element, because four keys define the uniqueness
        element, mpid, pspotid, ecutwfc = labels[itst] # get element, mpid, pspotid, and converged ecutwfc information
        if element not in elements: # if element is not in elements, append it
            elements.append(element)
            pspotids.append([]) # store a new element, thus append a new list for it
            decomposed_x.append([]) # store a new element, thus append a new list for it
            decomposed_y.append([]) # store a new element, thus append a new list for it
        idx = elements.index(element)
        pspotids[idx].append(pspotid)
        decomposed_x[idx].append(data[itst][0])
        decomposed_y[idx].append(data[itst][1])

    return elements, pspotids, decomposed_x, decomposed_y

import apns.analysis.postprocess.conv.ecutwfc_eks as amapeks
import apns.analysis.postprocess.conv.ecutwfc_press as amapprs
import apns.analysis.postprocess.conv.ecutwfc_istate as amapbs
def search(path: str, ethr: float = 1e-3, pthr: float = 0.1, bsthr: float = 1e-2):
    """this function will return three indexes data, they are, eks, prs and bs.
    For each tuple, the two variables have sequence correspondence, but not
    guaranteed across different tuples."""
    return amapeks.run(path, ethr),\
           amapprs.run(path, pthr),\
           amapbs.run(path, bsthr)

def summarize_conv(convs: list, result: dict):
    for conv in convs:
        for element, mpid, pnid, ecutwfc in conv:
            result.setdefault(element, {}).setdefault(pnid, ecutwfc)
            result[element][pnid] = max(result[element][pnid], ecutwfc)
    return result

import apns.outdated.identifier as amwi
def ecutwfc_convergence(convs: list):
    import os
    # summarize the convergence results
    fresult = amwi.TEMPORARY_FOLDER + "/apns_ecutwfc_db.json"
    if os.path.exists(fresult):
        with open(fresult, "r") as f:
            conv_result = json.load(f)
    else:
        conv_result = {}
    conv_result = summarize_conv(convs, conv_result)
    with open(fresult, "w") as f:
        json.dump(conv_result, f, indent=4)

def run():
    from apns.analysis.apns1_utils import ppfullname, apns1zval
    import json
    #import apns.analysis.postprocess.pseudopotential as amppspot
    # always make independent on path for run() function, instead, use entry()
    path = entry()
    # postprocess
    (conv_eks, data_eks), (conv_prs, data_prs), (conv_bs, data_bs) = search(path)
    # convergence
    # ecutwfc_convergence([conv_eks, conv_prs, conv_bs])
    # categorize by element
    # elements, pspotids, eks_x, eks_y = categorize_byelement(conv_eks, data_eks)
    # _, _, prs_x, prs_y = categorize_byelement(conv_prs, data_prs)
    # _, _, bs_x, bs_y = categorize_byelement(conv_bs, data_bs)
    # print(elements, file=open("elements.txt", "w")) # a single layer list of elements symbol
    # print(pspotids, file=open("pspotids.txt", "w")) # a two-layer list, the first is element and second is pspotid for each element
    # print(eks_x, file=open("eks_x.txt", "w")) # a three-layer list, the first is element, the second is pspotid, and the third is ecutwfc value
    # print(eks_y, file=open("eks_y.txt", "w")) # the same as eks_x, but the third layer is the data
    # print(prs_x, file=open("prs_x.txt", "w")) # the same as eks_x, but the third layer is the data
    # print(prs_y, file=open("prs_y.txt", "w")) # the same as eks_x, but the third layer is the data
    # print(bs_x, file=open("bs_x.txt", "w")) # the same as eks_x, but the third layer is the data
    # print(bs_y, file=open("bs_y.txt", "w")) # the same as eks_x, but the third layer is the data
    
    # print(conv_eks, file=open("conv_eks.txt", "w")) # a disordered list, contains data (elemm, mpid, pnid, ecut)
    # print(conv_prs, file=open("conv_prs.txt", "w")) # the same as conv_eks
    # print(conv_bs, file=open("conv_bs.txt", "w")) # the same as conv_eks
    # print(data_eks, file=open("data_eks.txt", "w")) # the list indexed like [elem][pspotid][i] -> list of [x, y]
    # print(data_prs, file=open("data_prs.txt", "w")) # the same as data_eks
    # print(data_bs, file=open("data_bs.txt", "w")) # the same as data_eks
    out = []
    ntests = len(conv_eks)
    assert ntests == len(conv_prs) and ntests == len(conv_bs), f"mismatch: {ntests} vs. {len(conv_prs)} vs. {len(conv_bs)}"
    # sort conv_eks and record the index mapping
    temp = sorted(conv_eks, key=lambda x: (x[0], x[1], x[2]))
    index_map_ = [conv_eks.index(t) for t in temp]
    conv_eks = temp
    data_eks = [data_eks[i] for i in index_map_]
    # sort conv_prs and record the index mapping
    temp = sorted(conv_prs, key=lambda x: (x[0], x[1], x[2]))
    index_map_ = [conv_prs.index(t) for t in temp]
    conv_prs = temp
    data_prs = [data_prs[i] for i in index_map_]
    # sort conv_bs and record the index mapping
    temp = sorted(conv_bs, key=lambda x: (x[0], x[1], x[2]))
    index_map_ = [conv_bs.index(t) for t in temp]
    conv_bs = temp
    data_bs = [data_bs[i] for i in index_map_]

    for i in range(ntests): # loop over all tests
        elem, mpid, pnid, ecut_eks = conv_eks[i]
        _e, _m, _p, ecut_prs = conv_prs[i]
        assert elem == _e and mpid == _m and pnid == _p, f"mismatch: {elem} vs. {_e}, {mpid} vs. {_m}, {pnid} vs. {_p}"
        _e, _m, _p, ecut_bs = conv_bs[i]
        assert elem == _e and mpid == _m and pnid == _p, f"mismatch: {elem} vs. {_e}, {mpid} vs. {_m}, {pnid} vs. {_p}"
        fcif = "mp-" + mpid + ".cif" if not mpid.startswith("mp-") and not mpid.endswith(".cif") else mpid + ".cif"
        pp = ppfullname(pnid)
        zvals = apns1zval(elem, pnid)
        ecuts = data_eks[i][0]
        eks = [float(e) for e in data_eks[i][1]]
        prs = [float(p) for p in data_prs[i][1]]
        bs = [float(b) for b in data_bs[i][1]]
        iconv = max(ecuts.index(ecut_eks), ecuts.index(ecut_prs), ecuts.index(ecut_bs))
        temp = {
            "name": elem, "fcif": fcif, "pp": pp, "zvals": [float(zvals)], "ecutwfc": ecuts,
            "de": eks, "dp": prs, "dbs": bs, "iconv": iconv
        }
        out.append(temp)
        print(temp)
    with open("out.json", "w") as f:
        json.dump(out, f, indent=4)
    exit()
    
    conv_results = {}
    for iconv in range(len(conv_eks)): # loop over all system-mpid-pnid tests...
        element = conv_eks[iconv][0]
        pnid = conv_eks[iconv][2]
        if not element in conv_results.keys():
            conv_results[element] = []
        conv_results[element].append((pnid, 
                                      max(conv_eks[iconv][3], 
                                          conv_prs[iconv][3],
                                          conv_bs[iconv][3])))
        
    print(conv_results)
    exit()
    # set Arial as default font
    plt.rcParams["font.family"] = "Arial"

    # remember: eks_*/prs_*/... is indiced by [element][pspotid][ecutwfc] to get a float value
    for i in range(len(elements)):
        element = elements[i]
        _pspotnames, _pspotids = zip(*[amppspot.testname_pspotid(element, id) for id in pspotids[i]])
        # DESIGN NOTE OF LOG_PLOTS AND STACK_LINEPLOTS
        # -------------------------------------------
        # LOG_PLOTS:
        # in some case, Bohrium crashes when calculating the pressure, this will lead to absent of a pair of (ecutwfc, pressure).
        # therefore, the length of eks_x and prs_x may not be the same, (the same for eks_y and prs_y).
        # In the following discrete_logplots() function, xs and ys support the input arranged in a way such that xs and ys are indiced
        # by [pspotid][ecutwfc], thus the eks_x[i][j] would be of meaning of ecutwfc list for j-th pseudopotential for element i.

        # what want to draw is, two subplots, one for energy convergence, one for pressure convergence.
        # For each subplot, there are several lines, each line corresponds to one pseudopotential. All x-y pairs are not needed
        # to be the same length, but for each pair, x and y should have identical length.
        # Therefore the most natural way to organize the data would be, first indiced by suplot, second by line and last by data.
        # the xs and ys should be instead indiced by [subplots][lines][data]
        
        # eks_x[i]/prs_x[i]/eks_y[i]/prs_y[i], indiced by [pspotid][ecutwfc]
        # or say shape = (npspotid, necutwfc)
        # let them be (nproperties, npspotid, necutwfc)
        xs = [eks_x[i], prs_x[i], bs_x[i]]
        ys = [eks_y[i], prs_y[i], bs_y[i]]

        logplot_style = {"highlight_ys": [1e-3, 0.1, 1e-2], "nrows": 1, 
                         "xtitle": "Planewave kinetic energy cutoff (ecutwfc, in Ry)", 
                         "ytitle": ["Absolute Kohn-Sham energy difference per atom (eV)", 
                                    "Absolute pressure difference (kbar)",
                                    "Band structure difference (eV)"], 
                         "ysymbols": ["$|\Delta E_{KS}|$", "$|\Delta P|$", "$|\eta_{all, 00}|$"],
                         "suptitle": element, 
                         "supcomment": "NOTE: Absence of data points result from SCF convergence failure or walltime limit.",
                         "labels": _pspotnames, "fontsize": 19}
        fig, ax = discrete_logplots(xs=xs, ys=ys, **logplot_style)
        plt.savefig(f"{element}_logplot.svg")
        plt.close()
        # STACK_LINEPLOTS
        # stack lineplots are organized in a way that, for each element, first indiced by pseudopotential,
        # then by property, and last by data. Therefore, the xs and ys should be indiced by [pspotid][property][data]
        # the shape of xs and ys should be (npspotid, nproperties, ndata)

        # eks_x[i]/prs_x[i]/eks_y[i]/prs_y[i], indiced by [pspotid][ecutwfc]
        # or say shape = (npspotid, necutwfc)
        # let them be (npspotid, nproperties, necutwfc)
        xs = [[eks_x[i][j], prs_x[i][j], bs_x[i][j]] for j in range(len(_pspotnames))]
        ys = [[eks_y[i][j], prs_y[i][j], bs_y[i][j]] for j in range(len(_pspotnames))]

        lineplot_style = {"highlight_xs": conv_results[element], "ncols": 1, 
                          "subtitles": _pspotnames, 
                          "z_vals": [z_valence(element, _pspotid) for _pspotid in _pspotids], 
                          "grid": True,
                          "xtitle": "Planewave kinetic energy cutoff (ecutwfc, in Ry)", 
                          "ytitle": ["Kohn-Sham energy difference per atom (eV)", 
                                     "Pressure difference (kbar)",
                                     "Band structure difference (eV)"],
                          "suptitle": element, 
                          "supcomment": "NOTE: The red circle indicates the converged ecutwfc wrt. ecutwfc$_{max}$\
 with precision threshold (1.0 meV/atom, 0.1 kbar, 10 meV) respectively.\n \
Absence of data points result from SCF convergence failure or walltime limit.",
                          "fontsize": 13, "alpha": 0.8}
        # fig, ax = stack_lineplots(xs=xs, ys=ys, **lineplot_style)
        # plt.savefig(f"{element}_stack.svg")
        # plt.close()

        shift_style = {"shifts": [5, 500, 10], "ld": "pseudopotential",
                       "ysymbols": ["$\Delta E_{KS}$", "$\Delta P$", "$\eta_{all, 00}$"],
                      }
        fig, ax = shift_lineplots(xs=xs, ys=ys, **shift_style, **lineplot_style)
        plt.savefig(f"{element}.svg")
        plt.close()

        import apns.analysis.external_frender.htmls as amaeh
        html = amaeh.pseudopotentials(element=element, 
                                      xc_functional="PBE", 
                                      software="ABACUS",
                                      fconv=f"{element}.svg",
                                      fconvlog=f"{element}_logplot.svg")
        with open(f"{element}.md", "w") as f:
            f.write(html)

import matplotlib.pyplot as plt
from matplotlib.legend import Legend
import numpy as np
import apns.analysis.external_frender.styles as amefs
def stack_lineplots(xs: list, ys: list, **kwargs):
    nsubplts = len(xs)
    nlines = len(ys[0])
    for i in range(1, nsubplts):
        assert len(xs[i]) == nlines
        assert len(ys[i]) == nlines

    highlight_xs = kwargs.get("highlight_xs", None)
    ncols = kwargs.get("ncols", 1)
    subtitles = kwargs.get("subtitles", None)
    labels = kwargs.get("labels", None)
    colors = kwargs.get("colors", None)
    markers = kwargs.get("markers", None)
    markersizes = kwargs.get("markersizes", None)
    linestyles = kwargs.get("linestyles", None)
    linewidths = kwargs.get("linewidths", None)
    alpha = kwargs.get("alpha", 1.0)
    grid = kwargs.get("grid", True)
    fontsize = kwargs.get("fontsize", 12)
    z_vals = kwargs.get("z_vals", None)
    subplotsize = kwargs.get("subplotsize", (20, 1))

    npspots = nsubplts
    subtitles = ["subplot " + str(i) for i in range(npspots)] if subtitles is None else subtitles
    assert len(subtitles) == npspots
    nprptys = nlines
    labels = ["data " + str(i) for i in range(nprptys)] if labels is None else labels
    assert len(labels) == nprptys

    colors = amefs.styles_factory(property="color", ndim=nprptys) if colors is None else colors
    assert len(colors) == nprptys
    markers = amefs.styles_factory(property="marker", ndim=nprptys) if markers is None else markers
    assert len(markers) == nprptys
    markersizes = amefs.styles_factory(property="markersize", ndim=nprptys) if markersizes is None else markersizes
    assert len(markersizes) == nprptys
    linestyles = amefs.styles_factory(property="linestyle", ndim=nprptys) if linestyles is None else linestyles
    assert len(linestyles) == nprptys
    linewidths = amefs.styles_factory(property="linewidth", ndim=nprptys) if linewidths is None else linewidths
    assert len(linewidths) == nprptys
    alpha = amefs.styles_factory(property="alpha", val=alpha, ndim=nprptys)
    assert len(alpha) == nprptys

    assert z_vals is None or len(z_vals) == npspots
    # create figure and axes
    nrows = nsubplts // ncols + (nsubplts % ncols > 0)
    fig, ax = plt.subplots(nrows, ncols, figsize=(subplotsize[0] * ncols + (nprptys - 1) * 0.02, # reserve space for multi-y axis
                                                  subplotsize[1] * nrows), 
                           squeeze=False)
    
    xtitle_styles = {"ha": "center", "va": "center", "transform": fig.transFigure, "fontsize": fontsize * 1.2}
    ytitle_styles = {"ha": "center", "va": "center", "rotation": "vertical", "transform": fig.transFigure, "fontsize": fontsize * 1.2}

    # for there may be data failed to converge or due to failure of Bohrium platform, the length of
    # xs between different pseudopotentials may not be the same, so we need to find the minimum and maximum
    # of ecutwfc for each pseudopotential

    ecutwfc_min = min([min([min(xs[i][j]) for j in range(nprptys)]) for i in range(npspots)])
    ecutwfc_max = max([max([max(xs[i][j]) for j in range(nprptys)]) for i in range(npspots)])
    ecutwfcs = [20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200]
    ecutwfcs.append(ecutwfc_max) if ecutwfc_max > ecutwfcs[-1] else None
    xticklabels = []
    for ecutwfc in np.arange(ecutwfc_min, ecutwfc_max + 10, 10):
        xticklabels.append(str(ecutwfc)) if ecutwfc in ecutwfcs else xticklabels.append("")
    xticks = np.linspace(ecutwfc_min, ecutwfc_max, int((ecutwfc_max - ecutwfc_min)/10) + 1)

    twinxs = []
    # plot
    for i in range(npspots): # loop over all pseudopotentials (although I am not willing to make it too specific)
                              # , but it is a rather easy way to explain the code
        row, col = i // ncols, i % ncols
        twinxs.append([])
        for j in range(nprptys):
            style = {"color": colors[j], "marker": markers[j], "markersize": markersizes[j],
                     "linestyle": linestyles[j], "linewidth": linewidths[j], "alpha": alpha[j]}
            twinxs[i].append(ax[row, col] if j == 0 else ax[row, col].twinx())
            twinxs[i][j].plot(xs[i][j], ys[i][j], label=labels[j], **style)
            # add yticks for each line
            ylim = np.max(np.abs(ys[i][j]))
            yticks = [-ylim*0.75, 0, ylim*0.75]
            twinxs[i][j].set_yticks(yticks)
            # set color for yaxis, yticks, yticklabels the same as the line, yticklabels .2f, fontsize = fontsize
            # thickness of yaxis = 2, yticklabels are left aligned
            twinxs[i][j].yaxis.label.set_color(colors[j])
            twinxs[i][j].tick_params(axis="y", colors=colors[j])
            twinxs[i][j].set_yticklabels(["%.2f"%ytick for ytick in yticks])
            twinxs[i][j].spines["right"].set_color(colors[j])
            twinxs[i][j].spines["right"].set_linewidth(2)

            # translate y axis with some distance
            twinxs[i][j].spines["right"].set_position(("axes", 1 + (j-1)*0.08)) if j >= 1 else None

            # set ylims
            ylims = [-ylim*1.25, ylim*1.25]
            twinxs[i][j].set_ylim(ylims[0], ylims[1])

            # set xlims
            twinxs[i][j].set_xlim(ecutwfc_min - 10, ecutwfc_max + 10)

        # add grid
        twinxs[i][0].grid(grid, alpha = 0.1)
        # add xticks
        twinxs[i][0].set_xticks(xticks)
        # add subtitle at right, text left aligned
        pspotid_style = {"horizontalalignment": "right", "verticalalignment": "top", 
                         "transform": twinxs[i][0].transAxes, "fontsize": fontsize,
                         "backgroundcolor": "white"}
        twinxs[i][0].text(0.995, 0.9, subtitles[i], **pspotid_style)
        twinxs[i][0].text(0.995, 0.3, "$Z$ = %d"%z_vals[i], **pspotid_style) if z_vals is not None else None
        twinxs[i][0].axhline(0, color="black", linewidth=2, alpha=0.1)

        conv_marker_style = {"markersize": 15, "markerfacecolor": "none", "markeredgecolor": "red", "markeredgewidth": 2,
                             "zorder": 10, "alpha": 0.5}
        if highlight_xs is not None:
            # add a circle at the x position, on the topmost layer
            twinxs[i][0].plot(highlight_xs[i][1], 0, "o", **conv_marker_style)
    # set fontsize for xticks
    twinxs[-1][0].set_xticklabels(xticklabels, fontsize=fontsize)

    # subplot size adjustment, leave no space between subplots vertically
    plt.subplots_adjust(hspace=0.0, wspace=0.2)

    # xtitle, only add xtitle to the last row
    xtitle = kwargs.get("xtitle", None)
    if xtitle is not None:
        plt.text(0.5, 0.05, xtitle, **xtitle_styles)
    # ytitle, only add ytitle to the left-middle
    ytitle = kwargs.get("ytitle", None)
    if ytitle is not None:
        ytitle = [ytitle] if isinstance(ytitle, str) else ytitle
        for i in range(len(ytitle)):
            x = 0.075 if i == 0 else 0.945 + (i - 1)*0.05
            plt.text(x, 0.5, ytitle[i], **ytitle_styles)
    # suptitle
    suptitle = kwargs.get("suptitle", None)
    if suptitle is not None:
        plt.suptitle(suptitle, fontsize=fontsize*1.5)
    # supcomment, below the suptitle
    supcomment = kwargs.get("supcomment", None)
    if supcomment is not None:
        supcomment_style = {"ha": "center", "va": "center", "transform": fig.transFigure, 
                            "fontsize": fontsize, "style": "italic"}
        plt.text(0.5, 0.91, supcomment, **supcomment_style)

    # set overall fontstyle
    plt.rcParams["font.family"] = "Arial"
    return fig, ax

def discrete_logplots(xs: list, ys: list, **kwargs):

    nsubplts = len(xs)
    nlines = len(ys[0])
    for i in range(1, nsubplts):
        assert len(xs[i]) == nlines, f"xs[{i}] has different length than xs[0]: {len(xs[i])} != {nlines}"
        assert len(ys[i]) == nlines, f"ys[{i}] has different length than ys[0]: {len(ys[i])} != {nlines}"
    
    xtitle, ytitle = kwargs.get("xtitle", None), kwargs.get("ytitle", None)
    ysymbols, labels = kwargs.get("ysymbols", None), kwargs.get("labels", None)
    colors = kwargs.get("colors", None)
    markers, markersizes = kwargs.get("markers", None), kwargs.get("markersizes", None)
    linestyles, linewidths = kwargs.get("linestyles", None), kwargs.get("linewidths", None)
    alpha = kwargs.get("alpha", 1.0)
    grid = kwargs.get("grid", True)
    fontsize = kwargs.get("fontsize", 12)
    nrows, subplotsize = kwargs.get("nrows", 1), kwargs.get("subplotsize", (10, 10))
    highlight_ys = kwargs.get("highlight_ys", None)

    npspots = nlines
    labels = ["data " + str(i) for i in range(npspots)] if labels is None else labels
    assert len(labels) == npspots
    nprptys = nsubplts
    ytitle = ["subplot " + str(i) for i in range(nprptys)] if ytitle is None else ytitle
    if isinstance(ytitle, str):
        ytitle = [ytitle] * nprptys
    assert len(ytitle) == nprptys

    colors = amefs.styles_factory(property="color", ndim=npspots) if colors is None else colors
    assert len(colors) == npspots
    markers = amefs.styles_factory(property="marker", ndim=npspots) if markers is None else markers
    assert len(markers) == npspots
    markersizes = amefs.styles_factory(property="markersize", ndim=npspots) if markersizes is None else markersizes
    assert len(markersizes) == npspots
    linestyles = amefs.styles_factory(property="linestyle", ndim=npspots) if linestyles is None else linestyles
    assert len(linestyles) == npspots
    linewidths = amefs.styles_factory(property="linewidth", ndim=npspots) if linewidths is None else linewidths
    assert len(linewidths) == npspots
    alpha = amefs.styles_factory(property="alpha", val=alpha, ndim=npspots)
    assert len(alpha) == npspots
    
    # create figure and axes
    ncols = nprptys // nrows + (nprptys % nrows > 0)
    fig, ax = plt.subplots(nrows, ncols, 
                           figsize=(subplotsize[0] * ncols, subplotsize[1] * nrows + 4), 
                           squeeze=False)

    # plot
    for i in range(nprptys): # for each property
        for j in range(npspots): # for each pseudopotential
            # because the last value of y is always to be 0, so it is not plotted
            styles = {"color": colors[j], "marker": markers[j], "markersize": markersizes[j], 
                      "linestyle": linestyles[j], "linewidth": linewidths[j], "alpha": alpha[j]}
            ax[0, i].plot(xs[i][j][:-1], np.abs(np.array(ys[i][j]))[:-1], label=labels[j], **styles)
        # add threshold line
        if highlight_ys is not None:
            thr_style = {"color": "red", "linewidth": 5, "alpha": 0.1, 
                         "label": ysymbols[i] + " < " + "%.2e"%highlight_ys[i]}
            ax[0, i].axhline(highlight_ys[i], **thr_style)

        ax[0, i].set_yscale("log")
        ax[0, i].grid(grid)
        # add two legends, the first is normal legend, the second is threshold legend
        lns = ax[0, i].get_lines()[:-1] if highlight_ys is not None else ax[0, i].get_lines()
        ax[0, i].legend(handles=lns, fontsize=fontsize, shadow=True, loc="upper right")
        if highlight_ys is not None:
            thr_legend = Legend(ax[0, i], [ax[0, i].get_lines()[-1]], [ax[0, i].get_lines()[-1].get_label()], 
                                fontsize=fontsize, frameon=True, shadow=False, loc="lower left")
            ax[0, i].add_artist(thr_legend)
        # set x/y title
        ax[0, i].set_xlabel(xtitle, fontsize=fontsize * 1.2)
        ax[0, i].set_ylabel(ytitle[i], fontsize=fontsize * 1.2)
        # set x/y label fontsize
        ax[0, i].tick_params(axis="both", labelsize=fontsize)

    # suptitle
    suptitle = kwargs.get("suptitle", None)
    if suptitle is not None:
        plt.suptitle(suptitle, fontsize=fontsize*1.5)
    # supcomment, below the suptitle
    supcomment = kwargs.get("supcomment", None)
    if supcomment is not None:
        # italicize the supcomment
        supcomment_style = {"ha": "center", "va": "center", "transform": fig.transFigure, 
                            "fontsize": fontsize, "style": "italic"}
        plt.text(0.5, 0.925, supcomment, **supcomment_style)

    return fig, ax

def shift_lineplots(xs: list, ys: list, **kwargs):
    """draw all lines in one figure, but add shift to distinguish between
    different pseudopotentials, therefore there are shifts needed to be
    defined as many as properties. In this mode, EVERY SINGLE LINE can have
    different length or number of data points, but x-y should still have
    the same length.
    
    The shift is of the same property between different pseudopotentials, 
    thus the ylim of the whole figure should be shift*npspots.
    
    Design to let all properties of the same pseudopotential share the same
    color, but with different marker and linestyle. The loop should still be
    2-dimensional, the one is pseudopotential and the other is property. If
    pseudopotential is the outer loop, it means for each pseudopotential, 
    firstly to draw all properties, then move to the next pseudopotential.
    Vice versa, if property is the outer loop, it means first draw one property
    for all pseudopotentials, then move to the next property.
    """
    
    ld = kwargs.get("ld", None)
    assert ld in ["pseudopotential", "property"] 
    # leading dimension is necessary, therefore if not given, raise an error
    if ld == "pseudopotential":
        npspots = len(xs)
        assert len(ys) == npspots
        nprptys = len(xs[0])
        for i in range(1, npspots):
            assert len(xs[i]) == nprptys, f"xs[{i}] has different length than xs[0]: {len(xs[i])} != {nprptys}"
            assert len(ys[i]) == nprptys, f"ys[{i}] has different length than ys[0]: {len(ys[i])} != {nprptys}"
    elif ld == "property":
        nprptys = len(xs)
        assert len(ys) == nprptys
        npspots = len(xs[0])
        for i in range(1, nprptys):
            assert len(xs[i]) == npspots, f"xs[{i}] has different length than xs[0]: {len(xs[i])} != {npspots}"
            assert len(ys[i]) == npspots, f"ys[{i}] has different length than ys[0]: {len(ys[i])} != {npspots}"
    else:
        raise ValueError("Unknown loop direction")
    
    colors = kwargs.get("colors", None)
    colors = amefs.styles_factory(property="color", ndim=npspots) if colors is None else colors
    assert len(colors) == npspots

    markers = kwargs.get("markers", None)
    markers = amefs.styles_factory(property="marker", ndim=nprptys) if markers is None else markers
    assert len(markers) == nprptys

    linestyles = kwargs.get("linestyles", None)
    linestyles = amefs.styles_factory(property="linestyle", ndim=nprptys) if linestyles is None else linestyles
    assert len(linestyles) == nprptys

    shifts = kwargs.get("shifts", None)
    assert len(shifts) == nprptys # if shifts is not given, then will cause a assertation failure

    # add pseudopotential name, z_valence information
    pspotnames = kwargs.get("subtitles", None)
    assert pspotnames is not None and len(pspotnames) == npspots
    z_vals = kwargs.get("z_vals", None)
    assert z_vals is None or len(z_vals) == npspots

    highlight_xs = kwargs.get("highlight_xs", None)
    assert highlight_xs is None or len(highlight_xs) == npspots

    ysymbols = kwargs.get("ysymbols", None)
    assert ysymbols is not None and len(ysymbols) == nprptys
    
    fontsize = kwargs.get("fontsize", 12)
    # create figure and axes
    fig, ax = plt.subplots(1, 1, figsize=(20, 10))
    twinxs = [[None for _ in range(nprptys)] for _ in range(npspots)]
    # record xticks, the basic one is [20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200]
    # but if the maximum ecutwfc is larger than 200, then add it to the list
    # then the xticklabels should be the same length as xticks
    xticks = [20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200]
    for i in range(nprptys):
        for j in range(npspots):
            invj = npspots - j - 1
            # create a new twinx for each pseudopotential
            twinxs[j][i] = ax.twinx()
            # only add pseudopotential name and z_valence information for the first property
            if i == 0:
                xpos = 0.995
                dypos = 1/npspots
                ypos1 = invj * dypos + dypos*0.9
                ypos2 = invj * dypos + dypos*0.3
                pspotid_style = {"horizontalalignment": "right", "verticalalignment": "top", 
                                 "transform": twinxs[j][i].transAxes, "fontsize": fontsize}
                twinxs[j][i].text(xpos, ypos1, pspotnames[j], **pspotid_style)
                twinxs[j][i].text(xpos, ypos2, "$Z$ = %d"%z_vals[j], **pspotid_style) if z_vals is not None else None
                # also add circle at the x position for the converged ecutwfc
                conv_marker_style = {"markersize": 15, "markerfacecolor": "none", "markeredgecolor": "red", "markeredgewidth": 2,
                                     "zorder": 10, "alpha": 0.5}
                twinxs[j][i].plot(highlight_xs[j][1], shifts[i] * invj, "o", **conv_marker_style) if highlight_xs is not None else None

            # prepare data and make shift
            x, y = np.array(xs[j][i]), np.array(ys[j][i]) + shifts[i] * invj
            # check if xticks should be updated
            xticks = np.unique(np.concatenate((xticks, x)))
            # all properties of the same pseudopotential share the same color
            # but with different marker and linestyle
            style = {"color": colors[j], "marker": markers[i], "linestyle": linestyles[i],
                     "linewidth": 1, "markersize": 5, "alpha": 0.8, "label": ysymbols[i]}
            twinxs[j][i].plot(x, y, **style)
            # set ylim
            twinxs[j][i].set_ylim(-shifts[i]/2, shifts[i] * (npspots - 1/2))
            # turn off yticks
            twinxs[j][i].set_yticks([]) if j != 0 else None
            # turn off yticklabels
            twinxs[j][i].set_yticklabels([]) if j != 0 else None
            
        # translate y axis with some distance for j == -1
        twinxs[0][i].spines["right"].set_position(("axes", 1 + i*0.04))
        # set thickness of yaxis to 2
        twinxs[0][i].spines["right"].set_linewidth(2)
        # set fontsize
        twinxs[0][i].tick_params(axis="y", labelsize=fontsize)

    # set xlims to be xmin - 10, xmax + 10
    xmin, xmax = np.min(xticks), np.max(xticks)
    for i in range(npspots):
        for j in range(nprptys):
            twinxs[i][j].set_xlim(xmin - 10, xmax + 10)
    # set xticklabels
    xticklabels = []
    for ecutwfc in xticks:
        xticklabels.append(str(ecutwfc)) if ecutwfc in xticks else xticklabels.append("")
    twinxs[-1][0].set_xticks(xticks)
    # set fontsize
    ax.tick_params(axis="x", labelsize=fontsize)

    # only add legends for the first pseudopotential for all properties
    lns = [twinxs[0][i].get_lines()[-1] for i in range(nprptys)]
    labels = [ysymbols[i] for i in range(nprptys)]
    legend_style = {"fontsize": fontsize, "shadow": True, "loc": "upper right", 
                    "bbox_to_anchor": (0, 0)}
    ax.legend(lns, labels, **legend_style)

    # add one "y = 0" for reference for each pseudopotential
    for j in range(npspots):
        twinxs[j][0].axhline(0 + shifts[0] * j, color="black", linewidth=1, alpha=0.1)
    # turn on xgrid
    ax.grid(True, alpha=0.1)
    # and another line to seperate different pseudopotentials
    for j in range(npspots - 1):
        twinxs[j][0].axhline(shifts[0] * (j + 1 - 1/2), color="black", linewidth=1, alpha=0.8)
    # turn off yticks
    ax.set_yticks([])
    # turn off yticklabels
    ax.set_yticklabels([])
    # set xtitle
    xtitle = kwargs.get("xtitle", None)
    if xtitle is not None:
        xtitle_style = {"ha": "center", "va": "center", "transform": fig.transFigure, 
                        "fontsize": fontsize * 1.2}
        plt.text(0.5, 0.05, xtitle, **xtitle_style)
    # set ytitles at right
    ytitles = kwargs.get("ytitle", None)
    if ytitles is not None:
        style = {"ha": "center", "va": "bottom", "transform": fig.transFigure, 
                 "fontsize": fontsize}
        ytitles = [ytitles] if isinstance(ytitles, str) else ytitles
        for i in range(len(ytitles)):
            xpos, ypos = 0.91 + i*0.03, 0.90
            plt.text(xpos, ypos, ysymbols[i], **style)
    # set suptitle
    suptitle = kwargs.get("suptitle", None)
    if suptitle is not None:
        plt.suptitle(suptitle, fontsize=fontsize * 1.5)
    # set supcomment
    supcomment = kwargs.get("supcomment", None)
    if supcomment is not None:
        supcomment_style = {"ha": "center", "va": "center", "transform": fig.transFigure, 
                            "fontsize": fontsize, "style": "italic"}
        plt.text(0.5, 0.925, supcomment, **supcomment_style)

    return fig, ax
            
if __name__ == "__main__":
    # this should not be changed no matter what kind of postprocess is!
    run()