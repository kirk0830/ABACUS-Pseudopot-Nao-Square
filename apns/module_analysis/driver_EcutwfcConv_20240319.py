"""
Driver Name: EcutwfcConv_20240319

Driver description:
This is driver for APNS pages pseudopotential convergence test
Two kinds of data are post-processed:
1. ecutwfc - energy per atom
2. ecutwfc - pressure of whole cell
"""

import json
def testname_pspotid(testname: str):
    """return pspotid from testname"""
    with open("download/pseudopotentials/description.json", "r") as f:
        description = json.load(f)
    for key in description.keys():
        if key.replace("_", "") == testname:
            val = description[key]
            with open(val + "/description.json", "r") as f:
                data = json.load(f)
            data = {key: data[key] for key in data.keys() if key != "files"}
            label = data["kind"].upper() + "-" + data["version"]
            label += " (" + data["appendix"] + ")" if data["appendix"] != "" else ""
            label = label.strip()
            return label, key

def categorize_byelement(labels, data):
    assert len(labels) == len(data)
    for label in labels:
        assert len(label) == 4
        element, mpid, pspotid, ecutwfc = label
        assert isinstance(element, str)
        assert isinstance(mpid, str)
        assert isinstance(pspotid, str)
        assert isinstance(ecutwfc, float)
    for piece in data:
        assert len(piece) == 2
        assert isinstance(piece[0], list)
        assert isinstance(piece[1], list)
        assert len(piece[0]) == len(piece[1])
    elements = []
    pspotids = []
    decomposed_x = []
    decomposed_y = []
    for itst in range(len(labels)):
        element, mpid, pspotid, ecutwfc = labels[itst]
        if element not in elements:
            elements.append(element)
            pspotids.append([pspotid])
            decomposed_x.append([])
            decomposed_y.append([])
        else:
            pspotids[elements.index(element)].append(pspotid)
        decomposed_x[elements.index(element)].append(data[itst][0])
        decomposed_y[elements.index(element)].append(data[itst][1])
    return elements, pspotids, decomposed_x, decomposed_y

import apns.module_analysis.postprocess.conv.ecutwfc_eks as amapeks
import apns.module_analysis.postprocess.conv.ecutwfc_press as amapprs
def run(path: str, ethr: float = 1e-3, pthr: float = 0.1):
    # the first is converged data, the second is list of two lists, all data
    first_eks, second_eks = amapeks.run(path, ethr)
    # the first is converged data, the second is list of two lists, all data
    first_prs, second_prs = amapprs.run(path, pthr)

    return first_eks, second_eks, first_prs, second_prs

import matplotlib.pyplot as plt
import numpy as np
def stack_lineplots(xs: list, ys: list,         # compulsory
                    logy: bool = False,         # optional
                    ncols: int = 1,
                    subtitles: list = None,     # optional,
                    labels: list = None,        # optional
                    colors: list = None,        # optional
                    markers: list = None,       # optional
                    markersizes: list = None,   # optional
                    linestyles: list = None,    # optional
                    linewidths: list = None,    # optional
                    alpha: float = 1.0,         # optional
                    grid: bool = True,          # optional
                    xtitle: str = None,         # optional
                    ytitle: str = None,         # optional
                    fontsize: int = 12,         # optional
                    ):
    """xs is list of x, ys is list of list of y, 
    allowing plot x with multiple y for each subplot"""
    # force xs and ys to match
    nsubplt = len(xs)
    assert nsubplt == len(ys)
    # force all subplots to have same number of lines
    nlines = len(ys[0])
    for i in range(1, nsubplt):
        assert len(ys[i]) == nlines
    print("Number of subplots:", nsubplt)
    print("Number of lines to draw in each subplot:", nlines)
    assert ncols > 0
    assert nsubplt >= ncols
    assert nlines > 0
    assert nsubplt <= 20

    # subtitles
    subtitles = ["subplot " + str(i) for i in range(nsubplt)] if subtitles is None else subtitles
    assert len(subtitles) == nsubplt
    # labels
    labels = ["data " + str(i) for i in range(len(ys[0]))] if labels is None else labels
    assert len(labels) == len(ys[0])
    # styles: all styles are universal for each subplots, for each column and each row, therefore
    # only nlines styles are needed to be provided
    def _get_styles_frompool(pool: int, i: int):
        return pool[i % len(pool)]
    # colors
    _colorpool = ["#2A306A", "#24B5A5", "#1DB8A8", "#015BBC", "#EACE4F"]
    if colors is not None:
        assert len(colors) == nlines or len(colors) == 1
        if len(colors) == 1:
            colors = colors * nlines
    else:
        colors = [_get_styles_frompool(_colorpool, i) for i in range(nlines)]
    # markers
    _markerpool = ["o", "s", "D", "v", "^", "<", ">", "p", "P", "*", "h", "H", "+", "x", "X", "|", "_"]
    if markers is not None:
        assert len(markers) == nlines or len(markers) == 1
        if len(markers) == 1:
            markers = markers * nlines
    else:
        markers = [_get_styles_frompool(_markerpool, i) for i in range(nlines)]
    # markersizes
    if markersizes is not None:
        assert len(markersizes) == nlines or len(markersizes) == 1
        if len(markersizes) == 1:
            markersizes = markersizes * nlines
    else:
        markersizes = [5] * nlines
    # linestyles
    _linestylepool = ["-", "--", "-.", ":"]
    if linestyles is not None:
        assert len(linestyles) == nlines or len(linestyles) == 1
        if len(linestyles) == 1:
            linestyles = linestyles * nlines
    else:
        linestyles = [_get_styles_frompool(_linestylepool, i) for i in range(nlines)]
    # linewidths
    if linewidths is not None:
        assert len(linewidths) == nlines or len(linewidths) == 1
        if len(linewidths) == 1:
            linewidths = linewidths * nlines
    else:
        linewidths = [1] * nlines
    # alpha
    if isinstance(alpha, float):
        alpha = [alpha] * nlines
    else:
        assert len(alpha) == nlines

    # create figure and axes
    nrows = nsubplt // ncols + (nsubplt % ncols > 0)
    fig, ax = plt.subplots(nrows, ncols, figsize=(60 * ncols, 10 * nrows), squeeze=False)
    
    twinxs = []
    # plot
    for i in range(nsubplt):
        row, col = i // ncols, i % ncols
        twinxs.append([])
        for j in range(len(ys[i])):
            if j == 0:
                twinxs[i].append(ax[row, col])
            else:
                twinxs[i].append(ax[row, col].twinx())
            twinxs[i][j].plot(xs[i], ys[i][j], 
                              label=labels[j], 
                              color=colors[j], 
                              marker=markers[j], 
                              markersize=markersizes[j], 
                              linestyle=linestyles[j], 
                              linewidth=linewidths[j], 
                              alpha=alpha[j])
            # add yticks for each line
            ylim = np.max(np.abs(ys[i][j]))
            twinxs[i][j].set_yticks([-ylim*0.75, 0, ylim*0.75])
            twinxs[i][j].set_ylim(-ylim*1.25, ylim*1.25)
            
        # add grid
        twinxs[i][0].grid(grid)
        # add subtitle at right, text left aligned
        twinxs[i][0].text(0.995, 0.9, subtitles[i], 
                          horizontalalignment="right", 
                          verticalalignment="top", 
                          transform=twinxs[i][0].transAxes,
                          fontsize=fontsize,
                          backgroundcolor="white")
        twinxs[i][0].axhline(0, color="black", linewidth=5, alpha=0.1)
    # subplot size adjustment, leave no space between subplots vertically
    plt.subplots_adjust(hspace=0.0, wspace=0.2)

    # set overall fontstyle
    plt.rcParams["font.family"] = "Arial"

    # xtitle, only add xtitle to the last row
    if xtitle is not None:
        plt.text(0.5, 0.05, xtitle, 
                 ha="center", va="center", 
                 transform=fig.transFigure,
                 fontsize=fontsize)
    # ytitle, only add ytitle to the left-middle
    if ytitle is not None:
        if isinstance(ytitle, str):
            ytitle = [ytitle]
        for i in range(len(ytitle)):
            x = 0.075 if i == 0 else 0.95 + (i - 1)*0.05
            plt.text(x, 0.5, ytitle[i], 
                     ha="center", va="center",
                     rotation="vertical",
                     transform=fig.transFigure,
                     fontsize=fontsize)

    return fig, ax

if __name__ == "__main__":
    #print(testname_pspotid("pd04"))
    result = run("../11549321")
    elements, pspotids, eks_x, eks_y = categorize_byelement(result[0], result[1])
    elements, pspotids, prs_x, prs_y = categorize_byelement(result[2], result[3])

    ncols = len(elements)
    element = elements[0] # only get the 1st element to test
    pspotid = pspotids[0] # only get the 1st element to test
    pspotid = [testname_pspotid(id)[0] for id in pspotid]
    ecutwfc = eks_x[0]  # only get the 1st element to test
    eks_y = eks_y[0]  # only get the 1st element to test
    prs_y = prs_y[0]  # only get the 1st element to test
    ys = [[eks_y[i], prs_y[i]] for i in range(len(ecutwfc))]

    fig, ax = stack_lineplots(xs=ecutwfc, 
                              ys=ys,
                              logy=False,
                              ncols=1,
                              subtitles=pspotid,
                              markers=["s", "o"],
                              markersizes=[7],
                              linestyles=["--", "-"],
                              linewidths=[1],
                              alpha=0.8,
                              grid=True,
                              xtitle="ecutwfc (Ry)",
                              ytitle=["DeltaE(KS, Kohn-Sham) (eV/atom)", "Pressure (kbar)"],
                              fontsize=12)
    plt.show()