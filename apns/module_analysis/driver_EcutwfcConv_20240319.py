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
                    ncols: int = 1,             # optional, number of columns
                    subtitles: list = None,     # optional, subtitles for each subplot
                    labels: list = None,        # optional, labels for each line
                    colors: list = None,        # optional, colors for each line
                    markers: list = None,       # optional, markers for each line
                    markersizes: list = None,   # optional, markersizes for each line
                    linestyles: list = None,    # optional, linestyles for each line
                    linewidths: list = None,    # optional, linewidths for each line
                    alpha: float = 1.0,         # optional, alpha for each line
                    grid: bool = True,          # optional, grid for each subplot
                    xtitle: str = None,         # optional, xtitle for the whole figure
                    ytitle: str|list = None,    # optional, ytitle for the whole figure
                    suptitle: str = None,       # optional, suptitle for the whole figure
                    fontsize: int = 12,         # optional, fontsize for the whole figure
                    ):
    # xs provide a flexible way to plot multiple subplots, each subplot has its own x
    # although for pseudopotential ecutwfc convergence test this is not necessary, 
    # for all pseudopotential for one element will be identical
    # force xs and ys to match
    nsubplts = len(xs)
    # for ys, the flexibility is that each subplot can have not only one lines, but
    # should all subplots have the same number of lines. Say for each pseudopotential
    #, not only energies but also pressure and/or other properties can be plotted
    # together. Therefore, ys is a 3-dimensional list, where the first dimension is
    # about subplot, the second is about line, the third is about data.
    assert nsubplts == len(ys)
    # force all subplots to have same number of lines
    nlines = len(ys[0])
    for i in range(1, nsubplts):
        assert len(ys[i]) == nlines
    print("Number of subplots:", nsubplts)
    print("Number of lines to draw in each subplot:", nlines)
    assert ncols > 0
    assert nsubplts >= ncols
    assert nlines > 0
    assert nsubplts <= 20

    # subtitles
    subtitles = ["subplot " + str(i) for i in range(nsubplts)] if subtitles is None else subtitles
    assert len(subtitles) == nsubplts
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
    nrows = nsubplts // ncols + (nsubplts % ncols > 0)
    fig, ax = plt.subplots(nrows, ncols, figsize=(20*ncols, nrows), squeeze=False)
    
    twinxs = []
    # plot
    for i in range(nsubplts):
        row, col = i // ncols, i % ncols
        twinxs.append([])
        for j in range(len(ys[i])):
            ys[i][j] = np.array(ys[i][j])
            twinxs[i].append(ax[row, col] if j == 0 else ax[row, col].twinx())
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
            yticks = [-ylim*0.75, 0, ylim*0.75]
            ylims = [-ylim*1.25, ylim*1.25]
            twinxs[i][j].set_yticks(yticks)
            twinxs[i][j].set_ylim(ylims[0], ylims[1])
            
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
    # suptitle
    if suptitle is not None:
        plt.suptitle(suptitle, fontsize=fontsize*1.5)
    return fig, ax

def log_plots(xs: list, ys: list,         # compulsory
              nrows: int = 1,             # optional, number of rows
              subtitles: list = None,     # optional, subtitles for each subplot
              labels: list = None,        # optional, labels for each line
              colors: list = None,        # optional, colors for each line
              markers: list = None,       # optional, markers for each line
              markersizes: list = None,   # optional, markersizes for each line
              linestyles: list = None,    # optional, linestyles for each line
              linewidths: list = None,    # optional, linewidths for each line
              alpha: float = 1.0,         # optional, alpha for each line
              grid: bool = True,          # optional, grid for each subplot
              xtitle: str = None,         # optional, xtitle for the whole figure
              ytitle: str|list = None,    # optional, ytitle for the whole figure
              suptitle: str = None,       # optional, suptitle for the whole figure
              fontsize: int = 12,         # optional, fontsize for the whole figure
              ):
    # has the same format of input as function stack_lineplots, where y has
    # nsubplot * nlines of data to plot, x has nsubplot of data to plot
    
    # what is different is the actuall number of subplots is nline, say the
    # number of kinds of data. Therefore, for each line, all "nsubplot" data 
    # will be plotted in a single subplot, and y in log scale

    # first to rearrange data
    nsubplts = len(xs) # the nsubplts-dimensional data is originally designed for
    # plot like n pseudopotentials in one shot. Therefore, the first dimension is
    # about pseudopotential, the second is about element, the third is about ecutwfc
    # for present Log plot, for one element, number of suplots dependes on number
    # of lines saved in ys
    nlines = len(ys[0])
    # nlines is the real number of subplots. Should be careful about the ys dimension
    #, the first is about subplot, the second is about line, the third is about data.
    # Therefore, for one certain property, the ydata should be ys[i][j] and should go
    # over all i for each j.
    for i in range(1, nsubplts):
        assert len(ys[i]) == nlines
    print("Number of subplots:", nlines)
    print("Number of lines to draw in each subplot:", nsubplts)
    # labels, this variable is designed for each pseudopotential, therefore, its
    # length should be nsubplts
    labels = ["data " + str(i) for i in range(nsubplts)] if labels is None else labels
    assert len(labels) == nsubplts
    # subtitles, this variable is designed for each line, therefore, its length should
    # be nlines. Each subtitle corresponds to one property to test for convergence
    subtitles = ["subplot " + str(i) for i in range(nlines)] if subtitles is None else subtitles
    assert len(subtitles) == nlines
    # ytitle
    ytitle = ["subplot " + str(i) for i in range(nlines)] if ytitle is None else ytitle
    if isinstance(ytitle, str):
        ytitle = [ytitle] * nlines
    assert len(ytitle) == nlines
    # styles
    # colors
    _colorpool = ["#2A306A", "#24B5A5", "#1DB8A8", "#015BBC", "#EACE4F"]
    if colors is not None:
        assert len(colors) == nsubplts or len(colors) == 1
        if len(colors) == 1:
            colors = colors * nsubplts
    else:
        colors = [_colorpool[i % len(_colorpool)] for i in range(nsubplts)]
    # markers
    _markerpool = ["o", "s", "D", "v", "^", "<", ">", "p", "P", "*", "h", "H", "+", "x", "X", "|", "_"]
    if markers is not None:
        assert len(markers) == nsubplts or len(markers) == 1
        if len(markers) == 1:
            markers = markers * nsubplts
    else:
        markers = [_markerpool[i % len(_markerpool)] for i in range(nsubplts)]
    # markersizes
    if markersizes is not None:
        assert len(markersizes) == nsubplts or len(markersizes) == 1
        if len(markersizes) == 1:
            markersizes = markersizes * nsubplts
    else:
        markersizes = [5] * nsubplts
    # linestyles
    _linestylepool = ["-", "--", "-.", ":"]
    if linestyles is not None:
        assert len(linestyles) == nsubplts or len(linestyles) == 1
        if len(linestyles) == 1:
            linestyles = linestyles * nsubplts
    else:
        linestyles = [_linestylepool[i % len(_linestylepool)] for i in range(nsubplts)]
    # linewidths
    if linewidths is not None:
        assert len(linewidths) == nsubplts or len(linewidths) == 1
        if len(linewidths) == 1:
            linewidths = linewidths * nsubplts
    else:
        linewidths = [1] * nsubplts
    # alpha
    if isinstance(alpha, float):
        alpha = [alpha] * nsubplts
    else:
        assert len(alpha) == nsubplts

    # create figure and axes
    fig, ax = plt.subplots(1, nlines, figsize=(10 * nlines, 10), squeeze=False)

    # plot
    for i in range(nlines):
        for j in range(nsubplts):
            ys[j][i] = np.array(ys[j][i])
            y = np.abs(ys[j][i])
            # y = np.log10(np.abs(ys[j][i]) + 1e-10)
            ax[0, i].plot(xs[j], y, 
                          label=labels[j], 
                          color=colors[j], 
                          marker=markers[j], 
                          markersize=markersizes[j], 
                          linestyle=linestyles[j], 
                          linewidth=linewidths[j], 
                          alpha=alpha[j])
            ax[0, i].set_yscale("log")
            ax[0, i].grid(grid)
            ax[0, i].legend(fontsize=fontsize)
            # set xtitle
            ax[0, i].set_xlabel(xtitle, fontsize=fontsize)
            ax[0, i].set_ylabel(ytitle[i], fontsize=fontsize)
    # suptitle
    if suptitle is not None:
        plt.suptitle(suptitle, fontsize=fontsize*1.5)
    
    # set overall fontstyle
    plt.rcParams["font.family"] = "Arial"

    return fig, ax

if __name__ == "__main__":
    #print(testname_pspotid("pd04"))
    result = run("../11549321")
    elements, pspotids, eks_x, eks_y = categorize_byelement(result[0], result[1])
    elements, pspotids, prs_x, prs_y = categorize_byelement(result[2], result[3])

    conv_eks = result[0]
    conv_prs = result[2]
    
    ncols = len(elements)
    element = elements[0] # only get the 1st element to test
    pspotid = pspotids[0] # only get the 1st element to test
    pspotid = [testname_pspotid(id)[0] for id in pspotid]
    ecutwfc = eks_x[0]  # only get the 1st element to test
    eks_y = eks_y[0]  # only get the 1st element to test
    prs_y = prs_y[0]  # only get the 1st element to test
    ys = [[eks_y[i], prs_y[i]] for i in range(len(ecutwfc))]

    fig, ax = log_plots(xs=ecutwfc, ys=ys, 
                        nrows=1,
                        subtitles=["Absolute DeltaE(KS, Kohn-Sham) (eV/atom)", "Absolute DeltaP (kbar)"],
                        labels=pspotid,
                        suptitle=element,
                        xtitle="ecutwfc (Ry)",
                        ytitle=["Absolute DeltaE(KS, Kohn-Sham) (eV/atom)", "Absolute DeltaP (kbar)"],
                        fontsize=12)
    plt.savefig(f"{element}_logplot.svg")


    fig, ax = stack_lineplots(xs=ecutwfc, ys=ys,
                              ncols=1,
                              subtitles=pspotid,
                              grid=True,
                              xtitle="ecutwfc (Ry)",
                              ytitle=["DeltaE(KS, Kohn-Sham) (eV/atom)", "Pressure (kbar)"],
                              suptitle=element,
                              fontsize=12)
    plt.savefig(f"{element}.svg")
    import apns.module_analysis.external_frender.htmls as amaeh
    import apns.module_analysis.external_frender.markdown as amaem
    html = amaeh.pseudopotentials(element=element, 
                                 xc_functional="PBE", 
                                 software="ABACUS",
                                 fconv=f"{element}.svg")
    with open("index.html", "w") as f:
        f.write(html)
