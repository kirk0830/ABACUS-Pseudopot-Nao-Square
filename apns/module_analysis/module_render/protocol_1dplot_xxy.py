import matplotlib.pyplot as plt
import numpy as np

def test():
    pass

"""
single element XX-Y plot

Parameters
----------
data : dict
    data to be plotted, in the format of python dictionary:
    data = {
        "pseudopotentials: {
            "SG15-1.0": {
                "Zval": 13,
                ...
            },
            ...
        },
        "calculation": {
            "pw": {
                "ecutwfc": [10, 15, 20, 25, 30, 35, 40, 45, 50],
                },
            "lcao": {
                "rcut": [6, 7, 8, 9, 10],
                "basis": ["DZP", "TZDP"]
                },
        }
        "results": {
            "SG15-1.0": {
                "pw": [...]
                "lcao": {
                    "DZP": [...],
                    "TZDP": [...]
                    }
                },
            ...
            }
        }

styles : dict
    styles to be used for plotting, in the format of python dictionary:
    styles = {
        "figure": {
            "size": (10, 10),
            "font": "Arial",
            "fontsize": 15,
            "color_shininess": 0.95
        },
        "pseudopotentials: {
            "label": {
                "font": "Arial",
                "size": 15
                }
        },
        "calculation": {
            "pw": {
                "line": {
                    "color": "blue",
                    "linestyle": "-",
                    "linewidth": 1.0
                    },
                "marker": {
                    "color": (142, 207, 201),
                    "marker": "*",
                    "ms": 10
                    }
                "xlim": [5, 100],
            }
            "lcao": {
                "line": {
                    "color": [(255, 190, 122), (130, 176, 210)],
                    "linestyle": "-",
                    "linewidth": 1.0
                    },
                "marker": {
                    "color": "red",
                    "marker": "o",
                    "ms": 10
                    }
                "xlim": [5.5, 15],
                }
            }
        }          
"""
def single_element_xxy_plot(data: dict, styles: dict) -> None:

    # Create the figure and the first axis
    fig, ax1 = plt.subplots(nrows=len(data["pseudopotentials"]), ncols=1, sharex=True, sharey=True, figsize=styles["figure"]["size"])
    for i in range(len(data["pseudopotentials"])):
        # Generate the data
        x1 = data["calculation"]["pw"]["ecutwfc"]
        y1 = data["results"][data["pseudopotentials"][i]]["pw"]
        # Plot the first data set on the first axis wit linestyle - and symbol triangle up
        line_pw = ax1[i].plot(
            x1, y1, 
            marker = styles["calculation"]["pw"]["marker"]["marker"], 
            ms = styles["calculation"]["pw"]["marker"]["ms"], 
            linestyle = styles["calculation"]["pw"]["line"]["linestyle"], 
            color = styles["calculation"]["pw"]["line"]["color"], 
            label = 'pw')
        #ax1[i].set_xlabel('ecutwfc (Ry)')
        # name this line as pw
        # disable xticks for the first axis
        if i != len(data["pseudopotentials"])-1:
            ax1[i].set_xticks([])
        # Create the second axis and plot the second data set on it
        x2 = data["calculation"]["lcao"]["rcut"]
        y2 = data["results"][data["pseudopotentials"][i]]["lcao"][data["calculation"]["lcao"]["basis"][0]]
        ax2 = ax1[i].twiny()
        line_dzp = ax2.plot(
            x2, y2, 
            marker = styles["calculation"]["lcao"]["marker"]["marker"], 
            linestyle = styles["calculation"]["lcao"]["line"]["linestyle"], 
            color = styles["calculation"]["lcao"]["line"]["color"][0], 
            label = data["calculation"]["lcao"]["basis"][0])
        ax2.set_xlim(styles["calculation"]["lcao"]["xlim"][0], styles["calculation"]["lcao"]["xlim"][1])
        #ax2.set_xlabel('Rcut (a.u.)')
        # disable xticks for the second axis
        if i == 0:
            # add xlabel as "cutoff radius (a.u.)"
            ax2.set_xlabel('cutoff radius (a.u.)', fontsize=styles["figure"]["fontsize"])
        if i != 0:
            ax2.set_xticks([])
        x3 = data["calculation"]["lcao"]["rcut"]
        y3 = data["results"][data["pseudopotentials"][i]]["lcao"][data["calculation"]["lcao"]["basis"][1]]
        ax3 = ax1[i].twiny()
        # use rgb color deep green
        line_tzdp = ax3.plot(
            x3, y3, 
            marker = styles["calculation"]["lcao"]["marker"]["marker"], 
            linestyle = styles["calculation"]["lcao"]["line"]["linestyle"], 
            color = styles["calculation"]["lcao"]["line"]["color"][1], 
            label = data["calculation"]["lcao"]["basis"][1])
        # set xlim for ax3
        ax3.set_xlim(styles["calculation"]["lcao"]["xlim"][0], styles["calculation"]["lcao"]["xlim"][1])
        #ax3.set_xlabel('Rcut (a.u.)')
        # disable xticks for the second axis
        ax3.set_xticks([])
        # add the subplot title at bottom-right corner of each subplot
        ax1[i].text(
            0.98, 0.05, 
            data["pseudopotentials"][i] + "\nZval = " + str(data["pseudopotentials"][i]["Zval"]), 
            horizontalalignment='right', 
            verticalalignment='bottom', 
            transform=ax1[i].transAxes, 
            fontsize=styles["pseudopotentials"]["label"]["size"])
        # set ylim for each subplot
        ax1[i].set_ylim(-0.6, 0.1)
        # add black dashed line at y=0
        ax1[i].axhline(y=0, color='k', linestyle='--', alpha = 0.5)
        # add vertical grid lines, semi-transparent
        ax1[i].grid(axis='x', alpha=0.5)
        # only add legend for the last subplot
        if i == len(data["pseudopotentials"])-1:
            lines = line_pw + line_dzp + line_tzdp
            labels = [l.get_label() for l in lines]
            ax1[i].legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, -0.65), ncol=3, fontsize=styles["pseudopotentials"]["label"]["size"])
    plt.rcParams["font.family"] = styles["figure"]["font"]
    # only add xticks on the last subplot
    ax1[len(data["pseudopotentials"])-1].set_xticks(data["calculation"]["pw"]["ecutwfc"])
    # only add xlabel on the last subplot
    ax1[len(data["pseudopotentials"])-1].set_xlabel('ecutwfc (Ry)', fontsize=styles["figure"]["fontsize"])
    # set all font style to be Arial
    plt.rcParams["font.family"] = styles["figure"]["font"]
    # adjust size of subplots
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.1, hspace=0.0)
    # Add title at left of the whole figure, set direction as vertical, centered at the whole figure
    plt.suptitle('cell parameter error %', fontsize=styles["figure"]["fontsize"], x=0.05, y=0.5, rotation=90, va='center')
    #plt.xlabel('ecutwfc (Ry)')
    plt.ylabel('cell parameter error %')
    # save the figure as png
    plt.savefig('protocol_1dplot_xxy.png', dpi=300)
    # after plotting, clear the figure
    plt.clf()

if __name__ == "__main__":
    
    y = [-0.5, -0.25, -0.125, -0.0625, -0.03125, -0.015625, -0.0078125, -0.00390625, -0.001953125, -0.0009765625, -0.00048828125, -0.000244140625, -0.0001220703125, -6.103515625e-05, -3.0517578125e-05, -1.52587890625e-05, -7.62939453125e-06, -3.814697265625e-06, -1.9073486328125e-06, -4.76837158203125e-07, -1.1920928955078125e-07, -2.9802322387695312e-08, -7.450580596923828e-09, -1.862645149230957e-09]
    y = np.array(y)
    rcut = [6, 7, 8, 9, 10]
    rcut = np.array(rcut)
    ecutwfc = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 120, 140, 160, 180, 200]
    ecutwfc = np.array(ecutwfc)
    y_for_rcut = np.array([y[0], y[1], y[2], y[5], y[10]])

    pseudopots = ["SG15-1.0", "SG15-1.1", "SG15-1.2", "PD03", "PD04", "031US", "100US", "031PAW", "100PAW", "GBRV-1.4"]
    zval = [13, 13, 13, 3, 3, 13, 13, 13, 13, 13]

    # Create the figure and the first axis
    fig, ax1 = plt.subplots(nrows=len(pseudopots), ncols=1, sharex=True, sharey=True, figsize=(10, 10))

    color_factor = 0.95
    color1 = np.array([142, 207, 201])/255 * color_factor
    color2 = np.array([255, 190, 122])/255 * color_factor
    color3 = np.array([130, 176, 210])/255 * color_factor

    rcut_xlims = [5.5, 15]
    marker_nao = 'o'
    marker_pw = '*'

    for i in range(len(pseudopots)):
        # Generate the data
        x1 = ecutwfc
        y1 = y*np.random.random(len(ecutwfc))
        # Plot the first data set on the first axis wit linestyle - and symbol triangle up
        line_pw = ax1[i].plot(x1, y1, marker = marker_pw, ms = 10, linestyle = '-', color = (color1[0], color1[1], color1[2]), label = 'pw')
        #ax1[i].set_xlabel('ecutwfc (Ry)')
        # name this line as pw
        # disable xticks for the first axis
        if i != len(pseudopots)-1:
            ax1[i].set_xticks([])
        # Create the second axis and plot the second data set on it
        x2 = rcut
        y2 = y_for_rcut*np.random.random(len(rcut))
        ax2 = ax1[i].twiny()
        line_dzp = ax2.plot(x2, y2, marker = marker_nao, linestyle = '-', color = (color2[0], color2[1], color2[2]), label = 'DZP')
        ax2.set_xlim(rcut_xlims[0], rcut_xlims[1])
        #ax2.set_xlabel('Rcut (a.u.)')
        # disable xticks for the second axis
        if i == 0:
            # add xlabel as "cutoff radius (a.u.)"
            ax2.set_xlabel('cutoff radius (a.u.)', fontsize=15)
        if i != 0:
            ax2.set_xticks([])
        x3 = rcut
        y3 = y_for_rcut*np.random.random(len(rcut))*0.95
        ax3 = ax1[i].twiny()
        # use rgb color deep green
        line_tzdp = ax3.plot(x3, y3, marker = marker_nao, linestyle = '-', color = (color3[0], color3[1], color3[2]), label = 'TZDP')
        # set xlim for ax3
        ax3.set_xlim(rcut_xlims[0], rcut_xlims[1])
        #ax3.set_xlabel('Rcut (a.u.)')
        # disable xticks for the second axis
        ax3.set_xticks([])
        # add the subplot title at bottom-right corner of each subplot
        ax1[i].text(0.98, 0.05, pseudopots[i] + "\nZval = " + str(zval[i]), horizontalalignment='right', verticalalignment='bottom', transform=ax1[i].transAxes, fontsize=10)
        # set ylim for each subplot
        ax1[i].set_ylim(-0.6, 0.1)
        # add black dashed line at y=0
        ax1[i].axhline(y=0, color='k', linestyle='--', alpha = 0.5)
        # add vertical grid lines, semi-transparent
        ax1[i].grid(axis='x', alpha=0.5, linestyle='-.')
        # only add legend for the last subplot
        if i == len(pseudopots)-1:
            lines = line_pw + line_dzp + line_tzdp
            labels = [l.get_label() for l in lines]
            ax1[i].legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, -0.65), ncol=3, fontsize=15)

    # only add xticks on the last subplot
    ax1[len(pseudopots)-1].set_xticks(ecutwfc)
    # only add xlabel on the last subplot
    ax1[len(pseudopots)-1].set_xlabel('ecutwfc (Ry)', fontsize=15)

    # set all font style to be Arial
    plt.rcParams["font.family"] = "Arial"
    # adjust size of subplots
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.1, hspace=0.0)
    # Add title at left of the whole figure, set direction as vertical, centered at the whole figure
    plt.suptitle('cell parameter error %', fontsize=15, x=0.05, y=0.5, rotation=90, va='center')
    #plt.xlabel('ecutwfc (Ry)')
    plt.ylabel('cell parameter error %')
    # save the figure as png
    plt.savefig('protocol_1dplot_xxy.png', dpi=300)
    # Show the plot
    plt.show()