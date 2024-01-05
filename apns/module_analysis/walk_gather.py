import os
import result_generation as rg

"""regulation: only use / to split path, not \\"""
result = {}
for root, folders, files in list(os.walk('./')):
    if "running_scf.log" in files:
        # first check if the calculation terminated abnormally
        normal_termination = False
        root = root.replace("\\", "/")
        parent_folder = "/".join(root.split("/")[:-1])
        with open(parent_folder + "/" + "out.log", "r") as f:
            lines = f.readlines()
            for line in lines:
                if "TIME STATISTICS" in line:
                    normal_termination = True
                    break
        if not normal_termination:
            print("warning: scf calculation not terminated normally in " + root)
            continue
        line = ""
        natom = 0
        with open(root + "/" + "running_scf.log", "r") as f:
            while not line.startswith("!"):
                line = f.readline().strip()
                if line.startswith("TIME STATISTICS"):
                    print("warning: scf calculation not converged in " + root)
                    break
                if line.startswith("TOTAL ATOM NUMBER"):
                    natom = int(line.split()[-1])
            energy = float(line.split()[-2])
        with open(root + "/" + "INPUT", "r") as f:
            while not line.startswith("ecutwfc"):
                line = f.readline().strip()
            ecutwfc = float(line.split()[1])

        layers = root.split("/")
        test = ""
        for layer in layers:
            if layer.startswith("t_"):
                test = layer
        if test not in result:
            result[test] = {
                "ecutwfc": [],
                "energy": []
            }
        result[test]["ecutwfc"].append(ecutwfc)
        result[test]["energy"].append(energy)
        result[test]["natom"] = natom
# sort
import numpy as np
for test in result:
    ecutwfc = np.array(result[test]["ecutwfc"])
    energy = np.array(result[test]["energy"])
    idx = np.argsort(ecutwfc)
    ecutwfc = ecutwfc[idx]
    energy = energy[idx]
    result[test]["ecutwfc"] = ecutwfc
    result[test]["energy"] = energy

# divide system by system
result_system = {}
for test in result:
    system = test.split("_")[2]
    pseudopotential = test.split("_")[-1]
    if system not in result_system:
        result_system[system] = {}
    if pseudopotential not in result_system[system]:
        result_system[system][pseudopotential] = {
            "ecutwfc": [],
            "energy": []
        }
    result_system[system][pseudopotential]["ecutwfc"] = result[test]["ecutwfc"]
    result_system[system][pseudopotential]["energy"] = result[test]["energy"]
    result_system[system][pseudopotential]["natom"] = result[test]["natom"]

# save to file
for system in result_system:
    with open(system + ".dat", "w") as f:
        for pseudopotential in result_system[system]:
            f.write("# pseudopotential: " + pseudopotential + "\n")
            f.write("# ecutwfc (Ry) | energy (Ry)\n")
            for i in range(len(result_system[system][pseudopotential]["ecutwfc"])):
                f.write(str(result_system[system][pseudopotential]["ecutwfc"][i]) + " " + str(result_system[system][pseudopotential]["energy"][i]) + "\n")
            f.write("\n")

import matplotlib.pyplot as plt
for system in result_system:
    # subplot for each pseudopotential
    fig, axs = plt.subplots(len(result_system[system]), 1, figsize=(20, 10))
    # set font as Arial
    plt.rcParams["font.family"] = "Arial"
    for i, pseudopotential in enumerate(result_system[system]):
        natom = result_system[system][pseudopotential]["natom"]
        i_highlight = -1
        ecutwfc = result_system[system][pseudopotential]["ecutwfc"]
        energies = result_system[system][pseudopotential]["energy"]
        for j in range(len(energies)):
            delta_e = np.abs(energies[j] - energies[-1])/natom
            if delta_e < 1e-3: # 1meV
                i_highlight = j
                break
        #energies /= np.abs(np.max(energies))
        axs[i].plot(ecutwfc, energies, "o-", color="#2B316B")
        #axs[i].set_xlabel("ecutwfc (Ry)")
        #axs[i].set_ylabel("energy (Ry)")
        # set 5 y ticks
        ymin = np.min(energies)
        ymax = np.max(energies)
        deltay = np.abs(ymax - ymin)
        ymin = ymin - np.abs(ymax - ymin)
        yticks = np.linspace(ymin, ymax, 3)
        ymax += np.abs(ymax - ymin)*0.1
        ymin -= np.abs(ymax - ymin)*0.1
        axs[i].set_ylim(ymin, ymax)
        # disable scientific notation for y ticks, and set precision to 4
        axs[i].yaxis.set_major_formatter('{:.4f}'.format)
        # set y ticks
        axs[i].set_yticks(yticks)
        # red-circle the point indiced by i_highlight
        if i_highlight != -1:
            axs[i].plot(ecutwfc[i_highlight], energies[i_highlight], "o", markersize=15, markerfacecolor="none", markeredgecolor="red", markeredgewidth=2)
            # add a text to indicate the point
            axs[i].text(ecutwfc[i_highlight], energies[i_highlight], str(ecutwfc[i_highlight])+", < 1 meV", horizontalalignment='left', verticalalignment='bottom')
        # add vertical and horizontal grid lines
        axs[i].grid(True, color="#c6c6c6")
        # at inside right top corner, add title
        axs[i].text(1.00, 0.95, pseudopotential, horizontalalignment='right', verticalalignment='top', transform=axs[i].transAxes)
        # resize subplots
        axs[i].figure.subplots_adjust(hspace=0.25, wspace=0.5, left=0.1, right=0.9, top=0.95, bottom=0.10)
        # disable xtitle except the last one
        if i != len(result_system[system]) - 1:
            axs[i].set_xticklabels([])

    # add title
    fig.suptitle(system+" basic energy convergence test", fontsize = 16)
    # add one x title
    fig.text(0.5, 0.04, "ecutwfc (Ry)", ha='center', fontsize = 16)
    # add one y title
    fig.text(0.04, 0.5, "Electronic energy (eV)", va='center', rotation='vertical', fontsize = 16)
    # save svg
    plt.savefig(system + ".svg")
    # save png
    plt.savefig(system + ".png")
    plt.close()

    # generate html
    fmarkdown = system + ".md"
    with open(fmarkdown, "w") as f:
        f.writelines(rg.generate_result_page(element=system))