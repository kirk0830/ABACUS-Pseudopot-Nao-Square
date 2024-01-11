import os
import result_generation as rg
import re

def wash_pseudopot_name(name: str):

    result = name.upper()
    result = result.replace("FR", " (Full Relativistic)")
    result = result.replace("SPD", "-spd").replace("-HIGH", "-high")

    # if has four consecutive number, assume the last two as version number
    _match = re.match(r".*([0-9]{4}).*", result)
    if _match:
        result = result.replace(_match.group(1), _match.group(1)[:-2] + " v" + ".".join(_match.group(1)[-2:]))
    return result


def collect_result_by_test(path_to_work: str = "./"):
    """this function will scan all the folders under current folder, and collect the result of each test"""

    test_pattern = r"([\w]+_[0-9]+_[\w\-\.]+)((_([A-Z][0-9])?)?)(_[\w]+)(.*)"
    """regulation: only use / to split path, not \\"""
    result = {}
    for root, folders, files in list(os.walk(path_to_work)):
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
            energy = 0.0
            pressure = 0.0
            with open(root + "/" + "running_scf.log", "r") as f:
                while not line.startswith("!"):
                    line = f.readline().strip()
                    if line.startswith("TIME STATISTICS"):
                        print("warning: scf calculation not converged in " + root)
                        break
                    if line.startswith("TOTAL ATOM NUMBER"):
                        natom = int(line.split()[-1])
                    if line.startswith("TOTAL-PRESSURE"):
                        pressure = float(line.split()[-2]) # in kbar
                energy = float(line.split()[-2])
            with open(root + "/" + "INPUT", "r") as f:
                while not line.startswith("ecutwfc"):
                    line = f.readline().strip()
                ecutwfc = float(line.split()[1])

            layers = root.split("/")
            test = ""
            for layer in layers:
                _match = re.match(test_pattern, layer)
                if _match:
                    test = _match.group(1)
                    break
            if test not in result:
                result[test] = {
                    "ecutwfc": [],
                    "energy": [],
                    "pressure": [],
                }
            result[test]["ecutwfc"].append(ecutwfc)
            result[test]["energy"].append(energy)
            result[test]["pressure"].append(pressure)
            result[test]["natom"] = natom
    return result

import numpy as np
def sort_by_ecutwfc(result):
    """sort the result by ecutwfc"""
    for test in result:
        ecutwfc = result[test]["ecutwfc"]
        energy = result[test]["energy"]
        pressure = result[test]["pressure"]
        idx = np.argsort(ecutwfc)
        ecutwfc = np.array(ecutwfc)[idx]
        energy = np.array(energy)[idx]
        pressure = np.array(pressure)[idx]

        result[test]["ecutwfc"] = ecutwfc
        result[test]["energy"] = energy
        result[test]["pressure"] = pressure
    return result

def distribute_result_to_system(result):
    """distribute the result to different systems"""
    result_system = {}
    for test in result.keys():
        system = test.split("_")[0]
        pseudopotential = test.split("_")[-1]
        if system not in result_system:
            result_system[system] = {}
        if pseudopotential not in result_system[system]:
            result_system[system][pseudopotential] = {
                "ecutwfc": [],
                "energy": [],
                "pressure": []
            }
        result_system[system][pseudopotential]["ecutwfc"] = result[test]["ecutwfc"]
        result_system[system][pseudopotential]["energy"] = result[test]["energy"]
        result_system[system][pseudopotential]["pressure"] = result[test]["pressure"]
        result_system[system][pseudopotential]["natom"] = result[test]["natom"]
    return result_system

def export_dat(result_system):
    """export the result to dat file"""
    for system in result_system:
        with open(system + ".dat", "w") as f:
            for pseudopotential in result_system[system]:
                f.write("# pseudopotential: " + pseudopotential + "\n")
                f.write("# ecutwfc (Ry) | energy (Ry) | pressure (kbar) \n")
                for i in range(len(result_system[system][pseudopotential]["ecutwfc"])):
                    line = str(result_system[system][pseudopotential]["ecutwfc"][i])
                    line += " " + str(result_system[system][pseudopotential]["energy"][i])
                    line += " " + str(result_system[system][pseudopotential]["pressure"][i])
                    f.write(line + "\n")
                f.write("\n")

if __name__ == "__main__":
    
    result = collect_result_by_test("../11549318/")

    result = sort_by_ecutwfc(result)
    result_system = distribute_result_to_system(result)
    export_dat(result_system)

    import matplotlib.pyplot as plt
    for system in result_system:
        # subplot for each pseudopotential
        npseudo = len(result_system[system])
        if npseudo <= 1:
            print("warning: not enough result of pseudopotential for " + system)
            continue
        fig, axs = plt.subplots(npseudo, 1, figsize=(20, 10))
        # set font as Arial
        plt.rcParams["font.family"] = "Arial"
        for i, pseudopotential in enumerate(result_system[system]):

            natom = result_system[system][pseudopotential]["natom"]
            ecutwfc = result_system[system][pseudopotential]["ecutwfc"]
            # two convergence criteria
            energies = result_system[system][pseudopotential]["energy"]
            energies = [energy - energies[-1] for energy in energies]
            pressures = result_system[system][pseudopotential]["pressure"]
            pressures = [pressure - pressures[-1] for pressure in pressures]

            #energies /= np.abs(np.max(energies))
            energy_axis = axs[i]
            pressure_axis = axs[i].twinx()
            
            # plot
            energy_axis.plot(ecutwfc, energies, label="Energy per atom",
                             linestyle="-", marker="o", color="#2B316B", alpha=0.8)
            pressure_axis.plot(ecutwfc, pressures, label="Pressure",
                               linestyle="--", marker="s", color="#24B5A5", alpha=0.8)
            #energy_axis.set_xlabel("ecutwfc (Ry)")
            #energy_axis.set_ylabel("energy (Ry)")
            
            yaxis_stretch = 1.5
            # set 3 y ticks
            single_side_error = np.abs(np.max(energies) - np.min(energies))
            ylimmin_e = -yaxis_stretch * single_side_error
            ylimmax_e = -ylimmin_e
            yticks_e = np.linspace(-single_side_error, single_side_error, 3)

            deltay_e = (ylimmax_e - ylimmin_e) / 3

            single_side_error = np.abs(np.max(pressures) - np.min(pressures))
            ylimmin_p = -yaxis_stretch * single_side_error
            ylimmax_p = -ylimmin_p
            yticks_p = np.linspace(-single_side_error, single_side_error, 3)

            # allow different y scale

            energy_axis.set_ylim(ylimmin_e, ylimmax_e)
            pressure_axis.set_ylim(ylimmin_p, ylimmax_p)

            # disable scientific notation for y ticks, and set precision to 4
            energy_axis.yaxis.set_major_formatter('{:.4f}'.format)
            pressure_axis.yaxis.set_major_formatter('{:.4f}'.format)
            # set y ticks
            energy_axis.set_yticks(yticks_e)
            pressure_axis.set_yticks(yticks_p)

            # red-circle the point indiced by i_highlight

            i_highlight_e_high = -1
            i_highlight_p_high = -1
            for j in range(len(energies)):
                delta_e = np.abs(energies[j] - energies[-1])/natom
                if delta_e < 1e-3: # 1meV
                    i_highlight_e_high = j
                    break
            for j in range(len(pressures)):
                delta_p = np.abs(pressures[j] - pressures[-1])
                if delta_p < 0.1:
                    i_highlight_p_high = j
                    break
            if i_highlight_e_high != -1 and i_highlight_p_high != -1:
                i_highlight = max(i_highlight_e_high, i_highlight_p_high)
                energy_axis.plot(ecutwfc[i_highlight], 
                                energies[i_highlight], 
                                "o", markersize=15, markerfacecolor="none", markeredgecolor="#D8006A", markeredgewidth=2,
                                zorder=10)
                # add a text to indicate the point
                # energy_axis.text(ecutwfc[i_highlight], 
                #                 energies[i_highlight] - deltay_e * 1.5, 
                #                 str(ecutwfc[i_highlight])+" Ry, < 1 meV, < 0.1 kbar", 
                #                 fontsize=12,
                #                 horizontalalignment='left', verticalalignment='bottom',
                #                 backgroundcolor="#FFFFFF",
                #                 zorder=9)

            # add vertical and horizontal grid lines
            energy_axis.grid(True, color="#c6c6c6")
            #pressure_axis.grid(True, color="#c6c6c6")
            # at inside right top corner, add title
            energy_axis.text(0.995, 0.90, wash_pseudopot_name(pseudopotential), 
                            fontsize=14, backgroundcolor="#FFFFFF", 
                            horizontalalignment='right', verticalalignment='top', transform=energy_axis.transAxes)
            # resize subplots
            energy_axis.figure.subplots_adjust(hspace=0.25, wspace=0.5, left=0.1, right=0.9, top=0.95, bottom=0.10)
            # disable xtitle except the last one
            if i != len(result_system[system]) - 1:
                energy_axis.set_xticklabels([])

        # add title
        fig.suptitle(system+" basic energy convergence test", fontsize = 16)
        # add one x title
        fig.text(0.5, 0.04, "ecutwfc (Ry)", ha='center', fontsize = 16)
        # add one y title
        fig.text(0.04, 0.5, "Relative electronic energy per atom (eV)", va='center', rotation='vertical', fontsize = 16)
        # add one y title
        fig.text(0.96, 0.5, "Relative pressure (kbar)", va='center', rotation='vertical', fontsize = 16)

        # add annotation
        annotation_title = "APNS Pseudopotential test note"
        fig.text(0.01, 0.05, annotation_title,
                 fontstyle="italic", fontsize=8, backgroundcolor="#FFFFFF",
                 horizontalalignment='left', verticalalignment='bottom')
        annotation =  "* The red circle indicates the point with energy difference per atom less than 1 meV and pressure difference less than 0.1 kbar.\n"
        annotation += "* The line with red circle at the end indicates the pseudopotential does not reach convergence at < 200 Ry ecutwfc."
        # add a box to show annotation, with black edge
        fig.text(0.01, 0.02, annotation, 
                 fontsize=8, backgroundcolor="#FFFFFF", 
                 horizontalalignment='left', verticalalignment='bottom')

        # save svg
        plt.savefig(system + ".svg")
        plt.savefig(system + ".png")
        plt.close()

        # generate html
        fmarkdown = system + ".md"
        with open(fmarkdown, "w") as f:
            f.writelines(rg.generate_result_page(element=system))
            