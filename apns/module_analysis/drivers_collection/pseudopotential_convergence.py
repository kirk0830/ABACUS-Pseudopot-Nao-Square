import apns.module_analysis.result_import as amari
import apns.module_analysis.module_render.pseudopotential_html_generator as amaphg
import numpy as np

if __name__ == "__main__":
    
    root_path = "../11548850/"
    axis_label = "ecutwfc"
    labels = ["energy", "pressure"]
    
    result = amari.collect_result_by_test(path_to_work=root_path, labels=labels)
    result = amari.sort(result=result, axis_label=axis_label)
    result_system = amari.distribute_result_to_system(result=result, labels=labels, to_markdown=True)

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
            ecutwfc = result_system[system][pseudopotential][axis_label]
            # two convergence criteria
            energies = result_system[system][pseudopotential]["energy"]
            energies = [energy - energies[-1] for energy in energies]
            pressures = result_system[system][pseudopotential]["pressure"]
            pressures = [pressure - pressures[-1] for pressure in pressures]

            #energies /= np.abs(np.max(energies))
            energy_axis = axs[i]
            pressure_axis = axs[i].twinx()
            
            # plot
            eline = energy_axis.plot(ecutwfc, energies, label="Energy per atom",
                             linestyle="-", marker="o", color="#2B316B", alpha=0.8)
            pline = pressure_axis.plot(ecutwfc, pressures, label="Pressure",
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
            energy_axis.text(0.995, 0.90, amari.wash_pseudopot_name(pseudopotential), 
                            fontsize=14, backgroundcolor="#FFFFFF", 
                            horizontalalignment='right', verticalalignment='top', transform=energy_axis.transAxes)
            # resize subplots
            energy_axis.figure.subplots_adjust(hspace=0.25, wspace=0.5, left=0.1, right=0.9, top=0.95, bottom=0.10)
            # disable xtitle except the last one
            if i != len(result_system[system]) - 1:
                energy_axis.set_xticklabels([])
            else:
                lines = eline + pline
                labels = [l.get_label() for l in lines]
                axs[i].legend(lines, labels, loc='upper center', bbox_to_anchor=(0.95, -0.3), ncol=1, fontsize=10, fancybox=True, shadow=True)

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
        annotation += "* The line with red circle at the end indicates the pseudopotential does not reach convergence till the last ecutwfc (200 or 300 Ry)."
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
            f.writelines(amaphg.generate_result_page(element=system))
 