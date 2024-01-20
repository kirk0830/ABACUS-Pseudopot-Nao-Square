"""
This module is used to generate input_analysis.json

Author: Kirk0830
Date: 2023-11-22

Description:
    This module is used to generate input_analysis.json

Functions:
    generate_dos_section
    generate_plot_band_gap_2d_section
    generate_input_analysis_json
"""
import json

def generate_dos_section(
        draw = False,
        window_width = 4,
        delta_e = 0.01,
        smear = 0.02,
        ):
    """
    generate dos section in input_analysis.json

    Args:

    :param draw: whether to draw dos
    :param window_width: dos window width, in eV
    :param delta_e: dos delta_e, in eV
    :param smear: dos smear, in eV

    Returns:

    :return: dos section in input_analysis.json
    """
    dos_section = {
        "draw": draw,
        "window_width": window_width,
        "delta_e": delta_e,
        "smear": smear,
    }
    return dos_section

def generate_plot_band_gap_2d_section(
        draw = False,
        mode = "imshow",
        color_map = "jet",
        ncol = 2
        ):
    """
    generate plot_band_gap_2d section in input_analysis.json

    Args:

    :param draw: whether to draw band gap 2d
    :param mode: imshow or contourf
    :param color_map: color map
    :param ncol: number of columns in legend

    Returns:

    :return: plot_band_gap_2d section in input_analysis.json
    """
    plot_band_gap_2d_section = {
        "draw": draw,
        "mode": mode,
        "color_map": color_map,
        "ncol": ncol
    }
    return plot_band_gap_2d_section

def generate_input_analysis_json(
        work_folder = "./pseudopotentials/results",
        dos_draw = False,
        dos_window_width = 4,
        dos_delta_e = 0.01,
        dos_smear = 0.02,
        draw_statistics = False,
        plot_band_gap_2d_draw = False,
        plot_band_gap_2d_mode = "imshow",
        plot_band_gap_2d_color_map = "jet",
        plot_band_gap_2d_ncol = 2,
        pseudopotentials = {},
        systems = [],
        experimental_values = {},
        appendix = ""
        ):
    """
    generate input_analysis.json

    Args:

    :param work_folder: work folder
    :param dos_draw: whether to draw dos
    :param dos_window_width: dos window width, in eV
    :param dos_delta_e: dos delta_e, in eV
    :param dos_smear: dos smear, in eV
    :param draw_statistics: whether to draw statistics
    :param plot_band_gap_2d_draw: whether to draw band gap 2d
    :param plot_band_gap_2d_mode: imshow or contourf
    :param plot_band_gap_2d_color_map: color map
    :param plot_band_gap_2d_ncol: number of columns in legend
    :param pseudopotentials: pseudopotentials
    :param systems: systems
    :param experimental_values: experimental values

    Returns:

    :return: input_analysis.json

    """
    input_analysis = {
        "work_folder": work_folder,
        "dos": generate_dos_section(
            draw = dos_draw,
            window_width = dos_window_width,
            delta_e = dos_delta_e,
            smear = dos_smear,
        ),
        "draw_statistics": draw_statistics,
        "plot_band_gap_2d": generate_plot_band_gap_2d_section(
            draw = plot_band_gap_2d_draw,
            mode = plot_band_gap_2d_mode,
            color_map = plot_band_gap_2d_color_map,
            ncol = plot_band_gap_2d_ncol,
        )
    }
    input_analysis["pseudopot_kinds"] = pseudopotentials["kinds"]
    input_analysis["versions"] = pseudopotentials["versions"]
    input_analysis["appendices"] = pseudopotentials["appendices"]

    input_analysis["systems"] = {}
    for system in systems:
        input_analysis["systems"][system] = {}

    for system in systems:
        input_analysis["systems"][system]["experimental_value"] = {}
        if system in experimental_values.keys():
            ev_sys = experimental_values[system]
            for key in ev_sys.keys():
                input_analysis["systems"][system]["experimental_value"][key] = ev_sys[key]
    
    with open("input_analysis_"+appendix+".json", "w") as f:
        json.dump(input_analysis, f, indent = 4)
    
    return input_analysis

if __name__ == "__main__":

    import time
    # timestamp
    ts = time.time()
    # time string
    time_str = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime(ts))

    generate_input_analysis_json(
        work_folder = "./pseudopotentials/results",
        dos_draw = False,
        dos_window_width = 4,
        dos_delta_e = 0.01,
        dos_smear = 0.02,
        draw_statistics = False,
        plot_band_gap_2d_draw = False,
        plot_band_gap_2d_mode = "imshow",
        plot_band_gap_2d_color_map = "jet",
        plot_band_gap_2d_ncol = 2,
        pseudopotentials = {
            "kinds": ["pd", "sg15"],
            "versions": ["04", "10"],
            "appendices": [""]
        },
        systems = ["Li2O", "Li2O2"],
        experimental_values = {
            "Li2O": {
                "band_gap": 5.37,
                "formation_energy": -2.83,
                "total_energy": -14.91,
            },
            "Li2O2": {
                "band_gap": 5.37,
                "formation_energy": -2.83,
                "total_energy": -14.91,
            }
        },
        appendix = time_str
    )