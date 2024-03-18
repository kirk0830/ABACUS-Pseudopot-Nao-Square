"""The design of functions in this file is, irrelevant to the specific software used

In practice, we use qespresso or abacus"""

import os
import re
def wash_pseudopot_name(name: str):
    """wash the pseudopotential name, to make it more readable
    Is it general enough? yes"""
    result = name.upper()
    result = result.replace("FR", " (Full Relativistic)")
    if result.startswith("PD04"):
        result = result.replace("PD04", "").lower()
        result = "PD04" + "-" + result.strip()
    # if has four consecutive number, assume the last two as version number
    _match = re.match(r".*([0-9]{4}).*", result)
    if _match:
        result = result.replace(_match.group(1), _match.group(1)[:-2] + " v" + ".".join(_match.group(1)[-2:]))
    return result

def initialize_test_result(labels: list):
    """return a dictionary with all labels initialized to 0
    Is it general enough? yes
    
    Usage: once extend the labels, just change the _datatypes"""
    _datatypes = {
        "energy": float, "pressure": float, "natom": int
    }
    result = {}
    for label in labels:
        result[label] = _datatypes[label](0)
    return result

import apns.module_analysis.__deprecated__.abacus as amaa
import apns.module_analysis.__deprecated__.qespresso as amaq

def initialize_greps(software: str = "abacus", labels: list = None):
    """gather grep methods for specific software
    Is it general enough? yes
    
    Implement software-by-software
    
    Return:
        grep_funcs (dict): a dictionary with all grep functions
        grep_ecutwfc (function): the grep function for ecutwfc
        grep_natom (function): the grep function for natom
    """
    if labels is None:
        labels = []

    grep_funcs = {}
    if software == "abacus":
        _greps = amaa.return_all_greps()
        for label in labels:
            if label != "ecutwfc" and label != "natom":
                grep_funcs[label] = _greps[label]
    elif software == "qespresso":
        raise NotImplementedError("module_analysis for qespresso is not implemented yet")
        for label in labels:
            pass
    else:
        raise NotImplementedError("module_analysis for " + software + " is not implemented yet")
    
    return grep_funcs, _greps["ecutwfc"], _greps["natom"]

def fname_setting(software: str = "abacus"):
    """return the file names (stdout, log and input script) for specific software
    Is it general enough? yes
    
    Implement software-by-software
    
    Return:
        stdout (str): the name of stdout file
        log (str): the name of log file
        input_script (str): the name of input script file
    """
    if software == "abacus":
        return amaa.fname_setting
    elif software == "qespresso":
        return amaq.fname_setting
    else:
        raise NotImplementedError("module_analysis for " + software + " is not implemented yet")

def is_normal_end(software: str = "abacus"):
    """check if the job is normally ended for specific software
    Is it general enough? yes
    
    Implement software-by-software
    
    Return:
        is_normal_end (function): the function to check if the job is normally ended
    """
    if software == "abacus":
        return amaa.is_normal_end
    elif software == "qespresso":
        return amaq.is_normal_end
    else:
        raise NotImplementedError("module_analysis for " + software + " is not implemented yet")

def collect_result_by_test(path_to_work: str = "./", 
                           software: str = "abacus", 
                           calculation: str = "scf",
                           suffix: str = "ABACUS",
                           labels: list = ["energy"]):
    """new version of this function will scan all the folders under current folder, and collect the result of each test"""

    test_pattern = r"^(([A-Za-z]+)_([0-9]+)_([A-Za-z0-9\-\+\.]+))_(.*)$"
    """regulation: only use / to split path, not \\"""
    
    fname_setting_fn = fname_setting(software = software)
    stdout, log, input_script = fname_setting_fn(job_dir = path_to_work, calculation = calculation, suffix = suffix, include_path = False)
    is_normal_end_fn = is_normal_end(software = software)
    result = {}
    for root, folders, files in list(os.walk(path_to_work)):
        root = root.replace("\\", "/")
        parent_folder = "/".join(root.split("/")[:-1])
        if log in files: # this is a test folder
            if not is_normal_end_fn(root + "/" + log):
                continue
        else:
            continue

        test_result = initialize_test_result(labels)
        grep_funcs, grep_ecutwfc, grep_natom = initialize_greps(software = software, labels = labels)

        _result = {}
        for label in labels:
            _result[label] = grep_funcs[label](root + "/" + log) # exception handle here, if convergence is not reached, will return False
            if _result is False:
                continue
        
        test_result.update(_result)
        test_result["ecutwfc"] = grep_ecutwfc(root + "/" + input_script)
        test_result["natom"] = grep_natom(root + "/" + log)

        layers = root.split("/")
        test = ""
        for layer in layers:
            _match = re.match(test_pattern, layer)
            if _match:
                test = _match.group(1)
                break

        if test not in result:
            result[test] = {label: [] for label in labels}
            result[test]["ecutwfc"] = []

        for label in labels:
            result[test][label].append(test_result[label])

        result[test]["ecutwfc"].append(test_result["ecutwfc"])
        result[test]["natom"] = test_result["natom"]

    return result

import numpy as np
def sort(result, axis_label: str = "ecutwfc", scalar_labels: list = None):
    """because tests are not guaranteed to arrange with ascending or discending order of something,
    here sort the result by ecutwfc by default
    
    Args:
        result (dict): the result dictionary, saves the result of each test, say the identifier of each test is test
        
    Development:
        add more in compulsory_scalar_labels, to make it more general
    """
    compulsory_scalar_labels = ["natom"]
    if scalar_labels is None:
        scalar_labels = []
    scalar_labels += compulsory_scalar_labels
    for test in result:
        axis = result[test][axis_label]
        indices = np.argsort(axis)
        for property in result[test]:
            if property == axis_label or property in scalar_labels:
                continue
            result[test][property] = np.array(result[test][property])[indices]
        result[test][axis_label] = np.array(axis)[indices]
    return result

def distribute_result_to_system(result: dict, labels: list = ["energy"], to_markdown: bool = True):
    """new version of distribute result to system, which is more general.
    The returned dictionary is organized as:
    {'Lu': {'dojo04': {'energy': array([-18386.35796251, -18730.22198428, -19367.69804648, -19842.61019759,
        -20046.86643719, -20109.08328875, -20120.35067214, -20121.44724116,
        -20121.48064955, -20121.48366003, -20121.48508414, -20121.48594797]), 'pressure': array([-6.59920873e+03, -1.69629178e+04, -3.15224246e+04, -1.91310207e+04,
        -8.24061514e+03, -2.16305993e+03, -3.22985857e+02, -3.42355710e+01,
        -1.78385600e+01, -1.77731700e+01, -1.78898750e+01, -1.77732310e+01]), 'ecutwfc': array([ 20.,  30.,  40.,  50.,  60.,  70.,  80.,  90., 100., 150., 200.,
        300.]), 'natom': 3}, 'dojo043plus': {'energy': array([-3856.12899366, -3859.90167784, -3861.57512149, -3861.864847  ,
        -3861.88331607, -3861.88388465, -3861.88416508, -3861.88425301,
        -3861.88432591, -3861.88445787, -3861.88451156, -3861.88452629]), 'pressure': array([-115.259069, -128.672841,  -55.574045,  -25.721725,  -21.978221,
            -21.930457,  -21.903452,  -21.901915,  -21.897788,  -21.89699 ,
            -21.896496,  -21.89596 ]), 'ecutwfc': array([ 20.,  30.,  40.,  50.,  60.,  70.,  80.,  90., 100., 150., 200.,
        300.]), 'natom': 3}, 
         ... 'pd04sp': {'energy': array([-3779.75031704, -3806.66997801, -3815.64010498, -3817.43002042,
        -3817.62851742, -3817.63653153, -3817.63750554, -3817.6378734 ,
        -3817.63821961, -3817.63956157, -3817.64010473, -3817.64022908]), 'pressure': array([-924.415683, -604.063875, -207.454182,  -52.876959,  -22.364928,
            -20.185806,  -20.190214,  -20.147959,  -20.180524,  -20.188209,
            -20.133733,  -20.127488]), 'ecutwfc': array([ 20.,  30.,  40.,  50.,  60.,  70.,  80.,  90., 100., 150., 200.,
        300.]), 'natom': 3}}}
    """
    result_system = {}
    for test in result.keys():
        system = test.split("_")[0]
        pseudopotential = test.split("_")[-1]
        if system not in result_system:
            result_system[system] = {}
        if pseudopotential not in result_system[system]:
            result_system[system][pseudopotential] = {
                label: [] for label in labels
            }
            result_system[system][pseudopotential].update({
                "ecutwfc": [],
                "natom": []
            })
        for label in labels:
            result_system[system][pseudopotential][label] = result[test][label]
        result_system[system][pseudopotential].update({
            "ecutwfc": result[test]["ecutwfc"],
            "natom": result[test]["natom"]
        })
    if to_markdown:
        export_data_as_markdown(result_system, labels)
    return result_system

def export_data_as_markdown(result_system: dict, labels: list = ["energy"]):
    """new version of export_dat, which is more general"""
    for system in result_system:
        with open(system + "_datatable.md", "w") as f:
            for pseudopotential in result_system[system]:
                f.write("### pseudopotential: " + pseudopotential + "\n")
                f.write("| ecutwfc (Ry) | ")
                for label in labels:
                    f.write(label + " | ")
                f.write("\n")
                f.write("| --- | ")
                for label in labels:
                    f.write("--- | ")
                f.write("\n")
                for i in range(len(result_system[system][pseudopotential]["ecutwfc"])):
                    line = "| " + str(result_system[system][pseudopotential]["ecutwfc"][i]) + " |"
                    for label in labels:
                        line += " %.4f|"%result_system[system][pseudopotential][label][i]
                    f.write(line + "\n")
                f.write("\n")

"""the following function is not complete, and not used yet"""
import matplotlib.pyplot as plt
def draw_stack_plot(result_of_system: dict, 
                    ndim_stack: int, 
                    axis_label: str = "ecutwfc", 
                    labels: list = ["energy"],
                    conv_thrs: list = [1e-3, 0.1],
                    supertitle: str = "Basic energy convergence test"):
    """draw stack plot for each system"""
    color_board = ["#2B316B", "#24B5A5"]
    linestyles = ["-", "--"]
    markers = ["o", "s"]
    units = {"ecutwfc": "Ry", "energy": "Ry", "pressure": "kbar"}
    if ndim_stack < 1:
        print("warning: ndim_stack should be at least 1")
        return
    fig, axs = plt.subplots(ndim_stack, 1, figsize=(20, 10))
    plt.rcParams["font.family"] = "Arial"
    for i, key in enumerate(list(result_of_system.keys())): # the exact identifier to stack, for pseudopotential, it is the name of pseudopotential
        natom = result_of_system[key]["natom"]
        x = result_of_system[key][axis_label]
        ys = []
        
        for label in labels:
            ys.append(result_of_system[key][label])
            # convert to relative error respect to the last one
            ys[-1] = [(ys[-1][j] - ys[-1][-1])/natom for j in range(len(ys[-1]))]
        
        yaxis_stretch = 1.5
        lines = []
        for ilbl, label in enumerate(labels):
            color = color_board[ilbl%len(color_board)]
            linestyle = linestyles[ilbl%len(linestyles)]
            marker = markers[ilbl%len(markers)]
            label = label + "/natom"
            if ilbl == 0:
                lines.append(
                    axs[i].plot(x, ys[ilbl], label=label, # data
                                linestyle=linestyle, marker=marker, color=color, alpha=0.8 # style
                                )
                )
                single_side_error = np.abs(np.max(ys[ilbl]) - np.min(ys[ilbl]))
                ylimmin = -yaxis_stretch * single_side_error
                ylimmax = -ylimmin
                yticks = np.linspace(-single_side_error, single_side_error, 3)

            else:
                lines.append(
                    axs[i].twinx().plot(x, ys[ilbl], label=label, # data
                                        linestyle=linestyle, marker=marker, color=color, alpha=0.8 # style
                                        )
                )
                single_side_error = np.abs(np.max(ys[ilbl]) - np.min(ys[ilbl]))
                ylimmin = -yaxis_stretch * single_side_error
                ylimmax = -ylimmin
                yticks = np.linspace(-single_side_error, single_side_error, 3)

            lines[-1].axes.set_ylim(ylimmin, ylimmax)
            lines[-1].axes.yaxis.set_major_formatter('{:.4f}'.format)
            lines[-1].axes.set_yticks(yticks)
            lines[-1].axes.grid(True, color="#c6c6c6")
            lines[-1].axes.text(0.995, 0.90, wash_pseudopot_name(key),
                                fontsize=14, backgroundcolor="#FFFFFF",
                                horizontalalignment='right', verticalalignment='top', transform=lines[-1].axes.transAxes)
            lines[-1].axes.figure.subplots_adjust(hspace=0.25, wspace=0.5, left=0.1, right=0.9, top=0.95, bottom=0.10)
            if i != len(result_of_system) - 1:
                lines[-1].axes.set_xticklabels([])
        
        fig.suptitle(supertitle, fontsize = 16)
        fig.text(0.5, 0.04, axis_label + "("+units[axis_label]+")", ha='center', fontsize = 16)

