"""this can also be used as a temporary platform for running tasks.
The only one important thing is to avoid hard-code as much as possible."""

import re
def input_parse(finp: str):
    """parse the ABACUS input file"""
    keyvalue_pattern = r"^([a-zA-Z0-9_]+)\s*([^#]+)(\s*)(#.*)?$"
    with open(finp, "r") as f:
        lines = f.readlines()
    parsed = {}
    for line in lines:
        match = re.match(keyvalue_pattern, line)
        if match:
            key = match.group(1)
            value = match.group(2)
            parsed[key] = value
    return parsed

import os
def postprocess_bysearch(search_domain: str,
                         filter_fn,
                         action_fn,
                         **kwargs):
    """postprocess output by search, once a folder meets the requirement defined by filter_fn,
    do action_fn on it"""
    return_values_collect = []
    for root, dirs, files in os.walk(search_domain):
        if filter_fn(root, dirs, files, **kwargs):
            return_values_collect.append(action_fn(root, dirs, files, **kwargs))
    return return_values_collect
        
def get_atomic_species(fstru: str):
    """get the first atomic species in the STRU file. This is used to determine the element of the system."""
    with open(fstru, "r") as f:
        lines = f.readlines()
    
    il = 0
    for line in lines:
        il += 1
        line = line.strip()
        if line.startswith("ATOMIC_SPECIES"):
            break
    # it is possible for symbol having a number at the end, like "H1", therefore we need to remove the number
    raw_symbol = lines[il].split()[0]
    symbol = ""
    for char in raw_symbol:
        if char.isalpha():
            symbol += char
    return symbol

import apns.module_workflow.identifier as amwi
import json
def is_ecutwfc_converged(root, dirs, files):
    """check if the ecutwfc is converged"""
    fecutwfc = amwi.TEMPORARY_FOLDER + "/ecutwfc_convergence.json"
    with open(fecutwfc, "r") as f:
        data = json.load(f)
    
    path_backup = os.path.abspath(os.getcwd())
    if "running_scf.log" in files:
        folder_layers = root.replace("\\", "/").split("/")
        
        os.chdir(root)
        os.chdir("..")
        if os.path.exists("STRU"):
            element = get_atomic_species("STRU")
        else:
            os.chdir(path_backup)
            return False
        pspot_id = ""
        for layer in folder_layers:
            if element in layer:
                pspot_id = layer.split("_")[2]
                break
        if pspot_id == "":
            os.chdir(path_backup)
            return False
        if os.path.exists("INPUT"):
            ecutwfc = float(input_parse("INPUT")["ecutwfc"])
            ecutwfc0 = data[element][pspot_id]
            print(f"""Parse folder {root}...
element: {element}
ecutwfc: {ecutwfc} Ry
converged value in record: {ecutwfc0} Ry""")
            if abs(ecutwfc - ecutwfc0) < 1e-6: # well-defined equal for float number

                os.chdir(path_backup)
                return True
        else:
            os.chdir(path_backup)
            return False
    #print(f"present folder does not meet the requirement: {root}")
    os.chdir(path_backup)
    return False

import apns.module_analysis.postprocess.dos_integrator as amados
import apns.module_analysis.__deprecated__.abacus as amamga
def calculate_dos(root, dirs, files):
    path_backup = os.path.abspath(os.getcwd())
    folder = root.replace("\\", "/")
    print("Calculating DOS for " + folder)
    os.chdir(folder)
    system = "running_scf.log"
    if not system in files:
        raise FileNotFoundError("running_scf.log not found in " + root)
    band_energies, nelec, nband, efermi = amamga.grep_band(system)
    window_width = 5 # unit eV
    emin = efermi - window_width/2
    emax = emin + window_width
    band_gap, band_energies, homo, lumo = amados.cal_bandgap(band_energies, nelec, efermi)
    e, dos = amados.smooth_1d(source=band_energies, 
                              lb=emin, 
                              ub=emax, 
                              dx=0.01, 
                              smear=0.05)
    os.chdir(path_backup)
    pseudopot = root.replace("\\", "/").split("/")[-3].split("_")[2]
    element = root.replace("\\", "/").split("/")[-3].split("_")[0]
    return {
        "element": element,
        "pspot_id": pseudopot,
        "x": e.flatten(),
        "y": dos.flatten(),
        "xmin": emin,
        "xmax": emax,
        "homo": homo,
        "lumo": lumo,
        "efermi": efermi
    }

def draw_stack_dos(search_domain: str):
    result = postprocess_bysearch(search_domain=search_domain,
                                  filter_fn=is_ecutwfc_converged,
                                  action_fn=calculate_dos)
    result_elementwise = {}
    for r in result:
        element = r["element"]
        if not element in result_elementwise.keys():
            result_elementwise[element] = []
        result_elementwise[element].append(r)
    for e in result_elementwise.keys():
        result = result_elementwise[e]
        npspot = len(result)
        import matplotlib.pyplot as plt
        import apns.module_analysis.result_import as amari
        # create subplots, n rows, 1 column
        fig, axs = plt.subplots(npspot, 1, figsize=(20, 10))
        plt.rcParams["font.family"] = "Arial"
        for ipspot in range(npspot):
            axs[ipspot].plot(result[ipspot]["x"], result[ipspot]["y"])
            axs[ipspot].set_xlim(result[ipspot]["xmin"], result[ipspot]["xmax"])
            axs[ipspot].set_ylim(0, 1.2*max(result[ipspot]["y"]))
            axs[ipspot].fill_between(result[ipspot]["x"], result[ipspot]["y"], y2 = -1,
                                    where=(
                                        (result[ipspot]["x"] > result[ipspot]["xmin"])&(result[ipspot]["x"] < result[ipspot]["efermi"])),
                                    alpha=0.5)
            axs[ipspot].text(0.995, 0.90, 
                            amari.wash_pseudopot_name(result[ipspot]["pspot_id"]),
                            fontsize=14, backgroundcolor="#FFFFFF",
                            horizontalalignment="right",
                            verticalalignment="top",
                            transform=axs[ipspot].transAxes)
            axs[ipspot].figure.subplots_adjust(hspace=0, wspace=0.5, left=0.1, right=0.9, top=0.95, bottom=0.10)
            if ipspot != npspot - 1:
                axs[ipspot].set_xticklabels([])
            else:
                axs[ipspot].set_xlabel("Energy (eV)", fontsize = 16)
            axs[ipspot].set_yticklabels([])
        fig.text(0.05, 0.5, "Density of States (DOS)", va="center", rotation="vertical", fontsize = 16)
        fig.suptitle(f"Converged DOS of element {e} single crystal calculated with different pseudopotentials", fontsize = 16)
        # add a transparent box contain two squares, the first is filled and the second in the second line is not filled
        # instruction follows as occupied and unoccupied states, respectively
        fig.add_artist(plt.Rectangle((0.05, 0.01), 0.05, 0.05, color="blue", alpha=0.5))
        plt.savefig(f"{e}_dos.svg")
        plt.savefig(f"{e}_dos.png")

if __name__ == "__main__":
    draw_stack_dos("../11588012")