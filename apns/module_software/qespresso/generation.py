"""
Quantum ESPRESSO input file generator

Author: Kirk0830
Date: 2023-11-22

Description:
    This file contains functions that can generate Quantum ESPRESSO input file from CIF file.

Useful cases:
    as it is
"""
from apns.module_structure.crystal_information_file import read_1 as read
from apns.module_structure.crystal_information_file import _method1_atomic_position as m1ap
from apns.module_structure.crystal_information_file import _method1_lattice as m1l
import apns.module_database as ai

basic_parameters = {
    "control": {
        "calculation": "scf",
        "restart_mode": "from_scratch",
        "pseudo_dir": "./",
        "outdir": "./",
        "tprnfor": ".true.",
        "tstress": ".true.",
        "wf_collect": ".true.",
        "nstep": 100,
        "verbosity": "high",
    },
    "system": {
        "ibrav": 0,
        "nat": 0,
        "ntyp": 0,
        "ecutwfc": "ecutwfc_to_test",
        "occupations": "smearing",
        "smearing": "gaussian",
        "degauss": 0.02,
        "nspin": 1,
        "starting_magnetization": 0.0,
        "nosym": ".false.",
        "noinv": ".false.",
        "noncolin": ".false.",
        "lspinorb": ".false.",
        "input_dft": "functional_to_test",
    },
    "electrons": {
        "electron_maxstep": 100,
        "conv_thr": 1e-7,
        "mixing_beta": 0.7,
        "mixing_mode": "plain",
        "mixing_ndim": 8,
        "diagonalization": "david"
    },
    "ions": {
        "ion_dynamics": "bfgs",
        "upscale": 100
    },
    "cell": {
        "cell_dynamics": "bfgs",
        "cell_dofree": "all",
        "cell_factor": 2,
        "press_conv_thr": 0.2
    },
    "k_points": {
        "k_points": [1, 1, 1],
        "k_points_shift": [0, 0, 0],
    }
}

def write_control(other_parameters: dict) -> str:
    """
    write control parameters
    """
    control = basic_parameters["control"]
    if "control" in other_parameters.keys():
        for key, value in other_parameters["control"].items():
            control[key] = value
    control_str = "&CONTROL\n"
    for key, value in control.items():
        # if value is string, add quotation marks
        if isinstance(value, str) and not (".true." in value or ".false." in value) and not value.endswith("to_test"):
            control_str += "%s = '%s',\n"%(key, value)
        else:
            control_str += "%s = %s,\n"%(key, str(value))
    control_str += "/\n"
    return control_str

def write_system(other_parameters: dict, atoms: dict) -> str:
    """
    write system parameters
    """
    system = basic_parameters["system"]
    ntyp = len(atoms.keys())
    nat = 0
    for atom in atoms.keys():
        nat += len(atoms[atom]["positions"])
    system["ntyp"] = ntyp
    system["nat"] = nat

    if "system" in other_parameters.keys():
        for key, value in other_parameters["system"].items():
            system[key] = value
    system_str = "&SYSTEM\n"
    for key, value in system.items():
        # if value is string, add quotation marks
        if isinstance(value, str) and not (".true." in value or ".false." in value) and not value.endswith("to_test"):
            system_str += "%s = '%s',\n"%(key, value)
        else:
            system_str += "%s = %s,\n"%(key, str(value))
    system_str += "/\n"
    return system_str

def write_electrons(other_parameters: dict) -> str:
    """
    write electrons parameters
    """
    electrons = basic_parameters["electrons"]
    if "electrons" in other_parameters.keys():
        for key, value in other_parameters["electrons"].items():
            electrons[key] = value
    electrons_str = "&ELECTRONS\n"
    for key, value in electrons.items():
        # if value is string, add quotation marks
        if isinstance(value, str) and not (".true." in value or ".false." in value) and not value.endswith("to_test"):
            electrons_str += "%s = '%s',\n"%(key, value)
        else:
            electrons_str += "%s = %s,\n"%(key, str(value))
    electrons_str += "/\n"
    return electrons_str

def write_ions(other_parameters: dict) -> str:
    """
    write ions parameters
    """
    ions = basic_parameters["ions"]
    if "ions" in other_parameters.keys():
        for key, value in other_parameters["ions"].items():
            ions[key] = value
    ions_str = "&IONS\n"
    for key, value in ions.items():
        # if value is string, add quotation marks
        if isinstance(value, str) and not (".true." in value or ".false." in value) and not value.endswith("to_test"):
            ions_str += "%s = '%s',\n"%(key, value)
        else:
            ions_str += "%s = %s,\n"%(key, str(value))
    ions_str += "/\n"
    return ions_str

def write_cell(other_parameters: dict) -> str:
    """
    write cell parameters
    """
    cell = basic_parameters["cell"]
    if "cell" in other_parameters.keys():
        for key, value in other_parameters["cell"].items():
            cell[key] = value
    cell_str = "&CELL\n"
    for key, value in cell.items():
        # if value is string, add quotation marks
        if isinstance(value, str) and not (".true." in value or ".false." in value) and not value.endswith("to_test"):
            cell_str += "%s = '%s',\n"%(key, value)
        else:
            cell_str += "%s = %s,\n"%(key, str(value))
    cell_str += "/\n"
    return cell_str

def write_atomic_species(atomic_species: dict) -> str:
    """
    write atomic species
    """
    atomic_species_str = "ATOMIC_SPECIES\n"
    for iatom in range(len(atomic_species["elements"])):
        atomic_species_str += atomic_species["elements"][iatom] + " " + str(atomic_species["mass"][iatom]) + " " + atomic_species["pseudopotentials"][iatom] + "\n"
    atomic_species_str += "\n"
    return atomic_species_str

def write_cell_parameters(lattice: dict, mode = "cif") -> str:
    """
    write cell parameters

    mode: "cif" or "nao"
    """
    cell_parameters_str = "CELL_PARAMETERS (angstrom)\n"
    if mode == "cif":
        for i in range(3):
            cell_parameters_str += "%12.8f %12.8f %12.8f\n"%(lattice["lattice_vectors"][i][0], lattice["lattice_vectors"][i][1], lattice["lattice_vectors"][i][2])
        cell_parameters_str += "\n"
    elif mode == "nao":
        lattice_constant = 1.889725989
        default_lattice = [
            [20.0, 0.0, 0.0],
            [0.0, 20.0, 0.0],
            [0.0, 0.0, 20.0],
        ]
        for i in range(3):
            cell_parameters_str += "%12.8f %12.8f %12.8f\n"%(default_lattice[i][0] * lattice_constant, default_lattice[i][1] * lattice_constant, default_lattice[i][2] * lattice_constant)
            
    return cell_parameters_str

def write_atomic_positions(atoms: dict) -> str:
    """
    write atomic positions
    """
    atomic_positions_str = "ATOMIC_POSITIONS (crystal)\n"
    for atom in atoms.keys():
        natom = len(atoms[atom]["positions"])
        for iatom in range(natom):
            atomic_positions_str += "%s %12.8f %12.8f %12.8f\n"%(atom, atoms[atom]["positions"][iatom][0], atoms[atom]["positions"][iatom][1], atoms[atom]["positions"][iatom][2])
    atomic_positions_str += "\n"
    return atomic_positions_str

def write_k_points(other_parameters = {}, lattice = {}, system_type = "metal") -> str:
    """
    write k points
    """
    k_points_str = "K_POINTS automatic\n"

    if system_type != "isolated":
        
        lattice_length_band_folding_maximum = 31.0 # angstrom, means nk*lattice should be no less than 25 angstrom
        llbfm = lattice_length_band_folding_maximum
        if system_type != "metal":
            llbfm = 20.0
        nkpts = [
            int(llbfm / lattice["lattice_parameters"]["a"]),
            int(llbfm / lattice["lattice_parameters"]["b"]),
            int(llbfm / lattice["lattice_parameters"]["c"])
        ]
    
    k_points = basic_parameters["k_points"]

    k_points["k_points"][0] = nkpts[0]
    k_points["k_points"][1] = nkpts[1]
    k_points["k_points"][2] = nkpts[2]

    if "k_points" in other_parameters.keys():
        for key, value in other_parameters["k_points"].items():
            k_points[key] = value

    k_points_str += "%d %d %d %d %d %d\n"%(
        k_points["k_points"][0], 
        k_points["k_points"][1], 
        k_points["k_points"][2], 
        k_points["k_points_shift"][0], 
        k_points["k_points_shift"][1], 
        k_points["k_points_shift"][2]
        )
    k_points_str += "\n"
    return k_points_str

def cif_to_qespresso(
        cif_file: str,
        other_parameters = {}) -> str:
    qe_fname = cif_file.split(".")[0] + ".in"
    cif = read(cif_file)
    atoms = m1ap(cif)
    lattice = m1l(cif)

    atomic_species = {
        "elements": [],
        "mass": [],
        "pseudopotentials": [],
    }
    for element in atoms.keys():
        atomic_species["elements"].append(element)
        atomic_species["mass"].append(ai.get_element_mass(element))
        atomic_species["pseudopotentials"].append(element + "_pseudopot")
    
    return_str = ""
    return_str += write_control(other_parameters)
    return_str += write_system(other_parameters, atoms)
    return_str += write_electrons(other_parameters)
    return_str += write_ions(other_parameters)
    return_str += write_cell(other_parameters)
    return_str += write_atomic_species(atomic_species)
    return_str += write_cell_parameters(lattice)
    return_str += write_atomic_positions(atoms)
    return_str += write_k_points(other_parameters, lattice)
    with open(qe_fname, 'w') as f:
        f.write(return_str)
    return qe_fname

def write_dimer_structure(symbol = "H", distance = 3.0) -> str:
    """
    write atomic positions, specifically for dimer
    """
    atomic_positions_str = "ATOMIC_POSITIONS (angstrom)\n"
    atomic_positions_str += "%s %12.8f %12.8f %12.8f"%(symbol, 0.0, 0.0, 0.0) + " 1 1 1\n"
    atomic_positions_str += "%s %12.8f %12.8f %12.8f"%(symbol, 0.0, 0.0, distance) + " 0 0 1\n"
    atomic_positions_str += "\n"
    return atomic_positions_str

def write_trimer_structure(symbol = "H", distance = 3.0) -> str:
    """
    write atomic positions, specifically for trimer, with one vertical angle
    """
    atomic_positions_str = "ATOMIC_POSITIONS (angstrom)\n"
    atomic_positions_str += "%s %12.8f %12.8f %12.8f"%(symbol, 0.0, 0.0, 0.0) + " 1 1 1\n"
    atomic_positions_str += "%s %12.8f %12.8f %12.8f"%(symbol, 0.0, 0.0, distance) + " 0 0 1\n"
    atomic_positions_str += "%s %12.8f %12.8f %12.8f"%(symbol, 0.0, distance, 0.0) + " 0 1 0\n"
    atomic_positions_str += "\n"
    return atomic_positions_str

def reference_structure_from_quantum_espresso(reference_structure = "dimer", symbol = "H", distance = 3.0) -> None:
    """quantum espresso input file generation for finding the bond length of reference structure

    Args:
        reference_structure (str, optional): reference structure type, dimer or trimer supported. Defaults to "dimer".
        symbol (str, optional): element symbol. Defaults to "H".
        distance (float, optional): characteristic distence between atoms. Defaults to 3.0 Angstrom.

    Raises:
        ValueError: reference_structure not properly set

    Returns:
        None

    Additional information:
        This function will be involved in the process of generating numerical orbitals. More specifically:
        element-wise pseudopotential tests:
            1. on quantum_espresso the pseudopotential validity test: preparation -> test -> analysis
            2. on PTG_dpsi-ABACUS the numerical orbitals generation: preparation -> generation
            3. on ABACUS the numerical orbitals accuracy test: preparation -> test -> analysis
    
    Generated files:
        filename: [symbol]_[reference_structure].in, e.g. H_dimer.in
        CONTROL, SYSTEM, ELECTRONS, IONS, CELL, ATOMIC_SPECIES sections are identical with cif_to_quantum_espresso
        ATOMIC_POSITIONS, CELL_PARAMETERS, K_POINTS sections are different.
        ATOMIC_POSITIONS: only two or three atoms are involved
        CELL_PARAMETERS: 20 Bohr for each lattice vector
        K_POINTS: 1 1 1 0 0 0 due to isolated systems calculation
    """
    if reference_structure == "dimer":
        atomic_positions_str = write_dimer_structure(symbol, distance)
    elif reference_structure == "trimer":
        atomic_positions_str = write_trimer_structure(symbol, distance)
    else:
        raise ValueError("reference_structure should be 'dimer' or 'trimer'")

    return_str = ""
    return_str += write_control({
        "control": {
            "calculation": "relax",
        }
    })
    return_str += write_system({}, {})
    return_str += write_electrons({})
    return_str += write_ions({})
    return_str += write_cell({})
    return_str += write_atomic_species({
        "elements": [symbol],
        "mass": [ai.get_element_mass(symbol)],
        "pseudopotentials": [symbol + "_pseudopot"],
    })
    return_str += write_cell_parameters({}, "nao")
    return_str += atomic_positions_str
    return_str += write_k_points({}, {}, "isolated")
    with open(symbol + "_" + reference_structure + ".in", 'w') as f:
        f.write(return_str)

if __name__ == "__main__":

    cif_to_qespresso("Ce.cif", {})