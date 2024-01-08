"""
Quantum ESPRESSO input file generator

Author: Kirk0830
Date: 2023-11-22
Refactor: 2024-01-06

Description:
    This file contains functions that can generate Quantum ESPRESSO input file from CIF file.

Useful cases:
    as it is
"""
import apns.module_structure.crystal_information_file as amscif
import apns.module_database.database as amdd
import seekpath as skps

"""call the following function in this way:
sections = ["&control", "&system", "&electrons", "&ions", "&cell"]
for section in sections:
    _section = _calculation(section, ntype, nat, calculation_settings_suite)
    with open(fname, 'a') as f:
        f.write(_section)
"""
def _calculation(section: str = "", ntype: int = 0, natom: int = 0, **kwargs) -> str:
    """new version of _control, _system, _electrons, _ions, _cell writes QE input file sections
    and will be iterated on.
    
    Args:
        section (str, optional): section name. Defaults to "".
        ntype (int, optional): number of atom types. Defaults to 0.
        natom (int, optional): number of atoms. Defaults to 0.
        
    Raises:
        ValueError: section not properly set
    
    Returns:
        str: section string
    """

    _section = section[1:] if section.startswith("&") else section
    if _section not in INPUT_TEMPLATE.keys():
        raise ValueError("section should be one of the following: {}".format(list(INPUT_TEMPLATE.keys())))
    section_content = INPUT_TEMPLATE[_section]
    if _section.lower() == "system":
        if ntype == 0 or natom == 0:
            raise ValueError("ntype and natom should be set")
        section_content["ntyp"] = ntype
        section_content["nat"] = natom
        ecut_warning = False
        if "ecutwfc" in kwargs.keys() and "ecutrho" not in kwargs.keys():
            kwargs["ecutrho"] = kwargs["ecutwfc"] * 4
            ecut_warning = True
        elif "ecutrho" in kwargs.keys() and "ecutwfc" not in kwargs.keys():
            kwargs["ecutwfc"] = kwargs["ecutrho"] / 4
            ecut_warning = True
        if ecut_warning:
            print("WARNING: ecutwfc and ecutrho should be set together, here dual = 4 is used by default.\n"
                + "         This value is suitable for norm-conserving pseudopotentials, but not for ultrasoft pseudopotentials.")
    for key, value in kwargs.items():
        if key in section_content.keys():
            section_content[key] = value
    result = "&{}\n".format(_section.upper())
    for key, value in section_content.items():
        # if value is string, add quotation marks
        if isinstance(value, str) and not (".true." in value or ".false." in value):
            result += "%s = '%s',\n"%(key, value)
        else:
            result += "%s = %s,\n"%(key, str(value))
    result += "/\n"
    return result

def _ATOMIC_SEPCIES(pseudopotential: dict = {}, **kwargs) -> str:
    """new version of ATOMIC_SPECIES writes QE input file section

    Args:
        pseudopotential (dict, optional): pseudopotential dictionary. Defaults to {}.

    Raises:
        ValueError: pseudopotential not properly set

    Returns:
        str: section string
    """
    result = "ATOMIC_SPECIES\n"
    if len(pseudopotential.keys()) == 0:
        raise ValueError("pseudopotential should be set")
    mass = kwargs.get("mass", {key: amdd.element_mass(key) for key in pseudopotential.keys()})
    for element in pseudopotential.keys():
        result += "{} {} {}\n".format(element, mass[element], pseudopotential[element])
    result += "\n"
    return result

def _K_POINTS(fname: str = "", nkpoints_in_line: int = 0) -> str:

    cif = amscif.read_1(fname)
    cell_parameters = cif["cell_parameters"]
    result = "K_POINTS"
    if nkpoints_in_line <= 0:
        result += " automatic\n"
        if nkpoints_in_line == 0:
            lattice_length_band_folding_maximum = 35.0 # angstrom, means nk*lattice should be no less than 31 angstrom
            llbfm = lattice_length_band_folding_maximum
            nkpts = [int(llbfm / cell_parameters[component]) for component in ["a", "b", "c"]]
        elif nkpoints_in_line == -1:
            nkpts = [1, 1, 1]
        result += "%d %d %d %d %d %d\n"%(nkpts[0], nkpts[1], nkpts[2], 0, 0, 0)
    else:
        result += "crystal\n"
        cell_vectors = amscif.cellparam_to_latvec(cell_parameters)
        positions = []
        numbers = []

        for element in cif["atomic_positions"].keys():
            for line in cif["atomic_positions"][element]:
                positions.append(line)
                numbers.append(amdd.get_element_index(element))
        
        _skps_result = skps.get_path(structure=(cell_vectors, positions, numbers))

        k_points = []
        for _path in _skps_result["path"]:
            start_coord = _skps_result["point_coords"][_path[0]]
            end_coord = _skps_result["point_coords"][_path[1]]
            dkx, dky, dkz = [(end_coord[i] - start_coord[i])/(nkpoints_in_line-1) for i in range(3)]
            for i in range(nkpoints_in_line - 1):
                k_points.append([start_coord[0]+i*dkx, start_coord[1]+i*dky, start_coord[2]+i*dkz])

        for k_point in k_points:
            result += "%12.8f %12.8f %12.8f 1\n"%(k_point[0], k_point[1], k_point[2])
    result += "\n"
    return result

def _CIF(fname: str = "", cell_scaling: float = 0.0, constraints: list = []) -> str:
    """new version of cif file to CELL_PARAMETERS, ATOMIC_POSITIONS, K_POINTS sections
    
    Args:
        fname (str, optional): cif file name. Defaults to "".
        cell_scaling (float, optional): cell scaling factor. Defaults to 0.0.
        constraints (list, optional): constraints. Defaults to [].

    Example of constraints:
    >>> constraints = [
    ...     [1, 1, 1] # in QE it means constraint on x, y, z, here for the first atom
    ...     [0, 0, 1] # in QE it means constraint on z, here for the second atom
    ... ]

    Raises:
    """
    if fname == "":
        raise ValueError("fname should be set")
    result = "CELL_PARAMETERS (angstrom)\n"
    cif = amscif.read_1(fname)
    cell_vectors = amscif.cellparam_to_latvec(cif["cell_parameters"])
    cell_vectors = [[(cell_scaling + 1.0) * cell_vectors[i][j] for j in range(3)] for i in range(3)]
    result += "\n".join(["%12.8f %12.8f %12.8f"%(cell_vectors[i][0], cell_vectors[i][1], cell_vectors[i][2]) for i in range(3)])
    
    natom = 0
    for element in cif["atomic_positions"].keys():
        natom += len(cif["atomic_positions"][element])
    result += "\n\nATOMIC_POSITIONS (crystal)\n"
    if len(constraints) == 0:
        constraints = [[0, 0, 0]] * natom
    elif len(constraints) < natom:
        print("WARNING: constraints length is less than natom, here constraints are automatically extended to natom")
        constraints = constraints + [[0, 0, 0]] * (natom - len(constraints))
    elif len(constraints) > natom:
        print("WARNING: constraints length is greater than natom, here constraints are automatically truncated to natom")
        constraints = constraints[:natom]
    for element in cif["atomic_positions"].keys():
        for iatom in range(len(cif["atomic_positions"][element])):
            result += "%2s %12.8f %12.8f %12.8f %d %d %d\n"%(element, cif["atomic_positions"][element][iatom][0], 
                                                                      cif["atomic_positions"][element][iatom][1], 
                                                                      cif["atomic_positions"][element][iatom][2], 
                                                            constraints[iatom][0], constraints[iatom][1], constraints[iatom][2])
    result += "\n"
    return result, natom

def _ISOLATED(element: str = "", shape: str = "", bond_length: float = 0.0, **kwargs) -> str:
    """new version of isolated system input file generation, support dimer, trimer, tetramer
    
    Args:
        element (str, optional): element symbol. Defaults to "".
        shape (str, optional): shape of the isolated system. Defaults to "".
        bond_length (float, optional): characteristic length of the isolated system. Defaults to 0.0.
    
    Raises:
        ValueError: element, shape, bond_length not properly set
    
    Returns:
        str: input file string
    """
    if element == "":
        raise ValueError("element should be set")
    if shape == "":
        raise ValueError("shape should be set")
    if bond_length == 0.0:
        raise ValueError("bond_length should be set")
    if "constraint" in kwargs.keys():
        print("WARNING: constraint is not supported in _ISOLATED")

    result = "CELL_PARAMETERS (angstrom)\n20.00000000 0.00000000 0.00000000\n0.00000000 20.00000000 0.00000000\n0.00000000 0.00000000 20.00000000\n\n"
    result += "\nATOMIC_POSITIONS (angstrom)\n"
    natom = 0
    if shape == "dimer":
        result += "%s %12.8f %12.8f %12.8f\n"%(element, 0.0, 0.0, 0.0)
        result += "%s %12.8f %12.8f %12.8f\n"%(element, bond_length, 0.0, 0.0)
        natom = 2
    elif shape == "trimer":
        result += "%s %12.8f %12.8f %12.8f\n"%(element, 0.0, 0.0, 0.0)
        result += "%s %12.8f %12.8f %12.8f\n"%(element, bond_length, 0.0, 0.0)
        result += "%s %12.8f %12.8f %12.8f\n"%(element, bond_length/2, bond_length/2*3**0.5, 0.0)
        natom = 3
    elif shape == "tetramer":
        result += "%s %12.8f %12.8f %12.8f\n"%(element, 0.0, 0.0, 0.0)
        result += "%s %12.8f %12.8f %12.8f\n"%(element, bond_length, 0.0, 0.0)
        result += "%s %12.8f %12.8f %12.8f\n"%(element, bond_length/2, bond_length/2*3**0.5, 0.0)
        result += "%s %12.8f %12.8f %12.8f\n"%(element, bond_length/2, bond_length/2*3**0.5/3, bond_length/2*3**0.5*2/3)
        natom = 4
    else:
        raise ValueError("shape should be one of the following: dimer, trimer, tetramer")
    result += "\n"
    return result, natom

"""deprecated old version of functions, but can still use in some way"""
# DEPRECATED
def section_control(other_parameters: dict) -> str:
    """
    write control parameters
    """
    control = INPUT_TEMPLATE["control"]
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
# DEPRECATED
def section_system(other_parameters: dict, atoms: dict) -> str:
    """
    write system parameters
    """
    system = INPUT_TEMPLATE["system"]
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
# DEPRECATED
def section_electrons(other_parameters: dict) -> str:
    """
    write electrons parameters
    """
    electrons = INPUT_TEMPLATE["electrons"]
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
# DEPRECATED
def section_ions(other_parameters: dict) -> str:
    """
    write ions parameters
    """
    ions = INPUT_TEMPLATE["ions"]
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
# DEPRECATED
def section_cell(other_parameters: dict) -> str:
    """
    write cell parameters
    """
    cell = INPUT_TEMPLATE["cell"]
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
# DEPRECATED
def section_atomic_species(atomic_species: dict) -> str:
    """
    write atomic species
    """
    atomic_species_str = "ATOMIC_SPECIES\n"
    for iatom in range(len(atomic_species["elements"])):
        atomic_species_str += atomic_species["elements"][iatom] + " " + str(atomic_species["mass"][iatom]) + " " + atomic_species["pseudopotentials"][iatom] + "\n"
    atomic_species_str += "\n"
    return atomic_species_str
# DEPRECATED
def section_cell_parameters(lattice: dict, mode = "cif") -> str:
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
# DEPRECATED
def section_atomic_positions(atoms: dict) -> str:
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
# DEPRECATED
def k_points(other_parameters = {}, lattice = {}, system_type = "metal") -> str:
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
    
    k_points = INPUT_TEMPLATE["k_points"]

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
# DEPRECATED
def cif_to_qespresso(
        cif_file: str,
        other_parameters = {}) -> str:
    qe_fname = cif_file.split(".")[0] + ".in"
    cif = amscif.read_1(cif_file)
    atoms = amscif._method1_atomic_position(cif)
    lattice = amscif._method1_lattice(cif)

    atomic_species = {
        "elements": [],
        "mass": [],
        "pseudopotentials": [],
    }
    for element in atoms.keys():
        atomic_species["elements"].append(element)
        atomic_species["mass"].append(amdd.element_mass(element))
        atomic_species["pseudopotentials"].append(element + "_pseudopot")
    
    return_str = ""
    return_str += section_control(other_parameters)
    return_str += section_system(other_parameters, atoms)
    return_str += section_electrons(other_parameters)
    return_str += section_ions(other_parameters)
    return_str += section_cell(other_parameters)
    return_str += section_atomic_species(atomic_species)
    return_str += section_cell_parameters(lattice)
    return_str += section_atomic_positions(atoms)
    return_str += k_points(other_parameters, lattice)
    with open(qe_fname, 'w') as f:
        f.write(return_str)
    return qe_fname
# DEPRECATED
def write_dimer_structure(symbol = "H", distance = 3.0) -> str:
    """
    write atomic positions, specifically for dimer
    """
    atomic_positions_str = "ATOMIC_POSITIONS (angstrom)\n"
    atomic_positions_str += "%s %12.8f %12.8f %12.8f"%(symbol, 0.0, 0.0, 0.0) + " 1 1 1\n"
    atomic_positions_str += "%s %12.8f %12.8f %12.8f"%(symbol, 0.0, 0.0, distance) + " 0 0 1\n"
    atomic_positions_str += "\n"
    return atomic_positions_str
# DEPRECATED
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
# DEPRECATED
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
    return_str += section_control({
        "control": {
            "calculation": "relax",
        }
    })
    return_str += section_system({}, {})
    return_str += section_electrons({})
    return_str += section_ions({})
    return_str += section_cell({})
    return_str += section_atomic_species({
        "elements": [symbol],
        "mass": [amdd.element_mass(symbol)],
        "pseudopotentials": [symbol + "_pseudopot"],
    })
    return_str += section_cell_parameters({}, "nao")
    return_str += atomic_positions_str
    return_str += k_points({}, {}, "isolated")
    with open(symbol + "_" + reference_structure + ".in", 'w') as f:
        f.write(return_str)

INPUT_TEMPLATE = {
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
        "ecutwfc": 100,
        "ecutrho": 400,
        "occupations": "smearing",
        "smearing": "gaussian",
        "degauss": 0.02,
        "nspin": 1,
        "starting_magnetization": 0.0,
        "nosym": ".false.",
        "noinv": ".false.",
        "noncolin": ".false.",
        "lspinorb": ".false.",
        "input_dft": "pbe",
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

if __name__ == "__main__":

    cif_to_qespresso("Ce.cif", {})