import re
import seekpath as skps
import apns.module_structure.cifparser as amscif
import apns.module_database.database as amdd

def STRU_cif(fname: str, 
             pseudopotentials: dict, 
             template: bool = False, 
             **kwargs):
    """convert cif file to STRU file.

    Args:
        fname (str): name of the cif file.
        pseudopotentials (dict): pseudopotentials information.
        template (bool, optional): whether to use template pseudopotential file name. Defaults to False.
        **kwargs: mass: dict, 
                  numerical_orbitals: dict, 
                  magnetism: dict,
                  constraints: dict,
                  cell_rescaling_factors: list, 
                  ...

    Returns:
        return_str (str): readily usable string for STRU file.
    """
    cif = amscif.read_1(fname)
    return_str = "ATOMIC_SPECIES\n"
    for _element in pseudopotentials:
        mass = kwargs.get("mass", {}).get(_element, 1.0)
        pseudopotential = _element + "_pseudopot"
        if not template:
            pseudopotential = pseudopotentials[_element]
        return_str += "%s %8.4f %s\n" % (_element, mass, pseudopotential)
    return_str += "\n"
    if "numerical_orbitals" in kwargs:
        return_str += "NUMERICAL_ORBITAL\n"
        for _element in kwargs["numerical_orbitals"]:
            numerical_orbital = _element + "_numerical_orbital"
            if not template:
                numerical_orbital = kwargs["numerical_orbitals"][_element]
            return_str += "%s\n" % (numerical_orbital)
        return_str += "\n"
    if template:
        return_str += "LATTICE_CONSTANT\nlattice_constant_to_test\n\n"
    else:
        return_str += "LATTICE_CONSTANT\n1.889726877\n\n"
    return_str += "LATTICE_VECTORS\n"
    cell_parameters = cif["cell_parameters"]
    lattice_vectors = amscif.to_latvec(**cell_parameters)

    for ilv, lattice_vector in enumerate(lattice_vectors):
        if "cell_rescaling_factors" in kwargs:
            lattice_vector[0] *= kwargs["cell_rescaling_factors"][ilv][0]
            lattice_vector[1] *= kwargs["cell_rescaling_factors"][ilv][1]
            lattice_vector[2] *= kwargs["cell_rescaling_factors"][ilv][2]

        return_str += "%12.8f %12.8f %12.8f\n" % (lattice_vector[0], lattice_vector[1], lattice_vector[2])
    return_str += "\n"
    return_str += "ATOMIC_POSITIONS\nDirect\n"
    atomic_positions = cif["atomic_positions"]
    for atom in atomic_positions.keys():
        return_str += "%s\n" % (atom)
        magnetism = 0.0
        if "magnetism" in kwargs:
            if atom in kwargs["magnetism"]:
                magnetism = kwargs["magnetism"][atom]
        return_str += "%4.2f\n"%(magnetism)
        return_str += "%d\n"%(len(atomic_positions[atom]))
        for ia, position in enumerate(atomic_positions[atom]):
            return_str += "%12.8f %12.8f %12.8f"%(position[0], position[1], position[2])
            if "constraints" in kwargs:
                return_str += " m %d %d %d"%(kwargs["constraints"][atom][ia][0], kwargs["constraints"][atom][ia][1], kwargs["constraints"][atom][ia][2])
            else:
                return_str += " m 1 1 1"
            return_str += "\n"

    return return_str, cell_parameters

def STRU_dimer(element: str, bond_length: float, pseudopotentials: dict, **kwargs):
    """_summary_

    Args:
        bond_length (float): dimer bond length, in Angstrom
        pseudopotentials (dict): key value pairs, key is element, value is pseudopotential file name.
        **kwargs: mass, numerical_orbitals, cell_parameters, magnetism, constraints, ...

    Returns:
        str: readily usable string for STRU file.
    """
    return_str = "ATOMIC_SPECIES\n"
    for _element in pseudopotentials:
        if "mass" in kwargs:
            return_str += "%s %8.4f %s\n" % (_element, kwargs["mass"][_element], pseudopotentials[_element])
        else:
            return_str += "%s %8.4f %s\n" % (_element, 1.0, pseudopotentials[_element])
    return_str += "\n"
    if "numerical_orbitals" in kwargs:
        return_str += "NUMERICAL_ORBITALS\n"
        for _element in kwargs["numerical_orbitals"]:
            return_str += "%s\n" % (kwargs["numerical_orbitals"][_element])
        return_str += "\n"
    return_str += "LATTICE_CONSTANT\n1.889726877\n\n"
    return_str += "LATTICE_VECTORS\n"
    if "cell_parameters" in kwargs:
        lattice_vectors = amscif.to_latvec(**kwargs["cell_parameters"])
    else:
        lattice_vectors = [
            [20.0, 0.0, 0.0],
            [0.0, 20.0, 0.0],
            [0.0, 0.0, 20.0]
        ]
    for lattice_vector in lattice_vectors:
        return_str += "%12.8f %12.8f %12.8f\n" % (lattice_vector[0], lattice_vector[1], lattice_vector[2])
    return_str += "\n"
    return_str += "ATOMIC_POSITIONS\nCartesian\n%s\n0.00\n2\n0.0 0.0 0.0 m 0 0 0\n0.0 0.0 %20.10f m 0 0 1"%(element, bond_length)

    return return_str

def STRU_Molecule(**kwargs):
    shape = kwargs.get("shape", "").lower()
    if shape == "":
        raise ValueError("shape must be specified.")
    pseudopotentials = kwargs.get("pseudopotentials", {})
    if len(pseudopotentials) == 0:
        raise ValueError("pseudopotentials must be specified.")
    if len(pseudopotentials) > 1:
        raise ValueError("pseudopotentials must be specified for only one element.")
    numerical_orbitals = kwargs.get("numerical_orbitals", None)
    if numerical_orbitals is not None and len(numerical_orbitals) > 1:
        raise ValueError("numerical_orbitals must be specified for only one element.")
    bond_length = kwargs.get("bond_length", 0.0)
    if bond_length == 0.0:
        raise ValueError("bond_length must be specified.")
    
    assert len(pseudopotentials) == 1, "Only one element is supported."
    element = list(pseudopotentials.keys())[0]

    return_str = "ATOMIC_SPECIES\n"
    for _element in pseudopotentials:
        mass = kwargs.get("mass", {}).get(_element, amdd.element_label_tomass(_element))
        pseudopotential = pseudopotentials[_element]
        return_str += "%s %8.4f %s\n" % (_element, mass, pseudopotential)
    return_str += "\n"
    if numerical_orbitals is not None:
        return_str += "NUMERICAL_ORBITAL\n"
        for _element in numerical_orbitals.keys():
            numerical_orbital = _element + "_numerical_orbital"
            numerical_orbital = numerical_orbitals[_element]
            return_str += "%s\n" % (numerical_orbital)
        return_str += "\n"

    lattice_constant = 1.889726877
    return_str += "LATTICE_CONSTANT\n%s\n\n"%(lattice_constant)
    return_str += "LATTICE_VECTORS\n"
    return_str += "%12.8f %12.8f %12.8f\n" %(20.0, 0.0, 0.0)
    return_str += "%12.8f %12.8f %12.8f\n" %(0.0, 20.0, 0.0)
    return_str += "%12.8f %12.8f %12.8f\n" %(0.0, 0.0, 20.0)

    return_str += "\n"
    return_str += "ATOMIC_POSITIONS\nDirect\n"

    return_str += "%s\n 0.00" % (element)
    
    if shape == "dimer":
        return_str += "\n2\n"
        return_str += "%12.8f %12.8f %12.8f m 0 0 0\n"%(0.0, 0.0, 0.0)
        return_str += "%12.8f %12.8f %12.8f m 1 0 0\n"%(bond_length, 0.0, 0.0)
    elif shape == "trimer":
        return_str += "\n3\n"
        return_str += "%12.8f %12.8f %12.8f m 0 0 0\n"%(0.0, 0.0, 0.0)
        return_str += "%12.8f %12.8f %12.8f m 1 0 0\n"%(bond_length, 0.0, 0.0)
        return_str += "%12.8f %12.8f %12.8f m 1 1 0\n"%(bond_length/2, bond_length/2*3**0.5, 0.0)
    elif shape == "tetramer":
        return_str += "\n4\n"
        return_str += "%12.8f %12.8f %12.8f m 0 0 0\n"%(0.0, 0.0, 0.0)
        return_str += "%12.8f %12.8f %12.8f m 1 0 0\n"%(bond_length, 0.0, 0.0)
        return_str += "%12.8f %12.8f %12.8f m 1 1 0\n"%(bond_length/2, bond_length/2*3**0.5, 0.0)
        return_str += "%12.8f %12.8f %12.8f m 1 1 1\n"%(bond_length/2, bond_length/2/3**0.5, bond_length/2*3**0.5)
    else:
        raise ValueError(f"shape must be dimer, trimer or tetramer: {shape} is not supported.")

    return return_str, [20.0, 20.0, 20.0, 90.0, 90.0, 90.0]

def STRU_ACWFRef(**kwargs):
    """for generating STRU file that can take ACWF all-electron calculation results
    as reference. Not realistic structures. May coincide with real structures in some cases,
    but not always."""
    pseudopotentials = kwargs.get("pseudopotentials", {})
    assert len(pseudopotentials) > 0, "Pseudopotentials must be specified."
    numerical_orbitals = kwargs.get("numerical_orbitals", {})
    cell_scaling = kwargs.get("cell_scaling", 0.0)
    starting_magnetization = kwargs.get("starting_magnetization", None)
    bravis = kwargs.get("shape", "sc").lower()
    assert bravis in ["sc", "bcc", "fcc", "diamond"], "bravis must be one of sc, bcc, fcc, diamond."

    assert len(pseudopotentials) == 1, "Only one element is supported."
    element = list(pseudopotentials.keys())[0]
    token = f"{element.capitalize()}-X/"
    token += bravis.upper() if bravis != "diamond" else "Diamond"
    # will look up the database of ACWF...
    volume = acwf_refvolume(token) * (1.0 + cell_scaling) # in Angstrom
    assert volume is not None, "Volume not found in ACWF database."
    celldm = amscif.volume_to_celldm(bravis, volume) # in Angstrom
    print(f"Generate with Volume: {volume:.4f} A^3, celldm: {celldm:.4f} A")

    # constant maintained inside
    ALPHA_BETA_GAMMA = {"sc": [90.0, 90.0, 90.0], "bcc": [109.4712206, 109.4712206, 109.4712206], 
                        "fcc": [60.0, 60.0, 60.0], "diamond": [60.0, 60.0, 60.0]}
    ATOMIC_POSITIONS = {"sc": [[0.0, 0.0, 0.0]], "bcc": [[0.0, 0.0, 0.0]],
                        "fcc": [[0.0, 0.0, 0.0]], "diamond": [[0.5, 0.5, 0.5], [0.75, 0.75, 0.75]]}
    if starting_magnetization is None:
        starting_magnetization = {element: 0.0 for element in pseudopotentials.keys()}
    return_str = "ATOMIC_SPECIES\n"
    species = {element: {"pseudopotential": pseudopotentials.get(element, None),
                         "numerical_orbital": numerical_orbitals.get(element, None),
                         "starting_magnetization": starting_magnetization.get(element, 0.0),
                         "atomic_positions": ATOMIC_POSITIONS[bravis]}}
    for element in species.keys():
        return_str += "%s %8.4f %s\n" % (element,
                                         amdd.element_label_tomass(element),
                                         species[element]["pseudopotential"])
    return_str += "\n"
    if numerical_orbitals is not None and len(numerical_orbitals) > 0:
        return_str += "NUMERICAL_ORBITAL\n"
        for element in species.keys():
            return_str += "%s\n" % (species[element]["numerical_orbital"])
        return_str += "\n"
    # scaling lattice volume by cell_scaling
    lattice_constant = 1.889726877
    return_str += "LATTICE_CONSTANT\n%s\n\n"%(lattice_constant)
    return_str += "LATTICE_VECTORS\n"
    cell_param = dict(zip(["a", "b", "c"], [celldm]*3))
    cell_param.update(dict(zip(["alpha", "beta", "gamma"], ALPHA_BETA_GAMMA[bravis])))
    lattice_vectors = amscif.to_latvec(**cell_param) # 3*3 list of lists
    for lattice_vector in lattice_vectors:
        return_str += "%12.8f %12.8f %12.8f\n" %(lattice_vector[0],
                                                 lattice_vector[1],
                                                 lattice_vector[2])
    return_str += "\n"
    return_str += "ATOMIC_POSITIONS\nDirect\n"
    for element in species.keys():
        return_str += "%s\n" % (element)
        return_str += "%4.2f\n"%(species[element]["starting_magnetization"])
        return_str += "%d\n"%(len(species[element]["atomic_positions"]))
        for ia, position in enumerate(species[element]["atomic_positions"]):
            return_str += "%12.8f %12.8f %12.8f"%(position[0], position[1], position[2])
            return_str += " m 1 1 1"
            return_str += "\n"
    
    return return_str, [celldm, celldm, celldm, ALPHA_BETA_GAMMA[bravis][0], ALPHA_BETA_GAMMA[bravis][1], ALPHA_BETA_GAMMA[bravis][2]]

import apns.module_structure.cifparser_pymatgen as amcp
import apns.module_structure.basic as ambs
def STRU_Pymatgen(**kwargs):
    """with Pymatgen's CifParser, convert cif file to STRU file."""
    fname = kwargs.get("fname", None)
    assert fname is not None, "Cif file name must be specified."
    pseudopotentials = kwargs.get("pseudopotentials", {})
    numerical_orbitals = kwargs.get("numerical_orbitals", None)
    cell_scaling = kwargs.get("cell_scaling", 0.0)
    starting_magnetization = kwargs.get("starting_magnetization", None)

    symbols, xyzs = amcp.structure(fname)
    a, b, c, alpha, beta, gamma, lattice_vectors = amcp.lattice(fname)
    if starting_magnetization is None:
        starting_magnetization = {element: 0.0 for element in pseudopotentials.keys()}

    return_str = "ATOMIC_SPECIES\n"
    species = ambs.expand_atomic_species(symbols=symbols,
                                         atomic_positions=xyzs,
                                         pseudopotentials=pseudopotentials,
                                         numerical_orbitals=numerical_orbitals,
                                         starting_magnetization=starting_magnetization)
    for element in species.keys():
        return_str += "%s %8.4f %s\n" % (element,
                                         amdd.element_label_tomass(element),
                                         species[element]["pseudopotential"])
    return_str += "\n"
    if numerical_orbitals is not None and len(numerical_orbitals) > 0:
        return_str += "NUMERICAL_ORBITAL\n"
        for element in species.keys():
            return_str += "%s\n" % (species[element]["numerical_orbital"])
        return_str += "\n"
    # scaling lattice volume by cell_scaling
    lattice_constant = 1.889726877 * (1.0 + cell_scaling)**(1/3)
    
    return_str += "LATTICE_CONSTANT\n%s\n\n"%(lattice_constant)
    return_str += "LATTICE_VECTORS\n"
    for lattice_vector in lattice_vectors:
        return_str += "%12.8f %12.8f %12.8f\n" %(lattice_vector[0],
                                                 lattice_vector[1],
                                                 lattice_vector[2])
    return_str += "\n"
    return_str += "ATOMIC_POSITIONS\nDirect\n"
    for element in species.keys():
        return_str += "%s\n" % (element)
        return_str += "%4.2f\n"%(species[element]["starting_magnetization"])
        return_str += "%d\n"%(len(species[element]["atomic_positions"]))
        for ia, position in enumerate(species[element]["atomic_positions"]):
            return_str += "%12.8f %12.8f %12.8f"%(position[0], position[1], position[2])
            return_str += " m 1 1 1"
            return_str += "\n"
    
    return return_str, [a, b, c, alpha, beta, gamma]

def KPT(isolated: bool = False, gamma_centered: bool = True, **kwargs):
    """new version of KPT file generation

    Args:
        mode (str): can be "isolated" or "crystal", "isolated" means the calculation is for isolated system, "crystal" means the calculation is for crystal system.
        gamma_centered (bool, optional): If kpoints sampling is Gamma point centered. Defaults to True.
        **kwargs: cell_parameters: list, metallic: bool, ...

    Raises:
        NotImplementedError: Non-Gamma centered kpoints sampling is not implemented yet.
        ValueError: mode must be "isolated" or "crystal".

    Returns:
        _type_: readily usable string for KPT file.
    """
    if "cell_parameters" in kwargs and not isolated:
        cell_parameters = kwargs["cell_parameters"]
        if isinstance(cell_parameters, list):
            a, b, c = cell_parameters[0], cell_parameters[1], cell_parameters[2]
        elif isinstance(cell_parameters, dict):
            a, b, c = cell_parameters["a"], cell_parameters["b"], cell_parameters["c"]
        else:
            raise TypeError("Not supported data format of cell parameters, must be list or dict.")
        _fold = 35 if kwargs.get("metallic", False) else 25
        return_str = "K_POINTS\n0\n%s\n%d %d %d 0 0 0\n" % ("Gamma" if gamma_centered else "MP", int(_fold/a) + 1, int(_fold/b) + 1, int(_fold/c) + 1)
    else:
        return_str = "K_POINTS automatic\n0\nGamma\n1 1 1 0 0 0\n"
        print("Warning: KPT file generation failed, use default KPT file instead.") if not isolated else None
    
    return return_str

import numpy as np
def KSPACING(kspacing: float|list[float], kspc_unit: str = "bohr", cellparam_unit: str = "Angstrom", **kwargs):
    """for high precision, set kspacing dynamically instead of a fixed fold value
    (as in KPT)."""
    kspacing = [kspacing]*3 if isinstance(kspacing, float) else kspacing
    assert isinstance(kspacing, list) and len(kspacing) == 3, "kspacing must be a list of 3 numbers."
    kspacing = [kspc/amdd.unit_conversion(1, kspc_unit, "Angstrom") for kspc in kspacing] # convert to Angstrom-1
    cellparams = kwargs.get("cell", None)
    if cellparams is not None:
        if isinstance(cellparams, list):
            a, b, c = cellparams[0], cellparams[1], cellparams[2]
            alpha, beta, gamma = cellparams[3], cellparams[4], cellparams[5]
        elif isinstance(cellparams, dict):
            a, b, c = cellparams.get("a", 0.0), cellparams.get("b", 0.0), cellparams.get("c", 0.0)
            assert a != 0.0 and b != 0.0 and c != 0.0, "Cell parameters a, b and c must be specified."
            alpha, beta, gamma = cellparams.get("alpha", 90.0), cellparams.get("beta", 90.0), cellparams.get("gamma", 90.0)
        else:
            raise TypeError("Not supported data format of cell parameters, must be list or dict.")
    # convert unit of cell parameters
    a, b, c = [amdd.unit_conversion(cellparam, cellparam_unit, "Angstrom") for cellparam in [a, b, c]]
    latvecs = amscif.to_latvec(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)
    recip_latvecs = np.linalg.inv(np.array(latvecs)).T
    # get norm of reciprocal lattice vectors
    a, b, c = [np.linalg.norm(vec)*2*np.pi for vec in recip_latvecs]
    # kspacing is in 1/Angstrom
    nkpts = [max(1, int(a/kspacing[0]) + 1),
             max(1, int(b/kspacing[1]) + 1),
             max(1, int(c/kspacing[2]) + 1)]
    return_str = "K_POINTS\n0\nGamma\n"
    return_str += f"{nkpts[0]} {nkpts[1]} {nkpts[2]} 0 0 0\n"
    return return_str

def KLINE(fname: str, nkpts_in_line: int = 10):

    result = amscif.read_1(fname)

    cell_vectors = amscif.to_latvec(**result['cell_parameters'])
    positions = []
    numbers = []

    for element in result['atomic_positions']:
        for line in result['atomic_positions'][element]:
            positions.append(line)
            numbers.append(amdd.element_label_toindex(element))

    result = skps.get_path(structure=(cell_vectors, positions, numbers))

    kpts = []
    for _path in result["path"]:
        start_coord = result["point_coords"][_path[0]] # fractional coordinates
        end_coord = result["point_coords"][_path[1]] # fractional coordinates
        dkx, dky, dkz = [(end_coord[i] - start_coord[i])/(nkpts_in_line-1) for i in range(3)]
        for i in range(nkpts_in_line - 1):
            kpts.append([start_coord[0]+i*dkx, start_coord[1]+i*dky, start_coord[2]+i*dkz])

    return_str = "K_POINTS\n0\n"
    return_str += "line\n%d\n" % (len(kpts))
    for kpt in kpts:
        return_str += "%12.8f %12.8f %12.8f 1\n" % (kpt[0], kpt[1], kpt[2])

    return return_str

def INPUT(calculation: dict,
          minimal: bool = True,
          template: bool = False):
    """new version of function INPUT for generating INPUT file of ABACUS,
     support enumeration on any number of keywords specified as list """
    minimal_inp = {
        "calculation": "scf",
        "ecutwfc": 100,
        "basis_type": "pw",
        "symmetry": 0,
        "pseudo_dir": "./",
        "scf_nmax": 100,
        "scf_thr": 1e-6,
        "dft_functional": "pbe"
    }
    values = {}
    comments = {}
    # parse input to key, value and comments
    pattern = r"^([\w_]*)([\s\w\.+-]*)([#]?)(.*)$"
    for line in INPUT_TEMPLATE.split("\n"):
        if line == "INPUT_PARAMETERS" or line.startswith("#"):
            continue
        _match = re.match(pattern, line)
        if _match and _match.group(1) != "":
            values[_match.group(1)] = _match.group(2)
            comments[_match.group(1)] = _match.group(4)
    # for template, change its value to xxx_to_test placeholder, then sed in other functions
    work_status_specified_valid_keys = []
    for key in calculation.keys():
        valid_key = key
        if key == "functionals":
            valid_key = "dft_functional"
        if isinstance(calculation[key], list):
            if template:
                values[valid_key] = valid_key + "_to_test"
            else:
                raise TypeError("For non-template mode input script generation, list-organized input is not supported")
        else:
            values[valid_key] = calculation[key]
        work_status_specified_valid_keys.append(valid_key)

    # for minimal, only keep the minimal keys, but add manually set keys
    for key in work_status_specified_valid_keys:
        minimal_inp[key] = values[key]

    # write.
    return_str = "INPUT_PARAMETERS\n"
    if minimal:
        for key in minimal_inp.keys():
            return_str += "%s %s # %s\n"%(key, values[key], comments[key])
    else:
        for key in values.keys():
            return_str += "%s %s # %s\n"%(key, values[key], comments[key])
    return return_str

INPUT_TEMPLATE = """INPUT_PARAMETERS
#Parameters (1.General)
suffix                         ABACUS #the name of main output directory
latname                        none #the name of lattice name
stru_file                      STRU #the filename of file containing atom positions
kpoint_file                    KPT #the name of file containing k points
pseudo_dir                     ./ #the directory containing pseudo files
orbital_dir                    ./ #the directory containing orbital files
pseudo_rcut                    15 #cut-off radius for radial integration
pseudo_mesh                    0 #0: use our own mesh to do radial renormalization; 1: use mesh as in QE
lmaxmax                        2 #maximum of l channels used
dft_functional                 default #exchange correlation functional
xc_temperature                 0 #temperature for finite temperature functionals
calculation                    scf #test; scf; relax; nscf; get_wf; get_pchg
esolver_type                   ksdft #the energy solver: ksdft, sdft, ofdft, tddft, lj, dp
ntype                          1 #atom species number
nspin                          1 #1: single spin; 2: up and down spin; 4: noncollinear spin
kspacing                       0 0 0  #unit in 1/bohr, should be > 0, default is 0 which means read KPT file
min_dist_coef                  0.2 #factor related to the allowed minimum distance between two atoms
nbands                         10 #number of bands
nbands_sto                     256 #number of stochastic bands
nbands_istate                  5 #number of bands around Fermi level for get_pchg calulation
symmetry                       0 #the control of symmetry
init_vel                       0 #read velocity from STRU or not
symmetry_prec                  1e-06 #accuracy for symmetry
symmetry_autoclose             1 #whether to close symmetry automatically when error occurs in symmetry analysis
nelec                          0 #input number of electrons
out_mul                        0 # mulliken  charge or not
noncolin                       0 #using non-collinear-spin
lspinorb                       0 #consider the spin-orbit interaction
kpar                           1 #devide all processors into kpar groups and k points will be distributed among each group
bndpar                         1 #devide all processors into bndpar groups and bands will be distributed among each group
out_freq_elec                  0 #the frequency ( >= 0) of electronic iter to output charge density and wavefunction. 0: output only when converged
dft_plus_dmft                  0 #true:DFT+DMFT; false: standard DFT calcullation(default)
rpa                            0 #true:generate output files used in rpa calculation; false:(default)
printe                         100 #Print out energy for each band for every printe steps
mem_saver                      0 #Only for nscf calculations. if set to 1, then a memory saving technique will be used for many k point calculations.
diago_proc                     1 #the number of procs used to do diagonalization
nbspline                       -1 #the order of B-spline basis
wannier_card                   none #input card for wannier functions
soc_lambda                     1 #The fraction of averaged SOC pseudopotential is given by (1-soc_lambda)
cal_force                      1 #if calculate the force at the end of the electronic iteration
out_freq_ion                   0 #the frequency ( >= 0 ) of ionic step to output charge density and wavefunction. 0: output only when ion steps are finished
device                         cpu #the computing device for ABACUS

#Parameters (2.PW)
ecutwfc                        100 ##energy cutoff for wave functions
ecutrho                        400 ##energy cutoff for charge density and potential
erf_ecut                       0 ##the value of the constant energy cutoff
erf_height                     0 ##the height of the energy step for reciprocal vectors
erf_sigma                      0.1 ##the width of the energy step for reciprocal vectors
fft_mode                       0 ##mode of FFTW
pw_diag_thr                    0.01 #threshold for eigenvalues is cg electron iterations
scf_thr                        1e-08 #charge density error
scf_thr_type                   2 #type of the criterion of scf_thr, 1: reci drho for pw, 2: real drho for lcao
init_wfc                       atomic #start wave functions are from 'atomic', 'atomic+random', 'random' or 'file'
init_chg                       atomic #start charge is from 'atomic' or file
chg_extrap                     atomic #atomic; first-order; second-order; dm:coefficients of SIA
out_chg                        0 #>0 output charge density for selected electron steps
out_pot                        0 #output realspace potential
out_wfc_pw                     0 #output wave functions
out_wfc_r                      0 #output wave functions in realspace
out_dos                        0 #output energy and dos
out_band                       0 #output energy and band structure
out_proj_band                  0 #output projected band structure
restart_save                   0 #print to disk every step for restart
restart_load                   0 #restart from disk
read_file_dir                  auto #directory of files for reading
nx                             0 #number of points along x axis for FFT grid
ny                             0 #number of points along y axis for FFT grid
nz                             0 #number of points along z axis for FFT grid
ndx                            0 #number of points along x axis for FFT smooth grid
ndy                            0 #number of points along y axis for FFT smooth grid
ndz                            0 #number of points along z axis for FFT smooth grid
cell_factor                    1.2 #used in the construction of the pseudopotential tables
pw_seed                        1 #random seed for initializing wave functions

#Parameters (3.Stochastic DFT)
method_sto                     2 #1: slow and save memory, 2: fast and waste memory
npart_sto                      1 #Reduce memory when calculating Stochastic DOS
nbands_sto                     256 #number of stochstic orbitals
nche_sto                       100 #Chebyshev expansion orders
emin_sto                       0 #trial energy to guess the lower bound of eigen energies of the Hamitonian operator
emax_sto                       0 #trial energy to guess the upper bound of eigen energies of the Hamitonian operator
seed_sto                       0 #the random seed to generate stochastic orbitals
initsto_ecut                   0 #maximum ecut to init stochastic bands
initsto_freq                   0 #frequency to generate new stochastic orbitals when running md
cal_cond                       0 #calculate electronic conductivities
cond_che_thr                   1e-08 #control the error of Chebyshev expansions for conductivities
cond_dw                        0.1 #frequency interval for conductivities
cond_wcut                      10 #cutoff frequency (omega) for conductivities
cond_dt                        0.02 #t interval to integrate Onsager coefficiencies
cond_dtbatch                   0 #exp(iH*dt*cond_dtbatch) is expanded with Chebyshev expansion.
cond_smear                     1 #Smearing method for conductivities
cond_fwhm                      0.4 #FWHM for conductivities
cond_nonlocal                  1 #Nonlocal effects for conductivities

#Parameters (4.Relaxation)
ks_solver                      genelpa #cg; dav; lapack; genelpa; scalapack_gvx; cusolver
scf_nmax                       100 ##number of electron iterations
relax_nmax                     1 #number of ion iteration steps
out_stru                       0 #output the structure files after each ion step
force_thr                      0.001 #force threshold, unit: Ry/Bohr
force_thr_ev                   0.0257112 #force threshold, unit: eV/Angstrom
force_thr_ev2                  0 #force invalid threshold, unit: eV/Angstrom
relax_cg_thr                   0.5 #threshold for switching from cg to bfgs, unit: eV/Angstrom
stress_thr                     0.5 #stress threshold
press1                         0 #target pressure, unit: KBar
press2                         0 #target pressure, unit: KBar
press3                         0 #target pressure, unit: KBar
relax_bfgs_w1                  0.01 #wolfe condition 1 for bfgs
relax_bfgs_w2                  0.5 #wolfe condition 2 for bfgs
relax_bfgs_rmax                0.8 #maximal trust radius, unit: Bohr
relax_bfgs_rmin                1e-05 #minimal trust radius, unit: Bohr
relax_bfgs_init                0.5 #initial trust radius, unit: Bohr
cal_stress                     1 #calculate the stress or not
fixed_axes                     None #which axes are fixed
fixed_ibrav                    0 #whether to preseve lattice type during relaxation
fixed_atoms                    0 #whether to preseve direct coordinates of atoms during relaxation
relax_method                   cg #bfgs; sd; cg; cg_bfgs;
relax_new                      1 #whether to use the new relaxation method
relax_scale_force              0.5 #controls the size of the first CG step if relax_new is true
out_level                      ie #ie(for electrons); i(for ions);
out_dm                         0 #>0 output density matrix
out_bandgap                    0 #if true, print out bandgap
use_paw                        0 #whether to use PAW in pw calculation
deepks_out_labels              0 #>0 compute descriptor for deepks
deepks_scf                     0 #>0 add V_delta to Hamiltonian
deepks_bandgap                 0 #>0 for bandgap label
deepks_out_unittest            0 #if set 1, prints intermediate quantities that shall be used for making unit test
deepks_model                    #file dir of traced pytorch model: 'model.ptg

#Parameters (5.LCAO)
basis_type                     lcao #PW; LCAO in pw; LCAO
nb2d                           0 #2d distribution of atoms
gamma_only                     0 #Only for localized orbitals set and gamma point. If set to 1, a fast algorithm is used
search_radius                  -1 #input search radius (Bohr)
search_pbc                     1 #input periodic boundary condition
lcao_ecut                      100 #energy cutoff for LCAO
lcao_dk                        0.01 #delta k for 1D integration in LCAO
lcao_dr                        0.01 #delta r for 1D integration in LCAO
lcao_rmax                      30 #max R for 1D two-center integration table
out_mat_hs                     0 #output H and S matrix
out_mat_hs2                    0 #output H(R) and S(R) matrix
out_mat_dh                     0 #output of derivative of H(R) matrix
out_interval                   1 #interval for printing H(R) and S(R) matrix during MD
out_app_flag                   1 #whether output r(R), H(R), S(R), T(R), and dH(R) matrices in an append manner during MD
out_mat_t                      0 #output T(R) matrix
out_element_info               0 #output (projected) wavefunction of each element
out_mat_r                      0 #output r(R) matrix
out_wfc_lcao                   0 #ouput LCAO wave functions, 0, no output 1: text, 2: binary
bx                             0 #division of an element grid in FFT grid along x
by                             0 #division of an element grid in FFT grid along y
bz                             0 #division of an element grid in FFT grid along z

#Parameters (6.Smearing)
smearing_method                gaussian #type of smearing_method: gauss; fd; fixed; mp; mp2; mv
smearing_sigma                 0.02 #energy range for smearing

#Parameters (7.Charge Mixing)
mixing_type                    broyden #plain; pulay; broyden
mixing_beta                    0.4 #mixing parameter: 0 means no new charge
mixing_ndim                    8 #mixing dimension in pulay or broyden
mixing_gg0                     1 #mixing parameter in kerker
mixing_beta_mag                -10 #mixing parameter for magnetic density
mixing_gg0_mag                 0 #mixing parameter in kerker
mixing_gg0_min                 0.1 #the minimum kerker coefficient
mixing_angle                   -10 #angle mixing parameter for non-colinear calculations
mixing_tau                     0 #whether to mix tau in mGGA calculation
mixing_dftu                    0 #whether to mix locale in DFT+U calculation

#Parameters (8.DOS)
dos_emin_ev                    -15 #minimal range for dos
dos_emax_ev                    15 #maximal range for dos
dos_edelta_ev                  0.01 #delta energy for dos
dos_scale                      0.01 #scale dos range by
dos_sigma                      0.07 #gauss b coefficeinet(default=0.07)
dos_nche                       100 #orders of Chebyshev expansions for dos

#Parameters (9.Molecular dynamics)
md_type                        nvt #choose ensemble
md_thermostat                  nhc #choose thermostat
md_nstep                       10 #md steps
md_dt                          1 #time step
md_tchain                      1 #number of Nose-Hoover chains
md_tfirst                      -1 #temperature first
md_tlast                       -1 #temperature last
md_dumpfreq                    1 #The period to dump MD information
md_restartfreq                 5 #The period to output MD restart information
md_seed                        -1 #random seed for MD
md_prec_level                  0 #precision level for vc-md
ref_cell_factor                1 #construct a reference cell bigger than the initial cell
md_restart                     0 #whether restart
lj_rcut                        8.5 #cutoff radius of LJ potential
lj_epsilon                     0.01032 #the value of epsilon for LJ potential
lj_sigma                       3.405 #the value of sigma for LJ potential
pot_file                       graph.pb #the filename of potential files for CMD such as DP
msst_direction                 2 #the direction of shock wave
msst_vel                       0 #the velocity of shock wave
msst_vis                       0 #artificial viscosity
msst_tscale                    0.01 #reduction in initial temperature
msst_qmass                     -1 #mass of thermostat
md_tfreq                       0 #oscillation frequency, used to determine qmass of NHC
md_damp                        1 #damping parameter (time units) used to add force in Langevin method
md_nraise                      1 #parameters used when md_type=nvt
cal_syns                       0 #calculate asynchronous overlap matrix to output for Hefei-NAMD
dmax                           0.01 #maximum displacement of all atoms in one step (bohr)
md_tolerance                   100 #tolerance for velocity rescaling (K)
md_pmode                       iso #NPT ensemble mode: iso, aniso, tri
md_pcouple                     none #whether couple different components: xyz, xy, yz, xz, none
md_pchain                      1 #num of thermostats coupled with barostat
md_pfirst                      -1 #initial target pressure
md_plast                       -1 #final target pressure
md_pfreq                       0 #oscillation frequency, used to determine qmass of thermostats coupled with barostat
dump_force                     1 #output atomic forces into the file MD_dump or not
dump_vel                       1 #output atomic velocities into the file MD_dump or not
dump_virial                    1 #output lattice virial into the file MD_dump or not

#Parameters (10.Electric field and dipole correction)
efield_flag                    0 #add electric field
dip_cor_flag                   0 #dipole correction
efield_dir                     2 #the direction of the electric field or dipole correction
efield_pos_max                 0.5 #position of the maximum of the saw-like potential along crystal axis efield_dir
efield_pos_dec                 0.1 #zone in the unit cell where the saw-like potential decreases
efield_amp                     0 #amplitude of the electric field

#Parameters (11.Gate field)
gate_flag                      0 #compensating charge or not
zgate                          0.5 #position of charged plate
relax                          0 #allow relaxation along the specific direction
block                          0 #add a block potential or not
block_down                     0.45 #low bound of the block
block_up                       0.55 #high bound of the block
block_height                   0.1 #height of the block

#Parameters (12.Test)
out_alllog                     0 #output information for each processor, when parallel
nurse                          0 #for coders
colour                         0 #for coders, make their live colourful
t_in_h                         1 #calculate the kinetic energy or not
vl_in_h                        1 #calculate the local potential or not
vnl_in_h                       1 #calculate the nonlocal potential or not
vh_in_h                        1 #calculate the hartree potential or not
vion_in_h                      1 #calculate the local ionic potential or not
test_force                     0 #test the force
test_stress                    0 #test the force
test_skip_ewald                0 #skip ewald energy

#Parameters (13.vdW Correction)
vdw_method                     none #the method of calculating vdw (none ; d2 ; d3_0 ; d3_bj
vdw_s6                         default #scale parameter of d2/d3_0/d3_bj
vdw_s8                         default #scale parameter of d3_0/d3_bj
vdw_a1                         default #damping parameter of d3_0/d3_bj
vdw_a2                         default #damping parameter of d3_bj
vdw_d                          20 #damping parameter of d2
vdw_abc                        0 #third-order term?
vdw_C6_file                    default #filename of C6
vdw_C6_unit                    Jnm6/mol #unit of C6, Jnm6/mol or eVA6
vdw_R0_file                    default #filename of R0
vdw_R0_unit                    A #unit of R0, A or Bohr
vdw_cutoff_type                radius #expression model of periodic structure, radius or period
vdw_cutoff_radius              default #radius cutoff for periodic structure
vdw_radius_unit                Bohr #unit of radius cutoff for periodic structure
vdw_cn_thr                     40 #radius cutoff for cn
vdw_cn_thr_unit                Bohr #unit of cn_thr, Bohr or Angstrom
vdw_cutoff_period   3 3 3 #periods of periodic structure

#Parameters (14.exx)
exx_hybrid_alpha               default #fraction of Fock exchange in hybrid functionals
exx_hse_omega                  0.11 #range-separation parameter in HSE functional
exx_separate_loop              1 #if 1, a two-step method is employed, else it will start with a GGA-Loop, and then Hybrid-Loop
exx_hybrid_step                100 #the maximal electronic iteration number in the evaluation of Fock exchange
exx_mixing_beta                1 #mixing_beta for outer-loop when exx_separate_loop=1
exx_lambda                     0.3 #used to compensate for divergence points at G=0 in the evaluation of Fock exchange using lcao_in_pw method
exx_real_number                0 #exx calculated in real or complex
exx_pca_threshold              0.0001 #threshold to screen on-site ABFs in exx
exx_c_threshold                0.0001 #threshold to screen C matrix in exx
exx_v_threshold                0.1 #threshold to screen C matrix in exx
exx_dm_threshold               0.0001 #threshold to screen density matrix in exx
exx_cauchy_threshold           1e-07 #threshold to screen exx using Cauchy-Schwartz inequality
exx_c_grad_threshold           0.0001 #threshold to screen nabla C matrix in exx
exx_v_grad_threshold           0.1 #threshold to screen nabla V matrix in exx
exx_cauchy_force_threshold     1e-07 #threshold to screen exx force using Cauchy-Schwartz inequality
exx_cauchy_stress_threshold    1e-07 #threshold to screen exx stress using Cauchy-Schwartz inequality
exx_ccp_rmesh_times            default #how many times larger the radial mesh required for calculating Columb potential is to that of atomic orbitals
exx_opt_orb_lmax               0 #the maximum l of the spherical Bessel functions for opt ABFs
exx_opt_orb_ecut               0 #the cut-off of plane wave expansion for opt ABFs
exx_opt_orb_tolerence          0 #the threshold when solving for the zeros of spherical Bessel functions for opt ABFs

#Parameters (16.tddft)
td_force_dt                    0.02 #time of force change
td_vext                        0 #add extern potential or not
td_vext_dire                   1 #extern potential direction
out_dipole                     0 #output dipole or not
out_efield                     0 #output dipole or not
ocp                            0 #change occupation or not
ocp_set                         #set occupation

#Parameters (17.berry_wannier)
berry_phase                    0 #calculate berry phase or not
gdir                           3 #calculate the polarization in the direction of the lattice vector
towannier90                    0 #use wannier90 code interface or not
nnkpfile                       seedname.nnkp #the wannier90 code nnkp file name
wannier_spin                   up #calculate spin in wannier90 code interface
wannier_method                 1 #different implementation methods under Lcao basis set
out_wannier_mmn                1 #output .mmn file or not
out_wannier_amn                1 #output .amn file or not
out_wannier_unk                1 #output UNK. file or not
out_wannier_eig                1 #output .eig file or not
out_wannier_wvfn_formatted     1 #output UNK. file in text format or in binary format

#Parameters (18.implicit_solvation)
imp_sol                        0 #calculate implicit solvation correction or not
eb_k                           80 #the relative permittivity of the bulk solvent
tau                            1.0798e-05 #the effective surface tension parameter
sigma_k                        0.6 # the width of the diffuse cavity
nc_k                           0.00037 # the cut-off charge density

#Parameters (19.orbital free density functional theory)
of_kinetic                     wt #kinetic energy functional, such as tf, vw, wt
of_method                      tn #optimization method used in OFDFT, including cg1, cg2, tn (default)
of_conv                        energy #the convergence criterion, potential, energy (default), or both
of_tole                        1e-06 #tolerance of the energy change (in Ry) for determining the convergence, default=2e-6 Ry
of_tolp                        1e-05 #tolerance of potential for determining the convergence, default=1e-5 in a.u.
of_tf_weight                   1 #weight of TF KEDF
of_vw_weight                   1 #weight of vW KEDF
of_wt_alpha                    0.833333 #parameter alpha of WT KEDF
of_wt_beta                     0.833333 #parameter beta of WT KEDF
of_wt_rho0                     0 #the average density of system, used in WT KEDF, in Bohr^-3
of_hold_rho0                   0 #If set to 1, the rho0 will be fixed even if the volume of system has changed, it will be set to 1 automaticly if of_wt_rho0 is not zero
of_lkt_a                       1.3 #parameter a of LKT KEDF
of_full_pw                     1 #If set to 1, ecut will be ignored when collect planewaves, so that all planewaves will be used
of_full_pw_dim                 0 #If of_full_pw = true, dimention of FFT is testricted to be (0) either odd or even; (1) odd only; (2) even only
of_read_kernel                 0 #If set to 1, the kernel of WT KEDF will be filled from file of_kernel_file, not from formula. Only usable for WT KEDF
of_kernel_file                 WTkernel.txt #The name of WT kernel file.

#Parameters (20.dft+u)
dft_plus_u                     0 #true:DFT+U correction; false: standard DFT calcullation(default)
yukawa_lambda                  -1 #default:0.0
yukawa_potential               0 #default: false
omc                            0 #the mode of occupation matrix control
hubbard_u           0 #Hubbard Coulomb interaction parameter U(ev)
orbital_corr        -1 #which correlated orbitals need corrected ; d:2 ,f:3, do not need correction:-1

#Parameters (21.spherical bessel)
bessel_nao_ecut                100.000000 #energy cutoff for spherical bessel functions(Ry)
bessel_nao_tolerence           1e-12 #tolerence for spherical bessel root
bessel_nao_rcut                6 #radial cutoff for spherical bessel functions(a.u.)
bessel_nao_smooth              1 #spherical bessel smooth or not
bessel_nao_sigma               0.1 #spherical bessel smearing_sigma
bessel_descriptor_lmax         2 #lmax used in generating spherical bessel functions
bessel_descriptor_ecut         100.000000 #energy cutoff for spherical bessel functions(Ry)
bessel_descriptor_tolerence    1e-12 #tolerence for spherical bessel root
bessel_descriptor_rcut         6 #radial cutoff for spherical bessel functions(a.u.)
bessel_descriptor_smooth       1 #spherical bessel smooth or not
bessel_descriptor_sigma        0.1 #spherical bessel smearing_sigma

#Parameters (22.non-collinear spin-constrained DFT)
sc_mag_switch                  0 #0: no spin-constrained DFT; 1: constrain atomic magnetization
decay_grad_switch              0 #switch to control gradient break condition
sc_thr                         1e-06 #Convergence criterion of spin-constrained iteration (RMS) in uB
nsc                            100 #Maximal number of spin-constrained iteration
nsc_min                        2 #Minimum number of spin-constrained iteration
sc_scf_nmin                    2 #Minimum number of outer scf loop before initializing lambda loop
alpha_trial                    0.01 #Initial trial step size for lambda in eV/uB^2
sccut                          3 #Maximal step size for lambda in eV/uB
sc_file                        none #file name for parameters used in non-collinear spin-constrained DFT (json format)
"""

def acwf_refvolume(token: str):
    refdata = ACWF_REFVOLUME.split("\n")
    for line in refdata:
        if token in line:
            return float(line.split()[1])
    return None
# CITATION
# Bosoni E, Beal L, Bercx M, et al. 
# How to verify the precision of density-functional-theory implementations via 
# reproducible and universal workflows[J]. 
# Nature Reviews Physics, 2024, 6(1): 45-58.
ACWF_REFVOLUME = """
            Ac-X/BCC       45.9436890652
        Ac-X/Diamond      129.7619376361
            Ac-X/FCC       45.5506539230
             Ac-X/SC       49.8290621407
            Ag-X/BCC       17.9815994341
        Ag-X/Diamond       60.1384477535
            Ag-X/FCC       17.8385603821
             Ag-X/SC       20.8078166537
            Al-X/BCC       16.9256568722
        Al-X/Diamond       55.2581279179
            Al-X/FCC       16.4953590598
             Al-X/SC       20.1551394162
            Am-X/BCC       16.1910646994
        Am-X/Diamond       38.9260608838
            Am-X/FCC       17.3637417757
             Am-X/SC       16.1166162914
            Ar-X/BCC       53.3545313857
        Ar-X/Diamond      197.2037257853
            Ar-X/FCC       52.2763833345
             Ar-X/SC       65.2406881909
            As-X/BCC       19.0521510903
        As-X/Diamond       57.0334563896
            As-X/FCC       19.3175775190
             As-X/SC       20.3677497329
            At-X/BCC       40.0072565617
        At-X/Diamond      133.8244131673
            At-X/FCC       39.0306665408
             At-X/SC       46.1440668021
            Au-X/BCC       18.0420836978
        Au-X/Diamond       58.5415097463
            Au-X/FCC       17.9789494007
             Au-X/SC       20.7696027903
             B-X/BCC        6.1393552113
         B-X/Diamond       16.6263169088
             B-X/FCC        5.8915500544
              B-X/SC        6.7004880836
            Ba-X/BCC       63.3054105555
        Ba-X/Diamond      113.1699823796
            Ba-X/FCC       64.1140410061
             Ba-X/SC       61.6056180281
            Be-X/BCC        7.8157560018
        Be-X/Diamond       29.3739872910
            Be-X/FCC        7.8716206739
             Be-X/SC       10.2672356016
            Bi-X/BCC       31.6346675774
        Bi-X/Diamond       96.9022394411
            Bi-X/FCC       31.8104681271
             Bi-X/SC       35.2050582284
            Br-X/BCC       26.7842124781
        Br-X/Diamond       86.1711243924
            Br-X/FCC       26.4175363357
             Br-X/SC       29.8278511781
             C-X/BCC        6.6857220257
         C-X/Diamond       11.3915229468
             C-X/FCC        7.3216347712
              C-X/SC        5.5822112996
            Ca-X/BCC       42.1504734668
        Ca-X/Diamond      159.9693954411
            Ca-X/FCC       42.1942968575
             Ca-X/SC       43.5849131941
            Cd-X/BCC       23.4195726366
        Cd-X/Diamond       74.9080478398
            Cd-X/FCC       22.8412781751
             Cd-X/SC       26.9236002457
            Ce-X/BCC       27.3239643916
        Ce-X/Diamond       60.3696280203
            Ce-X/FCC       26.5223557737
             Ce-X/SC       24.9210139881
            Cl-X/BCC       21.4547444716
        Cl-X/Diamond       67.5148577560
            Cl-X/FCC       21.2881895239
             Cl-X/SC       23.4640292893
            Cm-X/BCC       16.5207822537
        Cm-X/Diamond       38.6523626094
            Cm-X/FCC       17.4924664960
             Cm-X/SC       16.3990024768
            Co-X/BCC       10.5447926821
        Co-X/Diamond       29.7694335456
            Co-X/FCC       10.3083927827
             Co-X/SC       11.8932667993
            Cr-X/BCC       11.5481993396
        Cr-X/Diamond       33.0709117517
            Cr-X/FCC       11.8859164649
             Cr-X/SC       12.8069580127
            Cs-X/BCC      116.8417228518
        Cs-X/Diamond      377.5119512083
            Cs-X/FCC      117.3605789432
             Cs-X/SC      128.3590552257
            Cu-X/BCC       12.0045162940
        Cu-X/Diamond       38.3543187692
            Cu-X/FCC       11.9522246523
             Cu-X/SC       13.9349731242
            Dy-X/BCC       32.2885692044
        Dy-X/Diamond       46.0203764469
            Dy-X/FCC       32.4766359634
             Dy-X/SC       31.8571930831
            Er-X/BCC       33.9307067591
        Er-X/Diamond      160.7182113742
            Er-X/FCC       34.8228848521
             Er-X/SC       35.9430359073
            Eu-X/BCC       26.1340832200
        Eu-X/Diamond       41.3650702486
            Eu-X/FCC       24.9919769933
             Eu-X/SC       17.7951482963
             F-X/BCC       10.0841273303
         F-X/Diamond       29.0048013514
             F-X/FCC       10.1469872576
              F-X/SC       10.5209419861
            Fe-X/BCC       10.5004841964
        Fe-X/Diamond       28.9244448081
            Fe-X/FCC       10.2602117762
             Fe-X/SC       11.6513477719
            Fr-X/BCC      116.4923124868
        Fr-X/Diamond      384.0321596406
            Fr-X/FCC      117.1629536099
             Fr-X/SC      132.1740858506
            Ga-X/BCC       19.2055876101
        Ga-X/Diamond       50.8380451983
            Ga-X/FCC       18.9465000846
             Ga-X/SC       20.1173121537
            Gd-X/BCC       28.9474037903
        Gd-X/Diamond       41.8942701551
            Gd-X/FCC       27.9938704976
             Gd-X/SC       20.8095989335
            Ge-X/BCC       19.2694953095
        Ge-X/Diamond       47.8265898364
            Ge-X/FCC       19.5824508015
             Ge-X/SC       19.9417483114
             H-X/BCC        2.9667742929
         H-X/Diamond        6.8311116562
             H-X/FCC        2.9648091126
              H-X/SC        3.0866686939
            He-X/BCC       18.0303880567
        He-X/Diamond       64.2148789366
            He-X/FCC       17.7725807210
             He-X/SC       21.4852711323
            Hf-X/BCC       22.3047194703
        Hf-X/Diamond       70.1577406488
            Hf-X/FCC       22.5676034192
             Hf-X/SC       24.7811393086
            Hg-X/BCC       29.2371602125
        Hg-X/Diamond      112.5823446507
            Hg-X/FCC       32.3477976932
             Hg-X/SC       29.8543233498
            Ho-X/BCC       33.2672852865
        Ho-X/Diamond       50.8501121827
            Ho-X/FCC       33.8917393870
             Ho-X/SC       34.2940697226
             I-X/BCC       35.9867670342
         I-X/Diamond      121.1452775012
             I-X/FCC       35.1048845030
              I-X/SC       41.5631544369
            In-X/BCC       27.7805802898
        In-X/Diamond       76.4320586824
            In-X/FCC       27.5100988618
             In-X/SC       29.5502109822
            Ir-X/BCC       15.0556415096
        Ir-X/Diamond       43.1935547990
            Ir-X/FCC       14.5049941291
             Ir-X/SC       16.9945041257
             K-X/BCC       73.7795179893
         K-X/Diamond      223.8044852857
             K-X/FCC       74.0044361285
              K-X/SC       79.3544129228
            Kr-X/BCC       67.4634311123
        Kr-X/Diamond      249.7079835831
            Kr-X/FCC       66.0417653419
             Kr-X/SC       82.6645424965
            La-X/BCC       37.8175702599
        La-X/Diamond       74.6380534648
            La-X/FCC       36.9468915219
             La-X/SC       36.7450629470
            Li-X/BCC       20.2674652133
        Li-X/Diamond       51.3867180352
            Li-X/FCC       20.2244785929
             Li-X/SC       20.4088019741
            Lu-X/BCC       29.6257064032
        Lu-X/Diamond      101.2182490481
            Lu-X/FCC       28.9714042450
             Lu-X/SC       32.9400269740
            Mg-X/BCC       22.9172569645
        Mg-X/Diamond       80.8516729950
            Mg-X/FCC       23.1252461440
             Mg-X/SC       27.5805850687
            Mn-X/BCC       10.7808268397
        Mn-X/Diamond       30.3432739437
            Mn-X/FCC       10.7471443132
             Mn-X/SC       11.8988403344
            Mo-X/BCC       15.7925810470
        Mo-X/Diamond       46.0075591929
            Mo-X/FCC       16.0351302591
             Mo-X/SC       17.5962610036
             N-X/BCC        7.2347133844
         N-X/Diamond       18.3533073287
             N-X/FCC        7.6012909873
              N-X/SC        6.4799301718
            Na-X/BCC       37.0150641100
        Na-X/Diamond      109.1430788906
            Na-X/FCC       37.0989860953
             Na-X/SC       39.7513602453
            Nb-X/BCC       18.1415011497
        Nb-X/Diamond       51.5612380786
            Nb-X/FCC       18.7677989613
             Nb-X/SC       20.1149109148
            Nd-X/BCC       21.0678209810
        Nd-X/Diamond       47.0991055430
            Nd-X/FCC       22.7647561474
             Nd-X/SC       18.0765197733
            Ne-X/BCC       24.7112468179
        Ne-X/Diamond       89.1474187039
            Ne-X/FCC       24.3028219458
             Ne-X/SC       29.7451962147
            Ni-X/BCC       10.8950930618
        Ni-X/Diamond       33.0121380345
            Ni-X/FCC       10.8348969207
             Ni-X/SC       12.5590185628
            Np-X/BCC       17.8079227825
        Np-X/Diamond       42.9533016321
            Np-X/FCC       19.2945468597
             Np-X/SC       17.2734160816
             O-X/BCC        7.7862826653
         O-X/Diamond       21.3628964776
             O-X/FCC        7.9987519984
              O-X/SC        7.9535032234
            Os-X/BCC       14.7808432251
        Os-X/Diamond       42.9003271345
            Os-X/FCC       14.3408590455
             Os-X/SC       16.7303790128
             P-X/BCC       14.2301643213
         P-X/Diamond       41.3183089373
             P-X/FCC       14.5635813315
              P-X/SC       14.6564493780
            Pa-X/BCC       24.7971150304
        Pa-X/Diamond       61.0193764146
            Pa-X/FCC       25.2979382957
             Pa-X/SC       24.0207081049
            Pb-X/BCC       31.9704284094
        Pb-X/Diamond       88.0899664601
            Pb-X/FCC       32.0330693820
             Pb-X/SC       34.4758419297
            Pd-X/BCC       15.4442899164
        Pd-X/Diamond       49.0125255975
            Pd-X/FCC       15.3254220778
             Pd-X/SC       17.8613427336
            Pm-X/BCC       20.3578671427
        Pm-X/Diamond       43.3046763305
            Pm-X/FCC       22.2454982473
             Pm-X/SC       17.2982038540
            Po-X/BCC       32.8538744328
        Po-X/Diamond      104.9634814464
            Po-X/FCC       32.5633304447
             Po-X/SC       37.5954709759
            Pr-X/BCC       23.1413663276
        Pr-X/Diamond       52.4401030751
            Pr-X/FCC       24.0941010592
             Pr-X/SC       20.1508127695
            Pt-X/BCC       15.8389750612
        Pt-X/Diamond       48.2299658537
            Pt-X/FCC       15.6559513760
             Pt-X/SC       18.0860881007
            Pu-X/BCC       16.5643372636
        Pu-X/Diamond       40.4075680763
            Pu-X/FCC       17.8021327417
             Pu-X/SC       16.3671826034
            Ra-X/BCC       70.9668672681
        Ra-X/Diamond      339.3371556818
            Ra-X/FCC       71.6269845364
             Ra-X/SC       75.3455997294
            Rb-X/BCC       91.1440710492
        Rb-X/Diamond      282.7553207612
            Rb-X/FCC       91.4274679406
             Rb-X/SC       98.9794397816
            Re-X/BCC       15.1044598171
        Re-X/Diamond       45.1488092366
            Re-X/FCC       15.0163498073
             Re-X/SC       17.1611268525
            Rh-X/BCC       14.4742145663
        Rh-X/Diamond       41.9004577066
            Rh-X/FCC       14.0504550878
             Rh-X/SC       16.3052925485
            Rn-X/BCC       95.4471230782
        Rn-X/Diamond      353.9049285433
            Rn-X/FCC       93.1564256810
             Rn-X/SC      117.6672115419
            Ru-X/BCC       14.2356061964
        Ru-X/Diamond       40.5863643979
            Ru-X/FCC       13.8366151623
             Ru-X/SC       15.8381291955
             S-X/BCC       15.7618462400
         S-X/Diamond       48.5621194269
             S-X/FCC       15.8806885358
              S-X/SC       17.2197235453
            Sb-X/BCC       27.2257889203
        Sb-X/Diamond       85.3650954131
            Sb-X/FCC       27.4896291741
             Sb-X/SC       30.0634121990
            Sc-X/BCC       24.8858518758
        Sc-X/Diamond       68.9151380008
            Sc-X/FCC       24.6868520662
             Sc-X/SC       26.1480548969
            Se-X/BCC       20.3600286206
        Se-X/Diamond       63.5107554328
            Se-X/FCC       20.3778732432
             Se-X/SC       22.6834109000
            Si-X/BCC       14.6451673299
        Si-X/Diamond       40.9149469095
            Si-X/FCC       14.4822031565
             Si-X/SC       16.2292642666
            Sm-X/BCC       21.6458724878
        Sm-X/Diamond       41.8054697730
            Sm-X/FCC       22.8284101769
             Sm-X/SC       17.1965125100
            Sn-X/BCC       27.6472928166
        Sn-X/Diamond       73.6877592925
            Sn-X/FCC       28.0088425653
             Sn-X/SC       29.4530306864
            Sr-X/BCC       54.0128807068
        Sr-X/Diamond      223.7821246728
            Sr-X/FCC       54.8922650546
             Sr-X/SC       57.1290143708
            Ta-X/BCC       18.2919921205
        Ta-X/Diamond       56.8185715387
            Ta-X/FCC       18.8393711237
             Ta-X/SC       20.6669151965
            Tb-X/BCC       30.9008337902
        Tb-X/Diamond       43.3929318627
            Tb-X/FCC       30.5520741878
             Tb-X/SC       27.8199432306
            Tc-X/BCC       14.6195757062
        Tc-X/Diamond       42.4918607430
            Tc-X/FCC       14.5126740985
             Tc-X/SC       16.2350067792
            Te-X/BCC       28.5152914419
        Te-X/Diamond       92.8002629931
            Te-X/FCC       28.2787257642
             Te-X/SC       32.7934878779
            Th-X/BCC       32.5677469591
        Th-X/Diamond       91.4653155955
            Th-X/FCC       32.1839319361
             Th-X/SC       35.3321431170
            Ti-X/BCC       17.2668387055
        Ti-X/Diamond       45.8421534674
            Ti-X/FCC       17.3945872331
             Ti-X/SC       18.4131750886
            Tl-X/BCC       31.4143781862
        Tl-X/Diamond       90.2879174095
            Tl-X/FCC       31.1402049921
             Tl-X/SC       34.3940693000
            Tm-X/BCC       34.3583819751
        Tm-X/Diamond      163.3512533799
            Tm-X/FCC       35.3320901938
             Tm-X/SC       37.2077575159
             U-X/BCC       20.2661778805
         U-X/Diamond       49.5199753572
             U-X/FCC       21.7133355135
              U-X/SC       19.1223788669
             V-X/BCC       13.4607654817
         V-X/Diamond       37.2400616134
             V-X/FCC       13.9049548185
              V-X/SC       14.6801017804
             W-X/BCC       16.1454752416
         W-X/Diamond       49.5504586465
             W-X/FCC       16.4578846361
              W-X/SC       18.3694598023
            Xe-X/BCC       89.0348898020
        Xe-X/Diamond      331.8911738830
            Xe-X/FCC       87.0069651578
             Xe-X/SC      109.8051846381
             Y-X/BCC       33.0304565057
         Y-X/Diamond       87.5187101972
             Y-X/FCC       32.4715444515
              Y-X/SC       34.8169697264
            Yb-X/BCC       34.6395882212
        Yb-X/Diamond      164.1166906281
            Yb-X/FCC       35.7043529273
             Yb-X/SC       38.6416396501
            Zn-X/BCC       15.3751088074
        Zn-X/Diamond       49.3473967449
            Zn-X/FCC       15.1620084721
             Zn-X/SC       18.1790323898
            Zr-X/BCC       22.8447760184
        Zr-X/Diamond       61.9100700077
            Zr-X/FCC       23.2133785116
             Zr-X/SC       24.7431321708
"""

import unittest
class TestAbacusGeneration(unittest.TestCase):
    def test_kspacing(self):
        # provide cell parameter in Angstrom
        param = {"cell": [4.22798145, 4.22798145, 4.22798145, 60, 60, 60]}
        # but kspacing is in Bohr-1
        result = KSPACING(kspacing=0.03, kspc_unit="Bohr", **param)
        self.assertEqual(result, "K_POINTS\n0\nGamma\n33 33 33 0 0 0\n")

if __name__ == "__main__":
    unittest.main()