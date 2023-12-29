from apns.module_structure.crystal_information_file import cell_parameters_to_lattice_vectors as cptlv
from apns.module_structure.crystal_information_file import read_1 as read
import re

def STRU_cif(fname: str, 
             pseudopotentials: dict = {}, 
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
    cif = read(fname)
    return_str = "ATOMIC_SPECIES\n"
    for _element in pseudopotentials:
        mass = 1.0
        if "mass" in kwargs:
            if _element in kwargs["mass"]:
                mass = kwargs["mass"][_element]
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
    lattice_vectors = cptlv(cell_parameters)

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
        lattice_vectors = cptlv(kwargs["cell_parameters"])
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

def KPT_generation(mode: str, gamma_centered: bool = True, **kwargs):
    """KPT file generation

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
    if mode == "isolated":
        if gamma_centered:
            return_str = "K_POINTS automatic\n0\nGamma\n1 1 1 0 0 0\n"
        else:
            raise NotImplementedError
    elif mode == "crystal":
        _fold = 35
        if gamma_centered:
            if "cell_parameters" in kwargs:
                cell_parameters = kwargs["cell_parameters"]
                if isinstance(cell_parameters, list):
                    a = cell_parameters[0]
                    b = cell_parameters[1]
                    c = cell_parameters[2]
                elif isinstance(cell_parameters, dict):
                    a = cell_parameters["a"]
                    b = cell_parameters["b"]
                    c = cell_parameters["c"]
                else:
                    raise TypeError("Not supported data format of cell parameters, must be list or dict.")
                if "metallic" in kwargs:
                    if kwargs["metallic"]:
                        _fold = 35
                    else:
                        _fold = 25
                else:
                    _fold = 25
                return_str = "K_POINTS\n0\nGamma\n%d %d %d 0 0 0\n" % (int(_fold/a) + 1, int(_fold/b) + 1, int(_fold/c) + 1)
        else:
            raise NotImplementedError
    else:
        raise ValueError
    
    return return_str

def INPUT(work_status: dict,
          minimal: bool = True,
          template: bool = False):
    """new version of function INPUT for generating INPUT file of ABACUS, if in work_status
    find in additional_keywords there are parameter set as xxx_to_test """
    minimal_keys = {
        "calculation": "scf",
        "ecutwfc": 100,
        "basis_type": "pw",
        "symmetry": 0,
        "pseudo_dir": "./",
        "scf_nmax": 100,
        "scf_thr": 1e-6
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
    if template:
        values["dft_functional"] = "functional_to_test"
        values["ecutwfc"] = "ecutwfc_to_test"
        for key in work_status["additional_keywords"].keys():
            if isinstance(work_status["additional_keywords"][key], list):
                values[key] = key + "_to_test"
            else:
                values[key] = work_status["additional_keywords"][key]
    # for minimal, only keep the minimal keys, but add manually set keys
    for key in work_status["additional_keywords"].keys():
        if key not in minimal_keys:
            minimal_keys[key] = values[key]
    # write.
    return_str = "INPUT_PARAMETERS\n"
    if minimal:
        for key in minimal_keys:
            return_str += "%s %s # %s\n"%(key, values[key], comments[key])
    else:
        for key in values:
            return_str += "%s %s # %s\n"%(key, values[key], comments[key])
    return return_str

def INPUT_generation(minimal: bool = True,
                     template: bool = False,
                     **kwargs):
    
    minimal_keys = {
        "calculation": "scf",
        "ecutwfc": 100,
        "basis_type": "pw",
        "symmetry": 0,
        "pseudo_dir": "./",
        "scf_nmax": 100,
        "scf_thr": 1e-6
    }
    values = {}
    comments = {}
    for line in INPUT_TEMPLATE.split("\n"):
        if line == "INPUT_PARAMETERS" or line.startswith("#"):
            continue
        # delete contents all behind #
        words = line.split("#")
        line = words[0]
        try:
            comment = words[1]
        except IndexError:
            comment = ""
        words = line.split()
        if len(words) == 0:
            continue
        key = words[0]
        value = " ".join(words[1:])
        values[key] = value
        comments[key] = comment
    for key in kwargs:
        if key not in values:
            comments[key] = ""
        values[key] = kwargs[key]

    if template:
        values["dft_functional"] = "functional_to_test"
        values["ecutwfc"] = "ecutwfc_to_test"
        
    return_str = "INPUT_PARAMETERS\n"
    for key in kwargs:
        if key not in minimal_keys:
            minimal_keys[key] = values[key]

    if minimal:
        for key in minimal_keys:
            return_str += "%s %s # %s\n"%(key, values[key], comments[key])
    else:
        for key in values:
            return_str += "%s %s # %s\n"%(key, values[key], comments[key])
    return return_str

INPUT_TEMPLATE = """
INPUT_PARAMETERS
#Parameters (1.General)
suffix                         ABACUS #the name of main output directory
latname                        none #the name of lattice name
stru_file                      STRU #the filename of file containing atom positions
kpoint_file                    KPT #the name of file containing k points
pseudo_dir                     ./ #the directory containing pseudo files
orbital_dir                    ./ #the directory containing orbital files
pseudo_rcut                    15 #cut-off radius for radial integration
pseudo_mesh                    0 #0: use our own mesh to do radial renormalization; 1: use mesh as in QE
dft_functional                 default #exchange correlation functional
xc_temperature                 0 #temperature for finite temperature functionals
calculation                    scf #test; scf; relax; nscf; get_wf; get_pchg
esolver_type                   ksdft #the energy solver: ksdft, sdft, ofdft, tddft, lj, dp
ntype                          1 #atom species number
nspin                          1 #1: single spin; 2: up and down spin; 4: noncollinear spin
kspacing                       0 0 0  #unit in 1/bohr, should be > 0, default is 0 which means read KPT file
min_dist_coef                  0.2 #factor related to the allowed minimum distance between two atoms
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
diago_proc                     16 #the number of procs used to do diagonalization
nbspline                       -1 #the order of B-spline basis
wannier_card                   none #input card for wannier functions
soc_lambda                     1 #The fraction of averaged SOC pseudopotential is given by (1-soc_lambda)
cal_force                      1 #if calculate the force at the end of the electronic iteration
out_freq_ion                   0 #the frequency ( >= 0 ) of ionic step to output charge density and wavefunction. 0: output only when ion steps are finished
device                         cpu #the computing device for ABACUS

#Parameters (2.PW)
ecutwfc                        100 ##energy cutoff for wave functions
ecutrho                        800 ##energy cutoff for charge density and potential
erf_ecut                       0 ##the value of the constant energy cutoff
erf_height                     0 ##the height of the energy step for reciprocal vectors
erf_sigma                      0.1 ##the width of the energy step for reciprocal vectors
fft_mode                       0 ##mode of FFTW
pw_diag_ndim                   4 #max dimension for davidson
pw_diag_thr                    0.01 #threshold for eigenvalues is cg electron iterations
scf_thr                        1e-09 #charge density error
scf_thr_type                   1 #type of the criterion of scf_thr, 1: reci drho for pw, 2: real drho for lcao
init_wfc                       atomic #start wave functions are from 'atomic', 'atomic+random', 'random' or 'file'
init_chg                       atomic #start charge is from 'atomic' or file
chg_extrap                     first-order #atomic; first-order; second-order; dm:coefficients of SIA
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
scf_nmax                       200 ##number of electron iterations
relax_nmax                     50 #number of ion iteration steps
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
cal_stress                     0 #calculate the stress or not
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
basis_type                     pw #PW; LCAO in pw; LCAO
gamma_only                     0 #Only for localized orbitals set and gamma point. If set to 1, a fast algorithm is used
search_radius                  -1 #input search radius (Bohr)
search_pbc                     1 #input periodic boundary condition
lcao_ecut                      0 #energy cutoff for LCAO
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
bx                             1 #division of an element grid in FFT grid along x
by                             1 #division of an element grid in FFT grid along y
bz                             1 #division of an element grid in FFT grid along z

#Parameters (6.Smearing)
smearing_method                gauss #type of smearing_method: gauss; fd; fixed; mp; mp2; mv
smearing_sigma                 0.05 #energy range for smearing

#Parameters (7.Charge Mixing)
mixing_type                    broyden #plain; pulay; broyden
mixing_beta                    0.8 #mixing parameter: 0 means no new charge
mixing_ndim                    8 #mixing dimension in pulay or broyden
mixing_gg0                     1 #mixing parameter in kerker
mixing_beta_mag                -10 #mixing parameter for magnetic density
mixing_gg0_mag                 0 #mixing parameter in kerker
mixing_gg0_min                 0.1 #the minimum kerker coefficient
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
bessel_nao_ecut                200.000000 #energy cutoff for spherical bessel functions(Ry)
bessel_nao_tolerence           1e-12 #tolerence for spherical bessel root
bessel_nao_rcut                6 #radial cutoff for spherical bessel functions(a.u.)
bessel_nao_smooth              1 #spherical bessel smooth or not
bessel_nao_sigma               0.1 #spherical bessel smearing_sigma
bessel_descriptor_lmax         2 #lmax used in generating spherical bessel functions
bessel_descriptor_ecut         200.000000 #energy cutoff for spherical bessel functions(Ry)
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
