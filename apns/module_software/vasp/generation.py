import apns.module_structure.crystal_information_file as amscif
import seekpath as skps
import apns.module_database.database as amdd

def KPOINTS(fname: str = "", nkpoints_in_line: int = 0):

    if fname == "":
        raise ValueError("fname should not be empty")

    cif = amscif.read_1(fname)
    cell_parameters = amscif.cellparam_to_latvec(cif["cell_parameters"])

    result = ""
    if nkpoints_in_line > 0:
        cell_vectors = amscif.cellparam_to_latvec(cell_parameters)
        positions = []
        numbers = []

        for element in cif["atomic_positions"].keys():
            for line in cif["atomic_positions"][element]:
                positions.append(line)
                numbers.append(amdd.get_element_index(element))
        
        _skps_result = skps.get_path(structure=(cell_vectors, positions, numbers))
        result = _skps_result
    elif nkpoints_in_line == 0:
        
        result = "Automatic generation\n0\nMonkhorst-pack\n0 0 0\n0 0 0"
    else:
        result = "Automatic generation\n0\nMonkhorst-pack\n1 1 1\n0 0 0"
    
    return result

def POSCAR(fname: str = "", constraints: list = None):

    if fname == "":
        raise ValueError("fname should not be empty")
    
    cif = amscif.read_1(fname)
    cell_vectors = amscif.cellparam_to_latvec(cif["cell_parameters"])

    natom_per_type = {}
    for element in cif["atomic_positions"].keys():
        natom_per_type[element] = len(cif["atomic_positions"][element])
    
    result = "APNS Auto-generated POSCAR\n1.0\n"
    for vector in cell_vectors:
        result += "%10.8f %10.8f %10.8f\n" % (vector[0], vector[1], vector[2])
    
    _line1 = ""
    _line2 = ""
    for element in natom_per_type.keys():
        _line1 += element + " "
        _line2 += str(natom_per_type[element]) + " "
    result += _line1 + "\n" + _line2 + "\n"
    result += "direct\n"
    for element in cif["atomic_positions"].keys():
        for line in cif["atomic_positions"][element]:
            result += "%10.8f %10.8f %10.8f %s\n" % (line[0], line[1], line[2], element)

    return result

INPUT_TEMPLATE = """!APNS VASP INCAR auto-generation!
!-- always necessary --!
 ENCUT = 550.0 ! eV, ENMAX*1.3 for structural opt.
 ISMEAR = -5 ! 0: Gaussian, 1: MP, -5: tetrahedron
! SIGMA = 0.1 ! eV, specify this for ISMEAR=0,1. 0.1--0.2 eV.

!-- calculation modes --!
!! structural opt.
! NSW = 50
! IBRION = 2
! ISIF = 3 ! 2: fixed lattice & relax ions, 3: all opt
! EDIFFG = 0 ! +:Energy diff, -:Force, 0:stop after NSW iterations
! PSTRESS = 10.0 ! pressure in kB

!! band calc.
! ISTART = 0 ! SCF orbitals are not used
! ICHARG = 11 ! SCF charge density is used (fixed)
! ISMEAR = 0 ! tetrahedron does not work for band calc.
! SIGMA = 0.1

!! input tags for LDOS calc. are shown later

!-- functionals --!
 GGA = PE
! VOSKOWN = 1 ! when using PW91

!! HSE06. initial GGA (PBE) run is recommended before HSE calc.
! LHFCALC = .TRUE.
! HFSCREEN = 0.2
! ICHARG = 0 ! restart from PBE orbitals
! ISTART = 1 ! restart from PBE orbitals
! PRECFOCK = Fast ! (optional) accuracy and speed control. For old versions, ENCUTFOCK = 0.
! NKRED = 2 ! (optional) lower accuracy but fast. Instead, NKREDX,NKREDY,NKREDZ can be used.

!! +U parameter
! LDAU = .TRUE.
! LDAUTYPE = 2
! LDAUL = -1 2 -1 ! l-quantum number (-1:noU, 2:d, 3:f)
! LDAUU = 0.0  5.0  0.0 ! U values (eV)
! LDAUJ = 0.0  0.0  0.0 ! J values (For LDAUTYPE=2, Ueff=U-J)
! LMAXMIX = 4 ! 4 for d-electrons, 6 for f-electrons.
! LDAUPRINT = 2

!-- spin (incl. spin-orbit coupling) --!
!! no spin pol. & no SOC
 ISPIN = 1

!! spin pol. & no SOC
! ISPIN = 2
! ISYM = -1 
! LSORBIT = .FALSE.
! LNONCOLLINEAR = .FALSE.
! ICHARG = 2 ! set ICHARG=2 to use MAGMOM
! MAGMOM = 5*0.0 ! scalar, 5 atoms

!! SOC calc (vasp_ncl)
! ISPIN = 2
! ISYM = -1
! LSORBIT = .TRUE.
! LNONCOLLINEAR = .TRUE.
! SAXIS = 0 0 1 ! quantization axis
! LORBMOM = .TRUE.
!! (case 1) SOC calc using vasp_ncl after non-SOC (but spin-polarized) calc using vasp_std (i.e. use non-SOC orbitals as initial guess)
! ISTART = 1
! NBANDS = 100 ! (NBANDS in non-SOC run) x 2. see OUTCAR in non-SOC run. (not INCAR because NBANDS can be changed by vasp run)
!! (case 2) SOC calc without initial orbitals
! ICHARG = 2 ! set ICHARG=2 to use MAGMOM
! MAGMOM = 15*0.0 ! spinor, 5 atoms.

!-- accuracy settings (usually no need to change) --!
 EDIFF = 1E-8 ! can be a larger value
 LREAL = .FALSE. ! LREAL = Auto is recommended for a large system (to reduce cost)
 ADDGRID = .TRUE. ! improves accuracy but can result in (small) negative density in real space. ADDGRID = .FALSE. is default
 GGA_COMPAT = .FALSE.
 PREC = Accurate
 LASPH = .TRUE.

!! *** other tags *** !!

!! output: ldos, band weight, local magnetic moment, etc.
 LORBIT = 11 ! optional but recommended

!! parallelization control (optional)
! NPAR = 32 ! (default: no. of cores) band parallelization
! NCORE = 32 ! (default: 1) plane-wave parallelization
! KPAR = 1 ! (default: 1) k-point parallelization, not available in old versions

!! to improve convergence (e.g. for magnetic systems, SOCcalc,...)
! ALGO = A ! CG method (default: blocked Davidson). For metal: ALGO=D with an appropriate TIME is recommended.
! AMIX = 0.1
! BMIX = 0.00001
! AMIX_MAG = 0.2
! BMIX_MAG = 0.00001

!! wannier
!! First run: get orbitals with parallelization, Second run: restart with ISTART=1 and NPAR=N (NCORE should be 1)
!! First & Second runs
! ISYM = -1
! NBANDS = 200 ! sufficiently large value to get unoccupied bands for wannierization
!! Second run
! ISTART = 1
! NELMIN = 6 ! set a large value (5--10) to get well-converged unoccupied bands
! LWANNIER90 = .TRUE.
! LWANNIER90_RUN = .TRUE.
! LWRITE_MMN_AMN = .TRUE.

!! LDOS (LORBIT=11 necessary)
!! assuming LDOS calc using a fine k-mesh after an SCF run
! ISYM = -1
! NPAR = 1 ! version-dependent. see VASPwiki.
! ISTART = 0
! ICHARG = 11 ! use an SCF charge density
! ISMEAR = -5 ! tetrahedron recommended
! SIGMA = 0.05
! NEDOS = 3000
! EMIN = -0.5
! EMAX = 14.5

!! charged-state calc.
! NELECT = 1224

!! static dielectric properties
! ISMEAR = 0
! SIGMA = 0.01
! LEPSILON = .TRUE. !! alternatively LCALCEPS=.TRUE. can be used (e.g. for hybrid functionals)
! LPEAD = .TRUE.
! IBRION = 8
! LRPA = .FALSE. !! .TRUE.: the local field effect is calculated at the Hartree level only
"""

if __name__ == "__main__":
    print(KPOINTS("apns_cache/mp-2460.cif"))