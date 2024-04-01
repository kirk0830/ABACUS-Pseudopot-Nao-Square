      
"""
Portable ABACUS output files parser

Author: kirk0830

Presently support:
------------------------------------------------------------------------------------------
file                    ABACUS keyword              function to use             tested
------------------------------------------------------------------------------------------
data-H/S-0              out_mat_hs 1                read_mat_hs                 yes
LOWF_K_0.txt            out_wfc_lcao 1              read_lowf                   yes
QO_ovlp_0.dat           qo_switch 1                 read_ao_proj                yes
QO_supercells.dat       qo_swtich 1                 read_qo_supercells          no
kpoints                                             read_kpoints                yes
istate.info                                         read_istate                 yes
BANDS_1.dat             out_band                    read_bands                  no
time.json                                           read_time_statistics        no
running_*.scf                                       read_etraj_fromlog          no
------------------------------------------------------------------------------------------
with utility functions:
------------------------------------------------------------------------------------------
function                description                 tested
------------------------------------------------------------------------------------------
cxx_topycomplex         convert cxx complex to py   yes
unit_conversion         convert units               yes
Gauss_smoothing         Gaussian smoothing          no
zero_padding            zero padding                yes
------------------------------------------------------------------------------------------
directly run will trigger unittest on functions
"""

from itertools import groupby
import numpy as np

"""utils"""
def cxx_topycomplex(num: str):
    """cxx prints complex numbers in the form (a,b)"""
    num = num.replace('(', '').replace(')', '')
    num = num.split(',')
    return complex(float(num[0]), float(num[1]))

def unit_conversion(val: float, unitfrom: str = "eV", unitto: str = "eV") -> float:
    # energy unit takes eV as the base
    eunit = {"eV": 1.0, "Ry": 1/13.6, "Hartree": 1/27.2, "meV": 1e3, "J": 1.6e-19, "kcal/mol": 23.1}
    # length unit takes Angstrom as the base
    lunit = {"A": 1.0, "m": 1e-10, "nm": 0.1, "Bohr": 1.8897259886, "a.u.": 1.8897259886}
    # weight unit takes atomic unit as the base
    wtunit = {"a.u.": 1.0, "kg": 9.1e-31, "g": 9.1e-28, "me": 0.00054858}

    units = [eunit, lunit, wtunit]
    # both unit should belong to the same category
    for u in units:
        if unitfrom in u and unitto in u:
            return val * u[unitto] / u[unitfrom]
    raise ValueError("Units are not in the same category")

def Gauss_smoothing(x, y, sigma, normalize = True, normalize_to = 1.0):
    """Gaussian smoothing of the data y(x) with a standard deviation sigma.

    Args:
        x (np.ndarray): x-axis values.
        y (np.ndarray): y-axis values.
        sigma (float): standard deviation of the Gaussian.

    Returns:
        np.ndarray: smoothed y-axis values.
    """
    n = len(x)
    m = len(y)
    assert n == m
    y_smoothed = np.zeros(n)
    for i in range(n):
        for j in range(n):
            y_smoothed[i] += y[j]*np.exp(-(x[i] - x[j])**2/(2*sigma**2))
    if normalize:
        norm = np.trapz(y_smoothed, x)
        y_smoothed = y_smoothed * normalize_to / norm
            
    return y_smoothed

def zero_padding(xmin: float, xmax: float, dx: float, x: np.ndarray, y: np.ndarray):
    """Zero padding for the x and y data, so that the x data is within the range [xmin, xmax] with a step of dx.
    The y data is interpolated to the new x data.
    
    Args:
        xmin (float): minimum value of the new x data
        xmax (float): maximum value of the new x data, EXCLUSIVE!
        dx (float): step of the new x data
        x (np.ndarray): original x data
        y (np.ndarray): original y data
        interpolate (str): interpolation method, "none", "linear", "cubic"
    
    Returns:
        np.ndarray: new x data
        np.ndarray: new y data
    """
    xnew = np.arange(xmin, xmax, dx)
    ynew = np.zeros(len(xnew))

    def nearest_index(x, xnew):
        """return the nearest index of x in xnew, and the difference between x and xnew[index]"""
        i = np.abs(xnew - x).argmin()
        return i, xnew[i] - x
    
    for i in range(len(xnew)):
        if xnew[i] < x[0] or xnew[i] > x[-1]:
            continue
        _i, d = nearest_index(xnew[i], x)
        if abs(d) < dx:
            ynew[i] = y[_i]
    return xnew, ynew

import matplotlib.pyplot as plt
def draw_stackdos(e: list, dos: list,           # source data, stack into list
                  ncols: int = 1, 
                  suptitle: str = None,
                  tight_xrange: bool = True,
                  tight_yrange: bool = False,
                  interactive: bool = True,
                  ):
    assert len(e) == len(dos)
    ndos = len(dos)
    nrows = ndos // ncols
    nrows += 1 if abs(ndos - nrows*ncols) > 1e-3 else 0
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False)

    idos = 0
    for irow in range(nrows):
        for icol in range(ncols):
            axs[irow, icol].plot(e[idos], dos[idos])
            if tight_xrange:
                axs[irow, icol].set_xlim(np.min(e[idos]), np.max(e[idos]))
            if tight_yrange:
                axs[irow, icol].set_ylim(0, np.max(dos[idos]))
            idos += 1
    
    plt.axhline(0, linestyle = "-", lw=0.2, color="black")
    if suptitle is not None:
        plt.suptitle(suptitle)
    if interactive:
        plt.show()
    else:
        plt.savefig("draw_stackdos.png")

"""File: data-H/S-ik"""
def read_mat_hs(fhs):
    """read data-H/S-ik files from abacus for out_mat_hs 1"""
    with open(fhs, 'r') as f:
        data = "".join(f.readlines()).replace('\n', ' ').split()
    size = int(data[0])
    indexing = []
    for i in range(size):
        indexing += [(i, j) for j in range(i, size)]
    data = data[1:]
    mat = np.zeros((size, size), dtype=np.complex128)
    for inum, number in enumerate(data):
        mat[indexing[inum]] = cxx_topycomplex(number)
    mat = mat + mat.conj().T
    for i in range(size):
        mat[i, i] = mat[i, i] / 2
    return mat
"""File: LOWF_K_ik.txt"""
def read_lowf(flowf):
    """read the lowf matrix in k-space, and the k-vector.
    
    Args:
        flowf (str): path to the file containing the lowf matrix
    
    Returns:
        np.ndarray: lowf matrix in k-space
        np.ndarray: k-vector
    """
    nband, nlocal = 0, 0
    with open(flowf, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    for line in lines:
        if line.endswith("(number of bands)"):
            nband = int(line.split()[0])
        elif line.endswith("(number of orbitals)"):
            nlocal = int(line.split()[0])
        else:
            continue
    assert nband > 0
    assert nlocal > 0

    lowf_k = np.zeros((nlocal, nband), dtype=np.complex128)
    kvec_c = None

    ib, ilocal = 0, 0
    eband, occ = np.zeros(nband), np.zeros(nband)
    for line in lines:
        if line.endswith("(band)"):
            ib = int(line.split()[0]) - 1
            ilocal = 0
            continue
        if line.endswith("(Ry)"):
            eband[ib] = float(line.split()[0])
            continue
        if line.endswith("(Occupations)"):
            occ[ib] = float(line.split()[0])
            continue
        if not line.endswith(")"):
            if line.count(" ") != 2:
                # it means it is a band rather than coordinates of kpoint
                nums = line.split()
                c = [complex(float(nums[i]), float(nums[i+1])) for i in range(0, len(nums), 2)]
                for i in range(len(c)):
                    lowf_k[ilocal, ib] = c[i]
                    ilocal += 1
            else:
                kvec_c = np.array([float(x) for x in line.split()])
    return lowf_k, kvec_c, eband, occ
"""File: QO_ovlp_ik.dat"""
def read_ao_proj(fao_proj):
    """Read the atomic orbital projection matrix in k-space, and the k-vector.
    
    Args:
        fao_proj (str): path to the file containing the atomic orbital projection matrix
    
    Returns:
        np.ndarray: atomic orbital projection matrix in k-space
        np.ndarray: k-vector
    """
    with open(fao_proj, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    kvec_d = None
    if lines[0].startswith("KPOINT_COORDINATE"):
        kvec_d = np.array([float(x) for x in lines[0].split()[-3:]])
        lines = lines[1:]
    nlocal = len(lines[0].split())
    nao = len(lines)
    ao_proj = np.zeros((nao, nlocal), dtype=np.complex128)
    for i, line in enumerate(lines):
        nums = [cxx_topycomplex(num) for num in line.split()]
        for j in range(nlocal):
            ao_proj[i, j] = nums[j]
    # because in ao_proj, the row index is AO, we need to dagger it
    return ao_proj.conj().T, kvec_d
"""File: QO_supercells.dat"""
def read_qo_supercells(fqo_supercell):
    return np.loadtxt(fqo_supercell, dtype=np.int32)
"""File: kpoints"""
def read_kpoints(fkpoints, 
                 as_dict = False, 
                 titles: list = None,
                 table_titles: list = None
                ):
    """will read kpoints file and convert to plain np.ndarray if as_dict is not specified"""
    if titles is None:
        titles = [["index", "direct_x", "direct_y", "direct_z", "weight"],
                  ["index_0", "direct_x0", "direct_y0", "direct_z0", 
                   "index", "direct_x", "direct_y", "direct_z"]]
    if table_titles is None:
        table_titles = ["Reduced KPOINTS", 
                        "KPOINTS before-after correspondence of symmetry reduction"]
    with open(fkpoints, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    lines = [line for line in lines if len(line) > 0]
    lines = [line for line in lines if line[0].isdigit()]
    # sort lines into different lists according to number of columns
    lines = [line.split() for line in lines]
    ncols = [len(line) for line in lines]
    # make unique
    tables_range = []
    table_start = 0
    ncols_table = [ncols[0]]
    for i in range(1, len(ncols)):
        if ncols[i] != ncols[i-1]:
            tables_range.append((table_start, i))
            table_start = i
            ncols_table.append(ncols[i-1])
    tables_range.append((table_start, len(ncols)))
    # cut the lines into different tables
    tables = [lines[tables_range[i][0]:tables_range[i][1]] for i in range(len(tables_range))]
    # convert tables into np.ndarray
    tables = [np.array(table, dtype=np.float32) for table in tables]
    # makeup a list of dict
    if as_dict:
        result = []
        for table in tables:
            nrows, ncols = table.shape
            it = None
            for jt, title in enumerate(titles):
                if ncols == len(title):
                    it = jt
            if it is not None:
                result.append(dict(zip(titles[it], table.T.tolist())))
            else:
                raise ValueError("ill defined titles for kpoint source data.")
        tables = dict(zip(table_titles, result))
    return tables
"""File: istate.info"""
def read_istate(fistate):
    with open(fistate, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    lines = [line for line in lines if len(line) > 0]
    lines = [line for line in lines if line[0].isdigit()]
    data = np.array([line.split() for line in lines], dtype=np.float32)

    nrows, ncols = data.shape
    assert ncols == 3 or ncols == 5

    nks = nrows/np.max(data[:, 0])
    assert nks.is_integer()
    nks = int(nks)
    data_up = data.reshape((nks, -1, ncols))
    
    if ncols == 5:
        nband, ncols = data_up.shape[1], data_up.shape[2]
        data_dw = np.zeros((nks, nband, 3))
        data_dw[:, :, 0] = data_up[:, :, 0]
        data_dw[:, :, 1] = data_up[:, :, 3]
        data_dw[:, :, 2] = data_up[:, :, 4]
        data_up = data_up[:, :, :3]
        return [data_up, data_dw]
    else:
        return [data_up]
"""File: BANDS_*.dat"""
def read_bands(fband):
    return np.loadtxt(fband)
"""File: time.json"""
import json
def read_time_statistics(ftimejson):
    with open(ftimejson, "r") as f:
        timejson = json.load(f)
    return timejson
"""File: running_${}.log"""
def read_etraj_fromlog(flog, unit = "eV", term = "EKS"):
    harris = ["Harris", "harris", "HARRIS", "eharris", "EHARRIS", "eh"]
    fermi = ["Fermi", "fermi", "FERMI", "efermi", "EFERMI", "ef"]
    kohnsham = ["KohnSham", "kohnsham", "KOHN", "kohn", "KOHNSHAM", 
                "ekohnsham", "EKOHN", "EKOHNSHAM", "eks", "EKS", "e", "E", "energy"]
    if term in harris:
        header = "E_Harris"
    elif term in fermi:
        header = "E_Fermi"
    elif term in kohnsham:
        header = "E_KohnSham"
    else:
        raise ValueError("Unknown energy term")
    with open(flog, "r") as f:
        eners = [float(line.split()[-1]) for line in f.readlines() if line.strip().startswith(header)]
    eners = [unit_conversion(ener, "eV", unit) for ener in eners]
    return np.array(eners, dtype=np.float64)
"""File: running_${}.log"""
def read_natom_fromlog(flog):
    with open(flog, "r") as f:
        lines = [line.strip() for line in f.readlines()]
    for line in lines:
        if line.startswith("TOTAL ATOM NUMBER"):
            return int(line.split()[-1])
    return None
"""File: istate.info"""
def calculate_dos_fromistate(fistate, fkpoint = None, as_pair = False, with_occ = False, with_kwt = False):
    """Calculate the density of states from the istate.info file, optionally read kpoints weights from the kpoints file.

    Args:
    fistate (str): path to the istate.info file
    fkpoint (str): path to the kpoints file
    as_pair (bool): if True, return the dos as a pair of np.ndarray, otherwise return a list of np.ndarray
    with_occ (bool): if True, output the actual occpuation instead of the kpoint weight
    with_kwt (bool): if True, read kpoint weights from the kpoints file
    """
    istate = read_istate(fistate)
    kwts = None if (fkpoint is None or with_kwt) else read_kpoints(fkpoint)[0][:, -1]
    nspin = len(istate)
    assert nspin == 1 or nspin == 2

    dos = [[] for i in range(nspin)]
    for ispin in range(nspin):
        data = istate[ispin]
        nks, nband, ncols = data.shape
        assert ncols == 3
        # rearrange to list of tuple of (e, occ) pairs
        for ik in range(nks):
            wt = 1 if kwts is None else kwts[ik]
            for ib in range(nband):
                occ = wt
                occ *= data[ik, ib, 2] if with_occ else 1/nks
                dos[ispin].append((data[ik, ib, 1], occ))
        # sort the dos according to the energy
        dos[ispin].sort(key=lambda x: x[0])
        # accumulate occ with the same energy
        dos[ispin] = [(k, sum([x[1] for x in g])) for k, g in groupby(dos[ispin], key=lambda x: x[0])]
        # rearrange to e and occ list
        dos[ispin] = np.array(dos[ispin])
        if not as_pair:
            dos[ispin] = [[dos[ispin][i, 0], dos[ispin][i, 1]] for i in range(len(dos[ispin]))]
            dos[ispin] = np.array(dos[ispin]).T
    return dos
"""File: INPUT"""
import re
def read_keyvals_frominput(fin, keyword: str = None):
    kv = r"^([\w_-]+)(\s+)([^#]*)(.*)$"

    result = {}
    with open(fin, "r") as f:
        lines = [line.strip() for line in f.readlines()]
    lines = [line for line in lines if len(line) > 0]
    for line in lines:
        _match = re.match(kv, line)
        if _match is not None:
            result[_match.group(1)] = _match.group(3).strip()
    
    return result if keyword is None else result[keyword]
"""Path from Bohrium.dp.tech"""
def read_testconfig_fromBohriumpath(path: str):
    """    
    Bohrium will arrange the path in the following way:
    11548850/10430308/tmp/outputs/artifacts/outputs/Ar_23155_pd04_PBEecw100pwclfs1clstrs1scf/00010
    in this folder there will be OUT.ABACUS and INPUT, STRU, KPT, ...
    """
    path = path.replace("\\", "/")
    frags = path.split("/")

    _tp = r"^([\w]*)(_[\d]+)(_[^_]+)(_.*)$"
    test = ""
    for frag in frags:
        if "_" in frag and "OUT." not in frag and re.match(_tp, frag) is not None:
            _match = re.match(_tp, frag)
            system = _match.group(1)
            mpid = _match.group(2)[1:]
            pnid = _match.group(3)[1:]
            test = _match.group(4)[1:]
            return frag, system, mpid, pnid, test
    return None

def read_pressure_fromlog(flog):
    with open(flog, "r") as f:
        lines = [line.strip() for line in f.readlines()]
    for line in lines:
        if line.startswith("TOTAL-PRESSURE"):
            return float(line.split()[-2])
    return None

"""File: STRU/STRU_ION_D"""
def read_stru(fstruiond):
    """STRU_ION_D file is generated after relax, it for instance like:
    ```plaintext
    ATOMIC_SPECIES
    Si 28.0855 Si_ONCV_PBE-1.0.upf upf201

    LATTICE_CONSTANT
    1.92678975467

    LATTICE_VECTORS
        3.8492784000     0.0000000000     0.0000000000 #latvec1
        1.9246390700     3.3335741200     0.0000000000 #latvec2
        1.9246388000     1.1111903900     3.1429226100 #latvec3

    ATOMIC_POSITIONS
    Direct

    Si #label
    0 #magnetism
    2 #number of atoms
        0.8750000000     0.8750000000     0.8750000000 m  1  1  1
        0.1250000000     0.1250000000     0.1250000000 m  1  1  1
    ```
    """
    with open(fstruiond, "r") as f:
        lines = [line.strip() for line in f.readlines()]
    
    def read_lattice_constant(lines):
        for i, line in enumerate(lines):
            if line.startswith("LATTICE_CONSTANT"):
                return float(lines[i+1])
        return None

    def read_lattice_vectors(lines):
        for i, line in enumerate(lines):
            if line.startswith("LATTICE_VECTORS"):
                return np.array([line.split()[:3] for line in lines[i+1:i+4]], dtype=np.float32)
        return None

    return {"lattice_constant": read_lattice_constant(lines),
            "lattice_vectors": read_lattice_vectors(lines)}

"""File: STRU/STRU_ION_D"""
def read_volume_fromstru(fstruiond, unit = "Bohr"):

    result = read_stru(fstruiond)
    lattice_vectors = result["lattice_vectors"] # in angstrom
    lattice_constant = result["lattice_constant"] # in Bohr/angstrom * multiplier

    lv_inbohr = lattice_vectors * lattice_constant
    volume = np.linalg.det(lv_inbohr)*unit_conversion(1.0, "Bohr", unit)**3
    return volume

##########################
#  unittest starts here  #
##########################
import unittest
import os
class ReadAbacusOutTest(unittest.TestCase):

    def test_read_volume_fromstru(self):
        fstruiond = "STRU_ION_D"
        with open(fstruiond, "w") as f:
            f.write("""LATTICE_CONSTANT
1.0
                    
LATTICE_VECTORS
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
""")
        result = read_volume_fromstru(fstruiond)
        result_inAngstrom = read_volume_fromstru(fstruiond, "A")
        os.remove(fstruiond)
        self.assertAlmostEqual(result, 1.0, delta=1e-6)
        self.assertAlmostEqual(result_inAngstrom, result*unit_conversion(1.0, "Bohr", "A")**3, delta=1e-6)

    def test_read_stru(self):
        fstruiond = "STRU_ION_D"
        with open(fstruiond, "w") as f:
            f.write("""ATOMIC_SPECIES
Si 28.0855 Si_ONCV_PBE-1.0.upf upf201
                    
LATTICE_CONSTANT
1.92678975467
                    
LATTICE_VECTORS
    3.8492784000     0.0000000000     0.0000000000 #latvec1
    1.9246390700     3.3335741200     0.0000000000 #latvec2
    1.9246388000     1.1111903900     3.1429226100 #latvec3

ATOMIC_POSITIONS
Direct
                    
Si #label
0 #magnetism
2 #number of atoms
    0.8750000000     0.8750000000     0.8750000000 m  1  1  1
    0.1250000000     0.1250000000     0.1250000000 m  1  1  1
""")
        result = read_stru(fstruiond)
        os.remove(fstruiond)
        self.assertAlmostEqual(result["lattice_constant"], 1.92678975467, delta=1e-6)
        self.assertAlmostEqual(result["lattice_vectors"][0, 0], 3.8492784000, delta=1e-6)
        self.assertAlmostEqual(result["lattice_vectors"][1, 1], 3.3335741200, delta=1e-6)
        self.assertAlmostEqual(result["lattice_vectors"][2, 2], 3.1429226100, delta=1e-6)

    def test_unit_conversion(self):
        self.assertAlmostEqual(unit_conversion(1.0, "eV", "Ry"), 1/13.6)
        self.assertAlmostEqual(unit_conversion(1.0, "Ry", "eV"), 13.6)
        self.assertAlmostEqual(unit_conversion(1.0, "eV", "Hartree"), 1/27.2)
        self.assertAlmostEqual(unit_conversion(1.0, "Hartree", "eV"), 27.2)
        self.assertAlmostEqual(unit_conversion(1.0, "eV", "meV"), 1e3)
        self.assertAlmostEqual(unit_conversion(1.0, "meV", "eV"), 1e-3)
        self.assertAlmostEqual(unit_conversion(1.0, "eV", "J"), 1.6e-19)
        self.assertAlmostEqual(unit_conversion(1.0, "J", "eV"), 1/1.6e-19)
        self.assertAlmostEqual(unit_conversion(1.0, "eV", "kcal/mol"), 23.1)
        self.assertAlmostEqual(unit_conversion(1.0, "kcal/mol", "eV"), 1/23.1)

        self.assertAlmostEqual(unit_conversion(1.0, "A", "m"), 1e-10)
        self.assertAlmostEqual(unit_conversion(1.0, "m", "A"), 1e10)
        self.assertAlmostEqual(unit_conversion(1.0, "A", "nm"), 0.1)
        self.assertAlmostEqual(unit_conversion(1.0, "nm", "A"), 10)
        self.assertAlmostEqual(unit_conversion(1.0, "A", "Bohr"), 1.8897259886)
        self.assertAlmostEqual(unit_conversion(1.0, "Bohr", "A"), 1/1.8897259886)
        self.assertAlmostEqual(unit_conversion(1.0, "A", "a.u."), 1.8897259886)
        self.assertAlmostEqual(unit_conversion(1.0, "a.u.", "A"), 1/1.8897259886)

        self.assertAlmostEqual(unit_conversion(1.0, "a.u.", "kg"), 9.1e-31)
        self.assertAlmostEqual(unit_conversion(1.0, "kg", "a.u."), 1/9.1e-31)
        self.assertAlmostEqual(unit_conversion(1.0, "a.u.", "g"), 9.1e-28)

    def test_read_mat_hs(self):
        contents = """26 (1.8656418657e-01,1.0594053029e-17) (2.7657956362e-01,-5.1942556735e-03) (2.1388163461e-02,-1.2201371002e-01) (2.1388163461e-02,-1.2201371002e-01) (-2.1388163461e-02,1.2201371002e-01) (1.4269496956e-02,6.2019770534e-03) (1.4269496956e-02,6.2019770534e-03) (-1.4269496956e-02,-6.2019770534e-03) (-7.6301557547e-17,3.8954198060e-17) (2.2415242205e-02,-3.4366702659e-02) (-2.2415242205e-02,3.4366702659e-02) (-2.3880192233e-16,1.3978796491e-17) (-2.2415242205e-02,3.4366702659e-02) (-1.2800502202e-01,-9.0333484133e-02) (-1.9486369032e-01,-2.4335232173e-01) (6.4137465115e-02,-5.8874644112e-02) (6.4137465115e-02,-5.8874644112e-02) (-6.4137465115e-02,5.8874644112e-02) (2.3377806518e-02,2.0353989708e-01) (2.3377806518e-02,2.0353989708e-01) (-2.3377806518e-02,-2.0353989708e-01) (-3.1288258746e-17,-1.9558131271e-16) (6.5014513162e-02,-1.6050331725e-01) (-6.5014513162e-02,1.6050331725e-01) (-7.5036478095e-17,-1.2491084505e-16) (-6.5014513162e-02,1.6050331725e-01)
(1.3216196889e+00,2.8948192551e-17) (-2.8313057267e-02,-1.5333446699e-01) (-2.8313057267e-02,-1.5333446699e-01) (2.8313057267e-02,1.5333446699e-01) (-1.7435850730e-02,-3.7452130576e-01) (-1.7435850730e-02,-3.7452130576e-01) (1.7435850730e-02,3.7452130576e-01) (4.1961645570e-18,-9.7002177464e-18) (4.2282790690e-01,4.3160147291e-02) (-4.2282790690e-01,-4.3160147291e-02) (-1.1516357279e-16,-2.0291521803e-17) (-4.2282790690e-01,-4.3160147291e-02) (-1.9486369032e-01,-2.4335232173e-01) (1.4866576902e-01,-1.5925677954e-01) (7.3880771317e-02,-1.7322807871e-01) (7.3880771317e-02,-1.7322807871e-01) (-7.3880771317e-02,1.7322807871e-01) (2.0517309715e-01,-1.6195489409e-01) (2.0517309715e-01,-1.6195489409e-01) (-2.0517309715e-01,1.6195489409e-01) (-3.3044239734e-17,2.5799029706e-17) (1.8061791562e-01,-2.0869635148e-01) (-1.8061791562e-01,2.0869635148e-01) (1.0155077684e-16,7.9989241024e-17) (-1.8061791562e-01,2.0869635148e-01)
(7.3971299983e-01,-3.0222136962e-18) (-9.0322401358e-02,-1.9092114550e-18) (9.0322401358e-02,3.4842275745e-17) (-6.5256227698e-01,1.0870407103e-02) (1.8635417304e-01,7.6725492934e-03) (-1.8635417304e-01,-7.6725492934e-03) (-2.6890883495e-02,-9.7401524340e-04) (-7.8343347244e-02,-8.7611606346e-02) (7.8343347244e-02,8.7611606346e-02) (-7.7979747647e-17,-1.9270429674e-17) (1.8595260320e-01,6.2981194225e-02) (-6.4137465115e-02,5.8874644112e-02) (-7.3880771317e-02,1.7322807871e-01) (-6.3069526754e-02,-3.9026758903e-02) (-6.0888107411e-02,2.0602460950e-02) (6.0888107411e-02,-2.0602460950e-02) (6.6457754247e-03,1.2386635615e-01) (5.7376891448e-03,2.2828218389e-01) (-5.7376891448e-03,-2.2828218389e-01) (8.7836536072e-03,2.2517661827e-01) (6.7161780304e-02,-6.8029897100e-02) (-6.7161780304e-02,6.8029897100e-02) (3.0928666739e-17,-3.8433263943e-17) (-2.3761508049e-01,-2.7633789845e-01)
(7.3971299983e-01,7.5792516417e-18) (9.0322401358e-02,5.0539068593e-17) (1.8635417304e-01,7.6725492934e-03) (-6.5256227698e-01,1.0870407103e-02) (-1.8635417304e-01,-7.6725492934e-03) (1.3445441748e-02,4.8700762170e-04) (-7.8343347244e-02,-8.7611606346e-02) (1.8595260320e-01,6.2981194225e-02) (-2.3288188237e-02,-8.4352194446e-04) (7.8343347244e-02,8.7611606346e-02) (-6.4137465115e-02,5.8874644112e-02) (-7.3880771317e-02,1.7322807871e-01) (-6.0888107411e-02,2.0602460950e-02) (-6.3069526754e-02,-3.9026758903e-02) (6.0888107411e-02,-2.0602460950e-02) (5.7376891448e-03,2.2828218389e-01) (6.6457754247e-03,1.2386635615e-01) (-5.7376891448e-03,-2.2828218389e-01) (-4.3918268036e-03,-1.1258830914e-01) (6.7161780304e-02,-6.8029897100e-02) (-2.3761508049e-01,-2.7633789845e-01) (7.6068671619e-03,1.9500867176e-01) (-6.7161780304e-02,6.8029897100e-02)
(7.3971299983e-01,-9.6050360371e-18) (-1.8635417304e-01,-7.6725492934e-03) (-1.8635417304e-01,-7.6725492934e-03) (-6.5256227698e-01,1.0870407103e-02) (-1.3445441748e-02,-4.8700762170e-04) (1.8595260320e-01,6.2981194225e-02) (-7.8343347244e-02,-8.7611606346e-02) (-2.3288188237e-02,-8.4352194446e-04) (-7.8343347244e-02,-8.7611606346e-02) (6.4137465115e-02,-5.8874644112e-02) (7.3880771317e-02,-1.7322807871e-01) (6.0888107411e-02,-2.0602460950e-02) (6.0888107411e-02,-2.0602460950e-02) (-6.3069526754e-02,-3.9026758903e-02) (-5.7376891448e-03,-2.2828218389e-01) (-5.7376891448e-03,-2.2828218389e-01) (6.6457754247e-03,1.2386635615e-01) (4.3918268036e-03,1.1258830914e-01) (-2.3761508049e-01,-2.7633789845e-01) (6.7161780304e-02,-6.8029897100e-02) (7.6068671619e-03,1.9500867176e-01) (6.7161780304e-02,-6.8029897100e-02)
(1.7520218717e+00,-6.4130554457e-17) (5.3852842957e-01,3.6808708660e-17) (-5.3852842957e-01,1.3687670390e-18) (8.5593991303e-03,9.3151719182e-03) (-9.9468742746e-02,-5.5647140926e-02) (9.9468742746e-02,5.5647140926e-02) (2.1211821164e-17,1.1333274876e-17) (4.5930771095e-02,2.1107061039e-02) (-2.3377806518e-02,-2.0353989708e-01) (-2.0517309715e-01,1.6195489409e-01) (6.6457754247e-03,1.2386635615e-01) (5.7376891448e-03,2.2828218389e-01) (-5.7376891448e-03,-2.2828218389e-01) (2.7337024918e-01,3.7926301870e-01) (2.5168924767e-01,2.1329119083e-01) (-2.5168924767e-01,-2.1329119083e-01) (2.2981752577e-02,4.6301338252e-02) (9.9621785078e-02,3.5219380578e-01) (-9.9621785078e-02,-3.5219380578e-01) (-9.4782982403e-18,7.4090819498e-17) (-8.9479721787e-02,-2.5886080674e-01)
(1.7520218717e+00,-1.3498315304e-17) (-5.3852842957e-01,-5.0144118980e-18) (-4.2796995652e-03,-4.6575859591e-03) (-9.9468742746e-02,-5.5647140926e-02) (4.5930771095e-02,2.1107061039e-02) (7.4126570880e-03,8.0671755217e-03) (9.9468742746e-02,5.5647140926e-02) (-2.3377806518e-02,-2.0353989708e-01) (-2.0517309715e-01,1.6195489409e-01) (5.7376891448e-03,2.2828218389e-01) (6.6457754247e-03,1.2386635615e-01) (-5.7376891448e-03,-2.2828218389e-01) (2.5168924767e-01,2.1329119083e-01) (2.7337024918e-01,3.7926301870e-01) (-2.5168924767e-01,-2.1329119083e-01) (-1.1490876289e-02,-2.3150669126e-02) (9.9621785078e-02,3.5219380578e-01) (-8.9479721787e-02,-2.5886080674e-01) (1.9902781555e-02,4.0098135155e-02) (-9.9621785078e-02,-3.5219380578e-01)
(1.7520218717e+00,3.3068106531e-18) (4.2796995652e-03,4.6575859591e-03) (4.5930771095e-02,2.1107061039e-02) (-9.9468742746e-02,-5.5647140926e-02) (7.4126570880e-03,8.0671755217e-03) (-9.9468742746e-02,-5.5647140926e-02) (2.3377806518e-02,2.0353989708e-01) (2.0517309715e-01,-1.6195489409e-01) (-5.7376891448e-03,-2.2828218389e-01) (-5.7376891448e-03,-2.2828218389e-01) (6.6457754247e-03,1.2386635615e-01) (-2.5168924767e-01,-2.1329119083e-01) (-2.5168924767e-01,-2.1329119083e-01) (2.7337024918e-01,3.7926301870e-01) (1.1490876289e-02,2.3150669126e-02) (-8.9479721787e-02,-2.5886080674e-01) (9.9621785078e-02,3.5219380578e-01) (1.9902781555e-02,4.0098135155e-02) (9.9621785078e-02,3.5219380578e-01)
(2.4139923460e+00,7.6787772460e-18) (-1.4879110680e-01,1.3185469406e-02) (1.4879110680e-01,-1.3185469406e-02) (-1.1301121335e-16,1.1357441320e-17) (-2.9758221361e-01,2.6370938812e-02) (-5.9301014740e-17,-4.4345327139e-17) (2.5328385305e-18,1.4728981071e-16) (-8.7836536072e-03,-2.2517661827e-01) (4.3918268036e-03,1.1258830914e-01) (-4.3918268036e-03,-1.1258830914e-01) (-2.2981752577e-02,-4.6301338252e-02) (1.1490876289e-02,2.3150669126e-02) (-1.1490876289e-02,-2.3150669126e-02) (-4.1476995668e-01,-7.5124834078e-01) (-5.1394558523e-03,2.8467623007e-01) (5.1394558523e-03,-2.8467623007e-01) (1.1079842106e-16,5.2926584806e-17) (-1.0278911705e-02,5.6935246014e-01)
(2.4049917880e+00,6.2239981174e-18) (-2.7684284710e-01,2.6199152558e-17) (-2.5771375670e-01,-2.2837902932e-02) (-2.7684284710e-01,3.1191988247e-17) (6.5014513162e-02,-1.6050331725e-01) (1.8061791562e-01,-2.0869635148e-01) (-6.7161780304e-02,6.8029897100e-02) (-6.7161780304e-02,6.8029897100e-02) (2.3761508049e-01,2.7633789845e-01) (-9.9621785078e-02,-3.5219380578e-01) (-9.9621785078e-02,-3.5219380578e-01) (8.9479721787e-02,2.5886080674e-01) (-5.1394558523e-03,2.8467623007e-01) (5.4420852486e-03,-7.2381129508e-03) (7.8364208217e-02,2.1769780678e-01) (-8.9017986595e-03,4.9307369419e-01) (7.8364208217e-02,2.1769780678e-01)
(2.4049917880e+00,-1.0651439143e-17) (-2.5771375670e-01,-2.2837902932e-02) (2.7684284710e-01,-3.1851509209e-17) (-6.5014513162e-02,1.6050331725e-01) (-1.8061791562e-01,2.0869635148e-01) (6.7161780304e-02,-6.8029897100e-02) (2.3761508049e-01,2.7633789845e-01) (-6.7161780304e-02,6.8029897100e-02) (9.9621785078e-02,3.5219380578e-01) (8.9479721787e-02,2.5886080674e-01) (-9.9621785078e-02,-3.5219380578e-01) (5.1394558523e-03,-2.8467623007e-01) (7.8364208217e-02,2.1769780678e-01) (5.4420852486e-03,-7.2381129508e-03) (-8.9017986595e-03,4.9307369419e-01) (-7.8364208217e-02,-2.1769780678e-01)
(2.4139923460e+00,7.5470635514e-18) (3.4506503993e-16,1.7942063745e-17) (-3.1394905410e-17,1.6799067732e-17) (8.8650889450e-17,3.0878754156e-16) (-3.1617620965e-18,4.0010756360e-17) (-7.6068671619e-03,-1.9500867176e-01) (-7.6068671619e-03,-1.9500867176e-01) (-2.4143827671e-17,1.3050830379e-17) (-1.9902781555e-02,-4.0098135155e-02) (-1.9902781555e-02,-4.0098135155e-02) (1.2370460919e-16,7.3272519918e-17) (-8.9017986595e-03,4.9307369419e-01) (-8.9017986595e-03,4.9307369419e-01) (-4.1476995668e-01,-7.5124834078e-01) (6.5303063048e-18,3.1293434834e-16)
(2.4049917880e+00,6.7006661712e-18) (-6.5014513162e-02,1.6050331725e-01) (-1.8061791562e-01,2.0869635148e-01) (2.3761508049e-01,2.7633789845e-01) (6.7161780304e-02,-6.8029897100e-02) (-6.7161780304e-02,6.8029897100e-02) (8.9479721787e-02,2.5886080674e-01) (9.9621785078e-02,3.5219380578e-01) (-9.9621785078e-02,-3.5219380578e-01) (-1.0278911705e-02,5.6935246014e-01) (7.8364208217e-02,2.1769780678e-01) (-7.8364208217e-02,-2.1769780678e-01) (-1.4106380917e-17,2.0668479850e-16) (5.4420852486e-03,-7.2381129508e-03)
(1.8656418657e-01,-8.5533385074e-18) (2.7657956362e-01,5.1942556735e-03) (-2.1388163461e-02,-1.2201371002e-01) (-2.1388163461e-02,-1.2201371002e-01) (2.1388163461e-02,1.2201371002e-01) (-1.4269496956e-02,6.2019770534e-03) (-1.4269496956e-02,6.2019770534e-03) (1.4269496956e-02,-6.2019770534e-03) (-1.5775032225e-16,-3.8309818052e-17) (2.2415242205e-02,3.4366702659e-02) (-2.2415242205e-02,-3.4366702659e-02) (-2.7590890776e-16,9.8512050303e-18) (-2.2415242205e-02,-3.4366702659e-02)
(1.3216196889e+00,-5.8546928058e-17) (2.8313057267e-02,-1.5333446699e-01) (2.8313057267e-02,-1.5333446699e-01) (-2.8313057267e-02,1.5333446699e-01) (1.7435850730e-02,-3.7452130576e-01) (1.7435850730e-02,-3.7452130576e-01) (-1.7435850730e-02,3.7452130576e-01) (-2.3428176567e-16,8.6498997235e-17) (4.2282790690e-01,-4.3160147291e-02) (-4.2282790690e-01,4.3160147291e-02) (-3.1133883390e-17,4.4118564079e-17) (-4.2282790690e-01,4.3160147291e-02)
(7.3971299983e-01,-8.6973353376e-18) (-9.0322401358e-02,2.8677147428e-17) (9.0322401358e-02,3.3511165082e-17) (-6.5256227698e-01,-1.0870407103e-02) (1.8635417304e-01,-7.6725492934e-03) (-1.8635417304e-01,7.6725492934e-03) (2.6890883495e-02,-9.7401524340e-04) (7.8343347244e-02,-8.7611606346e-02) (-7.8343347244e-02,8.7611606346e-02) (-1.8303571849e-18,2.4146361618e-17) (-1.8595260320e-01,6.2981194225e-02)
(7.3971299983e-01,8.8887631800e-18) (9.0322401358e-02,1.4773947871e-17) (1.8635417304e-01,-7.6725492934e-03) (-6.5256227698e-01,-1.0870407103e-02) (-1.8635417304e-01,7.6725492934e-03) (-1.3445441748e-02,4.8700762170e-04) (7.8343347244e-02,-8.7611606346e-02) (-1.8595260320e-01,6.2981194225e-02) (2.3288188237e-02,-8.4352194446e-04) (-7.8343347244e-02,8.7611606346e-02)
(7.3971299983e-01,1.3979245040e-18) (-1.8635417304e-01,7.6725492934e-03) (-1.8635417304e-01,7.6725492934e-03) (-6.5256227698e-01,-1.0870407103e-02) (1.3445441748e-02,-4.8700762170e-04) (-1.8595260320e-01,6.2981194225e-02) (7.8343347244e-02,-8.7611606346e-02) (2.3288188237e-02,-8.4352194446e-04) (7.8343347244e-02,-8.7611606346e-02)
(1.7520218717e+00,2.7538695418e-17) (5.3852842957e-01,4.4831673896e-17) (-5.3852842957e-01,5.5457018359e-17) (-8.5593991303e-03,9.3151719182e-03) (9.9468742746e-02,-5.5647140926e-02) (-9.9468742746e-02,5.5647140926e-02) (-6.9997820844e-17,4.4042296901e-17) (-4.5930771095e-02,2.1107061039e-02)
(1.7520218717e+00,2.3906708428e-17) (-5.3852842957e-01,-1.0652324066e-17) (4.2796995652e-03,-4.6575859591e-03) (9.9468742746e-02,-5.5647140926e-02) (-4.5930771095e-02,2.1107061039e-02) (-7.4126570880e-03,8.0671755217e-03) (-9.9468742746e-02,5.5647140926e-02)
(1.7520218717e+00,-5.9359669572e-18) (-4.2796995652e-03,4.6575859591e-03) (-4.5930771095e-02,2.1107061039e-02) (9.9468742746e-02,-5.5647140926e-02) (-7.4126570880e-03,8.0671755217e-03) (9.9468742746e-02,-5.5647140926e-02)
(2.4139923460e+00,2.6685775856e-18) (-1.4879110680e-01,-1.3185469406e-02) (1.4879110680e-01,1.3185469406e-02) (-6.8867797932e-17,-1.1152458893e-17) (-2.9758221361e-01,-2.6370938812e-02)
(2.4049917880e+00,2.1775521270e-17) (-2.7684284710e-01,2.7818679575e-17) (-2.5771375670e-01,2.2837902932e-02) (-2.7684284710e-01,5.8113236464e-17)
(2.4049917880e+00,9.9145201054e-18) (-2.5771375670e-01,2.2837902932e-02) (2.7684284710e-01,-1.0719943177e-17)
(2.4139923460e+00,1.4541440928e-18) (3.7153746846e-16,-2.9663093932e-17)
(2.4049917880e+00,-9.8183814489e-18)
"""
        with open("unittest-data-H-0", "w") as f:
            f.write(contents)
        mat = read_mat_hs("unittest-data-H-0")
        os.remove("unittest-data-H-0")
        nrows, ncols = mat.shape
        self.assertEqual(nrows, ncols)
        self.assertListEqual(mat.tolist(), mat.conj().T.tolist())
        self.assertNotEqual(np.sum(np.abs(mat)), 0)

    def test_read_lowf_k(self):
        contents = """100 (index of k points)
-0.2000000000000000111022302 0.2000000000000000111022302 0.600000000000000088817842
4 (number of bands)
26 (number of orbitals)
1 (band)
-1.1860689608372088077459239e-01 (Ry)
1.5999999999999989924726052e-02 (Occupations)
-3.1977706914050274544791819e-01 1.4778402768247433929360568e-02 -9.9147889720000825786172527e-02 -4.1025830543366670344340719e-03 -4.0867853622418518627634754e-03 -1.9726774126268215026058783e-01 -3.5785065516329428447761529e-03 -5.1758015784750650845058573e-02 3.5785065516354130910059439e-03 5.1758015784753551302710406e-02 
-4.6569477447110148998765666e-03 1.3957924704586788244853324e-02 -8.5092652018143311032360998e-03 1.3662220254420518600335654e-02 8.5092652018162462379535782e-03 -1.3662220254418409176588867e-02 1.1450515384940938151681067e-02 -7.1349859911291023395474298e-04 7.6371822253800125679390476e-03 -4.9509551049409759870822967e-05 
-7.6371822253802068569683570e-03 4.9509551048958298086200003e-05 -6.4346264989109116859951282e-17 1.4609851835016357531593253e-17 -4.4196959698856565570967980e-03 -4.7708433192154703350595923e-03 1.7620240573500773467729985e-02 -3.1963307441283750343075098e-01 -3.2211191770134192657870642e-03 -9.9180439678873480446874566e-02 
1.9722362158670009302241510e-01 5.8400895287319381371515803e-03 5.1724162556428934789387597e-02 4.0384297629597575018944866e-03 -5.1724162556431481363450331e-02 -4.0384297629619588659855012e-03 -1.3998767776575644569692081e-02 4.5326951246921478796969573e-03 -1.3737317333977643041276195e-02 8.3874888387853001325833446e-03 
1.3737317333976014135932253e-02 -8.3874888387865977057433753e-03 -8.1525129551232534876742397e-04 1.1443720905233633819420014e-02 -1.1739267541587407144010324e-04 7.6364404338062717020330084e-03 1.1739267541522289961530845e-04 -7.6364404338062396096487028e-03 1.2318440084543081936088489e-16 -1.0212931888477439025837403e-16 
-4.7313692262129017351846016e-03 -4.4619282477396339167685646e-03 
2 (band)
6.2872916871888506751275827e-01 (Ry)
1.5999999999999989924726052e-02 (Occupations)
3.1260817161476051051494096e-01 3.9523051267608466652347943e-02 4.5350311236002971426728436e-02 -1.3018563646447239845471699e-02 -4.7060346417047149714107945e-03 -1.6264732290252792434870344e-01 7.9303110887917041549322050e-02 -2.9661651324172516597243998e-01 -7.9303110887925104544038390e-02 2.9661651324172666477352323e-01 
8.1223451744658059503478853e-03 4.0531141739573228743864775e-02 8.2593786618579551062069299e-03 3.3281571636863106766757170e-02 -8.2593786618616743533394242e-03 -3.3281571636860893259601824e-02 1.7671928469655628999440689e-03 9.4749649445254966600993285e-03 -1.8130792535928037800818657e-02 -1.1821746172252740880725952e-02 
1.8130792535928565156755354e-02 1.1821746172252563938931402e-02 -1.1673040885572585115361617e-16 -5.7350182674956975236210966e-18 2.7285521088207920215040758e-02 1.1035311447620044233985936e-02 -3.0189398196949296382385342e-01 -9.0254995387708467080578600e-02 -3.6213551274006730573695734e-02 -3.0244213210855561135836922e-02 
-7.0070524383306909133040108e-02 1.4685509872679483600599326e-01 -4.7405752814802276651207080e-02 3.0335298568538432428809415e-01 4.7405752814802797068249873e-02 -3.0335298568538998642551974e-01 2.3817480745071485659369870e-02 -3.3785700419033393160361811e-02 2.1011452149907520331462152e-02 -2.7099801199493577658072851e-02 
-2.1011452149905646830108097e-02 2.7099801199490420461346574e-02 -5.4474799553052991388235604e-03 7.9512825000055010332689065e-03 2.1362622492081966701915619e-02 -3.4810459351770437624873367e-03 -2.1362622492081258934737420e-02 3.4810459351761161191085581e-03 3.5562002555841115990996377e-16 -2.3916916217724609230636516e-18 
-2.9417584387218628783067231e-02 -9.3994075445877800674410185e-04 
3 (band)
9.8744543817293650711519604e-01 (Ry)
1.5999999999999989924726052e-02 (Occupations)
-8.3636292647237952623315733e-15 -1.3577240200388432423799089e-14 3.5572387634085824012124371e-15 5.2219313053405455284959666e-15 -4.5840919787120610610626527e-15 -2.6706490077950346691297862e-15 -7.8751487661521998862035332e-02 4.7687398465809860415021149e-01 -7.8751487661530700234990832e-02 4.7687398465810049152935335e-01 
-1.1213813484792038168774782e-15 -1.8589572847563108410883336e-15 2.2012576735092669738547144e-02 -2.0402005719270116523844649e-02 2.2012576735091177876357804e-02 -2.0402005719266758099195158e-02 -2.8627217008384311414880572e-16 1.5930459114340466411459415e-15 2.8455647149384297178720260e-02 -1.2874342576770836377164642e-02 
2.8455647149383964111812872e-02 -1.2874342576770317694845325e-02 -2.1175402412443172783795831e-02 1.0316413547243469261971427e-02 5.6479464920624599179231094e-17 4.5920182481040873357163753e-16 -3.9836088170890767092121206e-15 1.1259708453788227985349073e-14 8.5685132681265727678422681e-16 -4.3424682419959483095032841e-15 
-8.7588147963376015766579536e-16 -2.5259440664898476737129879e-15 -5.5164285999151907446957921e-03 4.8330131705601864222643371e-01 -5.5164285999093343182408944e-03 4.8330131705601869773758494e-01 -1.3173227355628467106665604e-15 -4.3622453229410693772900569e-15 -1.8122358639446951855678236e-02 -2.3924370195080882361970964e-02 
-1.8122358639445012434832094e-02 -2.3924370195081506862422316e-02 -6.6240975442659080759449097e-16 3.4005969263461408835908662e-16 2.5778136146996755151405978e-02 1.7634064995374659073901213e-02 2.5778136146996689231913891e-02 1.7634064995373288642355192e-02 -1.9054735759577438869438737e-02 -1.3847133311242322897816059e-02 
-7.1338993517695023144452640e-16 -1.0405777915451799600180700e-15 
4 (band)
9.9902696594316775513533457e-01 (Ry)
1.5999999999999989924726052e-02 (Occupations)
-2.9127789701333217653456131e-01 1.8193647520124978100852786e-02 -3.0777614071026970776756571e-02 -8.5553094794117434673808020e-03 -8.2735744738680394583241196e-02 5.2813834652818669734841706e-01 1.0393463249192083641680995e-01 -1.9055317264557419809634098e-01 -1.0393463249193035657924611e-01 1.9055317264557181111683803e-01 
4.4943443155054385584890042e-03 -3.8728531586091463623233722e-02 -1.0428752422613762895764467e-03 1.4997786016593572799848744e-02 1.0428752422590864545881573e-03 -1.4997786016591796443009343e-02 4.1704372554567853415008472e-02 2.5178418897788192032383137e-02 -6.8282815865452414771930023e-03 -1.6463734335311244749711079e-02 
6.8282815865454878079265910e-03 1.6463734335311827616799008e-02 4.2314612007130601624017984e-17 6.8151757356088314246584527e-17 1.5635201968096937297225679e-02 -2.6373239388834346375034556e-02 -1.7402489344011395910882811e-01 2.3428435406946329644384264e-01 -9.6038706340239640008027422e-03 3.0466711621679409827123663e-02 
4.8787861127001974637096282e-01 2.1851722402645196163817332e-01 -2.1640929383552776066323986e-01 -1.6731317286208545469516906e-02 2.1640929383552842679705464e-01 1.6731317286205797667530959e-02 -3.4919021166203996486832040e-02 -1.7342440728902463326077665e-02 1.3142296437640567374760536e-02 7.3007683510551528438536906e-03 
-1.3142296437636098727086420e-02 -7.3007683510556524442147719e-03 1.6234242929802460142507670e-03 -4.8688519816560978104025281e-02 1.0080854777147906994883009e-02 1.4698855215658366418773895e-02 -1.0080854777147214840216094e-02 -1.4698855215658702955128234e-02 4.6092781880857871623311441e-16 3.9939159460686655920808220e-17 
3.0633328662474972048368471e-02 1.2674665678546826255712610e-03 """
        with open("unittest-lowf", "w") as f:
            f.write(contents)
        result = read_lowf("unittest-lowf")
        os.remove("unittest-lowf")
        self.assertEqual(len(result), 4)
        psi, kvec_c, eband, occ = result
        nlocal, nbands = psi.shape
        self.assertEqual(kvec_c.shape, (3,))
        self.assertEqual(eband.shape, (nbands,))
        self.assertEqual(occ.shape, (nbands,))

    def test_read_ao_proj(self):
        contents = """KPOINT_COORDINATE:  -2.85714285714286e-01   0.00000000000000e+00  -4.28571428571429e-01
(-9.66773336319710e-03,4.04208265445529e-20) (-4.64885188287154e-02,-2.80278409684425e-19) (5.83346506713486e-20,-3.44101586037455e-04) (1.40866619761454e-20,-3.44101586037455e-04) (-3.73560157968366e-20,-4.01821902393477e-03) (-1.51789716628812e-19,1.21436746908368e-03) (1.94332398716500e-20,1.21436746908368e-03) (3.39558643832208e-20,1.41786296442356e-02) (2.48665724163603e-04,-1.96829131956368e-20) (-1.17824541765769e-04,2.74629123827173e-23) (-2.12312488011369e-04,6.01759589072413e-21) (4.30701668352267e-04,4.81199109123085e-21) (-2.12312488011369e-04,4.08899450759911e-21) (-2.23489627816395e-03,-3.71147572428485e-03) (-7.08055926387039e-03,-1.27309059614315e-02) (3.27583579473072e-03,1.97906603493135e-03) (3.27583579473072e-03,1.97906603493135e-03) (1.50303732608938e-02,-6.77610179403154e-03) (-5.91314922743143e-03,-2.34498807343220e-03) (-5.91314922743144e-03,-2.34498807343219e-03) (-2.70318811038195e-02,1.21805574370726e-02) (1.80396629197821e-04,5.08289741668323e-04) (1.04414973148315e-02,-4.79273738928242e-03) (1.97555516790357e-03,1.43656718674608e-03) (3.12456127284789e-04,8.80383657535595e-04) (1.97555516790357e-03,1.43656718674608e-03) 
(1.19528909651574e-01,6.84328194881948e-21) (1.01401587126593e-01,7.41441001437901e-19) (-9.31392048499040e-20,2.10018834695085e-03) (-7.97468775643641e-20,2.10018834695085e-03) (-1.90075185052089e-20,2.30616629955156e-02) (-2.41033720963447e-19,-7.35163662159392e-03) (-2.39380383648312e-18,-7.35163662159392e-03) (1.62099692397420e-19,-8.02885506000432e-02) (-1.28605165008403e-03,-1.94667693165630e-20) (5.72805793574380e-04,4.44815642732849e-20) (1.03216069858926e-03,-1.58357854707749e-20) (-2.22750679910332e-03,-1.83518404849661e-20) (1.03216069858926e-03,-1.72229790068017e-20) (1.32188222548314e-02,2.19285557738572e-02) (3.89088249896940e-02,7.01834438054466e-02) (-1.86963885935956e-02,-1.13032287756508e-02) (-1.86963885935956e-02,-1.13032287756508e-02) (-8.57856997240179e-02,3.86749159661674e-02) (3.24639262139521e-02,1.24569771571552e-02) (3.24639262139521e-02,1.24569771571552e-02) (1.48379368654628e-01,-6.68589073352137e-02) (-9.44098042451478e-04,-2.65922628609329e-03) (-6.00333776691774e-02,2.75128455141183e-02) (-1.15109816816545e-02,-8.45240793880238e-03) (-1.63522577685228e-03,-4.60591503633627e-03) (-1.15109816816545e-02,-8.45240793880239e-03) 
(-5.02689921880510e-01,4.80324694410518e-19) (2.46672729675783e-01,8.71673626240560e-19) (-1.88644539438013e-18,-7.09060083156112e-03) (7.50622104969982e-19,-7.09060083156112e-03) (-4.53455640479224e-19,-6.72002313275119e-02) (-3.87479067331042e-18,2.40757208430013e-02) (3.34646917032700e-18,2.40757208430013e-02) (-7.61840348944847e-19,2.21997872126945e-01) (2.03461766089111e-03,-1.73727663664651e-20) (-6.45190633076499e-04,4.09282420208846e-20) (-1.16591267977391e-03,5.07691905124611e-20) (3.52406116264035e-03,2.51243059630569e-20) (-1.16591267977391e-03,2.70877919992017e-20) (-4.18038420244833e-02,-6.91557733972685e-02) (-9.41804148257860e-02,-1.72420046087660e-01) (5.24309064272904e-02,3.15433823118449e-02) (5.24309064272904e-02,3.15433823118449e-02) (2.40719058246921e-01,-1.08567593274179e-01) (-7.82978617761724e-02,-2.49886112959492e-02) (-7.82978617761724e-02,-2.49886112959493e-02) (-3.58066978523747e-01,1.61487178578466e-01) (2.11361295146210e-03,5.85363994070313e-03) (1.73971016506751e-01,-7.94376033055252e-02) (3.45482773580546e-02,2.58733653754806e-02) (3.66088501946797e-03,1.01388017865123e-02) (3.45482773580546e-02,2.58733653754806e-02) 
(-2.48340518977958e-20,-1.02234217351818e-04) (-1.29256475159404e-19,-9.14989813333637e-04) (7.81146294060923e-02,9.65242472101171e-21) (5.32049801431030e-04,1.78421677343411e-20) (9.58720747451991e-04,3.23933099028097e-20) (9.54093337500756e-02,2.06750557628819e-20) (-1.67960809806152e-03,-7.53280646793652e-20) (-3.02654971487055e-03,1.70886390226571e-19) (-6.41832358184157e-21,9.31694948362761e-05) (-2.89455395064623e-22,-3.21102917973109e-04) (-3.83997312961536e-22,-7.73408859711252e-04) (8.47305457401495e-21,-2.49368845762530e-04) (-3.72420925206049e-27,-4.44869107936920e-11) (1.80620391809104e-03,1.27311966467630e-03) (2.27172799609056e-03,2.14173558366390e-04) (-1.25740754302192e-03,-2.52629792886868e-03) (1.27260730092084e-02,-5.63930379640086e-03) (3.11712763117411e-03,2.65399052555793e-03) (-7.10060132225724e-04,6.47663164341477e-04) (-9.49948039174898e-03,3.92457019390270e-03) (-3.32590854367089e-03,-3.25403313467156e-03) (-1.64325065008193e-03,-1.32140230821252e-03) (1.66958668674147e-03,1.32110665213213e-03) (7.35346158195902e-03,-3.29580517750373e-03) (-7.01392828854389e-07,7.86364276730358e-05) (-3.81686097082572e-03,-6.17059017724922e-03) 
(-4.78173138141991e-21,-1.02234217351818e-04) (-9.46313747369266e-20,-9.14989813333637e-04) (5.32049801431030e-04,1.78421677343411e-20) (7.81146294060923e-02,2.24073793290318e-20) (9.58720747451991e-04,2.60334337745676e-20) (-1.67960809806152e-03,-7.53280646793652e-20) (9.54093337500756e-02,8.23015943341441e-20) (-3.02654971487055e-03,1.65667162322279e-19) (4.71636505490217e-21,-2.62544502760893e-04) (1.44224971846458e-20,-3.21102917973109e-04) (-4.93032673159288e-28,-4.44869107936920e-11) (-2.03226413213726e-20,-4.39972734952866e-05) (-4.29705975257430e-21,-7.73408859711252e-04) (1.80620391809104e-03,1.27311966467630e-03) (2.27172799609056e-03,2.14173558366390e-04) (1.27260730092084e-02,-5.63930379640086e-03) (-1.25740754302192e-03,-2.52629792886868e-03) (3.11712763117411e-03,2.65399052555793e-03) (-9.49948039174898e-03,3.92457019390270e-03) (-7.10060132225722e-04,6.47663164341477e-04) (-3.32590854367089e-03,-3.25403313467156e-03) (8.21017901033145e-04,7.28802298133970e-04) (1.66958668674147e-03,1.32110665213213e-03) (-3.81686097082572e-03,-6.17059017724922e-03) (-1.42344750417067e-03,-1.10504975369492e-03) (7.35346158195902e-03,-3.29580517750373e-03) 
(7.87773371716313e-22,-9.61671816228982e-04) (-7.60085312959101e-20,-7.56402319396965e-03) (9.58720747451991e-04,3.23933099028097e-20) (9.58720747451991e-04,2.60334337745676e-20) (8.21228701949305e-02,-8.53878808333309e-21) (-3.02654971487055e-03,1.70886390226571e-19) (-3.02654971487055e-03,1.65667162322279e-19) (8.26758029929171e-02,1.01593971273167e-19) (1.50335923982333e-21,3.86597796060248e-04) (-4.93032673159288e-28,-4.44869107936920e-11) (7.62202317995500e-21,1.77634773551951e-04) (-7.28095388549642e-22,6.69607024870501e-04) (4.72004409820583e-21,1.77634773551951e-04) (8.30221652837401e-03,-3.74383099650222e-03) (1.03323086642701e-02,-4.65301108893189e-03) (3.11712763117411e-03,2.65399052555793e-03) (3.11712763117411e-03,2.65399052555793e-03) (-5.01653442769592e-04,-3.98251313914155e-04) (-3.32590854367089e-03,-3.25403313467156e-03) (-3.32590854367089e-03,-3.25403313467156e-03) (-3.41404353352137e-03,-6.96590514652101e-03) (3.49887239438645e-03,-1.56067009317330e-03) (-3.81686097082572e-03,-6.17059017724922e-03) (1.67098947239918e-03,1.16383379678605e-03) (6.06022475627751e-03,-2.70315989522942e-03) (1.67098947239918e-03,1.16383379678606e-03) 
(-9.28615504716790e-20,8.91023630880739e-04) (9.18458554205289e-20,7.69133743435419e-03) (-3.67244071226990e-01,-6.86292670692792e-19) (-3.77140194849150e-03,8.66025566557487e-20) (-6.79786412847233e-03,2.67062211434346e-19) (-4.18167443506277e-01,7.41810803911107e-19) (1.11458569399724e-02,4.60483928623010e-20) (2.00918472192600e-02,-1.61497023828207e-19) (-1.72993802726593e-20,-1.06086291419169e-03) (4.42934002790852e-20,3.04929478421846e-03) (-7.53755388013570e-20,6.06104565331798e-03) (-3.22644634878270e-20,2.41018605277799e-03) (4.15922448018390e-23,5.65374222994536e-07) (-1.29356948690621e-02,-9.15471965716643e-03) (-1.11815006349312e-02,2.07897803881372e-03) (7.69686057774022e-03,1.58688723965469e-02) (-8.36541031010350e-02,3.70040178888911e-02) (-2.06171224714912e-02,-1.76928200044234e-02) (8.76445920184040e-03,1.89110829763920e-03) (4.02731556071331e-02,-1.55933963668140e-02) (1.73724222780260e-02,1.82882443588506e-02) (1.15821045446951e-02,9.23772127036331e-03) (-1.17187974866053e-02,-9.34233520853913e-03) (-5.21236871496204e-02,2.33562287455120e-02) (1.00488195762600e-04,-4.92678854453438e-04) (2.70572535518723e-02,4.37780855464232e-02) 
(2.42304009910469e-20,8.91023630880739e-04) (3.74878922231158e-19,7.69133743435418e-03) (-3.77140194849150e-03,8.66025566557487e-20) (-3.67244071226990e-01,-2.85045285212852e-19) (-6.79786412847233e-03,1.40104102303459e-19) (1.11458569399724e-02,4.60483928623010e-20) (-4.18167443506277e-01,2.57203257045812e-19) (2.00918472192600e-02,6.84808745290357e-19) (4.23155606146378e-20,2.61771380664853e-03) (1.14313195990753e-19,3.04929478421847e-03) (-6.42868735995560e-23,5.65374222994537e-07) (-9.49195493291332e-22,2.86358792766198e-04) (1.75002705783758e-19,6.06104565331798e-03) (-1.29356948690621e-02,-9.15471965716642e-03) (-1.11815006349312e-02,2.07897803881372e-03) (-8.36541031010350e-02,3.70040178888911e-02) (7.69686057774022e-03,1.58688723965469e-02) (-2.06171224714912e-02,-1.76928200044234e-02) (4.02731556071331e-02,-1.55933963668140e-02) (8.76445920184039e-03,1.89110829763919e-03) (1.73724222780260e-02,1.82882443588506e-02) (-5.70402694203668e-03,-5.04553303904575e-03) (-1.17187974866053e-02,-9.34233520853914e-03) (2.70572535518723e-02,4.37780855464232e-02) (1.00806408628745e-02,7.75376186598776e-03) (-5.21236871496204e-02,2.33562287455120e-02) 
(5.63431481010176e-20,7.28392171043268e-03) (-4.63607427341875e-19,5.07091224303633e-02) (-6.79786412847233e-03,2.67062211434346e-19) (-6.79786412847233e-03,1.40104102303459e-19) (-3.96354455249854e-01,-9.75148235506928e-19) (2.00918472192600e-02,-1.61497023828207e-19) (2.00918472192600e-02,6.84808745290357e-19) (-3.30809802836842e-01,-1.57452375598033e-18) (-2.71456733927368e-20,-2.98262401841328e-03) (-6.42868735995560e-23,5.65374222994537e-07) (1.08732559097637e-19,-1.77107732133751e-03) (-2.08709563242369e-20,-5.16605633976705e-03) (6.84564410569095e-20,-1.77107732133751e-03) (-5.94936473800967e-02,2.68373101559808e-02) (-5.10200958339350e-02,2.30769543198857e-02) (-2.06171224714912e-02,-1.76928200044234e-02) (-2.06171224714912e-02,-1.76928200044234e-02) (2.21308565733664e-03,6.43141894374731e-04) (1.73724222780260e-02,1.82882443588506e-02) (1.73724222780260e-02,1.82882443588506e-02) (2.81624534401495e-02,5.56925776915159e-02) (-2.48224001565949e-02,1.10945541927754e-02) (2.70572535518723e-02,4.37780855464232e-02) (-1.19197738781305e-02,-8.35697749963225e-03) (-4.29936582370280e-02,1.92163315492133e-02) (-1.19197738781305e-02,-8.35697749963226e-03) 
(-2.23489627816395e-03,3.71147572428484e-03) (-7.08055926387039e-03,1.27309059614315e-02) (-3.27583579473072e-03,1.97906603493135e-03) (-3.27583579473073e-03,1.97906603493135e-03) (-1.50303732608938e-02,-6.77610179403154e-03) (5.91314922743144e-03,-2.34498807343220e-03) (5.91314922743144e-03,-2.34498807343220e-03) (2.70318811038195e-02,1.21805574370726e-02) (1.80396629197821e-04,-5.08289741668323e-04) (1.04414973148315e-02,4.79273738928243e-03) (1.97555516790357e-03,-1.43656718674608e-03) (3.12456127284789e-04,-8.80383657535595e-04) (1.97555516790357e-03,-1.43656718674608e-03) (-9.66773336319710e-03,4.04208265445529e-20) (-4.64885188287154e-02,-2.80278409684425e-19) (5.83346506713486e-20,-3.44101586037455e-04) (1.40866619761454e-20,-3.44101586037455e-04) (-3.73560157968366e-20,-4.01821902393477e-03) (-1.51789716628812e-19,1.21436746908368e-03) (1.94332398716500e-20,1.21436746908368e-03) (3.39558643832208e-20,1.41786296442356e-02) (2.48665724163603e-04,-1.96829131956368e-20) (-1.17824541765769e-04,2.74629123827173e-23) (-2.12312488011369e-04,6.01759589072413e-21) (4.30701668352267e-04,4.81199109123085e-21) (-2.12312488011369e-04,4.08899450759911e-21) 
(1.32188222548314e-02,-2.19285557738572e-02) (3.89088249896941e-02,-7.01834438054466e-02) (1.86963885935956e-02,-1.13032287756508e-02) (1.86963885935956e-02,-1.13032287756508e-02) (8.57856997240179e-02,3.86749159661674e-02) (-3.24639262139521e-02,1.24569771571552e-02) (-3.24639262139521e-02,1.24569771571551e-02) (-1.48379368654628e-01,-6.68589073352137e-02) (-9.44098042451478e-04,2.65922628609329e-03) (-6.00333776691774e-02,-2.75128455141183e-02) (-1.15109816816545e-02,8.45240793880238e-03) (-1.63522577685228e-03,4.60591503633627e-03) (-1.15109816816545e-02,8.45240793880238e-03) (1.19528909651574e-01,6.84328194881948e-21) (1.01401587126593e-01,7.41441001437901e-19) (-9.31392048499040e-20,2.10018834695085e-03) (-7.97468775643641e-20,2.10018834695085e-03) (-1.90075185052089e-20,2.30616629955156e-02) (-2.41033720963447e-19,-7.35163662159392e-03) (-2.39380383648312e-18,-7.35163662159392e-03) (1.62099692397420e-19,-8.02885506000432e-02) (-1.28605165008403e-03,-1.94667693165630e-20) (5.72805793574380e-04,4.44815642732849e-20) (1.03216069858926e-03,-1.58357854707749e-20) (-2.22750679910332e-03,-1.83518404849661e-20) (1.03216069858926e-03,-1.72229790068017e-20) 
(-4.18038420244833e-02,6.91557733972684e-02) (-9.41804148257860e-02,1.72420046087660e-01) (-5.24309064272904e-02,3.15433823118449e-02) (-5.24309064272904e-02,3.15433823118449e-02) (-2.40719058246921e-01,-1.08567593274179e-01) (7.82978617761724e-02,-2.49886112959492e-02) (7.82978617761724e-02,-2.49886112959492e-02) (3.58066978523747e-01,1.61487178578466e-01) (2.11361295146210e-03,-5.85363994070313e-03) (1.73971016506751e-01,7.94376033055253e-02) (3.45482773580546e-02,-2.58733653754806e-02) (3.66088501946797e-03,-1.01388017865123e-02) (3.45482773580546e-02,-2.58733653754806e-02) (-5.02689921880510e-01,4.80324694410518e-19) (2.46672729675783e-01,8.71673626240560e-19) (-1.88644539438013e-18,-7.09060083156112e-03) (7.50622104969982e-19,-7.09060083156112e-03) (-4.53455640479224e-19,-6.72002313275119e-02) (-3.87479067331042e-18,2.40757208430013e-02) (3.34646917032700e-18,2.40757208430013e-02) (-7.61840348944847e-19,2.21997872126945e-01) (2.03461766089111e-03,-1.73727663664651e-20) (-6.45190633076499e-04,4.09282420208846e-20) (-1.16591267977391e-03,5.07691905124611e-20) (3.52406116264035e-03,2.51243059630569e-20) (-1.16591267977391e-03,2.70877919992017e-20) 
(-1.80620391809103e-03,1.27311966467630e-03) (-2.27172799609057e-03,2.14173558366391e-04) (-1.25740754302192e-03,2.52629792886868e-03) (1.27260730092084e-02,5.63930379640086e-03) (3.11712763117412e-03,-2.65399052555793e-03) (-7.10060132225724e-04,-6.47663164341476e-04) (-9.49948039174898e-03,-3.92457019390270e-03) (-3.32590854367089e-03,3.25403313467156e-03) (1.64325065008193e-03,-1.32140230821253e-03) (-1.66958668674147e-03,1.32110665213213e-03) (-7.35346158195903e-03,-3.29580517750373e-03) (7.01392828854391e-07,7.86364276730358e-05) (3.81686097082571e-03,-6.17059017724922e-03) (-2.48340518977958e-20,-1.02234217351818e-04) (-1.29256475159404e-19,-9.14989813333637e-04) (7.81146294060923e-02,9.65242472101171e-21) (5.32049801431030e-04,1.78421677343411e-20) (9.58720747451991e-04,3.23933099028097e-20) (9.54093337500756e-02,2.06750557628819e-20) (-1.67960809806152e-03,-7.53280646793652e-20) (-3.02654971487055e-03,1.70886390226571e-19) (-6.41832358184157e-21,9.31694948362761e-05) (-2.89455395064623e-22,-3.21102917973109e-04) (-3.83997312961536e-22,-7.73408859711252e-04) (8.47305457401495e-21,-2.49368845762530e-04) (-3.72420925206049e-27,-4.44869107936920e-11) 
(-1.80620391809104e-03,1.27311966467630e-03) (-2.27172799609057e-03,2.14173558366390e-04) (1.27260730092084e-02,5.63930379640086e-03) (-1.25740754302192e-03,2.52629792886868e-03) (3.11712763117411e-03,-2.65399052555793e-03) (-9.49948039174898e-03,-3.92457019390270e-03) (-7.10060132225723e-04,-6.47663164341478e-04) (-3.32590854367089e-03,3.25403313467156e-03) (-8.21017901033143e-04,7.28802298133969e-04) (-1.66958668674147e-03,1.32110665213213e-03) (3.81686097082571e-03,-6.17059017724922e-03) (1.42344750417067e-03,-1.10504975369492e-03) (-7.35346158195903e-03,-3.29580517750373e-03) (-4.78173138141991e-21,-1.02234217351818e-04) (-9.46313747369266e-20,-9.14989813333637e-04) (5.32049801431030e-04,1.78421677343411e-20) (7.81146294060923e-02,2.24073793290318e-20) (9.58720747451991e-04,2.60334337745676e-20) (-1.67960809806152e-03,-7.53280646793652e-20) (9.54093337500756e-02,8.23015943341441e-20) (-3.02654971487055e-03,1.65667162322279e-19) (4.71636505490217e-21,-2.62544502760893e-04) (1.44224971846458e-20,-3.21102917973109e-04) (-4.93032673159288e-28,-4.44869107936920e-11) (-2.03226413213726e-20,-4.39972734952866e-05) (-4.29705975257430e-21,-7.73408859711252e-04) 
(-8.30221652837400e-03,-3.74383099650222e-03) (-1.03323086642701e-02,-4.65301108893188e-03) (3.11712763117412e-03,-2.65399052555793e-03) (3.11712763117411e-03,-2.65399052555793e-03) (-5.01653442769592e-04,3.98251313914155e-04) (-3.32590854367089e-03,3.25403313467156e-03) (-3.32590854367089e-03,3.25403313467156e-03) (-3.41404353352137e-03,6.96590514652101e-03) (-3.49887239438646e-03,-1.56067009317330e-03) (3.81686097082571e-03,-6.17059017724922e-03) (-1.67098947239918e-03,1.16383379678605e-03) (-6.06022475627751e-03,-2.70315989522942e-03) (-1.67098947239918e-03,1.16383379678605e-03) (7.87773371716313e-22,-9.61671816228982e-04) (-7.60085312959101e-20,-7.56402319396965e-03) (9.58720747451991e-04,3.23933099028097e-20) (9.58720747451991e-04,2.60334337745676e-20) (8.21228701949305e-02,-8.53878808333309e-21) (-3.02654971487055e-03,1.70886390226571e-19) (-3.02654971487055e-03,1.65667162322279e-19) (8.26758029929171e-02,1.01593971273167e-19) (1.50335923982333e-21,3.86597796060248e-04) (-4.93032673159288e-28,-4.44869107936920e-11) (7.62202317995500e-21,1.77634773551951e-04) (-7.28095388549642e-22,6.69607024870501e-04) (4.72004409820583e-21,1.77634773551951e-04) 
(1.29356948690621e-02,-9.15471965716642e-03) (1.11815006349312e-02,2.07897803881372e-03) (7.69686057774022e-03,-1.58688723965469e-02) (-8.36541031010349e-02,-3.70040178888911e-02) (-2.06171224714912e-02,1.76928200044234e-02) (8.76445920184040e-03,-1.89110829763920e-03) (4.02731556071331e-02,1.55933963668140e-02) (1.73724222780260e-02,-1.82882443588506e-02) (-1.15821045446951e-02,9.23772127036331e-03) (1.17187974866053e-02,-9.34233520853914e-03) (5.21236871496203e-02,2.33562287455120e-02) (-1.00488195762600e-04,-4.92678854453438e-04) (-2.70572535518723e-02,4.37780855464233e-02) (-9.28615504716790e-20,8.91023630880739e-04) (9.18458554205289e-20,7.69133743435419e-03) (-3.67244071226990e-01,-6.86292670692792e-19) (-3.77140194849150e-03,8.66025566557487e-20) (-6.79786412847233e-03,2.67062211434346e-19) (-4.18167443506277e-01,7.41810803911107e-19) (1.11458569399724e-02,4.60483928623010e-20) (2.00918472192600e-02,-1.61497023828207e-19) (-1.72993802726593e-20,-1.06086291419169e-03) (4.42934002790852e-20,3.04929478421846e-03) (-7.53755388013570e-20,6.06104565331798e-03) (-3.22644634878270e-20,2.41018605277799e-03) (4.15922448018390e-23,5.65374222994536e-07) 
(1.29356948690621e-02,-9.15471965716643e-03) (1.11815006349312e-02,2.07897803881372e-03) (-8.36541031010349e-02,-3.70040178888911e-02) (7.69686057774022e-03,-1.58688723965469e-02) (-2.06171224714912e-02,1.76928200044234e-02) (4.02731556071331e-02,1.55933963668140e-02) (8.76445920184040e-03,-1.89110829763919e-03) (1.73724222780260e-02,-1.82882443588506e-02) (5.70402694203668e-03,-5.04553303904575e-03) (1.17187974866053e-02,-9.34233520853914e-03) (-2.70572535518723e-02,4.37780855464233e-02) (-1.00806408628745e-02,7.75376186598776e-03) (5.21236871496203e-02,2.33562287455120e-02) (2.42304009910469e-20,8.91023630880739e-04) (3.74878922231158e-19,7.69133743435418e-03) (-3.77140194849150e-03,8.66025566557487e-20) (-3.67244071226990e-01,-2.85045285212852e-19) (-6.79786412847233e-03,1.40104102303459e-19) (1.11458569399724e-02,4.60483928623010e-20) (-4.18167443506277e-01,2.57203257045812e-19) (2.00918472192600e-02,6.84808745290357e-19) (4.23155606146378e-20,2.61771380664853e-03) (1.14313195990753e-19,3.04929478421847e-03) (-6.42868735995560e-23,5.65374222994537e-07) (-9.49195493291332e-22,2.86358792766198e-04) (1.75002705783758e-19,6.06104565331798e-03) 
(5.94936473800967e-02,2.68373101559808e-02) (5.10200958339349e-02,2.30769543198857e-02) (-2.06171224714912e-02,1.76928200044234e-02) (-2.06171224714912e-02,1.76928200044234e-02) (2.21308565733664e-03,-6.43141894374731e-04) (1.73724222780260e-02,-1.82882443588506e-02) (1.73724222780260e-02,-1.82882443588506e-02) (2.81624534401496e-02,-5.56925776915159e-02) (2.48224001565949e-02,1.10945541927754e-02) (-2.70572535518723e-02,4.37780855464233e-02) (1.19197738781305e-02,-8.35697749963226e-03) (4.29936582370280e-02,1.92163315492133e-02) (1.19197738781305e-02,-8.35697749963226e-03) (5.63431481010176e-20,7.28392171043268e-03) (-4.63607427341875e-19,5.07091224303633e-02) (-6.79786412847233e-03,2.67062211434346e-19) (-6.79786412847233e-03,1.40104102303459e-19) (-3.96354455249854e-01,-9.75148235506928e-19) (2.00918472192600e-02,-1.61497023828207e-19) (2.00918472192600e-02,6.84808745290357e-19) (-3.30809802836842e-01,-1.57452375598033e-18) (-2.71456733927368e-20,-2.98262401841328e-03) (-6.42868735995560e-23,5.65374222994537e-07) (1.08732559097637e-19,-1.77107732133751e-03) (-2.08709563242369e-20,-5.16605633976705e-03) (6.84564410569095e-20,-1.77107732133751e-03) 
"""
        with open("unittest-QO_ovlp.dat", "w") as f:
            f.write(contents)
        ao_proj, kvec_d = read_ao_proj("unittest-QO_ovlp.dat")
        os.remove("unittest-QO_ovlp.dat")
        nlocal, nao = ao_proj.shape
        self.assertEqual(nlocal, 26)
        self.assertEqual(nao, 18)
        self.assertTupleEqual(kvec_d.shape, (3,))

    def test_read_kpoints(self):
        contents = """                               nkstot now = 125
K-POINTS DIRECT COORDINATES
KPOINTS     DIRECT_X     DIRECT_Y     DIRECT_Z WEIGHT
    1  -0.40000000  -0.40000000  -0.40000000 0.0080
    2  -0.20000000  -0.40000000  -0.40000000 0.0080
    3   0.00000000  -0.40000000  -0.40000000 0.0080
    4   0.20000000  -0.40000000  -0.40000000 0.0080
    5   0.40000000  -0.40000000  -0.40000000 0.0080
    6  -0.40000000  -0.20000000  -0.40000000 0.0080
"""
        with open("unittest-kpoints", "w") as f:
            f.write(contents)
        kpoints = read_kpoints("unittest-kpoints")
        os.remove("unittest-kpoints")
        table0 = kpoints[0]
        self.assertEqual(table0.shape, (6, 5))
        contents = """                               nkstot now = 8
    KPT             DirectX             DirectY             DirectZ              Weight
    1                   0                   0                   0            0.015625
    2                0.25                0.25                0.25               0.125
    3                 0.5                 0.5                 0.5              0.0625
    4                   0                0.25                0.25             0.09375
    5                 0.5                 0.5                0.25               0.375
    6                 0.5                0.25                0.25              0.1875
    7                   0                 0.5                 0.5            0.046875
    8                 0.5                0.25               -0.25             0.09375
                                nkstot = 64                                                            ibzkpt
    KPT             DirectX             DirectY             DirectZ     IBZ             DirectX             DirectY             DirectZ
    1                   0                   0                   0       1                   0                   0                   0
    2                0.25                   0                   0       2                0.25                0.25                0.25
    3                 0.5                   0                   0       3                 0.5                 0.5                 0.5
    4               -0.25                   0                   0       2                0.25                0.25                0.25
    5                   0                0.25                   0       2                0.25                0.25                0.25
    6                0.25                0.25                   0       4                   0                0.25                0.25
    7                 0.5                0.25                   0       5                 0.5                 0.5                0.25
    8               -0.25                0.25                   0       6                 0.5                0.25                0.25
    9                   0                 0.5                   0       3                 0.5                 0.5                 0.5
    10                0.25                 0.5                   0       5                 0.5                 0.5                0.25
    11                 0.5                 0.5                   0       7                   0                 0.5                 0.5
    12               -0.25                 0.5                   0       5                 0.5                 0.5                0.25
    13                   0               -0.25                   0       2                0.25                0.25                0.25
    14                0.25               -0.25                   0       6                 0.5                0.25                0.25
    15                 0.5               -0.25                   0       5                 0.5                 0.5                0.25
    16               -0.25               -0.25                   0       4                   0                0.25                0.25
    17                   0                   0                0.25       2                0.25                0.25                0.25
    18                0.25                   0                0.25       4                   0                0.25                0.25
    19                 0.5                   0                0.25       5                 0.5                 0.5                0.25
    20               -0.25                   0                0.25       6                 0.5                0.25                0.25
    21                   0                0.25                0.25       4                   0                0.25                0.25
    22                0.25                0.25                0.25       2                0.25                0.25                0.25
    23                 0.5                0.25                0.25       6                 0.5                0.25                0.25
    24               -0.25                0.25                0.25       5                 0.5                 0.5                0.25
    25                   0                 0.5                0.25       5                 0.5                 0.5                0.25
    26                0.25                 0.5                0.25       6                 0.5                0.25                0.25
    27                 0.5                 0.5                0.25       5                 0.5                 0.5                0.25
    28               -0.25                 0.5                0.25       8                 0.5                0.25               -0.25
    29                   0               -0.25                0.25       6                 0.5                0.25                0.25
    30                0.25               -0.25                0.25       5                 0.5                 0.5                0.25
    31                 0.5               -0.25                0.25       8                 0.5                0.25               -0.25
    32               -0.25               -0.25                0.25       5                 0.5                 0.5                0.25
    33                   0                   0                 0.5       3                 0.5                 0.5                 0.5
    34                0.25                   0                 0.5       5                 0.5                 0.5                0.25
    35                 0.5                   0                 0.5       7                   0                 0.5                 0.5
    36               -0.25                   0                 0.5       5                 0.5                 0.5                0.25
    37                   0                0.25                 0.5       5                 0.5                 0.5                0.25
    38                0.25                0.25                 0.5       6                 0.5                0.25                0.25
    39                 0.5                0.25                 0.5       5                 0.5                 0.5                0.25
    40               -0.25                0.25                 0.5       8                 0.5                0.25               -0.25
    41                   0                 0.5                 0.5       7                   0                 0.5                 0.5
    42                0.25                 0.5                 0.5       5                 0.5                 0.5                0.25
    43                 0.5                 0.5                 0.5       3                 0.5                 0.5                 0.5
    44               -0.25                 0.5                 0.5       5                 0.5                 0.5                0.25
    45                   0               -0.25                 0.5       5                 0.5                 0.5                0.25
    46                0.25               -0.25                 0.5       8                 0.5                0.25               -0.25
    47                 0.5               -0.25                 0.5       5                 0.5                 0.5                0.25
    48               -0.25               -0.25                 0.5       6                 0.5                0.25                0.25
    49                   0                   0               -0.25       2                0.25                0.25                0.25
    50                0.25                   0               -0.25       6                 0.5                0.25                0.25
    51                 0.5                   0               -0.25       5                 0.5                 0.5                0.25
    52               -0.25                   0               -0.25       4                   0                0.25                0.25
    53                   0                0.25               -0.25       6                 0.5                0.25                0.25
    54                0.25                0.25               -0.25       5                 0.5                 0.5                0.25
    55                 0.5                0.25               -0.25       8                 0.5                0.25               -0.25
    56               -0.25                0.25               -0.25       5                 0.5                 0.5                0.25
    57                   0                 0.5               -0.25       5                 0.5                 0.5                0.25
    58                0.25                 0.5               -0.25       8                 0.5                0.25               -0.25
    59                 0.5                 0.5               -0.25       5                 0.5                 0.5                0.25
    60               -0.25                 0.5               -0.25       6                 0.5                0.25                0.25
    61                   0               -0.25               -0.25       4                   0                0.25                0.25
    62                0.25               -0.25               -0.25       5                 0.5                 0.5                0.25
    63                 0.5               -0.25               -0.25       6                 0.5                0.25                0.25
    64               -0.25               -0.25               -0.25       2                0.25                0.25                0.25
"""
        with open("unittest-kpoints", "w") as f:
            f.write(contents)
        kpoints = read_kpoints("unittest-kpoints")
        os.remove("unittest-kpoints")
        table0 = kpoints[0]
        table1 = kpoints[1]
        self.assertEqual(table0.shape, (8, 5))
        self.assertEqual(table1.shape, (64, 8))

        with open("unittest-kpoints", "w") as f:
            f.write(contents)
        kpoints = read_kpoints("unittest-kpoints", as_dict=True)
        os.remove("unittest-kpoints")
        self.assertEqual(len(kpoints), 2)
        for key in kpoints.keys():
            self.assertTrue(key in ["Reduced KPOINTS", "KPOINTS before-after correspondence of symmetry reduction"])
        reduced = kpoints["Reduced KPOINTS"]
        self.assertEqual(len(reduced), 5)
        self.assertEqual(len(reduced["index"]), 8)
        correspondence = kpoints["KPOINTS before-after correspondence of symmetry reduction"]
        self.assertEqual(len(correspondence), 8)
        self.assertEqual(len(correspondence["index"]), 64)

    def test_read_istate(self):
        # nspin = 1 case, from examples/scf/lcao_Si2
        contents = """BAND               Energy(ev)               Occupation                Kpoint = 1                        (0 0 0)
    1                 -5.33886                  0.03125
    2                  6.68544                  0.03125
    3                  6.68544                  0.03125
    4                  6.68544                  0.03125
    5                  9.41062                        0
    6                  9.41062                        0
    7                  9.41062                        0
    8                  10.2973                        0
    9                  15.4684                        0
10                   15.981                        0
11                   15.981                        0
12                  24.0591                        0
13                  24.0591                        0
14                  24.0591                        0


BAND               Energy(ev)               Occupation                Kpoint = 2                        (0.25 0.25 0.25)
    1                 -4.54612                     0.25
    2                  2.61208                     0.25
    3                  5.84516                     0.25
    4                  5.84516                     0.25
    5                  8.86596                        0
    6                  10.5133                        0
    7                  10.5133                        0
    8                  14.0081                        0
    9                  15.6813                        0
10                  15.6813                        0
11                  18.1386                        0
12                  23.7732                        0
13                  25.7553                        0
14                  25.7553                        0


BAND               Energy(ev)               Occupation                Kpoint = 3                        (0.5 0.5 0.5)
    1                 -2.94732                    0.125
    2                -0.568338                    0.125
    3                  5.35368                    0.125
    4                  5.35368                    0.125
    5                  8.34837                        0
    6                  10.4845                        0
    7                  10.4845                        0
    8                  17.3753                        0
    9                  18.6916                        0
10                  19.8038                        0
11                  19.8038                        0
12                   20.611                        0
13                   20.611                        0
14                  22.6587                        0


BAND               Energy(ev)               Occupation                Kpoint = 4                        (0 0.25 0.25)
    1                 -4.26166                   0.1875
    2                  3.15379                   0.1875
    3                  4.62975                   0.1875
    4                  4.62975                   0.1875
    5                  7.84453                        0
    6                  10.7321                        0
    7                  12.9122                        0
    8                  12.9122                        0
    9                  15.6898                        0
10                  17.5779                        0
11                  19.2999                        0
12                  24.1411                        0
13                  24.1411                        0
14                  28.8369                        0


BAND               Energy(ev)               Occupation                Kpoint = 5                        (0.5 0.5 0.25)
    1                 -2.58995                     0.75
    2                 0.194109                     0.75
    3                   2.8738                     0.75
    4                  4.28223                     0.75
    5                  8.45557                        0
    6                  11.8375                        0
    7                  12.9216                        0
    8                  13.2967                        0
    9                  18.1019                        0
10                  19.5933                        0
11                  20.6713                        0
12                  21.6336                        0
13                  26.3704                        0
14                  27.2454                        0


BAND               Energy(ev)               Occupation                Kpoint = 6                        (0.5 0.25 0.25)
    1                   -3.287                    0.375
    2                  1.13919                    0.375
    3                  2.60654                    0.375
    4                  5.22434                    0.375
    5                  9.46904                        0
    6                  11.8591                        0
    7                  12.2503                        0
    8                  13.7015                        0
    9                  14.8073                        0
10                  17.4078                        0
11                  21.8056                        0
12                  24.6549                        0
13                  25.7451                        0
14                  27.4814                        0


BAND               Energy(ev)               Occupation                Kpoint = 7                        (0 0.5 0.5)
    1                   -1.178                  0.09375
    2                   -1.178                  0.09375
    3                  3.58638                  0.09375
    4                  3.58638                  0.09375
    5                  7.56725                        0
    6                  7.56725                        0
    7                  17.2541                        0
    8                  17.2541                        0
    9                  20.3909                        0
10                  20.3909                        0
11                  20.7289                        0
12                  20.7289                        0
13                  22.7905                        0
14                  22.7905                        0


BAND               Energy(ev)               Occupation                Kpoint = 8                        (0.5 0.25 -0.25)
    1                -0.969272                   0.1875
    2                -0.969272                   0.1875
    3                  2.47473                   0.1875
    4                  2.47473                   0.1875
    5                   11.548                        0
    6                   11.548                        0
    7                  12.2037                        0
    8                  12.2037                        0
    9                  18.4438                        0
10                  18.4438                        0
11                  22.2824                        0
12                  22.2824                        0
13                  28.3376                        0
14                  28.3376                        0


"""
        with open("unittest-istate.info", "w") as f:
            f.write(contents)
        result = read_istate("unittest-istate.info")
        os.remove("unittest-istate.info")
        self.assertEqual(len(result), 1) # nspin = 1
        self.assertEqual(result[0].shape, (8, 14, 3))
        # nspin = 2, from examples/scf/lcao_Si2, change nspin from 1 to 2
        contents = """BAND       Spin up Energy(ev)               Occupation     Spin down Energy(ev)               Occupation                Kpoint = 1                        (0 0 0)
    1                 -5.33892                 0.015625                 -5.33892                 0.015625
    2                  6.68535                0.0156003                  6.68535                0.0156003
    3                  6.68535                0.0156003                  6.68535                0.0156003
    4                  6.68535                0.0156003                  6.68535                0.0156003
    5                  9.41057                        0                  9.41058                        0
    6                  9.41057                        0                  9.41058                        0
    7                  9.41057                        0                  9.41058                        0
    8                  10.2972                        0                  10.2972                        0
    9                  15.4684                        0                  15.4684                        0
10                   15.981                        0                   15.981                        0
11                   15.981                        0                   15.981                        0
12                   24.059                        0                   24.059                        0
13                   24.059                        0                   24.059                        0
14                   24.059                        0                   24.059                        0
15                  36.4983                        0                  36.4983                        0
16                  36.4983                        0                  36.4983                        0


BAND       Spin up Energy(ev)               Occupation     Spin down Energy(ev)               Occupation                Kpoint = 2                        (0.25 0.25 0.25)
    1                 -4.54618                    0.125                 -4.54618                    0.125
    2                  2.61201                    0.125                  2.61201                    0.125
    3                  5.84508                    0.125                  5.84508                    0.125
    4                  5.84508                    0.125                  5.84508                    0.125
    5                   8.8659                        0                   8.8659                        0
    6                  10.5132                        0                  10.5132                        0
    7                  10.5132                        0                  10.5132                        0
    8                   14.008                        0                   14.008                        0
    9                  15.6813                        0                  15.6813                        0
10                  15.6813                        0                  15.6813                        0
11                  18.1385                        0                  18.1385                        0
12                  23.7731                        0                  23.7731                        0
13                  25.7552                        0                  25.7552                        0
14                  25.7552                        0                  25.7552                        0
15                  30.6044                        0                  30.6044                        0
16                  33.0534                        0                  33.0534                        0


BAND       Spin up Energy(ev)               Occupation     Spin down Energy(ev)               Occupation                Kpoint = 3                        (0.5 0.5 0.5)
    1                 -2.94739                   0.0625                 -2.94739                   0.0625
    2                -0.568392                   0.0625                -0.568392                   0.0625
    3                   5.3536                   0.0625                   5.3536                   0.0625
    4                   5.3536                   0.0625                   5.3536                   0.0625
    5                  8.34831                        0                  8.34831                        0
    6                  10.4844                        0                  10.4844                        0
    7                  10.4844                        0                  10.4844                        0
    8                  17.3753                        0                  17.3753                        0
    9                  18.6915                        0                  18.6915                        0
10                  19.8038                        0                  19.8038                        0
11                  19.8038                        0                  19.8038                        0
12                  20.6109                        0                  20.6109                        0
13                  20.6109                        0                  20.6109                        0
14                  22.6587                        0                  22.6587                        0
15                   30.235                        0                   30.235                        0
16                  33.2213                        0                  33.2213                        0


BAND       Spin up Energy(ev)               Occupation     Spin down Energy(ev)               Occupation                Kpoint = 4                        (0 0.25 0.25)
    1                 -4.26172                  0.09375                 -4.26172                  0.09375
    2                  3.15371                  0.09375                  3.15371                  0.09375
    3                  4.62968                  0.09375                  4.62968                  0.09375
    4                  4.62968                  0.09375                  4.62968                  0.09375
    5                  7.84449              1.76191e-08                  7.84449              1.76185e-08
    6                   10.732                        0                   10.732                        0
    7                  12.9121                        0                  12.9121                        0
    8                  12.9121                        0                  12.9121                        0
    9                  15.6897                        0                  15.6897                        0
10                  17.5778                        0                  17.5778                        0
11                  19.2999                        0                  19.2999                        0
12                   24.141                        0                   24.141                        0
13                   24.141                        0                   24.141                        0
14                  28.8369                        0                  28.8369                        0
15                  31.3889                        0                  31.3889                        0
16                  31.9032                        0                  31.9032                        0


BAND       Spin up Energy(ev)               Occupation     Spin down Energy(ev)               Occupation                Kpoint = 5                        (0.5 0.5 0.25)
    1                 -2.59002                    0.375                 -2.59002                    0.375
    2                  0.19404                    0.375                 0.194039                    0.375
    3                  2.87374                    0.375                  2.87374                    0.375
    4                  4.28216                    0.375                  4.28216                    0.375
    5                  8.45553                        0                  8.45553                        0
    6                  11.8374                        0                  11.8374                        0
    7                  12.9215                        0                  12.9215                        0
    8                  13.2967                        0                  13.2967                        0
    9                  18.1018                        0                  18.1018                        0
10                  19.5932                        0                  19.5932                        0
11                  20.6712                        0                  20.6712                        0
12                  21.6335                        0                  21.6335                        0
13                  26.3704                        0                  26.3704                        0
14                  27.2454                        0                  27.2454                        0
15                  28.3509                        0                  28.3509                        0
16                  29.7558                        0                  29.7558                        0


BAND       Spin up Energy(ev)               Occupation     Spin down Energy(ev)               Occupation                Kpoint = 6                        (0.5 0.25 0.25)
    1                 -3.28707                   0.1875                 -3.28707                   0.1875
    2                  1.13912                   0.1875                  1.13912                   0.1875
    3                  2.60649                   0.1875                  2.60649                   0.1875
    4                  5.22427                   0.1875                  5.22427                   0.1875
    5                    9.469                        0                    9.469                        0
    6                   11.859                        0                   11.859                        0
    7                  12.2502                        0                  12.2502                        0
    8                  13.7014                        0                  13.7014                        0
    9                  14.8073                        0                  14.8073                        0
10                  17.4078                        0                  17.4078                        0
11                  21.8055                        0                  21.8055                        0
12                  24.6548                        0                  24.6548                        0
13                   25.745                        0                   25.745                        0
14                  27.4814                        0                  27.4814                        0
15                  28.0101                        0                  28.0101                        0
16                  29.3105                        0                  29.3105                        0


BAND       Spin up Energy(ev)               Occupation     Spin down Energy(ev)               Occupation                Kpoint = 7                        (0 0.5 0.5)
    1                 -1.17807                 0.046875                 -1.17807                 0.046875
    2                 -1.17807                 0.046875                 -1.17807                 0.046875
    3                  3.58632                 0.046875                  3.58632                 0.046875
    4                  3.58632                 0.046875                  3.58632                 0.046875
    5                  7.56722              3.70503e-05                  7.56722              3.70493e-05
    6                  7.56722              3.70503e-05                  7.56722              3.70493e-05
    7                   17.254                        0                   17.254                        0
    8                   17.254                        0                   17.254                        0
    9                  20.3909                        0                  20.3909                        0
10                  20.3909                        0                  20.3909                        0
11                  20.7289                        0                  20.7289                        0
12                  20.7289                        0                  20.7289                        0
13                  22.7904                        0                  22.7904                        0
14                  22.7904                        0                  22.7904                        0
15                  32.4774                        0                  32.4774                        0
16                  32.4774                        0                  32.4774                        0


BAND       Spin up Energy(ev)               Occupation     Spin down Energy(ev)               Occupation                Kpoint = 8                        (0.5 0.25 -0.25)
    1                -0.969349                  0.09375                 -0.96935                  0.09375
    2                -0.969349                  0.09375                 -0.96935                  0.09375
    3                  2.47468                  0.09375                  2.47468                  0.09375
    4                  2.47468                  0.09375                  2.47468                  0.09375
    5                   11.548                        0                   11.548                        0
    6                   11.548                        0                   11.548                        0
    7                  12.2037                        0                  12.2037                        0
    8                  12.2037                        0                  12.2037                        0
    9                  18.4437                        0                  18.4437                        0
10                  18.4437                        0                  18.4437                        0
11                  22.2824                        0                  22.2824                        0
12                  22.2824                        0                  22.2824                        0
13                  28.3375                        0                  28.3375                        0
14                  28.3375                        0                  28.3375                        0
15                  29.9083                        0                  29.9083                        0
16                  29.9083                        0                  29.9083                        0


"""
        with open("unittest-istate.info", "w") as f:
            f.write(contents)
        result = read_istate("unittest-istate.info")
        os.remove("unittest-istate.info")
        self.assertEqual(len(result), 2)
        result1 = result[0]
        result2 = result[1]
        self.assertTupleEqual(result1.shape, result2.shape)
        self.assertTupleEqual(result1.shape, (8, 16, 3))

    def test_zero_padding(self):
        x = np.array([1, 2, 3, 4], dtype=np.float64)
        y = np.array([1, 1, 1, 1], dtype=np.float64)
        x, y = zero_padding(xmin = 0, xmax = 5, dx = 1, x = x, y = y)
        self.assertEqual(len(x), 5)
        self.assertEqual(len(y), 5)
        self.assertEqual(y.tolist(), [0.0, 1.0, 1.0, 1.0, 1.0])

        x = np.array([1, 2, 3, 7], dtype=np.float64)
        y = np.array([1, 1, 1, 1], dtype=np.float64)
        x, y = zero_padding(xmin = 0, xmax = 10, dx = 1, x = x, y = y)
        self.assertEqual(len(x), 10)
        self.assertEqual(len(y), 10)
        self.assertEqual(y.tolist(), [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0])

        x = np.array([1, 2, 3, 7], dtype=np.float64)
        y = np.array([1, -1, 1, -1], dtype=np.float64)
        x, y = zero_padding(xmin = 0, xmax = 10, dx = 0.5, x = x, y = y)
        self.assertEqual(len(x), 20)
        self.assertEqual(len(y), 20)
        self.assertEqual(y.tolist(), 
                            [0.0, 0.0, 1.0, 0.0, -1.0,
                            0.0, 1.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0, -1.0, 
                            0.0, 0.0, 0.0, 0.0, 0.0])

    def test_read_input(self):
        content = """INPUT_PARAMETERS
#Parameters (1.General)
suffix                         autotest #the name of main output directory
latname                        none #the name of lattice name
stru_file                      STRU #the filename of file containing atom positions
kpoint_file                    KPT #the name of file containing k points
pseudo_dir                     ../../PP_ORB/ #the directory containing pseudo files
orbital_dir                    ../../PP_ORB/ #the directory containing orbital files
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
nbands                         6 #number of bands
nbands_sto                     256 #number of stochastic bands
nbands_istate                  5 #number of bands around Fermi level for get_pchg calulation
symmetry                       0 #the control of symmetry
init_vel                       0 #read velocity from STRU or not
symmetry_prec                  1e-06 #accuracy for symmetry
symmetry_autoclose             1 #whether to close symmetry automatically when error occurs in symmetry analysis
nelec                          0 #input number of electrons
nelec_delta                    0 #change in the number of total electrons
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
cal_force                      0 #if calculate the force at the end of the electronic iteration
out_freq_ion                   0 #the frequency ( >= 0 ) of ionic step to output charge density and wavefunction. 0: output only when ion steps are finished
device                         cpu #the computing device for ABACUS
precision                      double #the computing precision for ABACUS

#Parameters (2.PW)
ecutwfc                        50 ##energy cutoff for wave functions
ecutrho                        200 ##energy cutoff for charge density and potential
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
out_band                       0 #output energy and band structure (with precision 8)
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
basis_type                     lcao #PW; LCAO in pw; LCAO
nb2d                           0 #2d distribution of atoms
gamma_only                     0 #Only for localized orbitals set and gamma point. If set to 1, a fast algorithm is used
search_radius                  -1 #input search radius (Bohr)
search_pbc                     1 #input periodic boundary condition
lcao_ecut                      50 #energy cutoff for LCAO
lcao_dk                        0.01 #delta k for 1D integration in LCAO
lcao_dr                        0.01 #delta r for 1D integration in LCAO
lcao_rmax                      30 #max R for 1D two-center integration table
out_mat_hs                     0 #output H and S matrix (with precision 8)
out_mat_hs2                    0 #output H(R) and S(R) matrix
out_mat_dh                     0 #output of derivative of H(R) matrix
out_mat_xc                     0 #output exchange-correlation matrix in KS-orbital representation
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
smearing_method                gauss #type of smearing_method: gauss; fd; fixed; mp; mp2; mv
smearing_sigma                 0.002 #energy range for smearing

#Parameters (7.Charge Mixing)
mixing_type                    broyden #plain; pulay; broyden
mixing_beta                    0.7 #mixing parameter: 0 means no new charge
mixing_ndim                    8 #mixing dimension in pulay or broyden
mixing_restart                 0 #threshold to restart mixing during SCF
mixing_gg0                     1.5 #mixing parameter in kerker
mixing_beta_mag                -10 #mixing parameter for magnetic density
mixing_gg0_mag                 0 #mixing parameter in kerker
mixing_gg0_min                 0.1 #the minimum kerker coefficient
mixing_angle                   -10 #angle mixing parameter for non-colinear calculations
mixing_tau                     0 #whether to mix tau in mGGA calculation
mixing_dftu                    0 #whether to mix locale in DFT+U calculation
mixing_dmr                     0 #whether to mix real-space density matrix

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
efield_pos_max                 -1 #position of the maximum of the saw-like potential along crystal axis efield_dir
efield_pos_dec                 -1 #zone in the unit cell where the saw-like potential decreases
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
out_current                    0 #output current or not
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
out_wannier_unk                0 #output UNK. file or not
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
bessel_nao_ecut                50.000000 #energy cutoff for spherical bessel functions(Ry)
bessel_nao_tolerence           1e-12 #tolerence for spherical bessel root
bessel_nao_rcut                6 #radial cutoff for spherical bessel functions(a.u.)
bessel_nao_smooth              1 #spherical bessel smooth or not
bessel_nao_sigma               0.1 #spherical bessel smearing_sigma
bessel_descriptor_lmax         2 #lmax used in generating spherical bessel functions
bessel_descriptor_ecut         50.000000 #energy cutoff for spherical bessel functions(Ry)
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

#Parameters (23.Quasiatomic Orbital analysis)
qo_switch                      0 #0: no QO analysis; 1: QO analysis
qo_basis                       hydrogen #type of QO basis function: hydrogen: hydrogen-like basis, pswfc: read basis from pseudopotential
qo_thr                         1e-06 #accuracy for evaluating cutoff radius of QO basis function
"""
        with open("unittest-input.info", "w") as f:
            f.write(content)
        result = read_keyvals_frominput("unittest-input.info")
        os.remove("unittest-input.info")
        
        self.assertEqual(result["suffix"], "autotest")
        self.assertEqual(result["latname"], "none")
        self.assertEqual(result["stru_file"], "STRU")
        self.assertEqual(result["kpoint_file"], "KPT")
        self.assertEqual(result["pseudo_dir"], "../../PP_ORB/")
        self.assertEqual(result["orbital_dir"], "../../PP_ORB/")
        self.assertEqual(result["pseudo_rcut"], "15")
        self.assertEqual(result["pseudo_mesh"], "0")
        self.assertEqual(result["lmaxmax"], "2")
        self.assertEqual(result["dft_functional"], "default")
        self.assertEqual(result["xc_temperature"], "0")
        self.assertEqual(result["calculation"], "scf")
        self.assertEqual(result["esolver_type"], "ksdft")
        self.assertEqual(result["ntype"], "1")
        self.assertEqual(result["nspin"], "1")
        self.assertEqual(result["kspacing"], "0 0 0")
        self.assertEqual(result["min_dist_coef"], "0.2")
        self.assertEqual(result["nbands"], "6")
        self.assertEqual(result["nbands_sto"], "256")
        self.assertEqual(result["nbands_istate"], "5")
        self.assertEqual(result["symmetry"], "0")
        self.assertEqual(result["init_vel"], "0")
        self.assertEqual(result["symmetry_prec"], "1e-06")
        self.assertEqual(result["symmetry_autoclose"], "1")
        self.assertEqual(result["nelec"], "0")
        self.assertEqual(result["nelec_delta"], "0")
        self.assertEqual(result["out_mul"], "0")
        self.assertEqual(result["noncolin"], "0")
        self.assertEqual(result["lspinorb"], "0")
        self.assertEqual(result["kpar"], "1")
        self.assertEqual(result["bndpar"], "1")
        self.assertEqual(result["out_freq_elec"], "0")
        self.assertEqual(result["dft_plus_dmft"], "0")
        self.assertEqual(result["rpa"], "0")
        self.assertEqual(result["printe"], "100")
        self.assertEqual(result["mem_saver"], "0")
        self.assertEqual(result["diago_proc"], "1")
        self.assertEqual(result["nbspline"], "-1")
        self.assertEqual(result["wannier_card"], "none")
        self.assertEqual(result["soc_lambda"], "1")
        self.assertEqual(result["cal_force"], "0")
        self.assertEqual(result["out_freq_ion"], "0")
        self.assertEqual(result["device"], "cpu")
        self.assertEqual(result["precision"], "double")
        self.assertEqual(result["ecutwfc"], "50")
        self.assertEqual(result["ecutrho"], "200")
        self.assertEqual(result["erf_ecut"], "0")
        self.assertEqual(result["erf_height"], "0")
        self.assertEqual(result["erf_sigma"], "0.1")
        self.assertEqual(result["fft_mode"], "0")
        self.assertEqual(result["pw_diag_thr"], "0.01")
        self.assertEqual(result["scf_thr"], "1e-08")
        self.assertEqual(result["scf_thr_type"], "2")
        self.assertEqual(result["init_wfc"], "atomic")
        self.assertEqual(result["init_chg"], "atomic")
        self.assertEqual(result["chg_extrap"], "atomic")
        self.assertEqual(result["out_chg"], "0")
        self.assertEqual(result["out_pot"], "0")
        self.assertEqual(result["out_wfc_pw"], "0")
        # ...

    def test_read_testconfig_fromBohriumpath(self):
        result = read_testconfig_fromBohriumpath("D:/11548850/10430308/tmp/outputs/artifacts/outputs/Ar_23155_pd04_PBEecw100pwclfs1clstrs1scf/00010")
        self.assertNotEqual(result, None)
        frag, system, mpid, ppid, test = result
        self.assertEqual(frag, "Ar_23155_pd04_PBEecw100pwclfs1clstrs1scf")
        self.assertEqual(system, "Ar")
        self.assertEqual(mpid, "23155")
        self.assertEqual(ppid, "pd04")
        self.assertEqual(test, "PBEecw100pwclfs1clstrs1scf")

if __name__ == "__main__":

    unittest.main()

    