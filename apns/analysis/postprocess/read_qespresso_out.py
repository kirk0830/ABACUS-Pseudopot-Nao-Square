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

def read_bstraj(fout):
    """qe stdout is like (nspin 1):
    ```
     End of self-consistent calculation

          k =-0.0000 0.0000 0.0000 (   273 PWs)   bands (ev):

    -3.9000  18.6222  18.6520  18.6520  21.5375  24.8623  26.9277  26.9277
    31.6146  31.8351


          k = 0.0000-0.0000 0.5449 (   272 PWs)   bands (ev):

     1.3413   3.1273  20.9290  21.5726  22.0044  23.8594  25.8224  27.4089
    28.6819  29.4710
...
          k = 0.5226-0.4541-0.4541 (   265 PWs)   bands (ev):

     2.0535   9.6825  12.8988  14.3271  18.8888  20.7050  22.8824  28.0826
    31.7473  40.1283

     the Fermi energy is    10.3074 ev
    ```
    nspin 2:
    ```
     End of self-consistent calculation

 ------ SPIN UP ------------


          k =-0.0000 0.0000 0.0000 (   273 PWs)   bands (ev):

    -3.5782  20.7255  21.3223  21.3223  24.1388  25.0896  26.3857  26.3857
    31.3641  31.7074
...

          k = 0.5000-0.4782 0.0956 (   258 PWs)   bands (ev):

     6.0838   6.5519  10.6273  10.8943  26.5804  27.3612  30.7786  33.8417
    34.9311  35.7309

 ------ SPIN DOWN ----------


          k =-0.0000 0.0000 0.0000 (   273 PWs)   bands (ev):

    -3.5782  20.7255  21.3223  21.3223  24.1389  25.0896  26.3857  26.3857
    31.3641  31.7073
...
          k = 0.5000-0.4782 0.0956 (   258 PWs)   bands (ev):

     6.0838   6.5519  10.6273  10.8943  26.5804  27.3612  30.7786  33.8417
    34.9311  35.7310

     the Fermi energy is    11.1974 ev
    ```
    """
    import re
    # first read the whole file as one string
    with open(fout, 'r') as f:
        content = f.read()
    # then use re.split to split the string into blocks
    beg = r"End of self-consistent calculation"
    end = r"the Fermi energy is\s*[-|\s]\d+\.\d+\s*ev"
    # all blocks will start with the beginning pattern and end with the ending pattern
    blocks = re.split(beg, content, flags = re.MULTILINE)[1:]
    blocks = [re.split(end, b, flags = re.MULTILINE) for b in blocks]
    blocks = [b[0] for b in blocks]
    blocks = [[line for line in block.split("\n") if line.strip() != ""] for block in blocks]

    kstraj, bstraj = [], []
    for scf in blocks: # iterate on frames of trajectory
        buf, data = [], []
        for line in scf:
            if "SPIN" in line:
                data.append(buf)
                buf = []
            else:
                buf.append(line)
        data.append(buf)
        data = [d for d in data if len(d) > 0]
        ks_spin, bs_spin = [], []
        for spin in data:
            ks, bs = read_rks_bs(spin)
            ks_spin.append(ks)
            bs_spin.append(bs)
        kstraj.append(ks_spin)
        bstraj.append(bs_spin)
    return kstraj, bstraj

def read_rks_bs(lines: list):
    """read the band structure from the output of pw.x for nspin = 1 case"""
    import re
    import numpy as np
    pat = r"(\s*k =)([\s|-]?0\.\d+)([\s|-]?0\.\d+)([\s|-]?0\.\d+)(\s+\(\s+\d+\s+PWs\)\s+bands \(ev\):)(\s*)"
    buf, ks, bs = "", [], []
    i = 0
    while i < len(lines):
        if re.match(pat, lines[i]):
            buf = ""
            ks.append([float(n) for n in re.findall(r"[\s|-]0\.\d+", lines[i])])
            i += 1
            while not re.match(pat, lines[i]) and not "the Fermi energy is" in lines[i]:
                buf += "  " + lines[i]
                i += 1
                if i >= len(lines):
                    break
            bs.append([float(n) for n in buf.split() if re.match(r"[\s|-]?\d+\.\d+", n)])
        else:
            i += 1

    assert all([len(k) == 3 for k in ks]), "kpoints are not in 3D"
    return np.array(ks), np.array(bs)

def read_bs(fout):
    """only return the final band structure"""
    kstraj, bstraj = read_bstraj(fout)
    if len(kstraj) > 0 and len(bstraj) > 0:
        return kstraj[-1], bstraj[-1]
    return None, None

def to_abacus_istate(bs):
    """convert the band structure to the format of Abacus"""
    # Abacus format:
    # index ekb  occ
    # 1     -3.9 1.0
    # ...
    # however, can only be 1.0 for occ because qespresso does not provide occupation there
    if bs is None: return None
    return [[[[ib+1, e, 1.0] for ib, e in enumerate(bs[ispin][ik])] for ik in range(len(bs))] for ispin in range(len(bs))]

def read_istate(fout):
    import numpy as np
    bs = read_bs(fout)[1]
    if bs is None: return None
    return [np.array(b) for b in to_abacus_istate(bs)]

def read_etraj(fout, unit = "eV", term = "EKS"):
    """qe stdout is like:
    ```
     the Fermi energy is    11.1974 ev

!    total energy              =      -7.83826960 Ry
     estimated scf accuracy    <       0.00000001 Ry
     smearing contrib. (-TS)   =      -0.00003168 Ry
     internal energy E=F+TS    =      -7.83823792 Ry

     The total energy is F=E-TS. E is the sum of the following terms:
     one-electron contribution =       4.73249560 Ry
     hartree contribution      =       0.03475768 Ry
     xc contribution           =      -2.47216436 Ry
     ewald contribution        =     -10.13332683 Ry
    ```
    available terms: EKS, Fermi, Ehartree, Exc, Ewald
    """
    import re
    with open(fout, 'r') as f:
        lines = f.readlines()
    term = "fermi" if term in ["Fermi", "fermi", "FERMI", "efermi", "EFERMI", "ef"] else \
              "kohn-sham" if term in ["KohnSham", "kohnsham", "KOHN", "kohn", "KOHNSHAM", 
                "ekohnsham", "EKOHN", "EKOHNSHAM", "eks", "EKS", "e", "E", "energy"] else \
                "hartree" if term in ["hartree", "HARTREE", "ehartree", "EHARTREE", "eh", "EH"] else \
                "xc" if term in ["xc", "XC", "exc", "EXC"] else \
                "ewald" if term in ["ewald", "EWALD", "eewald", "EEWALD"] else None
    pat = r"the Fermi energy is\s+([\d\.\-]+)\s+ev" if term == "fermi" else \
            r"total energy\s+=\s+([\d\.\-]+)\s+Ry" if term == "kohn-sham" else \
            r"hartree contribution\s+=\s+([\d\.\-]+)\s+Ry" if term == "hartree" else \
            r"xc contribution\s+=\s+([\d\.\-]+)\s+Ry" if term == "xc" else \
            r"ewald contribution\s+=\s+([\d\.\-]+)\s+Ry" if term == "ewald" else None
    if pat is None:
        raise ValueError("term not recognized")
    etraj = []
    for line in lines:
        m = re.search(pat, line)
        if m:
            etraj.append(float(m.group(1)))
    return [unit_conversion(e, "Ry", unit) for e in etraj]

def read_e(fout, unit = "eV", term = "EKS"):
    """return the final one of the energy trajectory"""
    traj = read_etraj(fout, unit, term)
    if len(traj) == 0:
        return None
    return traj[-1]

def read_presstraj(fout):
    """qe stdout is like:
    ```
     Computing stress (Cartesian axis) and pressure

          total   stress  (Ry/bohr**3)                   (kbar)     P=       25.90
  -0.00017791   0.00000000   0.00000000          -26.17        0.00        0.00
   0.00000000   0.00035308  -0.00000000            0.00       51.94       -0.00
   0.00000000  -0.00000000   0.00035308            0.00       -0.00       51.94
    ```
    """
    import re
    with open(fout, 'r') as f:
        lines = f.readlines()
    presstraj = []
    for line in lines:
        m = re.search(r"P=\s+([\d\.\-]+)", line)
        if m:
            presstraj.append(float(m.group(1)))
    return presstraj

def read_press(fout):
    """return the final pressure"""
    presstraj = read_presstraj(fout)
    if len(presstraj) == 0:
        return None
    return presstraj[-1]

def read_natom(fout):
    """qe stdout is like:
    ```
    number of atoms/cell      =          2
    ```
    """
    import re
    with open(fout, 'r') as f:
        lines = f.readlines()
    for line in lines:
        m = re.search(r"number of atoms/cell\s+=\s+(\d+)", line)
        if m:
            return int(m.group(1))
    return None

def read_voltraj(fout, unit = "A"):
    """qe stdout is like:
    ```
     unit-cell volume          =      94.3341 (a.u.)^3
     ...
     new unit-cell volume =     95.57987 a.u.^3 (    14.16348 Ang^3 )
     new unit-cell volume =     96.69895 a.u.^3 (    14.32931 Ang^3 )
     ...
     new unit-cell volume =     10.0000 a.u.^3 (  100.0000 Ang^3 )
     ...
     new unit-cell volume =     99.93738 a.u.^3 (    14.80919 Ang^3 )
     unit-cell volume          =      99.9374 (a.u.)^3
    ```
    will grep both these two lines
    """
    import re
    with open(fout, 'r') as f:
        lines = f.readlines()
    voltraj = []
    for line in lines:
        m_init = re.match(r"unit-cell volume\s+=\s+([\d\.\-]+)\s+\(a\.u\.\)\^3", line)
        m_traj = re.match(r"new unit-cell volume\s+=\s+([\d\.\-]+)\s+a\.u\.\^3\s+\(\s+([\d\.\-]+)\s+Ang\^3\s+\)", line)
        if m_init or m_traj:
            voltraj.append(float(m_init.group(1) if m_init else m_traj.group(1)))
    fac_ = unit_conversion(1.0, "a.u.", unit)
    return [vol * fac_**3 for vol in voltraj]

def read_vol(fout, unit = "A"):
    """return the final volume"""
    voltraj = read_voltraj(fout, unit)
    if len(voltraj) == 0:
        return None
    return voltraj[-1]

import unittest
class TestReadQespressoOut(unittest.TestCase):
    def test_read_rks_bs(self):
        lines = """
     End of self-consistent calculation

          k =-0.0000 0.0000 0.0000 (   273 PWs)   bands (ev):

    -3.9000  18.6222  18.6520  18.6520  21.5375  24.8623  26.9277  26.9277
    31.6146  31.8351


          k = 0.0000-0.0000 0.5449 (   272 PWs)   bands (ev):

     1.3413   3.1273  20.9290  21.5726  22.0044  23.8594  25.8224  27.4089
    28.6819  29.4710
    """
        ks, bs = read_rks_bs(lines.split("\n"))
        self.assertEqual(ks.shape, (2, 3))
        self.assertEqual(bs.shape, (2, 10))

    def test_read_bstraj(self):
        import uuid, os
        # nspin = 1
        lines = """
     End of self-consistent calculation

          k =-0.0000 0.0000 0.0000 (   273 PWs)   bands (ev):

    -3.9000  18.6222  18.6520  18.6520  21.5375  24.8623  26.9277  26.9277
    31.6146  31.8351


          k = 0.0000-0.0000 0.5449 (   272 PWs)   bands (ev):

     1.3413   3.1273  20.9290  21.5726  22.0044  23.8594  25.8224  27.4089
    28.6819  29.4710


          k = 0.5226-0.4541-0.4541 (   265 PWs)   bands (ev):

     2.0535   9.6825  12.8988  14.3271  18.8888  20.7050  22.8824  28.0826
    31.7473  40.1283

     the Fermi energy is    10.3074 ev
    """
        fout = str(uuid.uuid4())
        with open(fout, 'w') as f:
            f.write(lines)
        
        kstraj, bstraj = read_bstraj(fout)
        os.remove(fout)
        ks, bs = kstraj[-1], bstraj[-1]
        self.assertEqual(len(ks), 1) # only one spin
        self.assertEqual(len(bs), 1)
        self.assertEqual(ks[0].shape, (3, 3))
        self.assertEqual(bs[0].shape, (3, 10))
        
        # nspin = 2
        lines = """
     End of self-consistent calculation

 ------ SPIN UP ------------


          k =-0.0000 0.0000 0.0000 (   273 PWs)   bands (ev):

    -3.5782  20.7255  21.3223  21.3223  24.1388  25.0896  26.3857  26.3857
    31.3641  31.7074

          k = 0.5000-0.4782 0.0956 (   258 PWs)   bands (ev):

     6.0838   6.5519  10.6273  10.8943  26.5804  27.3612  30.7786  33.8417
    34.9311  35.7309

 ------ SPIN DOWN ----------


          k =-0.0000 0.0000 0.0000 (   273 PWs)   bands (ev):

    -3.5782  20.7255  21.3223  21.3223  24.1389  25.0896  26.3857  26.3857
    31.3641  31.7073

          k = 0.5000-0.4782 0.0956 (   258 PWs)   bands (ev):

     6.0838   6.5519  10.6273  10.8943  26.5804  27.3612  30.7786  33.8417
    34.9311  35.7310

     the Fermi energy is    11.1974 ev
    """
        fout = str(uuid.uuid4())
        with open(fout, 'w') as f:
            f.write(lines)
        #read_bstraj(fout)
        kstraj, bstraj = read_bstraj(fout)
        os.remove(fout)
        ks, bs = kstraj[-1], bstraj[-1]
        self.assertEqual(len(ks), 2) # two spins
        self.assertEqual(len(bs), 2)
        self.assertEqual(ks[0].shape, (2, 3))
        self.assertEqual(bs[0].shape, (2, 10))
        self.assertEqual(ks[1].shape, (2, 3))
        self.assertEqual(bs[1].shape, (2, 10))

    def test_to_abacus_istate(self):
        lines = """
     End of self-consistent calculation

          k =-0.0000 0.0000 0.0000 (   273 PWs)   bands (ev):

    -3.9000  18.6222  18.6520  18.6520  21.5375  24.8623  26.9277  26.9277
    31.6146  31.8351


          k = 0.0000-0.0000 0.5449 (   272 PWs)   bands (ev):

     1.3413   3.1273  20.9290  21.5726  22.0044  23.8594  25.8224  27.4089
    28.6819  29.4710
    """
        ks, bs = read_rks_bs(lines.split("\n"))
        istate = to_abacus_istate(bs)
        print(istate)

if __name__ == "__main__":
    unittest.main()