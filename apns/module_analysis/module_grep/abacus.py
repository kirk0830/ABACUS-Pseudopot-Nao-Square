def is_normal_end(fname: str) -> bool:
    """check if the job is normally ended for Abacus"""
    with open(fname, "r") as f:
        for line in f:
            if line.startswith("TIME STATISTICS"):
                return True
    return False

def fname_setting(job_dir = "", calculation = "scf", suffix = "ABACUS", include_path = False):
    """return the file names (stdout, log and input script) for Abacus"""
    job_dir = job_dir if not job_dir.endswith("/") else job_dir[:-1]
    if include_path:
        return job_dir+"/out.log", job_dir+"/OUT.%s/running_%s.log"%(suffix, calculation), job_dir+"/INPUT"
    else:
        return "out.log", "running_%s.log"%calculation, "INPUT"

def return_all_greps():
    return {
        "energy": grep_energy,
        "pressure": grep_pressure,
        "natom": grep_natom,
        "ecutwfc": grep_ecutwfc
    }

def grep_energy(fname: str) -> float|bool:
    line = ""
    with open(fname, "r") as f:
        while True:
            line = f.readline()
            if not line:
                print("Abnormal job for test-case: " + fname + ", it is forced stopped.")
                return False
            line = line.strip()
            if line.startswith("TIME STATISTICS"):
                print("warning: scf calculation not converged in " + fname)
                return False
            if line.startswith("!"):
                break
    if "convergence has not been achieved" in line:
        print("warning: scf calculation not converged in " + fname)
        return False
    else:
        return float(line.split()[-2])

def grep_pressure(fname: str) -> float|bool:
    line = ""
    with open(fname, "r") as f:
        while True:
            line = f.readline()
            if not line:
                print("Abnormal job for test-case: " + fname + ", it is forced stopped.")
                return False
            line = line.strip()
            if line.startswith("TIME STATISTICS"):
                print("warning: scf calculation not converged in " + fname)
                return False
            if line.startswith("TOTAL-PRESSURE"):
                break
    return float(line.split()[-2])

def grep_natom(fname: str) -> int|bool:
    line = ""
    with open(fname, "r") as f:
        while True:
            line = f.readline()
            if not line:
                print("Abnormal job for test-case: " + fname + ", it is forced stopped.")
                return False
            line = line.strip()
            if line.startswith("TIME STATISTICS"):
                print("warning: scf calculation not converged in " + fname)
                return False
            if line.startswith("TOTAL ATOM NUMBER"):
                break
    return int(line.split()[-1])

def grep_ecutwfc(fname: str) -> float|bool:
    """grep ecutwfc from INPUT file"""
    with open(fname, "r") as f:
        lines = f.readlines()
    
    for line in lines:
        if line.startswith("ecutwfc"):
            return float(line.split()[1])
    
    return False

def AssertListEqual(list1: list, list2: list, nplaces: int = 7):
    """Assert two lists are equal"""
    if len(list1) != len(list2):
        raise AssertionError("Two lists have different lengths.")
    for i in range(len(list1)):
        if abs(list1[i] - list2[i]) > 10**(-nplaces):
            print("Two lists are not equal: ", list1, list2)
            raise AssertionError("Two lists are not equal.")

import re
def grep_band(fname: str) -> list|bool:
    """grep band energies from running_scf.log"""

    float_pattern = r"([0-9\-]+\.[0-9\-]+[eEdD][\+\-][0-9]+)|([0-9\-]+\.[0-9\-]+)"
    
    kpt_pattern = r"^(\s*)([0-9]+)(\s*)(%s)(\s+)(%s)(\s*)(%s)(\s*)(%s)(\s*)$"%(float_pattern, float_pattern, float_pattern, float_pattern)
    band_pattern = r"^(\s*)([0-9]+)(\s*)(%s)(\s+)(%s)(\s*)$"%(float_pattern, float_pattern)
    band_header_pattern = r"^(.*)(\=\s*)(%s)(\s+)(%s)(\s+)(%s)(.*)$"%(float_pattern, float_pattern, float_pattern)

    kpoints, kpt_wts, bands, occupations = [], [], [], []
    nelectrons, efermi = 0, 0
    with open(fname, "r") as f:
        read_band, read_kpoint, nkpts, nbands = False, False, 0, 0
        for line in f:
            line = line.strip()
            # not converged
            if line.startswith("TIME STATISTICS"):
                print("warning: scf calculation not converged in " + fname)
                return False
            
            if line.startswith("k-point number in this process"):
                nkpts = int(line.split()[-1])
                continue

            if line.startswith("NBANDS"):
                nbands = int(line.split()[-1])
                continue

            if line.startswith("K-POINTS CARTESIAN COORDINATES"):
                read_kpoint = True
                continue

            if line.startswith("EFERMI"):
                efermi = float(line.split()[-2])
                continue

            if line.startswith("AUTOSET number of electrons"):
                nelectrons = int(line.split()[-1])
                continue

            if read_kpoint:
                _match = re.match(kpt_pattern, line)
                if _match:
                    nkpts -= 1
                    kpoints.append([float(_match.group(4)), float(_match.group(6)), float(_match.group(8))])
                    kpt_wts.append(float(_match.group(16)))
                    continue
                else:
                    if nkpts > 0:
                        continue
                    else:
                        read_kpoint = False
                        #print("Read kpoints done: ", tuple(zip(kpoints, kpt_wts)))
                        nkpts = len(kpoints)
                        continue

            if line.startswith("STATE ENERGY(eV) AND OCCUPATIONS"):
                read_band = True
                if nbands == 0:
                    print("warning: nbands is 0 in " + fname)
                    return False
                else:
                    continue

            if read_band:
                if re.match(band_header_pattern, line):
                    _match = re.match(band_header_pattern, line)
                    #print(line)
                    index = int(line.split("/")[0]) - 1
                    AssertListEqual(
                        [float(_match.group(3)), float(_match.group(5)), float(_match.group(7))],
                        kpoints[index],
                        nplaces = 4
                    )
                    bands.append([])
                    occupations.append([])
                elif re.match(band_pattern, line):
                    _match = re.match(band_pattern, line)
                    bands[-1].append(float(_match.group(4)))
                    occupations[-1].append(float(_match.group(6)))
                else:
                    if len(bands) == nkpts:
                        read_band = False
                        #print("Read bands done.")
                        band_energies = []
                        for ik, band in enumerate(bands):
                            band_energies.extend(list(zip(band, [kpt_wts[ik]]*nbands)))
                        print(band_energies)
                        return band_energies, nelectrons, nbands, efermi
                    elif len(bands[-1]) == nbands:
                        #print("Read bands done for kpoint: ", kpoints[len(bands)-1], " with ", len(bands[-1]), " bands.")
                        continue
                    else:
                        #print("Unexpected error: len(bands[-1]) != nbands. Present line: ", line)
                        return False
                    
    return False

