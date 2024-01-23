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

def grep_band(fname: str) -> list|bool:
    """grep band energies from running_scf.log"""
    pass