"""qe2abacus provides conversion from Quantum ESPRESSO input script to that of ABACUS.

Following key points are especially considered:
1. basis_type: "pw" or "lcao". lcao is ABACUS specific basis_type, QE doesnot support it. Therefore
                if not specified, basis_type is set to "pw" by default.
2. NUMERICAL_ORBITALS: QE doesnot support numerical orbitals, therefore if not specified, this section
                        is not included in the STRU file.
"""
def control_to_INPUT(qe: dict) -> dict:
    """
    NOT COMPLETED
    """
    abacus = {
        "calculation": "scf",
        "pseudo_dir": "./",
        "orbitals_dir": "./",
        "symmetry": 0,
        "ecutwfc": 100,
        "scf_thr": 1e-6,
        "scf_nmax": 100,
        "basis_type": "pw",
        "smearing_method": "gauss",
        "smearing_sigma": 0.002,
        "cal_stress": False,
        "cal_force": False,
        "dft_functional": "pbe"
    }
    if "calculation" in qe["control"]:
        abacus["calculation"] = qe["control"]["calculation"]
    if "ecutwfc" in qe["control"]:
        abacus["ecutwfc"] = qe["control"]["ecutwfc"]
    if "tprnfor" in qe["control"]:
        abacus["cal_force"] = qe["control"]["tprnfor"]
    if "tstress" in qe["control"]:
        abacus["cal_stress"] = qe["control"]["tstress"]
    if "pseudo_dir" in qe["control"]:
        abacus["pseudo_dir"] = qe["control"]["pseudo_dir"]
    return abacus

def system_to_INPUT(qe: dict) -> dict:
    pass

def electrons_to_INPUT(qe: dict) -> dict:
    pass

def ions_to_INPUT(qe: dict) -> dict:
    pass

def cell_to_INPUT(qe: dict) -> dict:
    pass

def k_points_to_KPT(qe: dict) -> dict:
    pass

def atomic_species_to_STRU(qe: dict) -> dict:
    pass

def atomic_positions_to_STRU(qe: dict) -> dict:
    pass

def lattice_vectors_to_STRU(qe: dict) -> dict:
    pass

