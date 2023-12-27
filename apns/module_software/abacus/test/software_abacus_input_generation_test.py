from apns.module_software.abacus.generation import STRU_cif, INPUT_generation, KPT_generation
from apns.module_structure.crystal_information_file import read_1 as read_cif

if __name__ == "__main__":

    pseudopotentials = {
        "Er": "Er-sp.PD04.PBE.UPF",
        "O": "O_ONCV_PBE-1.0.upf"
    }

    numerical_orbitals = {
        "Er": "Er_numerical_orbital",
        "O": "O_numerical_orbital"
    }

    cif = read_cif("C:/Users/82108/Downloads/Er2O3_mp-13066_computed.cif")
    STRU = STRU_cif(
        "C:/Users/82108/Downloads/Er2O3_mp-13066_computed.cif", 
        pseudopotentials=pseudopotentials, 
        numerical_orbitals=numerical_orbitals
    )
    INPUT = INPUT_generation(
        basis_type="lcao",
        calculation="scf",
        ecutwfc=100,
        dft_functional="pbe",
        scf_thr=1.0e-8,
        ks_solver="genelpa",
        nspin=1,
        pseudo_dir="./",
        orbital_dir="./",
        symmetry = 0,
        smearing_method = "gaussian",
        smearing_sigma = 0.05
    )
    KPT = KPT_generation(
        mode = "crystal",
        gamma_centered = True,
        metallic = False,
        cell_parameters = [
            cif["cell_parameters"]["a"],
            cif["cell_parameters"]["b"],
            cif["cell_parameters"]["c"]
            ]
    )
    with open("INPUT", "w") as f:
        f.write(INPUT)
    with open("KPT", "w") as f:
        f.write(KPT)
    with open("STRU_ErO", "w") as f:
        f.write(STRU)