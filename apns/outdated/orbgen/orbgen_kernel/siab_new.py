"""Generate SIAB_INPUT.json for PTG_dpsi version 0.2.0

There are two ways to generate the SIAB_INPUT.json.
1. Manually set almost all information, this is not the one designed for
2. Read ecutwfc convergence test result and set the SIAB_INPUT.json automatically,
   if reference is given, then bond_lengths and nspin, smearing, will be set
   according to the reference.
Here only method 2 is implemented, method 1 is unnecessary to be implemented here.

Basic structure of SIAB_INPUT.json:
```json
{
    "environment": "module load intel/2019.5.281 openmpi/3.1.4 intel-mkl/2019.5.281 intel-mpi/2019.5.281",
    "mpi_command": "mpirun -np 1",
    "abacus_command": "abacus",

    "pseudo_dir": "./download/pseudopotentials/sg15_oncv_upf_2020-02-06/1.2",
    "pesudo_name": "Fe_ONCV_PBE-1.2.upf",
    "ecutwfc": 100,
    "bessel_nao_rcut": [6, 7, 8, 9, 10],
    "smearing_sigma": 0.015,

    "optimizer": "pytorch.SWAT",
    "spillage_coeff": [0.5, 0.5],
    "max_steps": 9000,

    "reference_systems": [
        {
            "shape": "dimer",
            "nbands": 8,
            "nspin": 1,
            "bond_lengths": "auto"
        },
        {
            "shape": "trimer",
            "nbands": 10,
            "nspin": 1,
            "bond_lengths": [1.9, 2.1, 2.6]
        }
    ],
    "orbitals": [
        {
            "zeta_notation": "SZ",
            "shape": "dimer",
            "nbands_ref": 4,
            "orb_ref": "none"
        },
        {
            "zeta_notation": "DZP",
            "shape": "dimer",
            "nbands_ref": 4,
            "orb_ref": "SZ"
        },
        {
            "zeta_notation": "TZDP",
            "shape": "trimer",
            "nbands_ref": 6,
            "orb_ref": "DZP"
        }
    ]
}
```
"""
import re
import os
def parse_reference(fname: str) -> dict:
    """Assume the reference is format in version earlier than 0.2.0"""
    """two pre-requisite function"""
    # --------------------------------------------------------------------
    def wash(inp: dict):
        """one pre-requisite function"""
        # ----------------------------------------------------------------
        def default(inp: dict):
            PTG_DPSI_DEFAULT = {
                "environment": "",
                "mpi_command": "mpirun -np 1",
                "abacus_command": "abacus",
                "optimizer": "pytorch.SWAT",
                "max_steps": 9000,
                "spillage_coeff": [0.5, 0.5]
            }
            # for version < 0.2.0
            if "EXE_opt" not in inp.keys():
                inp["EXE_opt"] = ""
            if inp["EXE_opt"] == "":
                inp["EXE_opt"] = "/opt_orb_pytorch_dpsi/main.py (default)" 
                optimizer_path = "/opt_orb_pytorch_dpsi" 
            else:
                optimizer_path =os.path.dirname(inp["EXE_opt"])
            if "EXE_env" not in inp.keys():
                inp["EXE_env"] = ""
            # for version >= 0.2.0
            for key in PTG_DPSI_DEFAULT.keys():
                if key not in inp.keys():
                    inp[key] = PTG_DPSI_DEFAULT[key]
            return inp, optimizer_path
        # ----------------------------------------------------------------
        inp, optimizer_path = default(inp)

        exe_pw = " ".join([str(word) for word in inp["EXE_pw"]]).replace("\\", "/")
        inp["EXE_pw"] = exe_pw

        exe_mpi = " ".join([str(word) for word in inp["EXE_mpi"]]).replace("\\", "/")
        inp["EXE_mpi"] = exe_mpi

        pseudo_dir = inp["Pseudo_dir"][0].strip().replace("\\", "/")
        if pseudo_dir.endswith("/"):
            pseudo_dir = pseudo_dir[:-1]
        inp["Pseudo_dir"] = pseudo_dir

        fpseudo = inp["Pseudo_name"][0].strip().replace("\\", "").replace("/", "")
        inp["Pseudo_name"] = fpseudo

        return inp
    # --------------------------------------------------------------------
    def compatibility_convert(inp: dict):
        """convert the old version input contents to new version"""
        result = {
            "environment": inp["EXE_env"],
            "mpi_command": inp["EXE_mpi"],
            "abacus_command": inp["EXE_pw"],
            "pseudo_dir": inp["Pseudo_dir"],
            "pseudo_name": inp["Pseudo_name"],
            "ecutwfc": inp["Ecut"],
            "bessel_nao_rcut": inp["Rcut"],
            "smearing_sigma": inp["sigma"],
            "optimizer": inp["optimizer"],
            "max_steps": inp["max_steps"],
            "spillage_coeff": inp["spillage_coeff"],
            "reference_systems": [],
            "orbitals": []
        }
        for key in inp.keys():
            if key.startswith("STRU"):
                result["reference_systems"].append({
                    "shape": inp[key][0],
                    "nbands": int(inp[key][1]), # maxL the useless parameter skipped
                    "nspin": int(inp[key][3]),
                    "bond_lengths": ["auto"] if inp[key][4] == "auto" else [float(v) for v in inp[key][4:]]
                })
            elif key.startswith("Save"):
                result["orbitals"].append({
                    "zeta_notation": inp[key][1],
                    "shape": inp[inp[inp[key][0]][0]][0],
                    "nbands_ref": "auto" if inp[inp[key][0]][1] == "auto" else int(inp[inp[key][0]][1]),
                    "orb_ref": "none" if inp[inp[key][0]][2] == "none" else inp[
                        "Save"+str(int(inp[key][0][5:]) - 1)][1]
                })

        return result
    # --------------------------------------------------------------------
    keyvalue_pattern = r"^(\w+)(\s+)([^#]*)(#.*)?"
    float_pattern = r"^\d+\.\d*$"
    int_pattern = r"^\d+$"
    scalar_keywords = ["Ecut", "sigma", "element", "max_steps"]
    result = {}
    if fname == "":
        raise ValueError("No filename provided")
    with open(fname, "r") as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        _match = re.match(keyvalue_pattern, line)
        if _match:
            key = _match.group(1).strip()
            value = _match.group(3).strip().split()
            value = [float(v) if re.match(float_pattern, v) else int(v) if re.match(int_pattern, v) else v for v in value]
            result[key] = value if key not in scalar_keywords else value[0]
    
    return compatibility_convert(wash(result))

def is_complete(setting: dict) -> bool:
    """Check if the setting is complete"""
    necessary_keys = [
        "environment",
        "mpi_command",
        "abacus_command",
        "pseudo_dir",
        "pseudo_name",
        "reference_systems",
        "orbitals"
    ]
    for key in necessary_keys:
        if key not in setting.keys():
            return False
    return True

def generate(rcuts: list,
             ecutwfc: str,
             pseudo_name: str,
             pseudo_dir: str,
             other_settings: dict,
             fref: str) -> dict:
    """generate one single SIAB_INPUT.json file based on existing reference file
    SIAB_INPUT in old version
    
    Parameters:
    rcuts: list, the cutoff radius for Sphbes, will be directly transferred to bessel_nao_rcut
    ecutwfc: float, the energy cutoff for wavefunctions
    pseudo_name: str, the pseudopotential file name, note that its path should be fed by pseudo_dir
    pseudo_dir: str, the directory for pseudopotential file
    other_settings: dict, other settings, can overwrite the reference settings
    fref: str, the reference file path
    """
    if not os.path.exists(fref):
        raise FileNotFoundError("The reference file does not exist")
    ref = parse_reference(fref)

    if ecutwfc < ref["ecutwfc"]:
        print("""Warning: ecutwfc is smaller than the reference, turn to use the reference value
However, you can also select to use the value from ecutwfc convergence test if
result is available.""")
        ecutwfc = ref["ecutwfc"]

    ref["bessel_nao_rcut"] = rcuts    
    ref["ecutwfc"] = ecutwfc
    ref["pseudo_dir"] = pseudo_dir
    ref["pseudo_name"] = pseudo_name

    # overwritten manner to update the ref
    if isinstance(other_settings, dict):
        for key in other_settings.keys():
            if key in ref.keys():
                ref[key] = other_settings[key]

    return ref