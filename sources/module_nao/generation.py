import os
import json

def get_valance_electronic_configuration(pspot_file: str) -> dict:
    """_summary_

    Args:
        pspot_file (str): filename of pseudopotential

    Returns:
        dict: _description_
    """

    clean_line = "start"
    with open(pspot_file, "r") as f:
        while clean_line != "</PP_INPUTFILE>":
            line = f.readline()
            clean_line = line.strip()
            if clean_line == "<PP_HEADER>":
                pass

def generate_reference_structure(reference_structure: dict) -> str:

    """
    Generate reference structure input script for SIAB calculation
    The reference_structure is a dict, with the following format:
    reference_structure = {
        "dimer": {
            "nbands": 9,
            "maxl": 2,
            "nspin": 1,
            "bond_length": [1.8, 2.1, 2.5, 3.0, 4.0]
        },
        "trimer": {
            "nbands": 12,
            "maxl": 2,
            "nspin": 1,
            "bond_length": [2.0, 2.3, 2.7]
        }
        ...
    }
    The output string is like:
    STRU1       dimer       9       2       1      1.8 2.1 2.5 3.0 4.0
    STRU2       trimer      12      2       1      2.0 2.3 2.7
    """
    return_str = ""
    # get the number of reference structures
    nref = len(reference_structure)
    # generate the return_str
    for i, (key, value) in enumerate(reference_structure.items()):
        return_str += "STRU{0}       {1}       {2}       {3}       {4}      ".format(
            i+1, key, value["nbands"], value["maxl"], value["nspin"]
        )
        for j, bond_length in enumerate(value["bond_length"]):
            if j == len(value["bond_length"]) - 1:
                return_str += "{0}\n".format(bond_length)
            else:
                return_str += "{0} ".format(bond_length)

    return return_str

def generate_orbital_configuration(orbital_configurations: dict) -> str:
    """
    Generate orbital configuration input script for SIAB calculation
    The orbital_configurations is a dict, with the following format:
    orbital_configurations = {
        "Level1": {
            "reference_structure": "STRU1",
            "input_orbital": "auto",
            "orbital_configuration": "none",
            "zeta": "1s1p"
        },
        "Level2": {
            "reference_structure": "STRU1",
            "input_orbital": "auto",
            "orbital_configuration": "fix",
            "zeta": "2s2p1d"
        },
        "Level3": {
            "reference_structure": "STRU2",
            "input_orbital": "auto",
            "orbital_configuration": "fix",
            "zeta": "3s3p2d"
        }
        ...
    }
    The output string is like:
     Level1      STRU1           auto        none        1s1p      
     Level2      STRU1           auto        fix         2s2p1d    
     Level3      STRU2           auto        fix         3s3p2d    
    """
    return_str = ""
    # get the number of orbital configurations
    nconf = len(orbital_configurations)
    # generate the return_str
    for i, (key, value) in enumerate(orbital_configurations.items()):
        return_str += "Level{0}      {1}           {2}        {3}         {4}\n".format(
            i+1, value["reference_structure"], value["input_orbital"], value["orbital_configuration"], value["zeta"]
        )

    return return_str

def generate_siab_input(
        exe_mpi: str, exe_pw: str, element: str, ecut: float, rcut: float,
        pseudo_dir: str, pseudo_name: str,
        reference_structure: dict,
        max_steps: int, 
        orbital_configurations: dict,
        save_orbitals: dict
):

    return_str = '''#--------------------------------------------------------------------------------
#1. CMD & ENV
 EXE_mpi      {0}
 EXE_pw       {1}

#-------------------------------------------------------------------------------- 
#2. Electronic calculatation
 element     {2}  # element name 
 Ecut        {3}  # cutoff energy (in Ry)
 Rcut        {4}  # cutoff radius (in a.u.)
 Pseudo_dir  {5}
 Pseudo_name {6}
 sigma       0.01 # energy range for gauss smearing (in Ry)

#--------------------------------------------------------------------------------
#3. Reference structure related parameters for PW calculation
#For the built-in structure types (including 'dimer', 'trimer' and 'tetramer'):
#STRU Name   #STRU Type  #nbands #MaxL   #nspin  #Bond Length list 
{7}

#-------------------------------------------------------------------------------- 
#4. SIAB calculatation
 max_steps    9000
#Orbital configure and reference target for each level
#LevelIndex  #Ref STRU name  #Ref Bands  #InputOrb    #OrbitalConf 
{8}

#--------------------------------------------------------------------------------
#5. Save Orbitals
#Index    #LevelNum   #OrbitalType 
 Save1    Level1      SZ
 Save2    Level2      DZP
 Save3    Level3      TZDP
'''
    return return_str

def automatical_ptg_dpsi(work_status: dict, element: str, pseudopotential: str):

    """ Generate numerical atomic orbitals in batch, use as-prepared SG15-type pseudopotential 
    corresponding orbital generation input script as reference """
    os.chdir(work_status["numerical_orbitals"]["resources"]) # check-in
    if not os.path.isdir("sg15_10"):
        print("warning: sg15_10 not found, will download from web")
        if not os.path.exists("resources.json"):
            print("error: resources.json not found, please check")
        else:
            with open(work_status["work_folder"] + "\\resources.json", "r") as f:
                resources = json.load(f)
                website = resources["numerical_orbitals"]["resources"]["sg15_10"]
                os.system("wget " + website)
                os.system("unzip SG15-Version1p0__AllOrbitals-Version2p0.zip")
                os.system("rm SG15-Version1p0__AllOrbitals-Version2p0.zip")
                os.system("mv SG15-Version1p0__AllOrbitals-Version2p0 sg15_10")
        # after download, check again
        if not os.path.isdir("sg15_10"):
            print("error: sg15_10 not found, please check")
            return
    # check-in
    os.chdir(work_status["work_folder"])

# unit test of this module
if __name__ == "__main__":

    out = generate_reference_structure({
        "dimer": {
            "nbands": 9,
            "maxl": 2,
            "nspin": 1,
            "bond_length": [1.8, 2.1, 2.5, 3.0, 4.0]
        },
        "trimer": {
            "nbands": 12,
            "maxl": 2,
            "nspin": 1,
            "bond_length": [2.0, 2.3, 2.7]
        }
    })
    print(out)

    out = generate_orbital_configuration({
        "Level1": {
            "reference_structure": "dimer",
            "input_orbital": "auto",
            "orbital_configuration": "none",
            "zeta": "1s1p"
        },
        "Level2": {
            "reference_structure": "dimer",
            "input_orbital": "auto",
            "orbital_configuration": "fix",
            "zeta": "2s2p1d"
        },
        "Level3": {
            "reference_structure": "trimer",
            "input_orbital": "auto",
            "orbital_configuration": "fix",
            "zeta": "3s3p2d"
        }
    })
    print(out)