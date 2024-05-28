      
"""It is the time to scan all available pseudopotentials"""
import json
import os
def load_ecutwfc(fecut: str):
    """load the ecutwfc from the file

    Args:
        fecut (str): file name including path of the ecutwfc

    Returns:
        _type_: _description_
    """
    with open(fecut, "r") as f:
        ecutwfc = json.load(f)
    return ecutwfc

def load_pseudo(fpspot: str):
    """load the pseudopotential database from the file

    Args:
        fpspot (str): file name including path of the pseudopotential database

    Returns:
        _type_: _description_
    """
    with open(fpspot, "r") as f:
        pspot = json.load(f)
    return pspot

def pspot_eq(val1: str, val2: str):
    """due to test name encoding, all underlines and dots are removed in the test name,
    however in pseudopotenital database, all underlines and dots are unchanged. It is
    needed to identify if the val1 and val2 indicate the same pseudopotential"""

    abbr = val1 if len(val1) < len(val2) else val2
    full = val1 if len(val1) > len(val2) else val2

    return full.replace("_", "").replace(".", "") == abbr

def attributes(val: str):
    """there are three attributes, kind, version and appendix"""
    words = val.split("_")
    words = [words[0], words[1], "_".join(words[2:])] if len(words) > 2 else words
    # map words to results with length 3, absent will be filled with empty string
    results = words + [""] * (3 - len(words))
    return results

def tranverse(ecut_db: dict, pspot_db: dict, included: list = None, excluded: list = None):
    included = included if included else ecut_db.keys()
    excluded = excluded if excluded else []
    for element in included:
        if element in excluded:
            continue
        ecuts = ecut_db[element] # many ecutwfcs indiced by abbr of pseudopotential name
        for abbr in ecuts.keys():
            for pspot in pspot_db[element].keys():
                if pspot_eq(abbr, pspot):
                    kind, version, appendix = attributes(pspot)
                    fpseudo = pspot_db[element][pspot]
                    assert os.path.exists(fpseudo), f"{fpseudo} does not exist"
                    ecut = ecuts[abbr]
                    print(f"Element: {element}, Pseudopotential: {pspot}, ecutwfc: {ecut}")
                    finp = write_apns_inp(element, kind, version, appendix, ecut)
                    os.system(f"python3 main.py -i {finp}")
                    os.remove(finp)
                    break

def write_apns_inp(element: str, kind: str, version: str, appendix: str, ecut: float):
    """Generate APNS input file specifically for EOS test on ideal Bravis lattice

    Args:
        element (str): element to test
        kind (str): pseudopotential kind, e.g. sg15
        version (str): pseudopotential version, e.g. 1.0
        appendix (str): pseudopotential appendix, e.g. fr
        ecut (float): ecutwfc value
    """
    result = {
        "global": {
            "test_mode": "pseudopotential",
            "software": "abacus",
            "work_dir": "./",
            "pseudo_dir": "./download/pseudopotentials",
            "orbital_dir": "./download/numerical_orbitals",
            "save_log": True
        },
        "calculation": {
            "basis_type": "pw",
            "ecutwfc": ecut,
            "cal_force": 1,
            "cal_stress": 1,
            "scf_nmax": 500,
            "nspin": 1,
            "symmetry": 1,
            "calculation": "relax"
        },
        "extensive": {
            "characteristic_lengths": [-0.06, -0.04, -0.02, 0.0, 0.02, 0.04, 0.06],
            "magnetism": "materials_project"
        },
        "systems": [f"{element}_bcc", f"{element}_fcc", f"{element}_diamond"],
        "materials_project": {
            "api_key": "__NOT_NEEDED_FOR_EOS__",
            "n_structures": 1,
            "theoretical": False,
            "most_stable": True
        },
        "pseudopotentials": {
            "kinds": [f"{kind}"],
            "versions": [f"{version}"],
            "appendices": [f"{appendix}"]
        },
        "numerical_orbitals": {
            "types": ["DZP"],
            "rcuts": [7, 8, 9, 10],
            "appendices": [""]
        }
    }
    finp = f"{element}_{kind}_{version}_{appendix}_{ecut}.json"
    with open(finp, "w") as f:
        json.dump(result, f, indent=4)
    return os.path.abspath("".join([os.getcwd(), "/", finp]))

import zipfile
import re
import time
import apns.test.abacustest as amia
import apns.test.compress as amic
def main(fecut: str = "./apns_cache/apns_ecutwfc_db.json", 
         fpspot: str = "./download/pseudopotentials/pseudo_db.json", 
         included: list = None, excluded: list = None):
    """main driver for iteratively generating EOS test input files

    Args:
        fecut (str, optional): _description_. Defaults to "./apns_cache/apns_ecutwfc_db.json".
        fpspot (str, optional): _description_. Defaults to "./download/pseudopotentials/pseudo_db.json".
        included (list, optional): _description_. Defaults to None.
        excluded (list, optional): _description_. Defaults to None.
    """
    ecut_db = load_ecutwfc(fecut=fecut)
    pspot_db = load_pseudo(fpspot=fpspot)
    tranverse(ecut_db, pspot_db, included=included, excluded=excluded)

    timestamp = time.strftime('%Y%m%d%H%M%S')
    job_folder = f"apns_eos_testsuite_{timestamp}"
    result = job_folder.replace("testsuite", "result")
    fzips = [f for f in os.listdir() if f.endswith(".zip") and f.startswith("apns")]
    # decompress all zip files to the job folder
    for fzip in fzips:
        amic.unpack(fzip, job_folder)
        os.system(f"rm -f {fzip}")
    abacustest_param = amia.write_abacustest_param(jobgroup_name=job_folder, 
                                                   bohrium_login={"username": "", "password": "", "project_id": ""},
                                                   save_dir=result, 
                                                   rundft=[{"ifrun": True, 
                                                            "command": amia.ABACUS_COMMAND, 
                                                            "ncores": 32, "memory": 64,
                                                            "job_folders": [f"{job_folder}/*"]}])
    fparam = job_folder.replace("testsuite", "param") + ".json"
    flog = job_folder.replace("testsuite", "submit") + ".log"
    with open(fparam, "w") as f:
        json.dump(abacustest_param, f, indent=4)
    print(f"abacustest parameter file is written to {fparam}")
    # run abacus test
    os.system(f"nohup abacustest mlops-submit -p {fparam} > {flog} 2>&1 &")
    # `mlops-submit` will not download files automatically
    # `submit` will.
    print(f"abacustest job is submitted, please check the log file {flog}")

    return None
    # to merge all produced zip files together into one zip file
    zipfiles = [f for f in os.listdir() if f.endswith(".zip") and f.startswith("apns")]
    zipfiles.sort(key=lambda x: int(re.findall(r"\d+", x)[0]))
    fzip = f"APNS_EOS_testsuite_{time.strftime('%Y%m%d%H%M%S')}.zip"
    # move all files and folders in the zipfiles to the new zip file
    with zipfile.ZipFile(fzip, "w") as z:
        for f in zipfiles:
            # get all files and folders in the zip file
            with zipfile.ZipFile(f, "r") as zf:
                for info in zf.infolist():
                    z.writestr(info, zf.read(info))
            os.remove(f)
    print(f"Files are merged into {fzip}")

import unittest
class TestEOSTranverse(unittest.TestCase):

    def test_attributes(self):
        result = attributes("sg15_1.0")
        self.assertEqual(result, ["sg15", "1.0", None])
        result = attributes("sg15")
        self.assertEqual(result, ["sg15", None, None])
        result = attributes("sg15_1.0_fr")
        self.assertEqual(result, ["sg15", "1.0", "fr"])

    def test_pspot_eq(self):
        self.assertTrue(pspot_eq("sg1510", "sg15_1.0"))
        self.assertTrue(pspot_eq("sg1510fr", "sg15_1.0_fr"))
        self.assertFalse(pspot_eq("sg1510fr", "sg15_1.0"))
        self.assertTrue(pspot_eq("sg15_1.0_fr", "sg1510fr"))

if __name__ == "__main__":
    #unittest.main()
    set0 = ["H" , "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", 
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
            "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", 
            "Cs", "Ba", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
            "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", 
            "Fr", "Ra", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", 
            "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]            
    excluded = ["Co", "W", "Al", "Cu", "Hf", "Mg", "Mo", "Ta", "Ti", "V", "Zr"] # elements that are already tested
    
    main(included=["Ac", "Th", "Pa", "U", "Np", "Pu"])

    exit()
    # seperate every five element into one group
    istart = set0.index("Ba")
    n_to_test = 1
    iend = istart + n_to_test
    for i in range(istart, iend, n_to_test):
        included = set0[i:i+n_to_test]
        main(included=included, excluded=excluded)
