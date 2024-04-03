"""APNS has three main functionalities: test, analysis and orbital generation. Each functionality has its own driver."""
import os
import apns.module_workflow.identifier as amwi
import apns.module_pseudo.manage as ampm
import apns.module_nao.manage as amnm
welcome = "\n"
welcome += "="*100
welcome += "\n"
welcome += """Welcome to ABACUS Pseudopot-Nao Square (APNS) workflow

This workflow collection is for performing tests on pseudopotential (pseudopot) and
pseudopot-nao (numerical atomic orbital) bundle. It also provides interface with ABACUS
numerical atomic orbital generation code, called SIAB, presently is also maintained by
AISI-ABACUS developers.

If it is the first run of APNS, please make sure your environment configured correctly:
1. store pseudopotentials in directory and specified in input file the `pseudo_dir`
, similarly for numerical atomic ortbials, set `orbital_dir`. The APNS workflow will
check and update the file `pseudo_db.json` everytime when starts. 

2. Also in `pseudo_dir`, please have a file `rules.json` in it. In this file, user should
explicitly define three domains to identify one kind of pseudopotential. Apart from
three domains needed to define for each pseudopotential, another two keys are needed 
to defined: re.folder and re.file. With these two keys, one can always get pseudopotentials
identified from mixture of pseudopotentials.

For more information about setting of input_*.json for different tasks, please refer to
online documentation.
"""
welcome += "="*100

class apns_driver:
    """APNS driver abstract class, unifying workflow drivers.
    For each kind of driver, must setup first, then run."""
    def __init__(self, finp: str):
        self.finp = finp
    def setup(self):
        """setup the driver"""
        print(welcome)
        # 1. initialize cache directory
        print("Current working directory: {}".format(os.getcwd()))
        """change id.TEMPORARY_FOLDER to absolute path"""
        amwi.TEMPORARY_FOLDER = os.path.join(os.getcwd(), amwi.TEMPORARY_FOLDER)
        """create cache directory if not exist"""
        if not os.path.exists(amwi.TEMPORARY_FOLDER):
            os.mkdir(amwi.TEMPORARY_FOLDER)
        else:
            print("Cache directory already exists.")
        # 2. refresh/reload pseudopotential and numerical orbital archive
        # because it is Pseudopot-Nao Square!
        print("Refreshing pseudopotential and numerical orbital archive...")
        with open(self.finp, "r") as f:
            inp = json.load(f)
        ampm.load(inp["global"]["pseudo_dir"])
        amnm.load(inp["global"]["orbital_dir"])
        
    def run(self):
        """run the driver"""
        pass

import apns.module_workflow.workflow_test.driver as amwtd
class test_driver(apns_driver):
    """test driver, for testing pseudopotentials and numerical orbitals"""
    def run(self):
        amwtd.run(self.finp)

import apns.module_workflow.workflow_analysis.driver as amwad
class analysis_driver(apns_driver):
    """analysis driver, for analyzing test results"""
    def run(self):
        amwad.run(self.finp)

import apns.module_workflow.workflow_orbgen.driver as amwod
class orbgen_driver(apns_driver):
    """orbgen driver, for generating numerical orbitals"""
    def run(self):
        amwod.run(self.finp)

import json
def spawn_driver(finp: str) -> apns_driver:
    """return corresponding driver according to detailed user settings"""
    with open(finp, "r") as f:
        inp = json.load(f)
    if inp["global"]["test_mode"] in ["pseudopotential", "numerical_orbital"]:
        print("Activate test mode: ", inp["global"]["test_mode"])
        return test_driver(finp)
    elif inp["global"]["test_mode"] == "analysis":
        print("Analysis mode activated.")
        return analysis_driver(finp)
    elif inp["global"]["test_mode"] == "orbgen":
        print("Orbgen mode activated.")
        return orbgen_driver(finp)
    else:
        raise ValueError("Invalid test mode.")