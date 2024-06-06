"""orbital generation workflow driver

There are several prerequisites for this workflow:
1. pseudopotentials have been properly configured. If not, run 
   `apns/module_pseudo/archive.py` to configure pseudopotentials.
2. pseudopotentials pw-convergence test has been done. If not, run
   convergence test first with input file `input_test.json` first
   by `apns/main.py -i input_test.json` first, then run analysis
   workflow to generate/grep converged values for pseudopotentials
"""
import apns.outdated.orbgen.initialize as amwwoi
def initialize(finp: str):
   """setup the orbgen workflow"""
   return amwwoi.initialize(finp)

import json
import os
import apns.outdated.orbgen.orbgen as amno
import time
def run(finp: str):
   """run orbital generation from this driver"""
   # initialize the workflow, read-in, check and load available pseudopotentials
   runtime_settings, abacus_settings, elements, vupfs = initialize(finp)
   ecutwfc = None if "ecutwfc" not in abacus_settings else abacus_settings["ecutwfc"]
   for element in elements:
      # this is just the reference mode. For general, ...
      siab_settings = {"element": element, "rcuts": runtime_settings["rcuts"], "ecutwfc": ecutwfc,
                       "fpseudos": vupfs[element], "other_settings": {
                       "environment": runtime_settings["environment"],
                       "mpi_command": runtime_settings["mpi_command"],
                       "abacus_command": runtime_settings["abacus_command"]}}
      siab_generator = amno.siab_generator(**siab_settings)
      finps_siab = []
      for siab_input in siab_generator:
         finp_siab = f"SIAB_INPUT_{time.strftime('%Y%m%d%H%M%S')}.json"
         with open(finp_siab, "w") as f:
            json.dump(siab_input, f, indent=4)
         if runtime_settings["generate_mode"] == "in-situ":
            insitu_siab(runtime_settings["generator"], finp_siab)
         else:
            finps_siab.append(finp_siab)
            time.sleep(1) # this is for updating time stamp
      if runtime_settings["generate_mode"] == "ex-situ":
         fscript = exsitu_siab(runtime_settings["generator"], finps_siab)
         print(f"Ex-situ SIAB generation script generated: {fscript}")
   print("Orbital generation workflow finished.")

def insitu_siab(siab: str, finp_siab: str):
   """this is driver for in-situ run siab generating numerical orbitals"""
   print(f"Start in-situ mode ABACUS numerical atomic orbital generation at {time.strftime('%Y-%m-%d %H:%M:%S')}")
   os.system("python " + siab + " -i " + finp_siab)

def exsitu_siab(siab: str, finps_siab: list):
   """this is driver for ex-situ run siab generating numerical orbitals
   will generate a Python script to iteratively run each siab input file"""
   script = f"""# ABACUS-Pseudopot-Nao-Square automated SIAB generation script
# for orbital generation ex-situ run
# generated by apns/module_workflow/workflow_orbgen/driver.py:exsitu_siab()
# Author: ABACUS-AISI developer team
import time
import os
def run(siab: str, fsiab: str):
   os.system("python " + siab + " -i " + fsiab)
   time.sleep(1)

if __name__ == "__main__":
   finps_siab = {finps_siab}
   for fsiab in finps_siab:
      run("{siab}", fsiab)
"""
   fscript = "run_siab_" + time.strftime("%Y%m%d%H%M%S") + ".py"
   with open(fscript, "w") as f:
      f.write(script)

   return fscript

import unittest
class TestDriver(unittest.TestCase):

    def test_exsitu_siab(self):
        
        fscript = exsitu_siab("./SIAB/SIAB_nouvelle.py", [
            "SIAB_INPUT_20210101000000.json", "SIAB_INPUT_20210101000001.json"])
        self.assertTrue(os.path.exists(fscript))
        #os.remove(fscript)

if __name__ == "__main__":
    unittest.main()