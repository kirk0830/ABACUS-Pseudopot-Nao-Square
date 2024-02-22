"""orbital generation workflow driver

There are several prerequisites for this workflow:
1. pseudopotentials have been properly configured. If not, run 
   `apns/module_pseudo/archive.py` to configure pseudopotentials.
2. pseudopotentials pw-convergence test has been done. If not, run
   convergence test first with input file `input_test.json` first
   by `apns/main.py -i input_test.json` first, then run analysis
   workflow to generate/grep converged values for pseudopotentials
"""
import apns.module_workflow.workflow_orbgen.initialize as amwwoi
def initialize(finp: str):
   """setup the orbgen workflow"""
   return amwwoi.initialize(finp)

import json
import os
import apns.module_nao.orbgen as amno
import time
def run(finp: str):
   """run orbital generation from this driver"""
   # initialize the workflow
   inp, pspot_ids = initialize(finp)
   ecutwfc = None if "ecutwfc" not in inp["abacus"] else inp["abacus"]["ecutwfc"]
   for element in inp["systems"]:
      # this is just the reference mode. For general, ...
      siab_generator = amno.siab_generator(element=element,
                                           rcuts=inp["numerical_orbitals"]["rcuts"],
                                           ecutwfc=ecutwfc,
                                           pspot_id=pspot_ids,
                                           other_settings={
                                              "environment": inp["environment"],
                                              "mpi_command": inp["mpi_command"],
                                              "abacus_command": inp["abacus_command"],
                                           })
      for siab_input in siab_generator:
         fsiab = f"SIAB_INPUT_{time.strftime('%Y%m%d%H%M%S')}.json"
         with open(fsiab, "w") as f:
            json.dump(siab_input, f, indent=4)
         siab = inp["orbgen"]["generator"]
         #os.system("python " + siab + " -i " + fsiab)

if __name__ == "__main__":
   run("input_orbgen.json")