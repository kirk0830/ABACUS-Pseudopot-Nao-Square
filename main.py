"""ABACUS Pseudopot-Nao Squared (APNS)

APNS is a workflow for testing the pseudopotential and pseudopotential-numerical atomic orbital
bundle. It also provides interface with ABACUS numerical atomic orbital generation code, called
SIAB, presently is also maintained by AISI-ABACUS developers.

Usage:
    main.py -v <version> -i <input_file>

Based on different task specified by keyword "global/test_mode", different workflow will be
executed. Currently available workflows are:
1. "global/test_mode" = "pseudopotential", which is the workflow for testing pseudopotential
    with ABACUS code, setting `basis_type = pw`, therefore only pseudopotential is variable.
    Several kinds of properties can be tested from a single SCF run, such as electronic 
    energies, forces and stress on present geometry, also Density of States (DOS).
2. "global/test_mode" = "nao", which is the workflow for testing pseudopotential-numerical
    atomic orbital bundle with ABACUS code, setting `basis_type = lcao`. In this case be sure
    to make pseudopotential and nao really consistent (nao strongly dependes on pseudopotenial).
    Similar kinds of properties can be tested from a single SCF run, such as electronic
    energies, forces and stress on present geometry, also Density of States (DOS).
3. "global/test_mode" = "orbgen", which is the workflow for testing the orbital generation
    with SIAB code. Two different execution modes are available, "ex-situ" and "in-situ".
    "ex-situ" mode is to generate the orbital generation configure file while "in-situ" mode
    pass configuration to SIAB code and execute the orbital generation directly.
4. "global/test_mode" = "analysis", which is the workflow for analyzing the results from
    either pseudopotential or nao test. It will generate a report for the test results, along
    with html test reports for posting on APNS Github Pages.
"""
import argparse
def initialize():
    """initialize the workflow by getting input file and version of the workflow
    Returns:
        input_file (str): input file specifying the workflow
        version (str): version of the workflow
    """
    parser = argparse.ArgumentParser(description="APNS")
    parser.add_argument("-v", "--version", help="Version of the workflow", default="v1")
    parser.add_argument("-i", "--input", help="input file specifying the workflow", default="input.json")

    input_file = parser.parse_args().input
    version = parser.parse_args().version
    return input_file, version

"""main"""
import apns.module_workflow.driver as amwd
import apns.module_workflow.workflow_test.driver as amwtd
def main():
    # initialize: get input file and version
    input_file, version = initialize()
    if version == "v1":
        """I use polymorphism here, because seems each task can really have similar interface
        But it might be heavy..."""
        driver = amwd.spawn_driver(input_file)
        driver.setup() # no matter which driver, run setup
        driver.run()   # no matter which driver, run the workflow
    else:
        """version 0 would be fully deprecated in the future"""
        amwtd.driver_v0(input_file)

if __name__ == "__main__":
    main()