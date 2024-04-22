"""
Usage:
    main.py -i <input_file>

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
def entry():
    """parse command line arguments"""
    parser = argparse.ArgumentParser(description="""
APNS is a workflow for testing the pseudopotential and pseudopotential-numerical atomic orbital
bundle. It also provides interface with ABACUS numerical atomic orbital generation code, called
SIAB, presently is also maintained by AISI-ABACUS developers.
""")
    parser.add_argument("-i", "--input", help="input file specifying the workflow", default="input.json")
    finp = parser.parse_args().input
    return finp

"""main"""
import apns.module_workflow.driver as amwd
def main():
    """main function"""
    finp = entry()
    driver = amwd.spawn_driver(finp)
    driver.setup()
    driver.run()

import unittest
class TestMain(unittest.TestCase):
    def test_entry(self):
        # test default
        self.assertEqual(entry(), "input.json")
        # mock command line arguments
        import sys
        sys.argv = ["main.py", "-i", "test.json"]
        self.assertEqual(entry(), "test.json")
    def test_main(self):
        print(f"Unittest on {__file__}, no test implemented for this workflow function main().")

if __name__ == "__main__":
    print(f"Unittest on {__file__} is skipped due to it is pure workflow function.")
    main()