"""APNS has three main functionalities: test, analysis and orbital generation. Each functionality has its own driver."""

class apns_driver:
    """APNS driver abstract class, unifying workflow drivers.
    For each kind of driver, must setup first, then run."""
    def __init__(self, finp: str):
        self.finp = finp
    def setup(self):
        """setup the driver"""
        pass
    def run(self):
        """run the driver"""
        pass

import apns.module_workflow.workflow_test.driver as amwtd
class test_driver(apns_driver):
    """test driver, for testing pseudopotentials and numerical orbitals"""
    def setup(self):
        """setup the driver"""
        pass
    def run(self):
        amwtd.driver_v1(self.finp)

import apns.module_workflow.workflow_analysis.driver as amwad
class analysis_driver(apns_driver):
    """analysis driver, for analyzing test results"""
    def setup(self):
        """setup the driver"""
        pass
    def run(self):
        pass

import apns.module_workflow.workflow_orbgen.driver as amwod
class orbgen_driver(apns_driver):
    """orbgen driver, for generating numerical orbitals"""
    def setup(self):
        """setup the driver"""
        pass
    def run(self):
        amwod.run(self.finp)

import json
def spawn_driver(finp: str) -> apns_driver:
    """return corresponding driver according to detailed user settings"""
    with open(finp, "r") as f:
        inp = json.load(f)
    if inp["global"]["test_mode"] == "pseudopotential" or inp["global"]["test_mode"] == "numerical_orbital":
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