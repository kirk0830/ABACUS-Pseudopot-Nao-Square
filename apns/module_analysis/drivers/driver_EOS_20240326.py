import os
import apns.module_analysis.postprocess.read_abacus_out as amapr
def is_outdir(files: list):
    return "running_relax.log" in files and\
           "STRU_ION_D" in files and\
           "INPUT" in files and\
           "STRU_SIMPLE.cif" in files and\
           "STRU_READIN_ADJUST.cif" in files

def search(path: str, fromcal: str = "relax"):
    
    result = {}
    for root, dirs, files in os.walk(path):
        if is_outdir(files):
            eks = amapr.read_etraj_fromlog(f"{root}/running_{fromcal}.log")[-1]
            v = amapr.read_volume_fromstru(f"{root}/STRU_ION_D", "A")
            _, sys, mpid, pnid, _ = amapr.read_testconfig_fromBohriumpath(root)
            result.setdefault((sys, mpid, pnid), []).append((v, eks))
    
    # sort by volume
    for key in result.keys():
        result[key] = sorted(result[key], key=lambda x: x[0])
        result[key] = list(zip(*result[key]))
        result[key][0] = list(result[key][0])
        result[key][1] = list(result[key][1])

    return result

import apns.module_analysis.postprocess.eos as amape
def calculate(sys_mpid_pnid_veks: dict):
    for key in sys_mpid_pnid_veks.keys():
        ve = sys_mpid_pnid_veks[key][0]
        eks = sys_mpid_pnid_veks[key][1]
        sys_mpid_pnid_veks[key] = amape.birch_murnaghan_eos(ve, eks, as_dict=True)
    return sys_mpid_pnid_veks

def run(path: str, fromcal: str = "relax"):
    ve_data = search(path, fromcal=fromcal)
    return calculate(ve_data)

import unittest
class TestEOS(unittest.TestCase):
    def test_is_outdir(self):
        files = ["running_relax.log", "STRU_ION_D", "INPUT", "STRU_SIMPLE.cif", "STRU_READIN_ADJUST.cif"]
        self.assertTrue(is_outdir(files))
        files = ["running_relax.log", "STRU_ION_D", "STRU_SIMPLE.cif", "STRU_READIN_ADJUST.cif"]
        self.assertFalse(is_outdir(files))
        files = ["running_relax.log", "STRU_ION_D", "STRU_SIMPLE.cif", "STRU_READIN_ADJUST.cif", "STRU_READIN_ADJUST.cif"]
        self.assertFalse(is_outdir(files))
        files = ["running_relax.log", "STRU_ION_D", "STRU_SIMPLE.cif", "STRU_READIN_ADJUST.cif", "STRU_READIN_ADJUST.cif", "STRU_READIN_ADJUST.cif"]
        self.assertFalse(is_outdir(files))
    def test_calculate(self):
        data = [[37.9098, -214.406018],
                [38.7163, -214.441355],
                [39.5229, -214.4649915],
                [40.3295, -214.4780729],
                [41.1361, -214.4816059],
                [41.9427, -214.4765565],
                [42.7493, -214.4637563]]
        smpv = {
            ("Si", "149", "dojo05"): ([x[0] for x in data], [x[1] for x in data])
        }
        result = calculate(smpv)
        self.assertAlmostEqual(
            result[("Si", "149", "dojo05")]["E0"], 
            -214.48165631975377, 
            delta=1e-4)
        self.assertAlmostEqual(
            result[("Si", "149", "dojo05")]["bulk_deriv"], 
            4.220117478624572, 
            delta=1e-4)
        self.assertAlmostEqual(
            result[("Si", "149", "dojo05")]["bulk_modulus_ev_ang3"], 
            0.5476332146301505, 
            delta=1e-4)
        # it is not values really fully reliable, just for testing
        # more reliable all-electron data can refer to https://doi.org/10.1038/s42254-023-00655-3
        """
        "Si-X/Diamond": {
            "E0": 0,
            "bulk_deriv": 4.31178461988603,
            "bulk_modulus_ev_ang3": 0.5524442002451444,
            "min_volume": 40.914946909495136,
            "residuals": 0
        },
        """
    
if __name__ == "__main__":
    # path = "../job-abacustest-v0.3.102-4edc4d"
    # result = run(path)
    # print(result)
    unittest.main()