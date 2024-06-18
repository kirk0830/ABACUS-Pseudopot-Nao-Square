ACWF_DATAPATH = "/root/documents/acwf-verification-scripts-main/acwf_paper_plots/code-data/"
ACWF_REFJSON = "results-unaries-verification-PBE-v1-AE-average.json"

def read_acwf_refdata(token: str, 
                      refdata_path: str = ACWF_DATAPATH + ACWF_REFJSON,
                      domain: str = "BM_fit_data"):
    """Read data from ACWF reference json file
    see Github repo:
    https://github.com/aiidateam/acwf-verification-scripts/tree/main
    
    and their paper:
    Bosoni, E., Beal, L., Bercx, M. et al. 
    How to verify the precision of density-functional-theory implementations via reproducible 
    and universal workflows. Nat Rev Phys 6, 45-58 (2024). 
    https://doi.org/10.1038/s42254-023-00655-3

    This project is referred in Standard Solid State Pseudopotential (SSSP) Library 2.0, which
    is maintained on website Materials Cloud:
    https://www.materialscloud.org/discover/sssp
    """
    import json
    with open(refdata_path, "r") as f:
        data = json.load(f)
    return data[domain][token]

def cal_delta_wrtacwf(element: str, bmfit: dict, vmin: float, vmax: float, bravis: str = "unknown",
                      refdata_path: str = ACWF_DATAPATH + ACWF_REFJSON) -> float:
    """Search the best fit data from ACWF reference json file,
    because the structural data provided by Materials Project
    does not explicitly include the conventional name of crystal
    phase like BCC, FCC, Diamond, ..."""
    import json, re
    import apns.analysis.postprocess.eos as amape
    print(f"""Calculate delta value wrt. ACWF all-electron calculation results.
element: {element}, 
bravis: {bravis}, 
vmin: {vmin}, 
vmax: {vmax}
""")
    delta = 1e10
    syspatn = r"([A-Z][a-z]*)(-X/)([(SC)(FCC)(BCC)(Diamond)])"
    with open(refdata_path, "r") as f:
        data = json.load(f)
    data = data["BM_fit_data"]

    phase = ""
    bravis = bravis.upper() if bravis.lower() != "diamond" else "Diamond"
    if bravis != "UNKNOWN" and bravis in ["SC", "FCC", "BCC", "Diamond"]:
        token = f"{element}-X/{bravis}"
        print("Get data from token", token)
        value = amape.delta_value(bm_fit1=bmfit, bm_fit2=data[token], vmin=vmin, vmax=vmax)
        return value, token
    # else...
    for key in data.keys():
        _match = re.match(syspatn, key)
        if _match is not None:
            if _match.group(1) == element:
                print("Scan crystal phase from ACWF reference data:", key)
                result = data[key]
                value = amape.delta_value(bm_fit1=bmfit, bm_fit2=result, vmin=vmin, vmax=vmax)
                if value < delta:
                    delta = value
                    phase = key
    return delta, phase

def birch_murnaghan(v, e0, b0, b0p, v0):
    return e0 + 9 * v0 * b0 / 16 * (b0p * ((v0 / v)**(2/3) - 1)**3 + ((v0 / v)**(2/3) - 1)**2 * (6 - 4 * (v0 / v)**(2/3)))

def fit_birch_murnaghan(volumes, energies, as_dict=False):
    """Fit the Birch-Murnaghan equation of state to the given volumes and energies.
    
    returns:
        v0: float, the minimum volume,
        e0: float, the minimum energy,
        b0: float, the bulk modulus,
        b0p: float, the pressure derivative of the bulk modulus.
    For aligning with the reference, their correspondances are:
        v0: min_volume,
        e0: E0,
        b0: bulk_modulus_ev_ang3,
        b0p: bulk_deriv."""
    import scipy.optimize as opt
    try:
        popt, _ = opt.curve_fit(birch_murnaghan, volumes, energies, p0=(energies[0], 1, 1, volumes[0]))
    except RuntimeError:
        data = "\n".join([f"{volumes[i]:.4f} {energies[i]:.4f}" for i in range(len(volumes))])
        print(f"Failed to fit the Birch-Murnaghan EOS. Source data:\n{data}")
        return None, None, None, None
    if as_dict:
        return {
            "min_volume": popt[3],
            "E0": popt[0],
            "bulk_modulus_ev_ang3": popt[1],
            "bulk_deriv": popt[2]
        }
    return popt[3], popt[0], popt[1], popt[2]

def delta_value(bm_fit1: dict, 
                bm_fit2: dict, 
                vmin: float, 
                vmax: float,
                natom: int = 1) -> float:
    """Calculate the delta value for two EOS fit results
    
    delta = sqrt(1/(vmax - vmin) * integral((E1 - E2)^2))
    
    Args:
        bm_fit1: dict, the first EOS fit result,
        bm_fit2: dict, the second EOS fit result,
        vmin: float, the minimum volume,
        vmax: float, the maximum volume,
        natom: int, the number of atoms in the cell.
    Returns:
        float, the delta value
    """
    from scipy.integrate import simpson
    v1, v2 = bm_fit1["min_volume"], bm_fit2["min_volume"]
    e1, e2 = bm_fit1["E0"], bm_fit2["E0"]
    b1, b2 = bm_fit1["bulk_modulus_ev_ang3"], bm_fit2["bulk_modulus_ev_ang3"]
    b1p, b2p = bm_fit1["bulk_deriv"], bm_fit2["bulk_deriv"]
    v = np.linspace(vmin, vmax, 100)
    e1 = birch_murnaghan(v, e1, b1, b1p, v1) - bm_fit1["E0"]
    e2 = birch_murnaghan(v, e2, b2, b2p, v2) - bm_fit2["E0"]
    delta_e_2 = simpson((e1 - e2)**2, x=v)/(vmax - vmin)
    return np.sqrt(delta_e_2)/natom

import unittest
import numpy as np
class TestEos(unittest.TestCase):
    def test_birch_murnaghan_eos(self):
        Ac_X2O = [
            [
                75.57259067144167,
                -3125.6610180711
            ],
            [
                77.18051813253214,
                -3125.6904674548
            ],
            [
                78.78844559362994,
                -3125.7069824516
            ],
            [
                80.39637305472405,
                -3125.7119481209
            ],
            [
                82.00430051581594,
                -3125.7065927846
            ],
            [
                83.61222797691082,
                -3125.6920087076
            ],
            [
                85.22015543800983,
                -3125.6691741391
            ]
        ]
        volumes = [x[0] for x in Ac_X2O]
        energies = [x[1] for x in Ac_X2O]
        reference = {
            "E0": -3125.7119562222942,
            "bulk_deriv": 4.608493219328219,
            "bulk_modulus_ev_ang3": 0.32175656481943105,
            "min_volume": 80.33593470120812,
        }
        result = fit_birch_murnaghan(volumes, energies, as_dict=True)
        for key in reference:
            self.assertAlmostEqual(reference[key], result[key], delta=1e-4)

        data = [[37.9098, -214.406018],
                [38.7163, -214.441355],
                [39.5229, -214.4649915],
                [40.3295, -214.4780729],
                [41.1361, -214.4816059],
                [41.9427, -214.4765565],
                [42.7493, -214.4637563]]
        energies_abacus = [x[1] for x in data]
        volumes_abacus = [x[0] for x in data]
        reference = {
            "E0": 0,
            "bulk_deriv": 4.31178461988603,
            "bulk_modulus_ev_ang3": 0.5524442002451444,
            "min_volume": 40.914946909495136,
            "residuals": 0
        }
        result_abacus = fit_birch_murnaghan(volumes_abacus, energies_abacus, as_dict=True)
        result_abacus = {
            "E0": -214.4820510584568,
            "bulk_deriv": 4.215165534680927,
            "bulk_modulus_ev_ang3": 0.5475682644107036,
            "min_volume": 41.05145826931731,
            "residuals": 0
        }

if __name__ == "__main__":
    unittest.main()