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

def cal_delta_wrtacwf(token: str, bmfit: dict, vmin: float, vmax: float,
                      refdata_path: str = ACWF_DATAPATH + ACWF_REFJSON) -> float:
    """Search the best fit data from ACWF reference json file,
    because the structural data provided by Materials Project
    does not explicitly include the conventional name of crystal
    phase like BCC, FCC, Diamond, ..."""
    import json
    print(f"""Calculate delta value wrt. ACWF all-electron calculation results.
query token: {token}, 
vmin: {vmin}, 
vmax: {vmax}
""")
    with open(refdata_path, "r") as f:
        data = json.load(f)["BM_fit_data"][token]
    return delta_value(bm_fit1=bmfit, bm_fit2=data, vmin=vmin, vmax=vmax)

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
    import numpy as np
    v1, v2 = bm_fit1["min_volume"], bm_fit2["min_volume"]
    e1, e2 = bm_fit1["E0"], bm_fit2["E0"]
    b1, b2 = bm_fit1["bulk_modulus_ev_ang3"], bm_fit2["bulk_modulus_ev_ang3"]
    b1p, b2p = bm_fit1["bulk_deriv"], bm_fit2["bulk_deriv"]
    v = np.linspace(vmin, vmax, 100)
    e1 = birch_murnaghan(v, e1, b1, b1p, v1) - bm_fit1["E0"]
    e2 = birch_murnaghan(v, e2, b2, b2p, v2) - bm_fit2["E0"]
    delta_e_2 = simpson((e1 - e2)**2, x=v)/(vmax - vmin)
    return np.sqrt(delta_e_2)/natom

class EOSSingleCase:
    """Definition: A single case in Equation-Of-States (EOS) test is:
    1. one system
    2. many volumes
    3. (therefore) many energies
    4. one (combination of) pseudopotential(s)"""
    system: str
    volumes: list
    energies: list
    natom: int
    pps: list

    def __init__(self, system: str, pps: list, cases: list):
        """init from a nested dict, data can be get like:
        ```python
        scratch = [{"volume": 10, "energy": -1.0, "natom": 1},
                   {"volume": 20, "energy": -2.0, "natom": 1},
                   {"volume": 30, "energy": -3.0, "natom": 1}]
        ```"""
        self.system = system
        self.pps = pps
        self.volumes = [float(c["volume"]) for c in cases]
        self.natom = cases[0]["natom"]
        assert all([c["natom"] == self.natom for c in cases]), "The number of atoms should be consistent for all volume tests"
        self.energies = [c["energy"] for c in cases]
    
    def sort(self):
        """sort the data according to volume"""
        import numpy as np
        idx = np.argsort(self.volumes)
        self.volumes = [self.volumes[i] for i in idx]
        self.energies = [self.energies[i] for i in idx]

    def calc_eos(self, as_dict=False):
        """calculate the EOS parameters, return the result as a dict"""
        return fit_birch_murnaghan(self.volumes, self.energies, as_dict=as_dict)
    
    def pp(self, as_list: bool = False):
        """return the pseudopotential string"""
        from apns.analysis.drivers.apns2_utils import convert_fpp_to_ppid
        return self.pps if as_list else "|".join([convert_fpp_to_ppid(pp) for pp in self.pps])
    
    def tokenize(system: str):
        """conver the system to token used by ACWF"""
        import re
        u_in = r"([A-Z][a-z]*)_(sc|fcc|bcc|diamond)"
        o_in = r"([A-Z][a-z]*O)_(x\dy\d)"
        assert re.match(u_in, system) or re.match(o_in, system), "The system name should be in the form of 'Element_Phase' or 'ElementO_xdyd'"
        if re.match(u_in, system):
            # expected format: unaries = r"([A-Z][a-z]*)(-X/)([SC|FCC|BCC|Diamond])"
            elem, phase = re.match(u_in, system).groups()
            return f"{elem}-X/{phase.upper()}" if phase != "diamond" else f"{elem}-X/Diamond"
        if re.match(o_in, system):
            # expected format: oxides = r"([A-Z][a-z]*)/(X\dO\d)"
            elem, phase = re.match(o_in, system).groups()
            return f"{elem}/X{phase}"
        return None
        
    def __call__(self):
        self.sort()
        bmfit = self.calc_eos(as_dict=True)
        token = self.tokenize(self.system)
        delta = cal_delta_wrtacwf(token, bmfit, min(self.volumes), max(self.volumes))
        return self.pp(), bmfit, delta
        # legend, line, scalarized data

