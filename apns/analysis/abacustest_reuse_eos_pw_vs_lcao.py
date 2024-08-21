"""with abacustest Reuse workflow, eos_pw_vs_lcao module,
can have following properties:
1. energies
2. volumes
3. band structures

with these quantities, can the following properties be calculated:
1. equation-of-state (EOS)
    by non-linearly fitting the Birch-Murnaghan equation of state,
    can get the equilibrium volume, bulk modulus, and its pressure derivative.
    With the fit data, calculate delta value
2. basis set completeness
    by directly comparing the energy minimum of the PW-LCAO pair,
    the energy difference is a good indicator of the basis set completeness.
    The lower the energy difference, the better the basis set completeness.
3. band structure similarity
    similarity is measured by the band structure difference.

For one shot of run, abacustest will produce following files:
abacustest file structure:
```bash
- task-main
- inputs
- outputs\
    - results\
        - test_system_folder_1
            - PW
                - eos-3 <- one single job folder, can be used for band structure similarity test
                - eos-2 <- EOS point "another"
                - ...
                - eos0
                - ...
                - eos3
            - LCAO1
                - eos-3
                - ...
            - ...
            - LCAOn
        - test_system_folder_2
        - ...
        - supermetrics.json
        - post.py
        - metrics.json <- it is where the results are stored (vols, eners, etc.)
        - test_system_result_pic_1
        - test_system_result_pic_2
        - ...
```
In metrics.json, the content is like:
```json
{
    "test_system_folder_1/PW/eos0": {...},
    "test_system_folder_1/PW/eos1": {...},
    //...
}
```
"""

def _sort_vols_energies(vols: list, eners: list):
    """sort the volumes and energies in ascending order of volumes.
    
    Args:
        vols (list): list of volumes
        eners (list): list of energies
    
    Returns:
        tuple: (sorted_vols, sorted_eners)
    """
    import numpy as np
    vols = np.array(vols)
    eners = np.array(eners)
    idx = np.argsort(vols)
    return vols[idx], eners[idx]

def _nested_matrices(mat: dict):
    """convert the matrices.json contents to two nested dict, the first stores
    result from PW and the second stores result from LCAO.
    
    Args:
        mat (dict): the matrices.json content
    
    Returns:
        pw (dict): the nested dict for PW results
        lcao (dict): the nested dict for LCAO results
    """
    pw, lcao = {}, {}
    for k, v in mat.items():
        assert k.count('/') == 2, f"key {k} is not in the correct format"
        system, basis, folder = k.split('/')
        # because the folder is hard coded to be name in range eos-3 to eos3,
        # one can always get the index of the folder by:
        index = int(folder[3:]) + 3 # move to start from 0
        if basis.upper() == "PW":
            pw.setdefault(system, [None] * 7)[index] = v
        else:
            iorb = int(basis[4:]) - 1 # LCAO1, LCAO2, ...
            lcao.setdefault(system, []).append([None] * 7)[iorb] = v
    return pw, lcao

def read_abacustest_metrices(fmat: str, nested: bool = True):
    """read the metrics.json file in abacustest output folder

    Args:
        fmat (str): the path to the metrics.json file
        nested (bool, optional): whether to return the nested dict. Defaults to True.

    Returns:
        dict: the metrics dict, or two nested dicts
    """
    import json
    with open(fmat, 'r') as f:
        metrices = json.load(f)
    return metrices if not nested else _nested_matrices(metrices)

def _cal_delta(dataset1: list, dataset2: list):
    """calculate the delta value between two datasets.
    
    Args:
        dataset1 (list): the first dataset, each element is a dict with keys: volume, energy_per_atom
        dataset2 (list): similar with dataset1
        
    Returns:
        float: the delta value (PER ATOM)
    """
    from apns.analysis.apns2_eos_utils import fit_birch_murnaghan, delta_value
    ener1 = [d["energy_per_atom"] for d in dataset1]
    vol1 = [d["volume"] for d in dataset1]
    ener2 = [d["energy_per_atom"] for d in dataset2]
    vol2 = [d["volume"] for d in dataset2]
    # sort the data
    vol1, ener1 = _sort_vols_energies(vol1, ener1)
    vol2, ener2 = _sort_vols_energies(vol2, ener2)
    # fit the data
    bm1 = fit_birch_murnaghan(vol1, ener1)
    bm2 = fit_birch_murnaghan(vol2, ener2)
    # calculate delta
    vmin = min(vol1[0], vol2[0])
    vmax = max(vol1[-1], vol2[-1])
    return delta_value(bm1, bm2, vmin, vmax)

def _cal_basis_complete(dataset1: list, dataset2: list):
    """calculate the basis completeness indicator (energy difference)
    
    Args:
        dataset1 (list): the first dataset, each element is a dict with keys: volume, energy_per_atom
        dataset2 (list): similar with dataset1
    
    Returns:
        float: the energy difference (PER ATOM)
    """
    ener1 = [d["energy_per_atom"] for d in dataset1]
    ener2 = [d["energy_per_atom"] for d in dataset2]
    return abs(min(ener1) - min(ener2))

def cal_delta_pw_vs_lcao(nested_pw: dict, nested_lcao: dict):
    """calculate the delta value between PW and LCAO results."""
    result = {}
    assert nested_pw.keys() == nested_lcao.keys(), "the systems in PW and LCAO are not the same"
    for system in nested_pw: # can loop over either keys of nested_pw or nested_lcao
        pw = nested_pw[system]
        for lcao in nested_lcao[system]:
            delta = _cal_delta(pw, lcao)
            result.setdefault(system, []).append(delta)
    return result

def _mat_key(system: str, basis: str, iorb = None, itest = None):
    """make up one key for the matrices.json"""
    iorb = "" if iorb is None else iorb
    itest = 0 if itest is None else itest
    assert basis.upper() in ["PW", "LCAO"], f"basis {basis} is not supported"
    assert iorb is None and basis.upper() == "PW" or iorb is not None and basis.upper() == "LCAO", \
        f"basis {basis} and iorb {iorb} are not matched"
    return f"{system}/{basis}{iorb}/eos{itest}"

def _get_test_feature(job_dir: str, mat_key: str):
    import os
    from apns.analysis.postprocess.read_abacus_out import read_stru
    parsed = read_stru(os.path.join(job_dir, "/".join(mat_key, "STRU")))

import unittest
class AbacustestReuseEOSPWvsLCAOTest(unittest.TestCase):
    def test_ignore(self):
        print("File abacustest_reuse_eos_pw_vs_lcao.py is not a unittest file.")