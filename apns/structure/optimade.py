"""Should you wish to cite the OPTIMADE specification, please use the following:

Evans et al, Developments and applications of the OPTIMADE API for materials discovery, design, and data exchange, Digital Discovery (2024) 10.1039/D4DD00039K (preprint: 10.48550/arXiv.2402.0057)
Andersen et al, OPTIMADE, an API for exchanging materials data, Sci. Data 8, 217 (2021) 10.1038/s41597-021-00974-z
Andersen et al, The OPTIMADE Specification, Zenodo, 10.5281/zenodo.4195050
If you use the optimade-python-tools to access or serve OPTIMADE APIs, please consider citing the following:

Evans et al, optimade-python-tools: A Python library for serving and consuming materials data via OPTIMADE APIs, Journal of Open Source Software, 6 (65), 3458 (2021), 10.21105/joss.03458"""

import pymatgen.ext.optimade as optimade
"""
structures = optimade.OptimadeRester().get_structures(elements=["Bk", "I", "O"])

for key, value in structures.items():
    print(key, value)
    for _k, _v in value.items():
        print(_k, _v)
"""
def download(formula: list[str]):
    """
    Download structures from the OPTIMADE API.

    Args:
        formula (list[str]): A list of chemical formulas.

    Returns:
        dict: A dictionary containing the chemical formulas as keys and a list of structures as values.
    """
    raise NotImplementedError("Download from Optimade is not implemented yet.")

import unittest
class TestOptimadeDatabase(unittest.TestCase):
    def test_download(self):
        self.assertRaises(NotImplementedError, download, ["TiO2", "Er2O3"])

if __name__ == "__main__":
    unittest.main()