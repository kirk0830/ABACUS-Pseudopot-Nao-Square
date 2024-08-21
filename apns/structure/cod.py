import pymatgen.ext.cod as cod
"""
cod.COD().get_structure_by_id(9009009).to(filename='CmO2.cif')

# other structure search method is not simply available because MySQL is needed.
"""
def download(formula: list[str]):
    """
    Download structures from the COD API.

    Args:
        formula (list[str]): A list of chemical formulas.

    Returns:
        dict: A dictionary containing the chemical formulas as keys and a list of structures as values.
    """
    raise NotImplementedError("Download from COD is not implemented yet.")

import unittest
class TestCodDatabase(unittest.TestCase):
    def test_download(self):
        self.assertRaises(NotImplementedError, download, ["TiO2", "Er2O3"])

if __name__ == "__main__":
    unittest.main()