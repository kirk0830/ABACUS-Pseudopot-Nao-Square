"""
import pymatgen.ext.optimade as optimade

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