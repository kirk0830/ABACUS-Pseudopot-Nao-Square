"""
import pymatgen.ext.optimade as optimade

structures = optimade.OptimadeRester().get_structures(elements=["Bk", "I", "O"])

for key, value in structures.items():
    print(key, value)
    for _k, _v in value.items():
        print(_k, _v)
"""