from pymatgen.io import cif

def export_primitive_cell(fcif: str):
    cifparser = cif.CifParser(fcif)
    structure = cifparser.parse_structures(primitive=True)[0]
    primitive = structure.to_primitive()
    fcif = fcif.replace(".cif", "_primitive.cif")
    primitive.to(filename=fcif, fmt="cif")
    return fcif

if __name__ == "__main__":
    fcif = "./apns/module_structure/POSCAR_Ac_XO3.cif"
    fcif = export_primitive_cell(fcif)
    print(f"Exported primitive cell to {fcif}")