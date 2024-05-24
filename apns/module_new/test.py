from pymatgen.io.cif import CifParser
parser = CifParser("apns_cache/Li3DyCl6-P321.cif")
print(parser.parse_structures(primitive=False)[0])