"""use function to get value from dict, to avoid unexpected change on dict
DO NOT GET VALUE FROM DICT DIRECTLY"""
PERIODIC_TABLE_TOFULLNAME = {
    'H': 'Hydrogen', 'He': 'Helium', 'Li': 'Lithium', 'Be': 'Beryllium', 'B': 'Boron', 'C': 'Carbon', 'N': 'Nitrogen',
    'O': 'Oxygen', 'F': 'Fluorine', 'Ne': 'Neon', 'Na': 'Sodium', 'Mg': 'Magnesium', 'Al': 'Aluminium',
    'Si': 'Silicon', 'P': 'Phosphorus', 'S': 'Sulfur', 'Cl': 'Chlorine', 'Ar': 'Argon', 'K': 'Potassium',
    'Ca': 'Calcium', 'Sc': 'Scandium', 'Ti': 'Titanium', 'V': 'Vanadium', 'Cr': 'Chromium', 'Mn': 'Manganese',
    'Fe': 'Iron', 'Co': 'Cobalt', 'Ni': 'Nickel', 'Cu': 'Copper', 'Zn': 'Zinc', 'Ga': 'Gallium',
    'Ge': 'Germanium', 'As': 'Arsenic', 'Se': 'Selenium', 'Br': 'Bromine', 'Kr': 'Krypton', 'Rb': 'Rubidium',
    'Sr': 'Strontium', 'Y': 'Yttrium', 'Zr': 'Zirconium', 'Nb': 'Niobium', 'Mo': 'Molybdenum', 'Tc': 'Technetium',
    'Ru': 'Ruthenium', 'Rh': 'Rhodium', 'Pd': 'Palladium', 'Ag': 'Silver', 'Cd': 'Cadmium', 'In': 'Indium',
    'Sn': 'Tin', 'Sb': 'Antimony', 'Te': 'Tellurium', 'I': 'Iodine', 'Xe': 'Xenon', 'Cs': 'Caesium',
    'Ba': 'Barium', 'La': 'Lanthanum', 'Ce': 'Cerium', 'Pr': 'Praseodymium', 'Nd': 'Neodymium', 'Pm': 'Promethium',
    'Sm': 'Samarium', 'Eu': 'Europium', 'Gd': 'Gadolinium', 'Tb': 'Terbium', 'Dy': 'Dysprosium', 'Ho': 'Holmium',
    'Er': 'Erbium', 'Tm': 'Thulium', 'Yb': 'Ytterbium', 'Lu': 'Lutetium', 'Hf': 'Hafnium', 'Ta': 'Tantalum',
    'W': 'Tungsten', 'Re': 'Rhenium', 'Os': 'Osmium', 'Ir': 'Iridium', 'Pt': 'Platinum', 'Au': 'Gold',
    'Hg': 'Mercury', 'Tl': 'Thallium', 'Pb': 'Lead', 'Bi': 'Bismuth', 'Po': 'Polonium', 'At': 'Astatine',
    'Rn': 'Radon', 'Fr': 'Francium', 'Ra': 'Radium', 'Ac': 'Actinium', 'Th': 'Thorium', 'Pa': 'Protactinium',
    'U': 'Uranium', 'Np': 'Neptunium', 'Pu': 'Plutonium', 'Am': 'Americium', 'Cm': 'Curium', 'Bk': 'Berkelium',
    'Cf': 'Californium', 'Es': 'Einsteinium', 'Fm': 'Fermium', 'Md': 'Mendelevium', 'No': 'Nobelium', 'Lr': 'Lawrencium',
    'Rf': 'Rutherfordium', 'Db': 'Dubnium', 'Sg': 'Seaborgium', 'Bh': 'Bohrium', 'Hs': 'Hassium', 'Mt': 'Meitnerium',
    'Ds': 'Darmstadtium', 'Rg': 'Roentgenium', 'Cn': 'Copernicium', 'Nh': 'Nihonium', 'Fl': 'Flerovium', 'Mc': 'Moscovium',
    'Lv': 'Livermorium', 'Ts': 'Tennessine', 'Og': 'Oganesson'
}
ELEMENT_FULLNAME_TOSYMBOL = {
    'Hydrogen': 'H', 'Helium': 'He', 'Lithium': 'Li', 'Beryllium': 'Be', 'Boron': 'B', 'Carbon': 'C', 'Nitrogen': 'N',
    'Oxygen': 'O', 'Fluorine': 'F', 'Neon': 'Ne', 'Sodium': 'Na', 'Magnesium': 'Mg', 'Aluminium': 'Al',
    'Silicon': 'Si', 'Phosphorus': 'P', 'Sulfur': 'S', 'Chlorine': 'Cl', 'Argon': 'Ar', 'Potassium': 'K',
    'Calcium': 'Ca', 'Scandium': 'Sc', 'Titanium': 'Ti', 'Vanadium': 'V', 'Chromium': 'Cr', 'Manganese': 'Mn',
    'Iron': 'Fe', 'Cobalt': 'Co', 'Nickel': 'Ni', 'Copper': 'Cu', 'Zinc': 'Zn', 'Gallium': 'Ga',
    'Germanium': 'Ge', 'Arsenic': 'As', 'Selenium': 'Se', 'Bromine': 'Br', 'Krypton': 'Kr', 'Rubidium': 'Rb',
    'Strontium': 'Sr', 'Yttrium': 'Y', 'Zirconium': 'Zr', 'Niobium': 'Nb', 'Molybdenum': 'Mo', 'Technetium': 'Tc',
    'Ruthenium': 'Ru', 'Rhodium': 'Rh', 'Palladium': 'Pd', 'Silver': 'Ag', 'Cadmium': 'Cd', 'Indium': 'In',
    'Tin': 'Sn', 'Antimony': 'Sb', 'Tellurium': 'Te', 'Iodine': 'I', 'Xenon': 'Xe', 'Caesium': 'Cs',
    'Barium': 'Ba', 'Lanthanum': 'La', 'Cerium': 'Ce', 'Praseodymium': 'Pr', 'Neodymium': 'Nd', 'Promethium': 'Pm',
    'Samarium': 'Sm', 'Europium': 'Eu', 'Gadolinium': 'Gd', 'Terbium': 'Tb', 'Dysprosium': 'Dy', 'Holmium': 'Ho',
    'Erbium': 'Er', 'Thulium': 'Tm', 'Ytterbium': 'Yb', 'Lutetium': 'Lu', 'Hafnium': 'Hf', 'Tantalum': 'Ta',
    'Tungsten': 'W', 'Rhenium': 'Re', 'Osmium': 'Os', 'Iridium': 'Ir', 'Platinum': 'Pt', 'Gold': 'Au',
    'Mercury': 'Hg', 'Thallium': 'Tl', 'Lead': 'Pb', 'Bismuth': 'Bi', 'Polonium': 'Po', 'Astatine': 'At',
    'Radon': 'Rn', 'Francium': 'Fr', 'Radium': 'Ra', 'Actinium': 'Ac', 'Thorium': 'Th', 'Protactinium': 'Pa',
    'Uranium': 'U', 'Neptunium': 'Np', 'Plutonium': 'Pu', 'Americium': 'Am', 'Curium': 'Cm', 'Berkelium': 'Bk',
    'Californium': 'Cf', 'Einsteinium': 'Es', 'Fermium': 'Fm', 'Mendelevium': 'Md', 'Nobelium': 'No', 'Lawrencium': 'Lr',
    'Rutherfordium': 'Rf', 'Dubnium': 'Db', 'Seaborgium': 'Sg', 'Bohrium': 'Bh', 'Hassium': 'Hs', 'Meitnerium': 'Mt',
    'Darmstadtium': 'Ds', 'Roentgenium': 'Rg', 'Copernicium': 'Cn', 'Nihonium': 'Nh', 'Flerovium': 'Fl', 'Moscovium': 'Mc',
    'Livermorium': 'Lv', 'Tennessine': 'Ts', 'Oganesson': 'Og'
}
PERIODIC_TABLE_TOINDEX = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 
    'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13,
    'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
    'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
    'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
    'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
    'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
    'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
    'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
    'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61,
    'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67,
    'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73,
    'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
    'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85,
    'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
    'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
    'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103,
    'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
    'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115,
    'Lv': 116, 'Ts': 117, 'Og': 118
}
PERIODIC_TABLE_TOSYMBOL = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N',
    8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al',
    14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K',
    20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn',
    26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga',
    32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb',
    38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
    44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In',
    50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs',
    56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm',
    62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho',
    68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta',
    74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au',
    80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
    86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa',
    92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk',
    98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr',
    104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt',
    110: 'Ds', 111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc',
    116: 'Lv', 117: 'Ts', 118: 'Og'
}
RVDW_ZERO_CHG = {
    "H": 1.96, "C": 1.85, "Li": 2.72, "Si": 2.25, "Na": 2.82, "Ge": 2.23, "K": 3.08, "Sn": 2.34, 
    "Rb": 3.22, "Pb": 2.34, "Cs": 3.38, "Nb": 2.5, "Cu": 2.3, "Ta": 2.44, "Ag": 2.34, "N": 1.7, 
    "Be": 2.32, "P": 2.09, "Mg": 2.45, "As": 2.16, "Ca": 2.77, "Sb": 2.33, "Sr": 2.9, "Bi": 2.4, 
    "Ba": 3.05, "Cr": 2.23, "Zn": 2.25, "Mo": 2.4, "Cd": 2.32, "W": 2.35, "Hg": 2.25, "O": 1.64, 
    "Sc": 2.64, "S": 2.0, "Y": 2.73, "Se": 2.1, "La": 2.86, "Te": 2.3, "B": 2.05, "Mn": 2.29, 
    "Al": 2.47, "Re": 2.38, "Ga": 2.38, "Br": 2.0, "In": 2.44, "I": 2.15, "Tl": 2.46, "Fe": 2.34, 
    "Ti": 2.52, "Co": 2.3, "Zr": 2.63, "Ni": 2.26, "Hf": 2.54, "Th": 2.78, "U": 2.8
}
RVDW_EQ = {
    "H": 1.56, "C": 1.97, "Li": 2.46, "Si": 2.27, "Na": 2.68, "Ge": 2.42, "K": 3.07, "Sn": 2.57, 
    "Rb": 3.23, "Pb": 2.72, "Cs": 3.42, "Nb": 2.41, "Cu": 2.24, "Ta": 2.41, "Ag": 2.41, "N": 1.88, 
    "Be": 2.14, "P": 2.2, "Mg": 2.41, "As": 2.34, "Ca": 2.79, "Sb": 2.5, "Sr": 2.98, "Bi": 2.64, 
    "Ba": 3.05, "Cr": 2.23, "Zn": 2.27, "Mo": 2.37, "Cd": 2.48, "W": 2.37, "Hg": 2.51, "O": 1.78, 
    "Sc": 2.59, "S": 2.13, "Y": 2.69, "Se": 2.27, "La": 2.76, "Te": 2.42, "B": 2.06, "Mn": 2.22, 
    "Al": 2.34, "Re": 2.35, "Ga": 2.44, "Br": 2.2, "In": 2.62, "I": 2.34, "Tl": 2.57, "Fe": 2.21, 
    "Ti": 2.37, "Co": 2.21, "Zr": 2.52, "Ni": 2.2, "Hf": 2.51, "Th": 2.72, "U": 2.5
}
RVDW_CRYSTAL = {
    "H": 1.2, "F": 1.47, "Cl": 1.75, "Br": 1.85, "I": 1.98, "O": 1.52, "S": 1.8, "N": 1.55, "C": 1.7, 
    "Li": 2.24, "Na": 2.57, "K": 3.0, "Rb": 3.12, "Cs": 3.31, "Cu": 2.0, "Ag": 2.13, "Au": 2.13, 
    "Be": 1.86, "Mg": 2.27, "Ca": 2.61, "Sr": 2.78, "Ba": 2.85, "Zn": 2.02, "Cd": 2.17, "Hg": 2.17, 
    "Sc": 2.28, "Y": 2.45, "La": 2.51, "B": 1.74, "Al": 2.11, "Ga": 2.08, "In": 2.24, "Tl": 2.25, 
    "Ti": 2.14, "Zr": 2.25, "Hf": 2.24, "Si": 2.06, "Ge": 2.13, "Sn": 2.29, "Pb": 2.36, "V": 2.03, 
    "Nb": 2.13, "Ta": 2.13, "As": 2.16, "Sb": 2.33, "Bi": 2.42, "Cr": 1.97, "Mo": 2.06, "W": 2.07, 
    "Mn": 1.96, "Tc": 2.04, "Re": 2.05, "Fe": 1.96, "Co": 1.95, "Ni": 1.94, "Ru": 2.02, "Rh": 2.02, 
    "Pd": 2.05, "Os": 2.03, "Ir": 2.03, "Pt": 2.06, "Th": 2.43, "U": 2.17
}
RVDW_CALCULATED = {
    "Li": 1.9, "B": 1.2, "P": 1.63, "Br": 1.73, "Na": 2.32, "Al": 1.75, "As": 1.81, "I": 1.98, 
    "K": 2.88, "Ga": 1.75, "Sb": 2.03, "Mn": 1.66, "Rb": 3.04, "In": 1.96, "Bi": 2.17, "Tc": 1.73, 
    "Cs": 3.27, "Tl": 1.98, "V": 1.72, "Re": 1.75, "Cu": 1.73, "Sc": 1.98, "Nb": 1.86, "Fe": 1.65, 
    "Ag": 1.77, "Y": 2.2, "Ta": 1.87, "Co": 1.64, "Au": 1.86, "La": 2.29, "S": 1.73, "Ni": 1.63, 
    "Be": 1.38, "Si": 1.68, "Se": 1.9, "Ru": 1.81, "Mg": 1.96, "Ge": 1.77, "Te": 2.14, "Rh": 1.75, 
    "Ca": 2.41, "Sn": 1.99, "Cr": 1.67, "Pd": 1.86, "Sr": 2.63, "Pb": 2.09, "Mo": 1.76, "Os": 1.83, 
    "Ba": 2.71, "Ti": 1.8, "W": 1.77, "Ir": 1.77, "Zn": 1.77, "Zr": 1.96, "Th": 2.2, "Pt": 1.87, 
    "Cd": 1.98, "Hf": 1.94, "U": 2.14, "Hg": 1.98
}
RVDW_WIKIPEDIA = {
    "He": 1.4, "Ne": 1.54, "Ar": 1.88, "Kr": 2.02, "Xe": 2.16, "Rn": 2.3, "H": 1.2,
    "Li": 1.82, "Na": 2.27, "K": 2.75, "Rb": 3.03, "Cs": 3.43, "Fr": 3.48, "Be": 1.53,
    "Mg": 1.73, "Ca": 2.31, "Sr": 2.49, "Ba": 2.68, "Ra": 2.83, "B": 1.92, "Al": 1.84,
    "Ga": 1.87, "In": 1.93, "Tl": 1.96, "C": 1.7, "Si": 2.1, "Ge": 2.11, "Sn": 2.17,
    "Pb": 2.2, "N": 1.55, "P": 1.8, "As": 1.85, "Sb": 2.06, "Bi": 2.07, "O": 1.52,
    "S": 1.8, "Se": 1.9, "Te": 2.06, "Po": 1.97, "F": 1.47, "Cl": 1.75, "Br": 1.85,
    "I": 1.98, "At": 2.02, "U": 1.86
}
RCOVALENT = {
    "Fr": 2.6, "Cs": 2.44, "Ra": 2.21, "Rb": 2.2, "Ac": 2.15, "Ba": 2.15, "La": 2.07, 
    "Th": 2.06, "Ce": 2.04, "Pr": 2.03, "K": 2.03, 
    "Nd": 2.01, "Pa": 2.0, "Pm": 1.99, "Eu": 1.98, "Sm": 1.98, "U": 1.96, 
    "Gd": 1.96, "Sr": 1.95, "Tb": 1.94, "Ho": 1.92, "Dy": 1.92, "Np": 1.90, 
    "Tm": 1.90, "Y": 1.90, "Er": 1.89, 
    "Pu": 1.87, "Lu": 1.87, "Yb": 1.87, "Am": 1.8, "Ca": 1.76, "Hf": 1.75, "Zr": 1.75, 
    "Ta": 1.7, "Sc": 1.7, "Cm": 1.69, "Na": 1.66, "Nb": 1.64, 
    "W": 1.62, "Ti": 1.6, "Mo": 1.54, "V": 1.53, "Re": 1.51, "Rn": 1.5, "At": 1.5, 
    "Bi": 1.48, "Tc": 1.47, "Pb": 1.46, "Ru": 1.46, "Tl": 1.45, "Ag": 1.45, "Os": 1.44, 
    "Cd": 1.44, "In": 1.42, "Rh": 1.42, "Ir": 1.41, "Mg": 1.41, "Po": 1.40, 
    "Xe": 1.40, "I": 1.39, "Sb": 1.39, 
    "Sn": 1.39, "Pd": 1.39, "Mn": 1.39, 
    "Cr": 1.39, "Te": 1.38, "Au": 1.36, "Pt": 1.36, "Hg": 1.32, 
    "Cu": 1.32, "Fe": 1.32, "Li": 1.28, "Co": 1.26, "Ni": 1.24, "Ga": 1.22, "Zn": 1.22, 
    "Al": 1.21, "Br": 1.2, "Se": 1.2, "Ge": 1.2, "As": 1.19, "Kr": 1.16, "Si": 1.11, "P": 1.07, 
    "Ar": 1.06, "S": 1.05, "Cl": 1.02, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, 
    "Ne": 0.58, "F": 0.57, "H": 0.31, "He": 0.28
}

def element_index_tolabel(label: int) -> str:

    return str(label) + "_" + PERIODIC_TABLE_TOSYMBOL[label]

def element_label_toindex(label: str) -> int:
    
    return PERIODIC_TABLE_TOINDEX[label]

def symbol_tol(label: str):

    data = ["s", "p", "d", "f", "g", "h", "i", "k", "l", "m", "n"]
    if label in data:
        return data.index(label)
    else:
        raise RuntimeError("too high angular momentum")

def l_tosymbol(l: int) -> str:

    data = ["s", "p", "d", "f", "g", "h", "i", "k", "l", "m", "n"]
    if l in range(1, 11):
        return data[l - 1]
    elif l == 0:
        return "s"
    else:
        raise RuntimeError("too high angular momentum")

def element_label_tomass(label: str) -> float:
    """
    element mass, partially from NIST, National Institute of Standards and Technology
    generated by Github.copilot
    """
    mass = {
        "H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182,
        "B": 10.811, "C": 12.0107, "N": 14.0067, "O": 15.9994,
        "F": 18.9984032, "Ne": 20.1797, "Na": 22.98976928,
        "Mg": 24.305, "Al": 26.9815386, "Si": 28.0855, "P": 30.973762,
        "S": 32.065, "Cl": 35.453, "Ar": 39.948, "K": 39.0983,
        "Ca": 40.078, "Sc": 44.955912, "Ti": 47.867, "V": 50.9415,
        "Cr": 51.9961, "Mn": 54.938045, "Fe": 55.845, "Co": 58.933195,
        "Ni": 58.6934, "Cu": 63.546, "Zn": 65.409, "Ga": 69.723,
        "Ge": 72.64, "As": 74.9216, "Se": 78.96, "Br": 79.904,
        "Kr": 83.798, "Rb": 85.4678, "Sr": 87.62, "Y": 88.90585,
        "Zr": 91.224, "Nb": 92.90638, "Mo": 95.94, "Tc": 98,
        "Ru": 101.07, "Rh": 102.9055, "Pd": 106.42, "Ag": 107.8682,
        "Cd": 112.411, "In": 114.818, "Sn": 118.71, "Sb": 121.76,
        "Te": 127.6, "I": 126.90447, "Xe": 131.293, "Cs": 132.9054519,
        "Ba": 137.327, "La": 138.90547, "Ce": 140.116, "Pr": 140.90765,
        "Nd": 144.242, "Pm": 145, "Sm": 150.36, "Eu": 151.964,
        "Gd": 157.25, "Tb": 158.92535, "Dy": 162.5, "Ho": 164.93032,
        "Er": 167.259, "Tm": 168.93421, "Yb": 173.054, "Lu": 174.9668,
        "Hf": 178.49, "Ta": 180.94788, "W": 183.84, "Re": 186.207,
        "Os": 190.23, "Ir": 192.217, "Pt": 195.084, "Au": 196.966569,
        "Hg": 200.59, "Tl": 204.3833, "Pb": 207.2, "Bi": 208.9804,
        "Po": 209, "At": 210, "Rn": 222, "Fr": 223, "Ra": 226,
        "Ac": 227, "Th": 232.03806, "Pa": 231.03588, "U": 238.02891,
        "Np": 237, "Pu": 244, "Am": 243, "Cm": 247, "Bk": 247,
        "Cf": 251, "Es": 252, "Fm": 257, "Md": 258, "No": 259,
        "Lr": 262, "Rf": 261, "Db": 262, "Sg": 266, "Bh": 264,
        "Hs": 277, "Mt": 278, "Ds": 281, "Rg": 282, "Cn": 285,
        "Nh": 286, "Fl": 289, "Mc": 289, "Lv": 293, "Ts": 294,
        "Og": 294
    }
    # remove all numbers in label
    label = ''.join([i for i in label if not i.isdigit()])
    return mass[label]

def number_tomultiplicity(n: int) -> str:

    data = ["s", "d", "t", "q", "p", "h"]
    if n in range(1, 7):
        return data[n - 1]
    else:
        raise ValueError("too high multiplicity")
    
def multiplicity_tonumber(m: str) -> int:

    data = ["s", "d", "t", "q", "p", "h"]
    if m in data:
        return data.index(m) + 1
    else:
        raise RuntimeError("too high multiplicity")

def radius_interpolation(element: str, radius_data: dict):
    """Interpolate rvdw data for elements not in the database."""
    if element not in radius_data.keys():
        left = PERIODIC_TABLE_TOINDEX[element]
        while True:
            left -= 1
            if PERIODIC_TABLE_TOSYMBOL[left] in radius_data.keys():
                left = PERIODIC_TABLE_TOSYMBOL[left]
                break
            if left < 1:
                raise ValueError("Element %s not in database, interpolation is also not possible" % element)
        right = PERIODIC_TABLE_TOINDEX[element]
        while True:
            right += 1
            if PERIODIC_TABLE_TOSYMBOL[right] in radius_data.keys():
                right = PERIODIC_TABLE_TOSYMBOL[right]
                break
            if right > 118:
                raise ValueError("Element %s not in database, interpolation is also not possible" % element)
        """linear interpolation"""
        distance_left = PERIODIC_TABLE_TOINDEX[element] - PERIODIC_TABLE_TOINDEX[left]
        distance_right = PERIODIC_TABLE_TOINDEX[right] - PERIODIC_TABLE_TOINDEX[element]
        return (radius_data[left] * distance_right + radius_data[right] * distance_left) / (distance_left + distance_right)
    else:
        return radius_data[element]

def element_label_torvdw(element: str, method: str = "0chg"):
    """Get Van der Waals radius of an element.
    
    Reference:
    0chg: Komissarov, A.V. and Heaven, M.C., J. Chem. Phys., 2000, vol. 113, no. 5, pp. 1775-1780
    eq: Komissarov, A.V. and Heaven, M.C., J. Chem. Phys., 2000, vol. 113, no. 5, pp. 1775-1780
    crystal: Bondi, A., J. Phys. Chem., 1964, vol. 68, no. 3, pp. 441-451.
             Batsanov, S.S., Zh. Fiz. Khim., 2000, vol. 74, no. 7, pp. 1273-1276.
    calculated: Batsanov, S.S., Inorganic Materials, Vol. 37, No. 9, 2001, pp. 871-885."""
    if method == "0chg":
        return radius_interpolation(element, RVDW_ZERO_CHG)
    elif method == "eq":
        return radius_interpolation(element, RVDW_EQ)
    elif method == "crystal":
        return radius_interpolation(element, RVDW_CRYSTAL)
    elif method == "calculated":
        return radius_interpolation(element, RVDW_CALCULATED)
    else:
        raise ValueError("Method %s is not supported." % method)

def element_label_toradius(element: str, radius_kind: str = "covalent", **kwargs):
    """Get radius of an element.
    
    Args:
        element (str): Element symbol.
        radius_kind (str): Radius kind, can be "covalent", "vdw", "ionic", "atomic".
        **kwargs: Parameters for radius calculation.
    
    Returns:
        float: Radius of the element.
    """
    if radius_kind == "covalent":
        return radius_interpolation(element, RCOVALENT)
    elif radius_kind == "vdw":
        if "method" in kwargs.keys():
            return element_label_torvdw(element, method=kwargs["method"])
        else:
            return element_label_torvdw(element, method="0chg")
    else:
        raise ValueError("Radius kind %s is not supported." % radius_kind)

def unit_conversion(val: float, unitfrom: str = "eV", unitto: str = "eV") -> float:
    # energy unit takes eV as the base
    eunit = {"eV": 1.0, "meV": 1e3, 
             "Ry": 1/13.6, "Rydberg": 1/13.6, 
             "Hartree": 1/27.21138602, "Ha": 1/27.21138602,
             "a.u.": 1/27.21138602, "au": 1/27.21138602,
             "J": 1.6e-19, "kJ": 1.6e-22, "kcal/mol": 23.1}
    # length unit takes Angstrom as the base
    lunit = {"A": 1.0, "Angstrom": 1.0, 
             "Bohr": 1.889726877, "bohr": 1.889726877, 
             "a.u.": 1.889726877, "au": 1.889726877,
             "fm": 1e5, "pm": 100.0, "nm": 0.1, "um": 1e-4, "mm": 1e-7, "cm": 1e-8, "m": 1e-10, "km": 1e-13}
    # weight unit takes atomic unit as the base
    wtunit = {"a.u.": 1.0, "kg": 9.1e-31, "g": 9.1e-28, "me": 0.00054858}

    units = [eunit, lunit, wtunit]
    # both unit should belong to the same category
    for u in units:
        if unitfrom in u and unitto in u:
            return val * u[unitto] / u[unitfrom]
    raise ValueError("Units are not in the same category")

import re
def orbital_configration2list(configuration, lmax):
    """Convert orbital configuration to list.
    
    Args:
        configuration (str): Orbital configuration, e.g. '1s1d'.
        lmax (int): Maximum angular momentum, this param will be used to determine the length of the list.
    
    Returns:
        list: Orbital configuration list, e.g. [1, 0, 0, 0, 1, 0, 0, 0]
    """
    l_label = ["s", "p", "d", "f", "g", "h", "i", "j", "k", "l"]
    orb_config_unit_pattern = r"(\d+)([spdfghijkl])"
    _match = re.findall(orb_config_unit_pattern, configuration)
    _result = [0] * (lmax + 1)
    for _m in _match:
        _result[l_label.index(_m[1])] = int(_m[0])
    return _result

import unittest
class TestDatabase(unittest.TestCase):
    
        def test_element_index_tolabel(self):
    
            self.assertEqual(element_index_tolabel(1), "1_H")
            self.assertEqual(element_index_tolabel(118), "118_Og")
    
        def test_element_label_toindex(self):
    
            self.assertEqual(element_label_toindex("H"), 1)
            self.assertEqual(element_label_toindex("Og"), 118)
    
        def test_symbol_tol(self):
    
            self.assertEqual(symbol_tol("s"), 0)
            self.assertEqual(symbol_tol("n"), 10)
    
        def test_l_tosymbol(self):
    
            self.assertEqual(l_tosymbol(0), "s")
    
        def test_element_label_tomass(self):
    
            self.assertAlmostEqual(element_label_tomass("H"), 1.00794, places=5)
            self.assertAlmostEqual(element_label_tomass("Og"), 294, places=5)
    
        def test_number_tomultiplicity(self):
    
            self.assertEqual(number_tomultiplicity(1), "s")
            self.assertEqual(number_tomultiplicity(6), "h")
    
        def test_multiplicity_tonumber(self):
    
            self.assertEqual(multiplicity_tonumber("s"), 1)
            self.assertEqual(multiplicity_tonumber("h"), 6)
    
        def test_radius_interpolation(self):
    
            self.assertAlmostEqual(radius_interpolation("H", RVDW_ZERO_CHG), 1.96, places=5)
    
        def test_element_label_torvdw(self):
    
            self.assertAlmostEqual(element_label_torvdw("H", method="0chg"), 1.96, places=5)
            
        def test_unit_conversion(self):
        
            self.assertEqual(unit_conversion(1, "Hartree", "eV"), 27.21138602)
            self.assertEqual(unit_conversion(1, "Hartree", "kcal/mol"), 627.509469)
            self.assertEqual(unit_conversion(1, "eV", "Hartree"), 1 / 27.21138602)
            self.assertEqual(unit_conversion(1, "eV", "kcal/mol"), 627.509469 / 27.21138602)
            self.assertEqual(unit_conversion(1, "kcal/mol", "Hartree"), 1 / 627.509469)
            self.assertEqual(unit_conversion(1, "kcal/mol", "eV"), 27.21138602 / 627.509469)
        
        def test_orbital_configration2list(self):

            self.assertListEqual(orbital_configration2list("1s1d", 3), [1, 0, 1, 0])
            self.assertListEqual(orbital_configration2list("1s1d", 4), [1, 0, 1, 0, 0])
            self.assertListEqual(orbital_configration2list("1s1d", 5), [1, 0, 1, 0, 0, 0])

if __name__ == "__main__":
    unittest.main()