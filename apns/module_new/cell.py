class AtomSpecies:
    """this class stores the information of an atom species, including
    basically the atomic number, the atomic mass, covalent radii, and
    its pseudopotentials, numerical atomic orbitals"""
    name, fullname, symbol, index, rcovalent, mass, magmom = \
        None, None, None, None, None, None, None
    pp, nao = None, None
    
    def __init__(self, symbol: str, name: str = None) -> None:
        assert isinstance(symbol, str), f'symbol should be a string: {symbol}'
        assert isinstance(name, str) or name is None, f'name should be a string: {name}'
        from apns.module_database import database
        self.name = name if name is not None else symbol
        self.fullname = database.PERIODIC_TABLE_TOFULLNAME[symbol]
        self.symbol = symbol
        self.index = database.PERIODIC_TABLE_TOINDEX[symbol]
        self.rcovalent = database.RCOVALENT[symbol]
        self.mass = database.ATOMIC_MASS[symbol]
    
    def set_pp(self, pseudo_dir: str, pptags: list):
        """set the pseudopotentials of this atom species by the tags"""
        from os.path import join as pjoin
        from apns.module_database.search import TagSearcher
        searcher = TagSearcher(pjoin(pseudo_dir, 'database.json'))
        self.pp = searcher(False, False, *pptags)
    
    def set_nao(self, orbital_dir: str, naotags: list):
        """set the numerical atomic orbitals of this atom species by the tags"""
        from os.path import join as pjoin
        from apns.module_database.search import TagSearcher
        searcher = TagSearcher(pjoin(orbital_dir, 'database.json'))
        self.nao = searcher(False, False, *naotags)

    # allow visit like a dict
    def __getitem__(self, key):
        return getattr(self, key)

class Cell:
    """this class stores one structure and its derivatives, such as
    varying the characteristic length, or the bond length, etc"""
    # lattice parameters
    a, b, c, alpha, beta, gamma = None, None, None, None, None, None
    # specific structures
    species, coords, species_map = None, None, None
    # whether the structure is periodic or not. For simple molecules, they are in principle not periodic
    periodic = None
    # scaling factors
    scale_fs = None
    # kpoints
    high_symmetry_ks, mpmesh_nks = None, None

    def __init__(self) -> None:
        pass

    def kernel_build(self, species: list, lat: list, coords: list, species_map: list):
        """build structure(s) from given species, lattice, and coordinates

        Args:
            species (list): list of AtomSpecies' or dicts
            lat (list): a, b, c and alpha, beta, gamma
            coords (list): coordinates of atoms, size = (n_atoms, 3)
            species_map (list): mapping the index of coords to species, size = n_atoms
        """
        import numpy as np
        self.species = species
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = lat
        self.coords = np.array(coords)
        self.species_map = species_map

    def standard_build(self, species: list, fname: str, **kwargs):
        """read structure file from external, and construct the structure.
        Currently only support CIF file.
        
        Args:
            species (list): list of AtomSpecies' or dicts
            fname (str): filename of the structure file
        Additional keywords will be assumed to be for generating kpath

        Examples:
        ```python
        from apns.module_new.cell import Cell
        species = [AtomSpecies('Si')]
        cell = Cell()
        cell.standard_build(species, 'Si.cif')
        ```
        """
        self.species = species
        self.lat0 = 1.8897259886 # Bohr to Angstrom
        from pymatgen.io.cif import CifParser
        parser = CifParser(fname)
        structure = parser.parse_structures(primitive=True)[0]
        # check consistency of species in the structure and the species list
        symbls = [spec.symbol for spec in species]
        symbls_ = [site.species.elements[0].symbol for site in structure]
        assert set(symbls) == set(symbls_), f'species not consistent: {symbls} {symbls_}'
        self.coords = structure.cart_coords
        self.species_map = [symbls.index(symbls_[i]) for i in range(len(structure))]
        self.a, self.b, self.c = structure.lattice.abc
        self.alpha, self.beta, self.gamma = structure.lattice.angles

        import seekpath as skps
        cell = self.cell_to_vec([self.a, self.b, self.c, self.alpha, self.beta, self.gamma], True)
        positions = self.coords.tolist()
        numbers = self.species_map
        self.high_symmetry_ks = skps.get_path((cell, positions, numbers), **kwargs)['point_coords']

    def ideal_build(self, species: list, config: str, characteristic: float):
        """build structures with simple shapes, such as ideal Bravis lattices,
        including SimpleCubic (SC), FaceCenteredCubic (FCC), BodyCenteredCubic (BCC),
        Diamond (diamond) and simple molecules such as dimer, trimer, tetramer, etc.
        
        Args:
            species (list): list of AtomSpecies' or dicts
            config (str): the configuration of the structure
            characteristic (float): the characteristic length of the structure
        
        Examples:
        ```python
        from apns.module_new.cell import Cell
        species = [AtomSpecies('Si')]
        cell = Cell()
        cell.ideal_build(species, 'sc', 5.43)
        ```
        """
        from apns.module_new.ideal import structure as ideal
        pure = config.lower() in ['sc', 'fcc', 'bcc', 'diamond', 'dimer', 'trimer', 'tetramer']
        periodic = config.lower() in ['sc', 'fcc', 'bcc', 'diamond', 'x2y', 'x2y3', 'x2y5', 'xy', 'xy2', 'xy3']
        composite = not pure and periodic
        config = config.lower() if pure else config.upper()
        assert (config == config.upper() and composite) or (config == config.lower() and pure)
        self.periodic = periodic
        # True, True: ideal sc, fcc, bcc, diamond
        # True, False: dimer, trimer, tetramer
        # False, True: X2Y, X2Y3, X2Y5, XY, XY2, XY3
        # False, False: not supported
        assert pure or periodic, f'config not supported for ideal_build: {config}'
        characteristic = self.vol_to_celldm(self.bravis_angles(config), -characteristic)\
            if characteristic < 0 and periodic else characteristic
        assert characteristic > 0, f'characteristic should be positive: {characteristic}'
        abc_angles, species_map, coords = ideal(config, characteristic)
        self.kernel_build(species, abc_angles, coords, species_map)

    def vol_to_celldm(self, angles: list, volume: float):
        """calculate the characteristic length from the volume of the cell"""
        import numpy as np
        assert len(angles) == 3, f'angles should be a list of 3 floats: {angles}'
        if np.linalg.norm(np.array(angles) - np.array([90, 90, 90])) < 1e-6:
            return volume**(1/3)
        elif np.linalg.norm(np.array(angles) - np.array([60, 60, 60])) < 1e-6:
            return np.sqrt(2)/2 * (4*volume)**(1/3)
        elif np.linalg.norm(np.array(angles) - np.array([109.4712206, 109.4712206, 109.4712206])) < 1e-6:
            return np.sqrt(3)/2 * (2*volume)**(1/3)
        else:
            raise NotImplementedError(f'angles not supported: {angles}')
    
    def cell_to_vec(self, lat: list, as_list: bool = False):
        """convert lattice parameters to vectors"""
        import numpy as np
        assert len(lat) == 6, f'lat should be a list of 6 floats: {lat}'
        a, b, c, alpha, beta, gamma = lat
        alpha, beta, gamma = np.deg2rad(alpha), np.deg2rad(beta), np.deg2rad(gamma)
        e11, e12, e13 = a, 0, 0
        e21, e22, e23 = b * np.cos(gamma), b * np.sin(gamma), 0
        e31 = c * np.cos(beta)
        e32 = (b * c * np.cos(alpha) - e21 * e31) / e22
        e33 = np.sqrt(c**2 - e31**2 - e32**2)
        return np.array([[e11, e12, e13], [e21, e22, e23], [e31, e32, e33]])\
            if not as_list else [[e11, e12, e13], [e21, e22, e23], [e31, e32, e33]]

    def bravis_angles(config: str):
        """get the angles of the simple perodic structures"""
        if config in ['sc', 'X2Y3', 'X2Y5']:
            return [90, 90, 90]
        elif config in ['fcc', 'diamond', 'X2Y', 'XY', 'XY2', 'XY3']:
            return [60, 60, 60]
        elif config in ['bcc']:
            return [109.4712206, 109.4712206, 109.4712206]
        else:
            raise NotImplementedError(f'easy_angle not supported: {config}')
    
    def rescale(self, scales: list):
        self.scale_fs = scales

    def to_abacus(self, suffix: str):
        """generate the STRU and KPT file(s) of ABACUS"""


import unittest
class TestAtomSpecies(unittest.TestCase):

    def test_constructor(self):
        Si = AtomSpecies('Si')
        self.assertEqual(Si.name, 'Si')
        self.assertEqual(Si.fullname, 'Silicon')
        self.assertEqual(Si.symbol, 'Si')
        self.assertEqual(Si.index, 14)
        self.assertEqual(Si.rcovalent, 1.11)
        self.assertEqual(Si.mass, 28.0855)
    
    def test_set_pp(self):
        import os
        import json
        # prepare a fake database
        database = {"the_file_name_with_path": ["tag1", "tag2", "Si"]}
        fdatabase = "./database.json"
        with open(fdatabase, "w") as f:
            json.dump(database, f)
        # test the search
        Si = AtomSpecies('Si')
        Si.set_pp(os.getcwd(), ["tag1", "tag2"])
        self.assertEqual(Si.pp, {"the_file_name_with_path"})
        os.remove(fdatabase)

class TestCell(unittest.TestCase):

    def test_constructor(self):
        cell = Cell()
        self.assertIsNone(cell.a)
        self.assertIsNone(cell.b)
        self.assertIsNone(cell.c)
        self.assertIsNone(cell.alpha)
        self.assertIsNone(cell.beta)
        self.assertIsNone(cell.gamma)
        self.assertIsNone(cell.species)
        self.assertIsNone(cell.coords)
        self.assertIsNone(cell.species_map)
        self.assertIsNone(cell.periodic)
        self.assertIsNone(cell.scale_fs)
        self.assertIsNone(cell.high_symmetry_ks)

    def test_standard_build(self):
        import os
        cif = """# generated using pymatgen
data_Ac2O
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a   4.84451701
_cell_length_b   4.84451701
_cell_length_c   4.84451701
_cell_angle_alpha   60.00000000
_cell_angle_beta   60.00000000
_cell_angle_gamma   60.00000000
_symmetry_Int_Tables_number   1
_chemical_formula_structural   Ac2O
_chemical_formula_sum   'Ac2 O1'
_cell_volume   80.39637305
_cell_formula_units_Z   1
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Ac  Ac0  1  0.75000000  0.75000000  0.75000000  1
  Ac  Ac1  1  0.25000000  0.25000000  0.25000000  1
  O  O2  1  0.00000000  0.00000000  0.00000000  1"""
        fcif = 'Ac2O.cif'
        species = [AtomSpecies('Ac'), AtomSpecies('O')]
        cell = Cell()
        with open(fcif, "w") as f:
            f.write(cif)
        cell.standard_build(species, fcif)
        self.assertAlmostEqual(cell.a, 4.84451701, delta=1e-5)
        self.assertAlmostEqual(cell.b, 4.84451701, delta=1e-5)
        self.assertAlmostEqual(cell.c, 4.84451701, delta=1e-5)
        self.assertAlmostEqual(cell.alpha, 60.00000000, delta=1e-5)
        self.assertAlmostEqual(cell.beta, 60.00000000, delta=1e-5)
        self.assertAlmostEqual(cell.gamma, 60.00000000, delta=1e-5)
        self.assertEqual(cell.species, species)
        self.assertEqual(cell.coords.shape, (3, 3))
        self.assertEqual(cell.species_map, [0, 0, 1])
        #print(cell.high_symmetry_ks)
        os.remove(fcif)

    def test_build(self):
        species = [AtomSpecies('Si')]
        cell = Cell()
        cell.kernel_build(species, [5, 5, 5, 90, 90, 90], [[0, 0, 0]], [0])
        self.assertEqual(cell.a, 5)
        self.assertEqual(cell.b, 5)
        self.assertEqual(cell.c, 5)
        self.assertEqual(cell.alpha, 90)
        self.assertEqual(cell.beta, 90)
        self.assertEqual(cell.gamma, 90)
        self.assertEqual(cell.species, species)
        self.assertEqual(cell.coords.shape, (1, 3))
        self.assertEqual(cell.species_map, [0])

    def test_rescale(self):
        species = [AtomSpecies('Si')]
        cell = Cell()
        cell.kernel_build(species, [5, 5, 5, 90, 90, 90], [[0, 0, 0]], [0])
        cell.rescale([0.98, 1, 1.02])
        self.assertEqual(cell.scale_fs, [0.98, 1, 1.02])
    
    def test_easy_build(self):
        species = [AtomSpecies('Si')]
        cell = Cell()
        cell.ideal_build(species, 'sc', 5.43)
        self.assertEqual(cell.a, 5.43)
        self.assertEqual(cell.b, 5.43)
        self.assertEqual(cell.c, 5.43)
        self.assertEqual(cell.alpha, 90)
        self.assertEqual(cell.beta, 90)
        self.assertEqual(cell.gamma, 90)
        self.assertEqual(cell.species, species)
        self.assertEqual(cell.coords.shape, (1, 3))
        self.assertEqual(cell.species_map, [0])


if __name__ == "__main__":
    unittest.main()
                         