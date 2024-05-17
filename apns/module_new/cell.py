"""
Refactor notes:
this is the file involved in refactor project of APNS for dealing with the more and more
complicated pseudopot-nao test cases.
"""

class AtomSpecies:
    """
    WHAT IS AN ATOMSPECIES?

    An atom species is such a object that know its name, symbol, index, covalent radius, mass,
    magnetic moment (in isolated state), pseudopotentials, and numerical atomic orbitals, (
    or other physical properties)

    WHAT IS NOT REALLY AN ATOMSPECIES KNOWS?

    the magmom on the fly, the position of one specific atom, or everything related to the
    simulation. Because it is not the property of one atom species, but the property of one
    object (but as coincidience) called "atom" in simulation.
    """
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
    
    def set_nao(self, pp: str, orbital_dir: str, naotags: list):
        """set the numerical atomic orbitals of this atom species by the tags"""
        from os.path import join as pjoin
        from apns.module_database.search import TagSearcher
        searcher = TagSearcher(pjoin(orbital_dir, 'database.json'))
        # there should be a subtag called pp, then in that pp, search those naotags corresponding orbitals
        self.nao = searcher(False, False, *naotags)

    # allow visit like a dict
    def __getitem__(self, key):
        return getattr(self, key)

class Cell:
    """
    WHAT IS A CELL?
    
    A cell is such a object that know its structure and lattice parameters, but this is just the
    half part, it also knows its reciprocal lattice information, such as the high symmetry k-points.
    From a kspacing, it can also know how its Monkhorst-Pack mesh should be.

    WHAT IS NOT REALLY A CELL KNOWS?
    
    the atom species. Actually the cell only knows there are many atomic sites and for each site, 
    there is a symbol. Or more abstractly, it knows symbols and kinds of sites that occupied by
    different atoms. The cell, does not own any AtomSpecies!
    """
    # basic
    periodic = None
    # real space lattice information
    a, b, c, alpha, beta, gamma = None, None, None, None, None, None
    lat0 = 1.8897259886
    # reciprocal space lattice information
    sym_ks, mpmesh_nks = None, None
    # "points" information
    coords, vels, magmoms, mobs = None, None, None, None
    symbols, species_map = None, None
    
    def __init__(self) -> None:
        """create one cell instance with the scaling factors"""
        pass

    def kernel_build(self, lat: list, coords: list, species_map: list):
        """build structure(s) from given species, lattice, and coordinates

        Args:
            lat (list): a, b, c and alpha, beta, gamma
            coords (list): coordinates of atoms, size = (n_atoms, 3)
            species_map (list): mapping the index of coords to species, size = n_atoms
        """
        import numpy as np
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = lat
        self.coords = np.array(coords)
        self.species_map = species_map

    def standard_build(self, fname: str, **kwargs):
        """read structure file from external, and construct the structure.
        Currently only support CIF file.
        """
        self.periodic = True
        from pymatgen.io.cif import CifParser
        parser = CifParser(fname)
        structure = parser.parse_structures(primitive=True)[0]
        symbols = [site.species.elements[0].symbol for site in structure]
        self.symbols = list(dict.fromkeys(symbols)) # find unique set but keep the sequence
        self.species_map = [self.symbols.index(symbl) for symbl in symbols]
        self.coords = structure.cart_coords
        self.magmoms = [0] * len(self.coords)
        self.a, self.b, self.c = structure.lattice.abc
        self.alpha, self.beta, self.gamma = structure.lattice.angles

        import seekpath as skps
        cell = Cell._abc_angles_to_vec([self.a, self.b, self.c, self.alpha, self.beta, self.gamma], True)
        self.sym_ks = skps.get_path((cell, self.coords.tolist(), self.species_map), **kwargs)['point_coords']

    def ideal_build(self, config: str, characteristic: float):
        """build structures with simple shapes, such as ideal Bravis lattices,
        including SimpleCubic (SC), FaceCenteredCubic (FCC), BodyCenteredCubic (BCC),
        Diamond (diamond) and simple molecules such as dimer, trimer, tetramer, etc.
        """
        import re
        symbol, config = config.split('_')
        pure = config.lower() in ['sc', 'fcc', 'bcc', 'diamond', 'dimer', 'trimer', 'tetramer']
        periodic = config.lower() in ['sc', 'fcc', 'bcc', 'diamond', 'x2y', 'x2y3', 'x2y5', 'xy', 'xy2', 'xy3']
        composite = not pure and periodic
        assert (re.match(r'^[A-Z][a-z]?$', symbol) and pure) or (re.match(r'^[A-Z][a-z]?[A-Z][a-z]?$', symbol) and composite), \
            f'symbol not supported for ideal_build: {symbol}'
        config = config.lower() if pure else config.upper()
        assert (config == config.upper() and composite) or (config == config.lower() and pure)
        self.periodic = periodic
        # True, True: ideal sc, fcc, bcc, diamond
        # True, False: dimer, trimer, tetramer
        # False, True: X2Y, X2Y3, X2Y5, XY, XY2, XY3
        # False, False: not supported
        assert pure or periodic, f'config not supported for ideal_build: {config}'
        characteristic = Cell._vcell_to_celldm(self.bravis_angles(config), -characteristic)\
            if characteristic < 0 and periodic else characteristic
        assert characteristic > 0, f'characteristic should be positive: {characteristic}'
        abc_angles, species_map, coords = Cell._structure(config, characteristic)
        self.kernel_build(abc_angles, coords, species_map)
        self.symbols = re.findall(r'([A-Z][a-z]?)', symbol)
        assert (len(self.symbols) == 2 and composite) or (len(self.symbols) == 1 and pure)

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

    def kmeshgen(self, kspacing: float|list[float]):
        """get the Monkhorst-Pack mesh"""
        import numpy as np
        kspacing = kspacing if isinstance(kspacing, list) else [kspacing] * 3
        vecs = np.array(Cell._abc_angles_to_vec([self.a, self.b, self.c, self.alpha, self.beta, self.gamma], True))
        recvecs = np.linalg.inv(vecs).T
        recvecs = 2*np.pi * recvecs
        norms = np.linalg.norm(recvecs, axis=1).tolist()
        assert len(norms) == len(kspacing), f'kspacing should be a list of 3 floats: {kspacing}'
        norms = [int(norm / kspac) for norm, kspac in zip(norms, kspacing)]
        self.mpmesh_nks = list(map(lambda x: max(1, x + 1), norms))

    #############################
    #  private static methods   #
    #############################
    
    def _vcell_to_celldm(angles: list, volume: float):
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
    
    def _abc_angles_to_vec(lat: list, as_list: bool = False):
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

    def _structure(config: str, characteristic: float):
        if config == 'sc':
            return Cell._simple_cubic(characteristic)
        elif config == 'fcc':
            return Cell._face_centered_cubic(characteristic)
        elif config == 'bcc':
            return Cell._body_centered_cubic(characteristic)
        elif config == 'diamond':
            return Cell._diamond(characteristic)
        elif config == 'dimer':
            return Cell._dimer(characteristic)
        elif config == 'trimer':
            return Cell._trimer(characteristic)
        elif config == 'tetramer':
            return Cell._tetramer(characteristic)
        elif config == 'X2Y':
            return Cell._X2Y(characteristic)
        elif config == 'X2Y3':
            return Cell._X2Y3(characteristic)
        elif config == 'X2Y5':
            return Cell._X2Y5(characteristic)
        elif config == 'XY':
            return Cell._XY(characteristic)
        elif config == 'XY2':
            return Cell._XY2(characteristic)
        elif config == 'XY3':
            return Cell._XY3(characteristic)
        else:
            raise NotImplementedError(f'config not supported: {config}')
        
    def _dimer(bl: float):
        """build dimer with bond length `bl`
        
        Args:
            bl (float): bond length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [20, 20, 20, 90, 90, 90], \
            [0, 0], np.array([[0, 0, 0], [bl, 0, 0]])

    def _trimer(bl: float):
        """build triangle with bond length `bl`
        
        Args:
            bl (float): bond length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [20, 20, 20, 90, 90, 90], \
            [0, 0, 0], np.array([[0, 0, 0], [bl, 0, 0], [bl/2, bl/2*3**0.5, 0]])

    def _tetramer(bl: float):
        """build tetrahedron with bond length `bl`
        
        Args:
            bl (float): bond length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [20, 20, 20, 90, 90, 90], \
            [0, 0, 0, 0], np.array([[0, 0, 0], [bl, 0, 0], [bl/2, bl/2*3**0.5, 0], [bl/2, bl/2*3**0.5/3, bl/2*3**0.5]])

    def _simple_cubic(celldm: float):
        """build simple cubic with characteristic length `celldm`
        
        Args:
            celldm (float): characteristic length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [celldm, celldm, celldm, 90, 90, 90], \
            [0], np.array([[0, 0, 0]])

    def _body_centered_cubic(celldm: float):
        """build body centered cubic with characteristic length `celldm`
        
        Args:
            celldm (float): characteristic length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [celldm, celldm, celldm, 109.4712206, 109.4712206, 109.4712206], \
            [0], np.array([[0, 0, 0]])

    def _face_centered_cubic(celldm: float):
        """build face centered cubic with characteristic length `celldm`
        
        Args:
            celldm (float): characteristic length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [celldm, celldm, celldm, 60.0, 60.0, 60.0], \
            [0], np.array([[0, 0, 0]])

    def _diamond(celldm: float):
        """build diamond with characteristic length `celldm`
        
        Args:
            celldm (float): characteristic length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [celldm, celldm, celldm, 60.0, 60.0, 60.0], \
            [0, 0], np.array([[0, 0, 0], [0.25, 0.25, 0.25]])

    def _X2Y(celldm: float):
        """build X2Y structure with characteristic length `celldm`
        
        Args:
            celldm (float): characteristic length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [celldm, celldm, celldm, 60, 60, 60], \
            [0]*2 + [1]*1, \
            np.array([[0.75, 0.75, 0.75], [0.25, 0.25, 0.25], [0, 0, 0]])

    def _X2Y3(celldm: float):
        """build X2Y3 structure with characteristic length `celldm`
        
        Args:
            celldm (float): characteristic length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [celldm, celldm, celldm, 90, 90, 90], \
            [0]*4 + [1]*6, \
            np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], 
                        [0.25, 0.75, 0.75], 
                        [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0], 
                        [0.5, 0, 0.5], [0, 0.5, 0.5]])

    def _X2Y5(celldm: float):
        """build X2Y5 structure with characteristic length `celldm`
        
        Args:
            celldm (float): characteristic length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [celldm, celldm, celldm, 90, 90, 90], \
            [0]*4 + [1]*10, \
            np.array([[0.75, 0.75, 0.75], [0.25, 0.25, 0.75], [0.25, 0.75, 0.25],
                        [0.75, 0.25, 0.25],
                        [0.25, 0.75, 0.75], [0.75, 0.25, 0.75], [0.75, 0.75, 0.25],
                        [0.25, 0.25, 0.25], [0.5, 0.5, 0], [0.5, 0, 0.5],
                        [0, 0.5, 0.5], [0.5, 0, 0], [0, 0.5, 0],
                        [0, 0, 0.5]])

    def _XY(celldm: float):
        """build XO structure with characteristic length `celldm`
        
        Args:
            celldm (float): characteristic length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [celldm, celldm, celldm, 60, 60, 60], \
            [0]*1 + [1]*1, \
            np.array([[0, 0, 0], [0.5, 0.5, 0.5]])

    def _XY2(celldm: float):
        """build XO2 structure with characteristic length `celldm`
        
        Args:
            celldm (float): characteristic length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [celldm, celldm, celldm, 60, 60, 60], \
            [0]*1 + [1]*2, \
            np.array([[0, 0, 0], [0.75, 0.75, 0.75], [0.25, 0.25, 0.25]])

    def _XY3(celldm: float):
        """build XO3 structure with characteristic length `celldm`
        
        Args:
            celldm (float): characteristic length in Angstrom
        Returns:
            lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
            species_map (list): mapping the index of coords to species, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [celldm, celldm, celldm, 60, 60, 60], \
            [0]*1 + [1]*3, \
            np.array([[0, 0, 0], 
                        [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]])

class StructureManager:
    """structure manager is for managing structures, including downloading, generating, etc."""
    atom_species: list[AtomSpecies] = []
    structures: list[Cell] = []
    desc = []
    
    def __init__(self, atom_species) -> None:
        """initialize the structure manager with a list of atom species
        Example:
        ```python
        elements = ['Si', 'O', 'Na', 'Rh']
        # create a structure manager with knowledge of Si, O, Na, Rh in-one-shot
        species = list(map(AtomSpecies, elements))
        # search for pseudopotentials of Si, O, Na, Rh in-one-shot
        # in which the pptags is a dict of the pseudopotential tags of each species
        map(lambda s: s.set_pp('./download/pseudopotentials/', pptags[s.symbol]), species)
        sm = StructureManager(species)
        # complete! the StructureManager is ready to use
        # AtomSpecies.set_nao function will be called later, after the pp's have been combined
        ```
        """
        assert isinstance(atom_species, list), f'atom_species should be a list: {atom_species}'
        self.atom_species = atom_species
    
    #######################
    # structure prototype #
    #######################
    """fucntions in this part is for importing structures without additional scalling, say
    not for EOS, not for varying bond lengths, ..."""
    def describe_structure(self, desc):
        """
        For lazy build. First describe, then build explicitly and output instantly.
        : if the cell can be built with several parameters, then those parameters are
        representation of those exact structures. Calculate on those parameters means
        to vary and generate structures.
        
        NEED A UNIFIED INTERFACE TO CREATE CELL:
        1. standard_build: fcif
        2. ideal_build: config, characteristic_length
        ```python
        # descriptor
        #         str                str                    list[float]                list[str] list[str]
        # (cif/bravis/molecule, fcif/sc/dimer, EOS_scalings/EOS_scalings/bond_lengths, pptags,   naotags)
        # pptags and naotags are associated.
        # once from pptags determines a list of pseudopotentials, 
        # for each pseudopotential, search numerical atomic orbital according to naotags.
        # example
        desc = [('cif', '/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/structures/Si.cif',
                [0.95, 0.97, 0.99, 1.01, 1.03], ['Si', 'PBE', 'NC', 'sg15', 'sr'], 
                ['6au', '100Ry' , '2s2p1d']),
                ('bravis', 'sc', [0.95, 0.97, 0.99, 1.01, 1.03], ['Si', 'PBE', 'NC', 'sg15', 'sr'], 
                None),
                ('molecule', 'dimer', [1.0, 1.02, 1.04, 1.06, 1.08], ['Si', 'PBE', 'NC', 'sg15', 'sr'], 
                ['6au', '100Ry' , '2s2p1d'])]
        sm = StructureManager(species)
        sm.describe_structure(desc)
        ```
        """
        import re
        import os
        assert isinstance(desc, list), f'desc should be a list: {desc}'
        assert all([isinstance(d, tuple) for d in desc]), f'each element in desc should be a tuple: {desc}'
        assert all([len(d) == 5 for d in desc]), f'each element in desc should be a tuple of 5 elements: {desc}'
        assert all([isinstance(d[0], str) and d[0] in ['cif', 'bravis', 'molecule'] for d in desc])
        assert all([isinstance(d[1], str) and \
               (re.match(r'([A-Z][a-z]?_[dimer|trimer|tetramer|sc|bcc|fcc|diamond])|([A-Z][a-z]?[A-Z][a-z]?_[xy2|xy3|x2y|x2y3|x2y5])'\
                or os.path.exists(d[1]))) for d in desc]), f'descriptor is neither a cif file nor a bravis/molecule type: {desc}'
        assert all([isinstance(d[2], list) and all([isinstance(x, float) for x in d[2]]) for d in desc]), \
            f'characteristic scaling should be list of floats: {desc}'
        assert all([isinstance(d[3], list) and all([isinstance(x, str) for x in d[3]]) for d in desc]), \
            f'pptags should be list of strings: {desc}'
        assert all([\
            (isinstance(d[4], list) and all([isinstance(x, str) for x in d[4]])) or (d[4] is None) \
                for d in desc]), f'naotags should be list of strings: {desc} or None in basis_type pw calculation'
        self.desc = desc
        """will add replace, pop and append functionalities later"""
    
    def inject_atomspecies(self, icell: int):
        """This function estabilishes the actual connection between the abstract Cell object whose
        "atom(s)" are just names of points rather than physically meaningful atom (in real world),
        while the AtomSpecies' are the "real" one: inject specific AtomSpecies instances to Cell, 
        which reads both symbols and species_map. Create a composited tuple, the first is desc and 
        the second is list of AtomSpecies."""
        # first find the AtomSpecies of the structure
        _cell = self.structures[icell]
        symbols = _cell.symbols
        assert len(symbols) == len(_cell.species_map), f'symbols and species_map should have the same length: {symbols}, {_cell.species_map}'
        # then find the AtomSpecies of the structure
        atom_species = [s for s in self.atom_species if s.symbol in symbols]
        assert len(atom_species) == len(symbols), f'atom_species should have the same length as symbols: {atom_species}, {symbols}'
        # then inject the AtomSpecies to the Cell
        self.structures[icell] = (_cell, atom_species)

    def build(self, pseudo_dir: str, orbital_dir: str):
        """expand all tests, this function is unique in APNS, but should not appear in
        ABACUS newly-refactored UnitCell module, becuase StructureManager is not exactly
        the StructureIterator or something, while it is less necessary to implement such
        a class or function."""
        # check all elements in self.desc are tuples of three elements, the first
        # is str, the second is str and the third is list of float
        # build should be based on each Cell. For one cell, iterate to generate set of AtomSpecies
        # then save and inject them after the exact cell is built
        import itertools as itools
        lcao = [d[4] for d in self.desc if d[4] is not None]
        assert len(lcao) == 0 or len(lcao) == len(self.desc), f'naotags should be all None or all not None: {lcao}'
        lcao = (len(lcao) == len(self.desc))
        for d in self.desc: # for each structure prototype
            # build AtomSpecies: search pseudopotential first
            map(lambda s: s.set_pp(pseudo_dir, d[3]), self.atom_species)
            pps_comb = [list(s.pp) for s in self.atom_species]
            pps_comb = list(itools.product(*pps)) # convert to combinations of pps of different species
            for pps in pps_comb: # for each AtomSpecies, one pseudopotential
                assert len(pps) == len(self.atom_species), f'pps should be a list of pseudopotentials: {pps}'
                if lcao:
                    for i, s in enumerate(self.atom_species):
                        s.set_nao(pps[i], orbital_dir, d[4])
                    naos_comb = [list(s.nao) for s in self.atom_species]
                    naos_comb = list(itools.product(*naos_comb)) # for each AtomSpecies, one nao
                else:
                    naos_comb = [[None] * len(pps)] # only one possibility of naos: all is none
                assert all([len(naos) == len(pps) for naos in naos_comb]), f'naos should be a list of naos: {naos_comb}'
                # then for each pps-naos, inject AtomSpecies information into Cell.



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

        Si.set_pp("./download/pseudopotentials/", ["Si", "PBE", "NC", "sg15", "sr"])
        self.assertEqual(Si.pp, {'/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.2.upf', 
                                 '/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.0.upf', 
                                 '/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.1.upf'})

class TestCell(unittest.TestCase):

    def test_constructor(self):
        cell = Cell()
        self.assertIsNone(cell.a)
        self.assertIsNone(cell.b)
        self.assertIsNone(cell.c)
        self.assertIsNone(cell.alpha)
        self.assertIsNone(cell.beta)
        self.assertIsNone(cell.gamma)
        self.assertIsNone(cell.coords)
        self.assertIsNone(cell.species_map)
        self.assertIsNone(cell.periodic)
        self.assertIsNone(cell.sym_ks)

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
        cell = Cell()
        with open(fcif, "w") as f:
            f.write(cif)
        cell.standard_build(fcif)
        self.assertAlmostEqual(cell.a, 4.84451701, delta=1e-5)
        self.assertAlmostEqual(cell.b, 4.84451701, delta=1e-5)
        self.assertAlmostEqual(cell.c, 4.84451701, delta=1e-5)
        self.assertAlmostEqual(cell.alpha, 60.00000000, delta=1e-5)
        self.assertAlmostEqual(cell.beta, 60.00000000, delta=1e-5)
        self.assertAlmostEqual(cell.gamma, 60.00000000, delta=1e-5)
        self.assertEqual(cell.coords.shape, (3, 3))
        self.assertEqual(cell.species_map, [0, 0, 1])
        #print(cell.sym_ks)
        os.remove(fcif)

    def test_build(self):
        cell = Cell()
        cell.kernel_build([5, 5, 5, 90, 90, 90], [[0, 0, 0]], [0])
        self.assertEqual(cell.a, 5)
        self.assertEqual(cell.b, 5)
        self.assertEqual(cell.c, 5)
        self.assertEqual(cell.alpha, 90)
        self.assertEqual(cell.beta, 90)
        self.assertEqual(cell.gamma, 90)
        self.assertEqual(cell.coords.shape, (1, 3))
        self.assertEqual(cell.species_map, [0])
    
    def test_easy_build(self):
        cell = Cell()
        cell.ideal_build('Si_sc', 5.43)
        self.assertEqual(cell.a, 5.43)
        self.assertEqual(cell.b, 5.43)
        self.assertEqual(cell.c, 5.43)
        self.assertEqual(cell.alpha, 90)
        self.assertEqual(cell.beta, 90)
        self.assertEqual(cell.gamma, 90)
        self.assertEqual(cell.coords.shape, (1, 3))
        self.assertEqual(cell.species_map, [0])
        self.assertEqual(cell.symbols, ['Si'])

        with self.assertRaises(AssertionError):
            cell.ideal_build('SiO_diamond', 5.43)
        with self.assertRaises(AssertionError):
            cell.ideal_build('Si_x2y3', 5.43)
        with self.assertRaises(AssertionError):
            cell.ideal_build('Si_imagination', 5.43)
        cell.ideal_build('SiO_xy2', 5.43) # SiO2
        self.assertEqual(cell.a, 5.43)
        self.assertEqual(cell.b, 5.43)
        self.assertEqual(cell.c, 5.43)
        self.assertEqual(cell.alpha, 60)
        self.assertEqual(cell.beta, 60)
        self.assertEqual(cell.gamma, 60)
        self.assertEqual(cell.coords.shape, (3, 3))
        self.assertEqual(cell.species_map, [0, 1, 1])
        self.assertEqual(cell.symbols, ['Si', 'O'])

    def test_kspacing(self):
        cell = Cell()
        cell.kernel_build([4.22798145, 4.22798145, 4.22798145, 60, 60, 60], [[0, 0, 0]], [0])
        cell.kmeshgen(0.03*cell.lat0)
        self.assertEqual(cell.mpmesh_nks, [33, 33, 33])

class TestStructureManager(unittest.TestCase):

    def test_describe_structure_re(self):
        """regular expression in StructureManager.build function"""
        import re
        _re = r'([A-Z][a-z]?_[dimer|trimer|tetramer|sc|bcc|fcc|diamond])|([A-Z][a-z]?[A-Z][a-z]?_[xy2|xy3|x2y|x2y3|x2y5])'
        _match = re.match(_re, 'Si_dimer')
        self.assertIsNotNone(_match)
        _match = re.match(_re, 'O_sc')
        self.assertIsNotNone(_match)
        _match = re.match(_re, 'SiO_xy2')
        self.assertIsNotNone(_match)
        _match = re.match(_re, 'SiO_x2y5')
        self.assertIsNotNone(_match)
        _match = re.match(_re, 'SiO_x2y')
        self.assertIsNotNone(_match)
        _match = re.match(_re, 'CO_x2y3')
        self.assertIsNotNone(_match)
        _match = re.match(_re, 'Si_x2y3')
        self.assertIsNone(_match)
        _match = re.match(_re, 'SiO_x2y3')
        self.assertIsNotNone(_match) # indicating preset Si2O3 structure
        _match = re.match(_re, 'SiO_diamond')
        self.assertIsNone(_match) # SiO2 in diamond-like Bravis lattice is not defined in this way

if __name__ == "__main__":
    unittest.main()
                         