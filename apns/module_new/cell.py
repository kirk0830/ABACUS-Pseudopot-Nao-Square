"""
Refactor notes:
this is the file involved in refactor project of APNS for dealing with the more and more
complicated pseudopot-nao test cases.
"""

class AtomSpecies:
    """# AtomSpecies
    AtomSpecies is a prototype decoupling static AtomSpecies information with the dynamic
    structure information from practical simulation. 
    For APNS, it can be used to describe the properties of one atom species, such as the name,
    symbol, index, covalent radius, mass, magnetic moment, pseudopotentials, and numerical atomic
    orbitals. With one AtomSpecies, the Cartesian direct product can be made in unit of AtomSpecies,
    instead of the unit of Atom.

    Hope would be helpful for ABACUS `UnitCell` class refactor plan.

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
    
    def __init__(self, symbol: str, **kwargs) -> None:
        assert isinstance(symbol, str), f'symbol should be a string: {symbol}'
        from apns.module_database import database
        self.name = kwargs.get('name', symbol)
        self.fullname = database.PERIODIC_TABLE_TOFULLNAME[symbol]
        self.symbol = symbol
        self.index = database.PERIODIC_TABLE_TOINDEX[symbol]
        self.rcovalent = database.RCOVALENT[symbol]
        self.mass = database.ATOMIC_MASS[symbol]
        self.magmom = kwargs.get('magmom', 0.0)
    
    def set_pp(self, pseudo_dir: str, pptags: list):
        """set the pseudopotentials of this atom species by the tags"""
        from os.path import join as pjoin
        from apns.module_database.search import TagSearcher
        searcher = TagSearcher(pjoin(pseudo_dir, 'database.json'))
        pptags = [] if pptags is None else pptags
        self.pp = searcher(False, False, *(list(set(pptags + [self.symbol]))))
    
    def set_nao(self, pp: str, orbital_dir: str, naotags: list):
        """set the numerical atomic orbitals of this atom species by the tags"""
        import os
        from apns.module_database.search import TagSearcher
        searcher = TagSearcher(os.path.join(orbital_dir, 'database.json'))
        naotags = [] if naotags is None else naotags
        self.nao = searcher(False, False, *(naotags + [pp]))

    # allow visit like a dict
    def __getitem__(self, key):
        return getattr(self, key)

class Cell:
    """# Cell
    A prototype decoupling static AtomSpecies information with the dynamic structure information
    from practical simulation. Hope would be helpful for ABACUS `UnitCell` class refactor plan

    ### WHAT IS A CELL?
    
    A cell is such a object that know its structure and lattice parameters, but this is just the
    half part, it also knows its reciprocal lattice information, such as the high symmetry k-points.
    From a kspacing, it can also know how its Monkhorst-Pack mesh should be.

    ### WHAT IS NOT REALLY A CELL KNOWS?
    
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
    sym_ks, possible_kpath, mpmesh_nks = None, None, None
    # "points" information
    coords, vels, magmoms, mobs = None, None, None, None
    labels, kinds, labels_kinds_map = None, None, None
    
    def __init__(self) -> None:
        """create one cell instance with the scaling factors"""
        pass

    def build(self, identifier: str, config: str, scale: float, kspacing: float = -1.0):
        assert identifier in ["cif", "bravis", "molecule"], f'config should be cif, bravis or molecule: {config}'
        self.standard_build(config, scale, kspacing) if identifier == "cif" else self.ideal_build(config, scale, kspacing)
        self.vels = [[0.0] * 3] * len(self.coords)
        self.magmoms = [0.0] * len(self.coords)
        self.mobs = [[1] * 3] * len(self.coords)

    def kernel_build(self, lat: list, coords: list, labels_kinds_map: list):
        """build structure(s) from given species, lattice, and coordinates

        Args:
            lat (list): a, b, c and alpha, beta, gamma
            coords (list): coordinates of atoms, size = (n_atoms, 3)
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
        """
        import numpy as np
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = lat
        self.coords = np.array(coords)
        self.labels_kinds_map = labels_kinds_map

    def standard_build(self, fname: str, scale: float, kspacing: float = -1.0, **kwargs):
        """read structure file from external, and construct the structure.
        Currently only support CIF file.
        """
        self.periodic = True
        from pymatgen.io.cif import CifParser
        parser = CifParser(fname)
        structure = parser.parse_structures(primitive=True)[0]
        self.labels = [site.species.elements[0].symbol for site in structure]
        self.kinds = list(dict.fromkeys(self.labels)) # find unique set but keep the sequence
        self.labels_kinds_map = [self.kinds.index(lable) for lable in self.labels]
        self.coords = structure.cart_coords
        self.magmoms = [0] * len(self.coords)
        self.a, self.b, self.c = [i*scale**(1/3) for i in structure.lattice.abc]
        self.alpha, self.beta, self.gamma = structure.lattice.angles

        import seekpath
        cell = Cell._abc_angles_to_vec([self.a, self.b, self.c, self.alpha, self.beta, self.gamma], True)
        seekpath_result = seekpath.get_path((cell, self.coords.tolist(), self.labels_kinds_map), **kwargs)
        self.sym_ks = seekpath_result['point_coords']
        self.possible_kpath = seekpath_result['path']
        if kspacing > 0: self.kmeshgen(kspacing)
        else: self.mpmesh_nks = [1] * 3
        
    def ideal_build(self, config: str, characteristic: float, kspacing: float = -1.0):
        """build structures with simple shapes, such as ideal Bravis lattices,
        including SimpleCubic (SC), FaceCenteredCubic (FCC), BodyCenteredCubic (BCC),
        Diamond (diamond) and simple molecules such as dimer, trimer, tetramer, etc.
        """
        import re
        match_ = re.match(r'^([A-Z][a-z]?)_(dimer|trimer|tetramer|sc|bcc|fcc|diamond)$', config) or \
            re.match(r'^([A-Z][a-z]?[A-Z][a-z]?)_(xy2|xy3|x2y|x2y3|x2y5)$', config)
        assert match_, f'config not supported for ideal_build: {config}'
        kind, config = match_.groups()
        pure = config.lower() in ['sc', 'fcc', 'bcc', 'diamond', 'dimer', 'trimer', 'tetramer']
        periodic = config.lower() in ['sc', 'fcc', 'bcc', 'diamond', 'x2y', 'x2y3', 'x2y5', 'xy', 'xy2', 'xy3']
        composite = not pure and periodic
        assert (re.match(r'^[A-Z][a-z]?$', kind) and pure) or (re.match(r'^[A-Z][a-z]?[A-Z][a-z]?$', kind) and composite), \
            f'atomic symbol not supported for ideal_build: {kind}'
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
        abc_angles, labels_kinds_map, coords = Cell._structure(config, characteristic)
        self.kernel_build(abc_angles, coords, labels_kinds_map)
        self.kinds = re.findall(r'([A-Z][a-z]?)', kind)
        self.labels = [self.kinds[i] for i in self.labels_kinds_map]
        assert (len(self.kinds) == 2 and composite) or (len(self.kinds) == 1 and pure)
        if kspacing > 0: self.kmeshgen(kspacing)
        else: self.mpmesh_nks = [1] * 3

    def divide_subset(self, identifiers: list, by: str = "magmom"):
        """divide the labels_kinds_map into more piece, by the identifiers"""
        if identifiers is None: return
        if identifiers == []: return
        assert isinstance(identifiers, list), f'identifiers should be a list: {identifiers}'
        assert len(identifiers) == len(self.coords), f"everytime identifiers should be defined for each atom: {identifiers}"
        idx_sp_, unique_, kind_count_, labels_ = 0, [], [0]*len(self.kinds), []
        for idx_sp, idtf in zip(self.labels_kinds_map, identifiers):
            if (idx_sp, idtf) not in unique_:
                unique_.append((idx_sp, idtf))
                idx_sp_ += 1
                kind_count_[idx_sp] += 1
            labels_.append(self.kinds[idx_sp] + str(kind_count_[idx_sp]))
        self.labels = labels_
        if by == "magmom":
            self.magmoms = identifiers

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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
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
            kinds_map (list): mapping the index of coords to kinds, size = n_atoms
            coords (np.ndarray): coordinates of atoms
        """
        import numpy as np
        return [celldm, celldm, celldm, 60, 60, 60], \
            [0]*1 + [1]*3, \
            np.array([[0, 0, 0], 
                        [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]])

class StructureManager:
    """# StructureManager
    structure manager is for managing structures, including downloading, generating, etc.
    
    See unittest for its usage. Most of time it is easy:
    ```python
    sm = StructureManager() # but can also describe structures at this step
    ```
    The second step is to describe the structures, which is to provide the information of the
    structures, including the type of the structure, the parameters of the structure, the
    pseudopotentials and numerical atomic orbitals used in the structure.
    - The type of the structure can be 'cif', 'bravis', 'molecule'.
    - The second arg is a tuple composed of fname of cif file or the name of the structure,
    then if provided, the magnetic moment of EACH atom in the structure.
    - The third arg is a list of scalings, which will be used to scale the structure, will be
    useful when doing EOS and varying bond lengths.
    - Pseudopotentials and numerical atomic orbitals will be searched by tags specified.

    ```python
    desc = [(
                'cif',
                ('/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/structures/Si.cif',
                [0.00, 1.00]),
                [0.95, 0.97, 0.99, 1.01, 1.03], 
                {"Si": ['PBE', 'NC', 'sg15', 'sr']}, 
                {"Si": ['6au', '100Ry' , '2s2p1d']}
            ),
            (
                'bravis',
                ('SiO_xy2', None), 
                [0.95, 0.97, 0.99, 1.01, 1.03], 
                {"Si": ['PBE', 'NC', 'sg15', 'sr'], "O": ['PBE', 'US', 'GBRV', '1.4']}, 
                None
            ),
            (
                'molecule',
                ('Ba_dimer', None), 
                [1.0, 1.02, 1.04, 1.06, 1.08], 
                {"Ba": ['PBE', 'NC', 'sg15', 'sr']}, 
                {"Ba": ['6au', '100Ry' , '2s2p1d']}
            )]
    sm.describe(desc) # enough information has been provided to build the structures
    ```
    then to build the structures, according to the description
    ```python
    sm.build(pseudo_dir, orbital_dir)
    ```
    at this step, all structures have been built, along with their AtomSpecies, can be
    accessed by `sm.structures`

    Then export the structures to the abacus input file
    ```python
    result = sm.export(fmt='abacus')
    ```
    if the `save` is not None, then the result will be saved to the file, otherwise
    result will contain the content of files.

    `StructureManager` allows decoupling between the DFT calculation setting and the structure
    definition.
    """
    # descriptors of structures
    desc: list[tuple[str, str, list[float], list[str], list[str]|None]] = []
    # structures generated
    structures: list[tuple[Cell, list[AtomSpecies]]] = []

    def idtfr_gen(the_second: str) -> str:
        """generate the identifier of the structure from the second arg of the descriptor"""
        import os
        if os.path.exists(the_second) and os.path.isfile(the_second) and the_second.endswith(".cif"):
            return "cif", the_second
        import re
        match_ = re.match(r'^([A-Z][a-z]?)_(sc|bcc|fcc|diamond)$', the_second)
        if match_:
            return "bravis", *match_.groups()
        match_ = re.match(r'^([A-Z][a-z]?[A-Z][a-z]?)_(xy2|xy3|x2y|x2y3|x2y5)$', the_second)
        if match_:
            return "bravis", *match_.groups()
        match_ = re.match(r'^([A-Z][a-z]?)_(dimer|trimer|tetramer)$', the_second)
        if match_:
            return "molecule", *match_.groups()
        return None

    def atomsets_transpose(atomsets: list):
        # transpose atomsets from {atom: [[pptags...], [naotags...]]} to [{atom: [pptags...]}, {atom: [naotags...]}]
        atomsets = [[{key: value[i] for key, value in atomset.items()} for i in range(len(next(iter(atomset.values()))))]
                    for atomset in atomsets]
        # for each atomset, there are two dicts, if all values of dict is None, then set the dict as None
        atomsets = [[None if all([value is None for value in tags.values()]) else tags for tags in atomset] for atomset in atomsets]
        return atomsets

    def make_desc(atomsets: list, structures: list):
        """
        iterate on the `structures` section in input, in-on-shot download all structures
        needed, and return the descriptors of the structures for calling `describe()`.
        This function is called like:
        ```python
        atomsets = inp["atomsets"]
        structures = inp["structures"]
        desc = StructureManager.make_desc(atomsets, structures)
        ```
        """
        from apns.module_structure.api import download
        api_keys = {}
        formula = {}
        db_formula_isid_map = {}
        for is_, s in enumerate(structures):
            db = s.get("database", "mp") # default structure database is Materials Project
            for id_, d in enumerate(s["desc"]):
                if d[0] == "search":
                    api_keys[db] = s.get("api_key", "")
                    formula.setdefault(db, []).append(d[1])
                    db_formula_isid_map.setdefault((db, d[1]), []).append((is_, id_))
        # should add local check to avoid download the same structure!
        log = download(api_keys=api_keys, formula=formula) # will be nested dict of list of tuples, [db][formula][icif] = (fname, magmoms)
        atomsets = StructureManager.atomsets_transpose(atomsets)
        # first convert not cif structures to descriptors
        desc_ = [[("cif" if d[0] == "local" else StructureManager.idtfr_gen(d[1])[0], (d[1], None), d[2], 
             atomsets[s["atomset"]][0], atomsets[s["atomset"]][1]) for d in s["desc"] if d[0] != "search"] for s in structures] 
        # then add those downloaded structures to the descriptors
        for db, search_results in log.items(): # for each database, their formulas would also be dict, keys are formula, values are list of tuples
            for formula, result in search_results.items(): # result is a list of tuples, each tuple is (fname, magmoms)
                for is_, id_ in db_formula_isid_map[(db, formula)]:
                    desc_[is_].extend([("cif", (fname, magmoms), structures[is_]["desc"][id_][2], atomsets[s["atomset"]][0], atomsets[s["atomset"]][1])
                        for fname, magmoms in result])
        return desc_

    def __init__(self, desc: list = None) -> None:
        """structure can be imported at the initialization of the structure manager,
        but can also be imported/overwritten later if call the describe() method."""
        if desc is not None:
            self.describe(desc)
    
    """fucntions in this part is for importing structures without additional scalling, say
    not for EOS, not for varying bond lengths, ..."""
    def describe(self, desc, overwrite: bool = True):
        """
        For lazy build. First describe, then build explicitly and output instantly.
        : if the cell can be built with several parameters, then those parameters are
        representation of those exact structures. Calculate on those parameters means
        to vary and generate structures.
        
        ```python
        # descriptor
        #        identifier               str                         list[float]                    list[str] list[str]
        # (cif/bravis/molecule, (fcif/X_sc/X_dimer, magmoms), EOS_scalings/EOS_scalings/bond_lengths, pptags,   naotags)
        # pptags and naotags are associated.
        # once from pptags determines a list of pseudopotentials, 
        # for each pseudopotential, search numerical atomic orbital according to naotags.
        # example
        desc = [(
                    'cif',
                    ('/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/structures/Si.cif',
                    [0.00, 1.00]),
                    [0.95, 0.97, 0.99, 1.01, 1.03], 
                    {"Si": ['PBE', 'NC', 'sg15', 'sr']}, 
                    {"Si": ['6au', '100Ry' , '2s2p1d']}
                ),
                (
                    'bravis',
                    ('SiO_xy2', None), 
                    [0.95, 0.97, 0.99, 1.01, 1.03], 
                    {"Si": ['PBE', 'NC', 'sg15', 'sr'], "O": ['PBE', 'US', 'GBRV', '1.4']}, 
                    None
                ),
                (
                    'molecule',
                    ('Ba_dimer', None), 
                    [1.0, 1.02, 1.04, 1.06, 1.08], 
                    {"Ba": ['PBE', 'NC', 'sg15', 'sr']}, 
                    {"Ba": ['6au', '100Ry' , '2s2p1d']}
                )]
        sm = StructureManager()
        sm.describe(desc)
        ```
        """
        assert isinstance(desc, list), f'desc should be a list: {desc}'
        assert all([isinstance(d, tuple) for d in desc]), f'each element in desc should be a tuple: {desc}'
        assert all([len(d) == 5 for d in desc]), f'each element in desc should be a tuple of 5 elements: {desc}'
        assert all([isinstance(d[0], str) and d[0] in ['cif', 'bravis', 'molecule'] for d in desc]), \
            f'identifier should be cif, bravis or molecule: {desc}'
        assert all([isinstance(d[1], tuple) and len(d[1]) == 2 for d in desc])
        assert all([isinstance(d[1][0], str) and (isinstance(d[1][1], list) or d[1][1] is None) for d in desc])
        assert all([isinstance(d[1][0], str) and StructureManager.idtfr_gen(d[1][0])[0] in ["cif", "bravis", "molecule"]
                    for d in desc]), f'descriptor is neither a cif file nor a bravis/molecule type: {desc}'
        assert all([isinstance(d[2], list) and all([isinstance(x, float) for x in d[2]]) for d in desc]), \
            f'characteristic scaling should be list of floats: {desc}'
        assert all([isinstance(d[3], dict) and all([all([isinstance(tag, str) for tag in tags]) for tags in d[3].values()]) for d in desc]), \
            f'pptags should be dict of lists of strings: {desc}'
        assert all([\
            (isinstance(d[4], dict) and all([all([isinstance(tag, str) for tag in tags]) for tags in d[4].values])) or (d[4] is None) \
                for d in desc]), f'naotags should be dict of lists of strings: {desc} or None in basis_type pw calculation'
        self.desc = desc if overwrite else self.desc + desc
        """will add replace, pop and append functionalities later"""
    
    def build(self, pseudo_dir: str, orbital_dir: str):
        """expand all tests, this function is unique in APNS, but should not appear in
        ABACUS newly-refactored UnitCell module, becuase StructureManager is not exactly
        the StructureIterator or something, while it is less necessary to implement such
        a class or function."""
        import itertools as itools
        lcao = [d[4] for d in self.desc if d[4] is not None]
        assert len(lcao) == 0 or len(lcao) == len(self.desc), f'naotags should be all None or all not None: {lcao}'
        lcao = (len(lcao) != 0)
        # print(f"Build {len(self.desc)} structures for ABACUS-LCAO calculation: {lcao}")
        self.structures = [] # everytime refresh all structures?
        for idtfir, stru, scales, pptags, naotags in self.desc: # for each structure prototype
            proto = Cell()
            config, magmom = stru
            proto.build(idtfir, config, 1.0) # leaves d[3] the pptags and d[4] the naotags apart
            on_the_fly_species = [AtomSpecies(symbol) for symbol in proto.kinds] # it is the way how AtomSpecies bind with Cell
            pptags = [pptags[k] for k in proto.kinds]
            for i, s in enumerate(on_the_fly_species): s.set_pp(pseudo_dir, pptags[i])
            pps_scale_comb = [s.pp for s in on_the_fly_species] + [scales]
            assert pps_scale_comb is not None, f'pps_scale_comb should not be None: {pps_scale_comb}'
            pps_scale_comb = list(itools.product(*pps_scale_comb))
            for pps_scale in pps_scale_comb:
                pps, scale = pps_scale[:-1], pps_scale[-1]
                specific = Cell()
                specific.build(idtfir, config, scale)
                specific.divide_subset(magmom)
                for i, s in enumerate(on_the_fly_species): s.pp = pps[i] # inject the pseudopotentials
                naotags = None if not lcao else [naotags[k] for k in proto.kinds]
                for s in on_the_fly_species: s.set_nao(pps, orbital_dir, naotags)
                naos_comb = [s.nao for s in on_the_fly_species]
                assert naos_comb is not None, f'naos_comb should not be None: {naos_comb}'
                naos_comb = list(itools.product(*naos_comb)) if lcao else [None]
                for naos in naos_comb:
                    for i, s in enumerate(on_the_fly_species): s.nao = naos[i] if lcao else None # inject the numerical atomic orbitals
                    self.structures.append((specific, on_the_fly_species))

    def export(self, fmt: str = "abacus", save: str = None):
        """export the structures to files. If save is given, the the value given to save would be the prefix.
        Return list of tuple of three strings, if save is specified, then return the file name, otherwise the
        content of the files. The first string is the structure file, the second is the kline file, and the third
        is the kpt file. The content of the files are in the order of the structures."""
        import uuid
        result = []
        assert fmt == "abacus", f'fmt should be abacus: {fmt}'
        for s in self.structures:
            assert isinstance(s[0], Cell) and isinstance(s[1], list) and all([isinstance(x, AtomSpecies) for x in s[1]]), \
                f'structure should be a tuple of Cell and list of AtomSpecies.'
            stru = StructureManager.write_abacus_stru(s[0], s[1])
            kline, kpt = StructureManager.write_abacus_kpt(s[0])
            if save is not None:
                stamp = uuid.uuid4().hex
                fstru, fkline, fkpt = f"{stamp}.STRU", f"{stamp}.KLINE", f"{stamp}.KPT"
                with open(fstru, "w") as f: f.write(stru)
                with open(fkline, "w") as f: f.write(kline if kline is not None else "")
                with open(fkpt, "w") as f: f.write(kpt)
                result.append((fstru, fkline, fkpt))
            else:
                result.append((stru, kline, kpt))
        return result
    #######################
    #      utilities      #
    #######################
    def write_abacus_stru(cell: Cell, species: list[AtomSpecies]):
        """write a single ABACUS STRU file with Cell and list of AtomSpecies instances"""
        assert max(cell.labels_kinds_map) + 1 == len(species), f'cell.labels_kinds_map should be consistent with species. {cell.labels_kinds_map} vs {len(species)}'
        result = "ATOMIC_SPECIES\n"
        for label in list(dict.fromkeys(cell.labels)):
            s = species[cell.labels_kinds_map[cell.labels.index(label)]]
            fpp = s.pp.replace("\\", "/").split("/")[-1]
            result += f"{label:4s} {s.mass:8.4f} {fpp}\n"
        if all([s.nao is not None for s in species]):
            result += "\nNUMERICAL_ORBITAL\n"
            for label in list(dict.fromkeys(cell.labels)):
                s = species[cell.labels_kinds_map[cell.labels.index(label)]]
                result += s.nao.replace("\\", "/").split("/")[-1] + "\n"
        # else:
        #     print(f"Warning: not all species have numerical atomic orbitals, ignore if PW calculation")
        result += f"\nLATTICE_CONSTANT\n{cell.lat0:<20.10f}\n"
        result += f"\nLATTICE_VECTOR\n"
        latv = Cell._abc_angles_to_vec([cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma], True)
        for i in range(3):
            result += f"{latv[i][0]:<20.10f} {latv[i][1]:<20.10f} {latv[i][2]:<20.10f}\n"
        result += f"\nATOMIC_POSITIONS\n"
        result += "Cartesian_angstrom\n" if not cell.periodic else "Direct\n"
        for label in list(dict.fromkeys(cell.labels)):
            ias = [il for il, lbl in enumerate(cell.labels) if lbl == label] # indices of atoms with the same label
            result += f"{label}\n{cell.magmoms[ias[0]]:<4.2f}\n{len(ias)}\n"
            for j in ias:
                result += f"{cell.coords[j][0]:<20.10f} {cell.coords[j][1]:<20.10f} {cell.coords[j][2]:<20.10f} m "
                result += f"{cell.mobs[j][0]:<1d} {cell.mobs[j][1]:<1d} {cell.mobs[j][2]:<1d}\n"
        return result

    def write_abacus_kpt(cell: Cell):
        """write both KLINE and KPT file with Cell"""
        if cell.possible_kpath is not None and cell.sym_ks is not None:
            nks_each = 20
            nks_tot = len(cell.possible_kpath) * nks_each + 1
            kline = f"KPOINTS\n{nks_tot}\nLine\n"
            # merge all consecutive segments connecting high symmetrical kpoints
            segs = []
            for seg in cell.possible_kpath:
                b, e = seg
                if not segs or b != segs[-1][-1]:
                    segs.append([b, e])
                else:
                    segs[-1].append(e)

            for seg in segs:
                for ik, k in enumerate(seg):
                    nks = nks_each if ik != len(seg) - 1 else 1
                    kline += f"{cell.sym_ks[k][0]:13.10f} {cell.sym_ks[k][1]:13.10f} {cell.sym_ks[k][2]:13.10f} {nks}\n"
        else:
            kline = None
        kpt = f"KPOINTS\n0\nGAMMA\n"
        kpt += f"{cell.mpmesh_nks[0]:<d} {cell.mpmesh_nks[1]:<d} {cell.mpmesh_nks[2]:<d} 0 0 0\n"
        return kline, kpt

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
        self.assertIsNone(cell.labels_kinds_map)
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
        cell.standard_build(fcif, 1.0)
        self.assertAlmostEqual(cell.a, 4.84451701, delta=1e-5)
        self.assertAlmostEqual(cell.b, 4.84451701, delta=1e-5)
        self.assertAlmostEqual(cell.c, 4.84451701, delta=1e-5)
        self.assertAlmostEqual(cell.alpha, 60.00000000, delta=1e-5)
        self.assertAlmostEqual(cell.beta, 60.00000000, delta=1e-5)
        self.assertAlmostEqual(cell.gamma, 60.00000000, delta=1e-5)
        self.assertEqual(cell.coords.shape, (3, 3))
        self.assertEqual(cell.labels_kinds_map, [0, 0, 1])
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
        self.assertEqual(cell.labels_kinds_map, [0])
    
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
        self.assertEqual(cell.labels_kinds_map, [0])
        self.assertEqual(cell.kinds, ['Si'])

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
        self.assertEqual(cell.labels_kinds_map, [0, 1, 1])
        self.assertEqual(cell.kinds, ['Si', 'O'])

    def test_kspacing(self):
        cell = Cell()
        cell.kernel_build([4.22798145, 4.22798145, 4.22798145, 60, 60, 60], [[0, 0, 0]], [0])
        cell.kmeshgen(0.03*cell.lat0)
        self.assertEqual(cell.mpmesh_nks, [33, 33, 33])

    def test_divide_subset(self):
        import numpy as np
        cell = Cell()
        cell.labels = ["Si", "Si", "O", "O", "O", "Na", "Ca", "Ca"]
        cell.kinds = ["Si", "O", "Na", "Ca"]
        cell.labels_kinds_map = [0, 0, 1, 1, 1, 2, 3, 3]
        cell.coords = np.random.random((8, 3)).tolist()
        magmom = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0]
        cell.divide_subset(identifiers=magmom)
        self.assertEqual(cell.labels_kinds_map, [0, 0, 1, 1, 1, 2, 3, 3])
        self.assertEqual(cell.labels, ["Si1", "Si1", "O1", "O1", "O1", "Na1", "Ca1", "Ca1"])
        # recover
        cell.labels = ["Si", "Si", "O", "O", "O", "Na", "Ca", "Ca"]
        cell.kinds = ["Si", "O", "Na", "Ca"]
        cell.labels_kinds_map = [0, 0, 1, 1, 1, 2, 3, 3]
        # AFM on Si
        magmom = [1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0]
        cell.divide_subset(identifiers=magmom)
        self.assertEqual(cell.labels_kinds_map, [0, 0, 1, 1, 1, 2, 3, 3])
        self.assertEqual(cell.labels, ["Si1", "Si2", "O1", "O1", "O1", "Na1", "Ca1", "Ca1"])
        # recover
        cell.labels = ["Si", "Si", "O", "O", "O", "Na", "Ca", "Ca"]
        cell.kinds = ["Si", "O", "Na", "Ca"]
        cell.labels_kinds_map = [0, 0, 1, 1, 1, 2, 3, 3]
        # AFM on Si, FM on O, FM on Na, AFM on Ca
        magmom = [1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0]
        cell.divide_subset(identifiers=magmom)
        self.assertEqual(cell.labels_kinds_map, [0, 0, 1, 1, 1, 2, 3, 3])
        self.assertEqual(cell.labels, ["Si1", "Si2", "O1", "O1", "O1", "Na1", "Ca1", "Ca2"])

class TestStructureManager(unittest.TestCase):

    def test_atomsets_transpose(self):
        atomsets = [
            {"H": [["pptagH1", "pptagH2"], None],
             "O": [["pptagO1", "pptagO2"], None]},
        ]
        result = StructureManager.atomsets_transpose(atomsets)
        self.assertEqual(result, [[{'H': ['pptagH1', 'pptagH2'], 'O': ['pptagO1', 'pptagO2']}, None]])

    def est_make_desc(self):
        atomsets = [
            {"H": [["pptagH1", "pptagH2"], None],
             "O": [["pptagO1", "pptagO2"], None]},
        ]
        structures = [
            {"database": "mp",
            "atomset": 0, 
            "api_key": "wV1HUdmgESPVgSmQj5cc8WvttCO8NTHp", "desc": [
            ["search", "BaTiO3", [1.0]]
            ]}
        ]
        result = StructureManager.make_desc(
            atomsets, structures
        )
        ref = [[('cif', ('apns_cache/mp-5777.cif', [0.0, 0.0, 0.0, 0.0, 0.0]), [1.0], 
                 {'H': ['pptagH1', 'pptagH2'], 'O': ["pptagO1", "pptagO2"]}, None)]]
        self.assertEqual(result, ref)

    def test_describe_re(self):
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

    def test_idtfr_gen(self):
        self.assertEqual(StructureManager.idtfr_gen('Si_dimer'), ('molecule', 'Si', 'dimer'))
        self.assertEqual(StructureManager.idtfr_gen('SiO_xy2'), ('bravis', 'SiO', 'xy2'))
        self.assertEqual(StructureManager.idtfr_gen('SiO_x2y3'), ('bravis', 'SiO', 'x2y3'))
        self.assertIsNone(StructureManager.idtfr_gen('SiO_diamond'))
        fcif = "Si.cif"
        with open(fcif, "w") as f:
            f.write("fake cif")
        self.assertEqual(StructureManager.idtfr_gen(fcif), ('cif', fcif))
        import os
        os.remove(fcif)

    def test_build(self):
        import os
        sm = StructureManager()
        # test Bravis lattice
        sm.describe([('bravis', ('Si_diamond', None), [0.95, 0.97, 0.99, 1.01, 1.03], {"Si": ['PBE', 'NC', 'sg15', 'sr']}, None)])
        sm.build("./download/pseudopotentials/", "./download/numerical_orbitals/")
        self.assertEqual(len(sm.structures), 5 * 3) # 3 is the number of pseudopotentials searched
        result = sm.export(fmt="abacus", save="unittest")
        for f in result:
            for ff in f:
                self.assertTrue(os.path.exists(ff))
                os.remove(ff)
        # test molecule
        sm.describe([('molecule', ('Si_dimer', None), [0.95, 0.97, 0.99, 1.01, 1.03], {"Si": ['PBE', 'NC', 'sg15', 'sr']}, None)])
        sm.build("./download/pseudopotentials/", "./download/numerical_orbitals/")
        self.assertEqual(len(sm.structures), 5 * 3)
        result = sm.export(fmt="abacus", save="unittest")
        for f in result:
            for ff in f:
                self.assertTrue(os.path.exists(ff))
                os.remove(ff)
        # test cif with magmom
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
        with open(fcif, "w") as f:
            f.write(cif)
        sm.describe([('cif', (fcif, [1.0, -1.0, 0.0]), [0.95, 0.97, 0.99, 1.01, 1.03], 
                    {"Ac": ['PBE', 'NC', 'sr'], "O": ['PBE', 'GBRV', '1.5']}, None)])
        sm.build("./download/pseudopotentials/", "./download/numerical_orbitals/")
        self.assertEqual(len(sm.structures), 2 * 1 * 5) # there are 2 available pseudopotentials for Ac, 1 for O and 5 scalings
        os.remove(fcif)
        result = sm.export(fmt="abacus", save="unittest")
        for f in result:
            for ff in f:
                self.assertTrue(os.path.exists(ff))
                os.remove(ff)

if __name__ == "__main__":
    unittest.main(exit=False)
