
class AtomSpeciesGeneartor:
    name, fullname, symbol, index, rcovalent, mass, magmom, u_minus_j = \
        None, None, None, None, None, None, None, None
    pptags, naotags = None, None

    def __init__(self, symbol: str, pptags: list, naotags: list = None, **kwargs) -> None:
        assert isinstance(symbol, str), f'symbol should be a string: {symbol}'
        from apns.module_database import database
        self.name = kwargs.get('name', symbol)
        assert isinstance(self.name, str), f'name should be a string: {self.name}'
        self.fullname = database.PERIODIC_TABLE_TOFULLNAME[symbol]
        self.symbol = symbol
        assert self.symbol in database.PERIODIC_TABLE_TOINDEX, f'symbol not supported: {self.symbol}'
        self.index = database.PERIODIC_TABLE_TOINDEX[symbol]
        self.rcovalent = database.RCOVALENT[symbol]
        self.mass = database.ATOMIC_MASS[symbol]
        self.magmom = kwargs.get('magmom', 0.0)
        assert isinstance(pptags, list), f'pseudopotential tags should be a list: {pptags}'
        assert all([isinstance(tag, str) for tag in pptags]), f'pseudopotential tags should be a list of strings: {pptags}'
        self.pptags = pptags
        assert isinstance(naotags, list) or naotags is None, f'nao tags should be a list or None: {naotags}'
        assert naotags is None or all([isinstance(tag, str) for tag in naotags]), f'nao tags should be a list of strings: {naotags}'
        self.naotags = naotags
        self.u_minus_j = kwargs.get('u_minus_j', None)
        assert self.u_minus_j is None or isinstance(self.u_minus_j, list), f'u_minus_j should be a list of floats: {self.u_minus_j}'
        assert self.u_minus_j is None or all([isinstance(u, float) for u in self.u_minus_j]), f'u_minus_j should be a list of floats: {self.u_minus_j}'

    def __call__(self, pseudo_dir, orbital_dir = None):
        """iteratively create AtomSpecies instances"""
        import itertools as it
        from os.path import join as pjoin
        from apns.module_database.search import TagSearcher
        searcher = TagSearcher(pjoin(pseudo_dir, 'database.json'))
        pps = searcher(False, False, *(list(set(self.pptags + [self.symbol]))))
        
        for pp in pps:
            if orbital_dir is not None:
                searcher = TagSearcher(pjoin(orbital_dir, 'database.json'))
                naos = searcher(False, False, *(list(set(self.naotags + [pp]))))
                for pp, nao in it.product(pps, naos):
                    yield AtomSpecies(name=self.name, fullname=self.fullname, symbol=self.symbol, index=self.index,
                                      rcovalent=self.rcovalent, mass=self.mass, magmom=self.magmom, pp=pp, nao=nao)
            else:
                yield AtomSpecies(name=self.name, fullname=self.fullname, symbol=self.symbol, index=self.index,
                                  rcovalent=self.rcovalent, mass=self.mass, magmom=self.magmom, pp=pp, nao=None)

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
    
    def __init__(self, name, fullname, symbol, index, rcovalent, mass, magmom, pp, nao) -> None:
        self.name = name
        self.fullname = fullname
        self.symbol = symbol
        self.index = index
        self.rcovalent = rcovalent
        self.mass = mass
        self.magmom = magmom
        self.pp = pp
        self.nao = nao

class CellGenerator:

    # basic
    identifier, config = None, None
    # reciprocal space lattice information
    kspacing = None
    # "points" information
    magmoms = None
    # generator specific
    scales = None

    def __init__(self, identifier: str, config: str, scales: list, **kwargs) -> None:
        """create one cell instance with the scaling factors"""
        assert identifier in ["cif", "bravis", "molecule"], f'config should be cif, bravis or molecule: {config}'
        self.identifier = identifier
        import re
        import os
        _match = re.match(r'^([A-Z][a-z]?)_(dimer|trimer|tetramer|sc|bcc|fcc|diamond)$', config) or \
            re.match(r'^([A-Z][a-z]?[A-Z][a-z]?)_(xy2|xy3|x2y|x2y3|x2y5)$', config)
        assert _match or os.path.exists(config), f'config should be a file or a simple structure: {config}'
        self.config = config
        assert isinstance(scales, list), f'scales should be a list of floats: {scales}'
        assert all([isinstance(scale, float) for scale in scales]), f'scales should be a list of floats: {scales}'
        self.scales = scales
        self.kspacing = kwargs.get('kspacing', -1.0)
        assert isinstance(self.kspacing, float), f'kspacing should be a float: {self.kspacing}'
        self.magmoms = kwargs.get('magmoms', None)
        assert self.magmoms is None or isinstance(self.magmoms, list), f'magmoms should be a list of floats: {self.magmoms}'

    def __call__(self):
        """iteratively create Cell instances"""
        for scale in self.scales:
            params = self.build(self.identifier, self.config, scale, self.kspacing)
            params["labels"] = CellGenerator.divide_subset(
                params["kinds"], params["labels"], self.magmoms, params["coords"])\
                    if self.magmoms is not None else params["labels"]
            cell = Cell(**params)
            yield cell

    def build(self, identifier: str, config: str, scale: float, kspacing: float) -> dict:
        assert identifier in ["cif", "bravis", "molecule"], f'config should be cif, bravis or molecule: {config}'
        build_func = CellGenerator.standard_build if identifier == "cif" else CellGenerator.ideal_build
        build_result = build_func(config, scale, kspacing)
        vels = [[0.0] * 3] * len(build_result["coords"])
        magmoms = [0.0] * len(build_result["coords"]) if self.magmoms is None else self.magmoms
        mobs = [[1] * 3] * len(build_result["coords"])

        keys = ["vels", "magmoms", "mobs"]
        vals = [vels, magmoms, mobs]
        build_result.update(dict(zip(keys, vals)))
        return build_result

    def standard_build(fname: str, scale: float, kspacing: float = -1.0):
        """read structure file from external, and construct the structure.
        Currently only support CIF file.
        """

        from pymatgen.io.cif import CifParser
        parser = CifParser(fname)
        structure = parser.parse_structures(primitive=True)[0]
        labels = [site.species.elements[0].symbol for site in structure]
        kinds = list(dict.fromkeys(labels)) # find unique set but keep the sequence
        labels_kinds_map = [kinds.index(lable) for lable in labels]
        coords = structure.cart_coords
        magmoms = [0] * len(coords)
        a, b, c = [i*scale**(1/3) for i in structure.lattice.abc]
        alpha, beta, gamma = structure.lattice.angles

        import seekpath
        cell = CellGenerator._abc_angles_to_vec([a, b, c, alpha, beta, gamma], True)
        seekpath_result = seekpath.get_path((cell, coords.tolist(), labels_kinds_map))
        sym_ks = seekpath_result['point_coords']
        possible_kpath = seekpath_result['path']
        mpmesh_nks = CellGenerator.kmeshgen(a, b, c, alpha, beta, gamma, kspacing) if kspacing > 0 else [1] * 3

        keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma', 'labels', 'kinds', 'labels_kinds_map', 'coords', 'magmoms', 'sym_ks', 'possible_kpath', 'mpmesh_nks']
        vals = [a, b, c, alpha, beta, gamma, labels, kinds, labels_kinds_map, coords, magmoms, sym_ks, possible_kpath, mpmesh_nks]
        return dict(zip(keys, vals))
        
    def ideal_build(config: str, characteristic: float, kspacing: float = -1.0):
        """build structures with simple shapes, such as ideal Bravis lattices,
        including SimpleCubic (SC), FaceCenteredCubic (FCC), BodyCenteredCubic (BCC),
        Diamond (diamond) and simple molecules such as dimer, trimer, tetramer, etc.
        """
        import re
        import numpy as np
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
        # True, True: ideal sc, fcc, bcc, diamond
        # True, False: dimer, trimer, tetramer
        # False, True: X2Y, X2Y3, X2Y5, XY, XY2, XY3
        # False, False: not supported
        assert pure or periodic, f'config not supported for ideal_build: {config}'
        characteristic = CellGenerator._vcell_to_celldm(CellGenerator.bravis_angles(config), -characteristic)\
            if characteristic < 0 and periodic else characteristic
        assert characteristic > 0, f'characteristic should be positive: {characteristic}'
        abc_angles, labels_kinds_map, coords = CellGenerator._structure(config, characteristic)
        a, b, c, alpha, beta, gamma = abc_angles
        coords = np.array(coords)
        kinds = re.findall(r'([A-Z][a-z]?)', kind)
        labels = [kinds[i] for i in labels_kinds_map]
        assert (len(kinds) == 2 and composite) or (len(kinds) == 1 and pure)
        mpmesh_nks = CellGenerator.kmeshgen(a, b, c, alpha, beta, gamma, kspacing) if kspacing > 0 else [1] * 3

        keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma', 'labels', 'kinds', 'labels_kinds_map', 'coords', 'mpmesh_nks']
        vals = [a, b, c, alpha, beta, gamma, labels, kinds, labels_kinds_map, coords, mpmesh_nks]
        return dict(zip(keys, vals))

    def divide_subset(fullset: list, dividee: list, divider: list, fullset_dividee_map: list):
        """divide the labels into more piece"""
        if divider is None: return
        if divider == []: return
        assert isinstance(divider, list), f'identifiers should be a list: {divider}'
        assert len(divider) == len(fullset), f"everytime identifiers should be defined for each atom: {divider}"
        idx_sp_, unique_, kind_count_, labels_ = 0, [], [0]*len(dividee), []
        for idx_sp, idtf in zip(fullset_dividee_map, divider):
            if (idx_sp, idtf) not in unique_:
                unique_.append((idx_sp, idtf))
                idx_sp_ += 1
                kind_count_[idx_sp] += 1
            labels_.append(dividee[idx_sp] + str(kind_count_[idx_sp]))
        return labels_

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

    def kmeshgen(a, b, c, alpha, beta, gamma, kspacing: float|list[float]):
        """get the Monkhorst-Pack mesh"""
        import numpy as np
        kspacing = kspacing if isinstance(kspacing, list) else [kspacing] * 3
        vecs = np.array(CellGenerator._abc_angles_to_vec([a, b, c, alpha, beta, gamma], True))
        recvecs = np.linalg.inv(vecs).T
        recvecs = 2*np.pi * recvecs
        norms = np.linalg.norm(recvecs, axis=1).tolist()
        assert len(norms) == len(kspacing), f'kspacing should be a list of 3 floats: {kspacing}'
        norms = [int(norm / kspac) for norm, kspac in zip(norms, kspacing)]
        return list(map(lambda x: max(1, x + 1), norms))

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
            return CellGenerator._simple_cubic(characteristic)
        elif config == 'fcc':
            return CellGenerator._face_centered_cubic(characteristic)
        elif config == 'bcc':
            return CellGenerator._body_centered_cubic(characteristic)
        elif config == 'diamond':
            return CellGenerator._diamond(characteristic)
        elif config == 'dimer':
            return CellGenerator._dimer(characteristic)
        elif config == 'trimer':
            return CellGenerator._trimer(characteristic)
        elif config == 'tetramer':
            return CellGenerator._tetramer(characteristic)
        elif config == 'X2Y':
            return CellGenerator._X2Y(characteristic)
        elif config == 'X2Y3':
            return CellGenerator._X2Y3(characteristic)
        elif config == 'X2Y5':
            return CellGenerator._X2Y5(characteristic)
        elif config == 'XY':
            return CellGenerator._XY(characteristic)
        elif config == 'XY2':
            return CellGenerator._XY2(characteristic)
        elif config == 'XY3':
            return CellGenerator._XY3(characteristic)
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
    # real space lattice information
    a, b, c, alpha, beta, gamma = None, None, None, None, None, None
    lat0 = 1.8897259886
    # reciprocal space lattice information
    sym_ks, possible_kpath, mpmesh_nks = None, None, None
    # "points" information
    coords, vels, magmoms, mobs = None, None, None, None
    labels, kinds, labels_kinds_map = None, None, None

    def __init__(self, **kwargs) -> None:
        self.a = kwargs.get('a', None)
        self.b = kwargs.get('b', None)
        self.c = kwargs.get('c', None)
        self.alpha = kwargs.get('alpha', None)
        self.beta = kwargs.get('beta', None)
        self.gamma = kwargs.get('gamma', None)
        self.sym_ks = kwargs.get('sym_ks', None)
        self.possible_kpath = kwargs.get('possible_kpath', None)
        self.mpmesh_nks = kwargs.get('mpmesh_nks', None)
        self.coords = kwargs.get('coords', None)
        self.vels = kwargs.get('vels', None)
        self.magmoms = kwargs.get('magmoms', None)
        self.mobs = kwargs.get('mobs', None)
        self.labels = kwargs.get('labels', None)
        self.kinds = kwargs.get('kinds', None)
        self.labels_kinds_map = kwargs.get('labels_kinds_map', None)


import unittest
class TestAtomSpeciesGenerator(unittest.TestCase):

    def test_constructor(self):
        asg = AtomSpeciesGeneartor('Si', None)
        self.assertEqual(asg.name, 'Si')
        self.assertEqual(asg.fullname, 'Silicon')
        self.assertEqual(asg.symbol, 'Si')
        self.assertEqual(asg.index, 14)
        self.assertEqual(asg.rcovalent, 1.11)
        self.assertEqual(asg.mass, 28.0855)
    
    def test_oncall(self):
        import os
        import json
        # prepare a fake database
        database = {"the_file_name_with_path": ["tag1", "tag2", "Si"]}
        fdatabase = "./database.json"
        with open(fdatabase, "w") as f:
            json.dump(database, f)
        asg = AtomSpeciesGeneartor('Si', ["tag1", "tag2"])
        times = 0
        for as_ in asg("./"):
            times += 1
            self.assertEqual(as_.pp, {"the_file_name_with_path"})
        self.assertEqual(times, 1)
        os.remove(fdatabase)

        asg = AtomSpeciesGeneartor('Si', ["Si", "PBE", "NC", "sg15", "sr"])
        times = 0
        ref = {'/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.2.upf', 
               '/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.0.upf', 
               '/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.1.upf'}
        for as_ in asg("./download/pseudopotentials/"):
            self.assertEqual(as_.pp, ref)
            times += 1
        self.assertEqual(times, 3)

class TestCellGeneartor(unittest.TestCase):

    def test_constructor(self):
        with self.assertRaises(AssertionError):
            CellGenerator('cif', 'Si', 1.0) # raise error because scales must be list

        cg = CellGenerator('cif', 'Si', [1.0])
        self.assertEqual(cg.identifier, 'cif')
        self.assertEqual(cg.config, 'Si')
        self.assertEqual(cg.scales, [1.0])

    def test_kspacing(self):
        result = CellGenerator.kmeshgen(4.22798145, 4.22798145, 4.22798145, 60, 60, 60, 0.03*1.889725989)
        self.assertEqual(result, [33, 33, 33])

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
        cg = CellGenerator('cif', fcif, [1.0])
        with open(fcif, "w") as f:
            f.write(cif)
        times = 0
        for cell in cg():
            self.assertEqual(cell.a, 4.84451701)
            self.assertEqual(cell.b, 4.84451701)
            self.assertEqual(cell.c, 4.84451701)
            self.assertEqual(cell.alpha, 60.0)
            self.assertEqual(cell.beta, 60.0)
            self.assertEqual(cell.gamma, 60.0)
            self.assertEqual(cell.labels, ["Ac", "Ac", "O"])
            self.assertEqual(cell.kinds, ["Ac", "O"])
            self.assertEqual(cell.labels_kinds_map, [0, 0, 1])
            self.assertEqual(cell.coords, [[0, 0, 0], [0.75, 0.75, 0.75], [0.25, 0.25, 0.25]])
            times += 1
        self.assertEqual(times, 1)

        cg = CellGenerator('cif', fcif, [1.0, 1.02, 1.04, 1.06, 1.08])
        times = 0
        for cell in cg():
            self.assertEqual(cell.a, 4.84451701 * cg.scales[times]**(1/3))
            self.assertEqual(cell.b, 4.84451701 * cg.scales[times]**(1/3))
            self.assertEqual(cell.c, 4.84451701 * cg.scales[times]**(1/3))
            self.assertEqual(cell.alpha, 60.0)
            self.assertEqual(cell.beta, 60.0)
            self.assertEqual(cell.gamma, 60.0)
            self.assertEqual(cell.labels, ["Ac", "Ac", "O"])
            self.assertEqual(cell.kinds, ["Ac", "O"])
            self.assertEqual(cell.labels_kinds_map, [0, 0, 1])
            self.assertEqual(cell.coords, [[0, 0, 0], [0.75, 0.75, 0.75], [0.25, 0.25, 0.25]])
            times += 1
        self.assertEqual(times, 5)

        # test with magmom
        cg = CellGenerator('cif', fcif, [1.0], magmoms = [1.0, -1.0, 0.0])
        times = 0
        for cell in cg():
            self.assertEqual(cell.a, 4.84451701)
            self.assertEqual(cell.b, 4.84451701)
            self.assertEqual(cell.c, 4.84451701)
            self.assertEqual(cell.alpha, 60.0)
            self.assertEqual(cell.beta, 60.0)
            self.assertEqual(cell.gamma, 60.0)
            self.assertEqual(cell.labels, ["Ac1", "Ac2", "O1"])
            self.assertEqual(cell.kinds, ["Ac", "O"])
            self.assertEqual(cell.labels_kinds_map, [0, 0, 1])
            self.assertEqual(cell.coords, [[0, 0, 0], [0.75, 0.75, 0.75], [0.25, 0.25, 0.25]])
            self.assertEqual(cell.magmoms, [1.0, -1.0, 0.0])
            times += 1

        os.remove(fcif)

    def test_ideal_build(self):
        cg = CellGenerator('bravis', 'Si_sc', [1.0]) # for bravis, the scale ought to be scaling factor but presently
        # is directly the volume of the cell
        times = 0
        for cell in cg():
            self.assertEqual(cell.a, 1.0)
            self.assertEqual(cell.b, 1.0)
            self.assertEqual(cell.c, 1.0)
            self.assertEqual(cell.alpha, 90.0)
            self.assertEqual(cell.beta, 90.0)
            self.assertEqual(cell.gamma, 90.0)
            self.assertEqual(cell.labels, ["Si"])
            self.assertEqual(cell.kinds, ["Si"])
            self.assertEqual(cell.labels_kinds_map, [0])
            self.assertEqual(cell.coords, [[0, 0, 0]])
            times += 1
        self.assertEqual(times, 1)

        cg = CellGenerator('bravis', 'Si_sc', [1.0, 1.02, 1.04, 1.06, 1.08])
        times = 0
        for cell in cg():
            self.assertEqual(cell.a, 1.0 * cg.scales[times]**(1/3))
            self.assertEqual(cell.b, 1.0 * cg.scales[times]**(1/3))
            self.assertEqual(cell.c, 1.0 * cg.scales[times]**(1/3))
            self.assertEqual(cell.alpha, 90.0)
            self.assertEqual(cell.beta, 90.0)
            self.assertEqual(cell.gamma, 90.0)
            self.assertEqual(cell.labels, ["Si"])
            self.assertEqual(cell.kinds, ["Si"])
            self.assertEqual(cell.labels_kinds_map, [0])
            self.assertEqual(cell.coords, [[0, 0, 0]])
            times += 1
        self.assertEqual(times, 5)

        cg = CellGenerator('molecule', 'Ba_dimer', [1.0])
        times = 0
        for cell in cg():
            self.assertEqual(cell.a, 20.0)
            self.assertEqual(cell.b, 20.0)
            self.assertEqual(cell.c, 20.0)
            self.assertEqual(cell.alpha, 90.0)
            self.assertEqual(cell.beta, 90.0)
            self.assertEqual(cell.gamma, 90.0)
            self.assertEqual(cell.labels, ["Ba", "Ba"])
            self.assertEqual(cell.kinds, ["Ba"])
            self.assertEqual(cell.labels_kinds_map, [0, 0])
            self.assertEqual(cell.coords, [[0, 0, 0], [1.0, 0, 0]])
            times += 1
        self.assertEqual(times, 1)

        cg = CellGenerator('molecule', 'Ba_dimer', [1.0, 1.02, 1.04, 1.06, 1.08])
        times = 0
        for cell in cg():
            self.assertEqual(cell.a, 20.0)
            self.assertEqual(cell.b, 20.0)
            self.assertEqual(cell.c, 20.0)
            self.assertEqual(cell.alpha, 90.0)
            self.assertEqual(cell.beta, 90.0)
            self.assertEqual(cell.gamma, 90.0)
            self.assertEqual(cell.labels, ["Ba", "Ba"])
            self.assertEqual(cell.kinds, ["Ba"])
            self.assertEqual(cell.labels_kinds_map, [0, 0])
            self.assertEqual(cell.coords, [[0, 0, 0], [cg.scales[times], 0, 0]])
            times += 1
        self.assertEqual(times, 5)


if __name__ == "__main__":
    unittest.main(exit=False)
