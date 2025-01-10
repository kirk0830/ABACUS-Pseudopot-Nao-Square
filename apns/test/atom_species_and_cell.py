class AtomSpeciesGeneartor:
    """An AtomSpeciesGenerator is for generating AtomSpecies instances based on the combination of
    pseudopotentials and their corresponding numerical atomic orbitals (if any, for LCAO calculation)
    . There is also a interface to ABACUS-ORBGEN package to generate the jY basis (truncated 
    spherical Bessel function multiplied by spherical harmonics) for the LCAO calculation. The full
    implementation will be in ABACUS-ORBGEN."""
    symbol = None
    pseudo_dir, orbital_dir = None, None
    pptags, naotags = None, None
    ecut_db = None
    def __init__(self, 
                 symbol: str, 
                 pseudo_dir: str, 
                 pptags: list, 
                 orbital_dir: str = None, 
                 naotags: list = None, 
                 **kwargs) -> None:
        """create one AtomSpeciesGenerator instance
        
        Parameters
        ----------
        symbol: str
            the symbol of the atom species
        pseudo_dir: str
            the path to the pseudopotential directory where the pseudopotentials are stored.
            For pseudopotential file management, the database.json file should present. For
            more information about the database.json file, please see example provided in
            the folder `/download/pseudopotentials/`
        pptags: list
            the tags for searching the pseudopotential file. Use tags to define/search all
            available pseudopotential for present element, this is supported by a well-indexed
            pseudopotential system.
        orbital_dir: str
            the path to the numerical atomic orbital directory where the numerical atomic orbitals
            are stored. Orbitals are managed in the same way as pseudopotentials.
        naotags: list
            similar to pptags

        Raises
        ------
        AssertionError
            if the symbol is not a string
            if the pseudo_dir is not a string
            if the pptags is not a list
            if the orbital_dir is not a string or None
            if the naotags is not a list or None  
        """
        assert isinstance(symbol, str), f'symbol should be a string: {symbol}'
        self.symbol = symbol
        assert isinstance(pseudo_dir, str), f'pseudo_dir should be a string: {pseudo_dir}'
        self.pseudo_dir = pseudo_dir
        assert isinstance(pptags, list), f'pseudopotential tags should be a list: {pptags}'
        assert all([isinstance(tag, str) for tag in pptags]), f'pseudopotential tags should be a list of strings: {pptags}'
        self.pptags = pptags
        assert isinstance(orbital_dir, str) or orbital_dir is None, f'orbital_dir should be a string or None: {orbital_dir}'
        self.orbital_dir = orbital_dir
        assert isinstance(naotags, list) or naotags is None, f'nao tags should be a list or None: {naotags}'
        assert naotags is None or all([isinstance(tag, str) for tag in naotags]), f'nao tags should be a list of strings: {naotags}'
        self.naotags = naotags
        print(f"""AtomSpeciesGenerator setup
Symbol: {self.symbol}
Pseudopotential directory: {self.pseudo_dir}
Pseudopotential tags: {self.pptags}
Orbital directory: {self.orbital_dir}
Numerical atomic orbital tags: {self.naotags}
""")

    def connect_ecut_db(self, fname: str):
        """connect with the ecut database"""
        import json
        import os
        if not os.path.exists(fname): return
        with open(fname, "r") as f:
            self.ecut_db = json.load(f)
        
    def __call__(self):
        """iteratively create AtomSpecies instances"""
        import itertools as it
        from os.path import join as pjoin
        from apns.test.tag_search import TagSearcher
        psplib = TagSearcher(pjoin(self.pseudo_dir, 'database.json'))
        pps = psplib.search(False, False, *(list(set(self.pptags + [self.symbol]))))
        assert len(pps) > 0, f'no pseudopotential found for {self.symbol} with tags {self.pptags}'

        from apns.test import database
        name = self.symbol
        fullname = database.PERIODIC_TABLE_TOFULLNAME[self.symbol]
        index = database.PERIODIC_TABLE_TOINDEX[self.symbol]
        rcovalent = database.RCOVALENT[self.symbol]
        mass = database.ATOMIC_MASS[self.symbol]
        magmom = 0.0
        
        for pp in pps:
            init_set = {"name": name, "fullname": fullname, "symbol": self.symbol, "index": index,
                        "rcovalent": rcovalent, "mass": mass, "magmom": magmom, "pp": pp, "nao": None}
            if self.ecut_db is not None:
                init_set.update({"ecutwfc": self.ecut_db.get(pp, None)})
            if self.orbital_dir is not None and self.naotags is not None:
                orblib = TagSearcher(pjoin(self.orbital_dir, 'database.json'))
                naos = orblib.search(False, False, *(list(set(self.naotags + [pp]))))
                for pp, nao in it.product(pps, naos):
                    init_set.update({"nao": nao})
                    yield AtomSpecies(**init_set)
            else:
                yield AtomSpecies(**init_set)

    def inflate(folded: dict):
        """expand the folded dict to the full dict if there are __element__ and __tags__
        two lists exist in folded dict. Expand to let every element in __element__ list
        has such value. If there is multiple definition, always let the one outside 
        __element__ overwrite the one inside __element__."""
        plain = {k: v for k, v in folded.items() if k not in ["__element__", "__tags__"]}
        for element in folded.get("__element__", []):
            plain.update({element: folded.get("__tags__", [])}) if element not in plain else None
        return plain

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
    label, fullname, symbol, index, rcovalent, mass, magmom = \
        None, None, None, None, None, None, None
    pp, nao = None, None
    ecutwfc = None
    
    def __init__(self, **kwargs) -> None:
        self.label = kwargs.get('label', None)
        self.fullname = kwargs.get('fullname', None)
        self.symbol = kwargs.get('symbol', None)
        self.index = kwargs.get('index', None)
        self.rcovalent = kwargs.get('rcovalent', None)
        self.mass = kwargs.get('mass', None)
        self.magmom = kwargs.get('magmom', None)
        self.pp = kwargs.get('pp', None)
        self.nao = kwargs.get('nao', None)
        self.ecutwfc = kwargs.get('ecutwfc', None)

    def __str__(self) -> str:
        """overload the print function to print the AtomSpecies information"""
        return (f"AtomSpecies: name: {self.label}, fullname: {self.fullname}, "
                f"symbol: {self.symbol}, index: {self.index}, covalent radius: {self.rcovalent}, "
                f"mass: {self.mass}, magmom: {self.magmom}, pseudopotential: {self.pp}, "
                f"Recommended ecutwfc: {self.ecutwfc} Ry, numerical atomic orbital: {self.nao}")

    def as_dict(self) -> dict:
        return {"name": self.label, "fullname": self.fullname, "symbol": self.symbol,
                "index": self.index, "rcovalent": self.rcovalent, "mass": self.mass,
                "magmom": self.magmom, "pp": self.pp, "ecutwfc": self.ecutwfc, "nao": self.nao}

class CellGenerator:
    """Python-Generator enpowered generator for generating Cell instances based on the
    simple imported proto (Bravis lattice, molecule, or CIF file), and several perturbation
    covering both the isolated/molecule and solid/periodic systems. It provides the flexibility
    to handling structure and poses to be a good abstraction. With it, developer can pay less
    attention to the details of the structure and focus other parts such as AtomSpecies,
    DFTParamSet, and so on."""
    # basic
    identifier, config = None, None
    # reciprocal space lattice information
    kspacing = None
    # "points" information
    magmoms = None
    # generator specific
    pertkind = None
    pertmags = None

    def __init__(self, config: str, pertmags: list, pertkind: str = "scale-abc", **kwargs) -> None:
        """create one cell instance with the scaling factors
        
        Parameters
        ----------
        config: str
            the configuration file, can be a CIF file, a Bravis lattice, or a molecule
        pertmags: list
            perturbation magnitudes, can be a list of floats
        pertkind: str
            the kind of perturbation, can be shear, twist, or scale. For molecule, do not use it
        kspacing: float
            the Monkhorst-Pack mesh spacing, only useful for Bravis lattice and CIF file
        magmoms: list
            magnetic moments for each atom, only useful for Bravis lattice and CIF file
        
        Raises
        ------
        AssertionError
            if the config is not a valid Bravis lattice, molecule or CIF file
            if the pertkind is not shear, twist, or scale
            if the pertmags is not a list of floats
            if the scales is not a list of floats
        """
        import re
        import os

        identifier = "bravis" if re.match(r"^([A-Z][a-z]?)_(sc|bcc|fcc|diamond)$", config) else "illegal"
        identifier = "bravis" if re.match(r"^([A-Z][a-z]?[A-Z][a-z]?)_(xy|xy2|xy3|x2y|x2y3|x2y5)$", config) else identifier
        identifier = "molecule" if re.match(r"^([A-Z][a-z]?)_(atom|monomer|dimer|trimer|tetramer)$", config) else identifier
        identifier = "cif" if os.path.exists(config) else identifier
        assert identifier != "illegal", f'config should be a valid Bravis lattice, molecule or CIF file: {config}'
        self.identifier = identifier
        self.config = config
        assert re.match(r"^(shear|twist|scale)(-[abc]+)?$", pertkind), f'pertkind should be shear, twist or scale: {pertkind}'
        self.pertkind = pertkind
        self.pertmags = pertmags

        assert isinstance(pertmags, list), f'scales should be a list of floats: {pertmags}'
        assert all([isinstance(scale, float) for scale in pertmags]), f'scales should be a list of floats: {pertmags}'

        self.kspacing = kwargs.get('kspacing', [-1.0])
        assert isinstance(self.kspacing, list) and all([isinstance(ks, float) for ks in self.kspacing])
        
        self.magmoms = kwargs.get('magmoms', None)
        assert self.magmoms is None or isinstance(self.magmoms, list), f'magmoms should be a list of floats: {self.magmoms}'
        print(f"""CellGenerator setup
Identifier (type of structure): {self.identifier}
Config (structure configuration): {self.config}
Scales: {self.pertmags}
K-spacing: {" ".join(map(str, self.kspacing))} in Bohr-1
Magnetic moments: {self.magmoms}
""")

    def __call__(self):
        """iteratively create Cell instances"""
        import itertools as it
        for s, k in it.product(*[self.pertmags, self.kspacing]):
            yield self.build(s, k, self.pertkind)

    def build(self, scale: float, kspacing: float, pertkind: str) -> dict:
        """called by __call__ to build the Cell instance
        
        Parameters
        ----------
        scale: float
            the scaling factor
        kspacing: float
            the Monkhorst-Pack mesh spacing
        pertkind: str
            the kind of perturbation
        
        Returns
        -------
        Cell
            the generated Cell instance"""
        build_func = CellGenerator.build_from_cif if self.identifier == "cif" \
            else CellGenerator.build_bravis
        build_func = CellGenerator.build_molecule if self.identifier == "molecule" else build_func
        unpertmag = scale if self.identifier == "molecule" else 1.0
        # for cif and bravis, can perform more kinds of perturbation, we will do it in the following
        param = build_func(self.config, unpertmag, kspacing)
        param["coords"] = param["coords"].tolist()
        if self.magmoms is not None:
            param["labels"] = CellGenerator.divide_subset(param["coords"], param["kinds"], self.magmoms, param["labels_kinds_map"])
            param["magmoms"] = self.magmoms
        else:
            param["magmoms"] = [0] * len(param["coords"])
        # hard code for now
        param["mobs"] = [[1, 1, 1]] * len(param["coords"])
        param["vels"] = [[0, 0, 0]] * len(param["coords"])

        proto = Cell(**param)
        if self.identifier == "molecule":
            return proto
        return CellTransformer._kernel(pertkind, proto, scale)

    def build_from_cif(fname: str, scale: float, kspacing: float = -1.0):
        """read structure file from external, and construct the structure.
        Currently only support CIF file.
        """
        from pymatgen.io.cif import CifParser
        parser = CifParser(fname)
        structure = parser.parse_structures(primitive=False)[0]
        labels = [site.species.elements[0].symbol for site in structure]
        kinds = list(dict.fromkeys(labels)) # find unique set but keep the sequence
        labels_kinds_map = [kinds.index(lable) for lable in labels]
        coords = structure.frac_coords
        magmoms = [0] * len(coords)
        a, b, c = [i*scale**(1/3) for i in structure.lattice.abc]
        alpha, beta, gamma = structure.lattice.angles

        import seekpath
        cell = CellGenerator.abc_angles_to_vec([a, b, c, alpha, beta, gamma], True)
        seekpath_result = seekpath.get_path((cell, coords, labels_kinds_map))
        sym_ks = seekpath_result['point_coords']
        possible_kpath = seekpath_result['path']
        mpmesh_nks = CellGenerator.kmeshgen(a, b, c, alpha, beta, gamma, kspacing) if kspacing > 0 else [1] * 3

        keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma', 'labels', 'kinds', 'labels_kinds_map', \
                'coords', 'magmoms', 'sym_ks', 'possible_kpath', 'mpmesh_nks', 'periodic']
        vals = [a, b, c, alpha, beta, gamma, labels, kinds, labels_kinds_map, \
                coords, magmoms, sym_ks, possible_kpath, mpmesh_nks, True]
        return dict(zip(keys, vals))
    
    def build_bravis(bravis: str, scale: float, kspacing: float = -1.0):
        from apns.test.bravis_and_molecule import lookup_acwf_db, lookup_bravis_angles, \
            vol_to_abc_angles, lookup_bravis_lattice
        vol = lookup_acwf_db(bravis)
        alpha, beta, gamma = lookup_bravis_angles(bravis)
        celldm = vol_to_abc_angles(vol, [alpha, beta, gamma])
        abc_angles, labels, kinds, labels_kinds_map, coords = lookup_bravis_lattice(bravis, celldm*scale**(1/3))
        a, b, c, alpha, beta, gamma = abc_angles
        magmoms = [0] * len(coords)

        import seekpath
        cell = CellGenerator.abc_angles_to_vec([a, b, c, alpha, beta, gamma], True)
        seekpath_result = seekpath.get_path((cell, coords, labels_kinds_map))
        sym_ks = seekpath_result['point_coords']
        possible_kpath = seekpath_result['path']
        mpmesh_nks = CellGenerator.kmeshgen(a, b, c, alpha, beta, gamma, kspacing) if kspacing > 0 else [1] * 3

        keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma', 'labels', 'kinds', 'labels_kinds_map', \
                'coords', 'magmoms', 'sym_ks', 'possible_kpath', 'mpmesh_nks', 'periodic']
        vals = [a, b, c, alpha, beta, gamma, labels, kinds, labels_kinds_map, \
                coords, magmoms, sym_ks, possible_kpath, mpmesh_nks, True]
        return dict(zip(keys, vals))

    def build_molecule(molecule: int, bond_length: float, kspacing: float = -1.0, celldm = 40.0):
        from apns.test.bravis_and_molecule import lookup_molecule
        print(f"kspacing: {kspacing} (bohr-1) setting is discarded due to isolated system generation.")
        abc_angles, labels, kinds, labels_kinds_map, coords = lookup_molecule(
            molecule, bond_length, 
            celldm = celldm, center = True)
        a, b, c, alpha, beta, gamma = abc_angles
        magmoms = [0] * len(coords)
        sym_ks, possible_kpath = None, None
        mpmesh_nks = [1, 1, 1]

        keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma', 'labels', 'kinds', 'labels_kinds_map', \
                'coords', 'magmoms', 'sym_ks', 'possible_kpath', 'mpmesh_nks', 'periodic']
        vals = [a, b, c, alpha, beta, gamma, labels, kinds, labels_kinds_map, \
                coords, magmoms, sym_ks, possible_kpath, mpmesh_nks, False]
        return dict(zip(keys, vals))

    def divide_subset(fullset: list, dividee: list, divider: list, dividee_fullset_map: list):
        """divide the labels into more piece"""
        if divider is None: return
        if divider == []: return
        assert isinstance(divider, list), f'identifiers should be a list: {divider}'
        assert len(divider) == len(fullset), f"everytime identifiers should be defined for each atom: {len(divider)} vs {len(fullset)}"
        idx_sp_, unique_, kind_count_, labels_ = 0, [], [0]*len(dividee), []
        for idx_sp, idtf in zip(dividee_fullset_map, divider):
            if (idx_sp, idtf) not in unique_:
                unique_.append((idx_sp, idtf))
                idx_sp_ += 1
                kind_count_[idx_sp] += 1
            labels_.append(dividee[idx_sp] + str(kind_count_[idx_sp]))
        return labels_

    def kmeshgen(a, b, c, alpha, beta, gamma, kspacing: float|list[float]):
        """get the Monkhorst-Pack mesh"""
        import numpy as np
        kspacing = kspacing if isinstance(kspacing, list) else [kspacing] * 3
        vecs = np.array(CellGenerator.abc_angles_to_vec([a, b, c, alpha, beta, gamma], True))
        recvecs = np.linalg.inv(vecs).T
        recvecs = 2*np.pi * recvecs
        norms = np.linalg.norm(recvecs, axis=1).tolist()
        assert len(norms) == len(kspacing), f'kspacing should be a list of 3 floats: {kspacing}'
        norms = [int(norm / kspac) for norm, kspac in zip(norms, kspacing)]
        return list(map(lambda x: max(1, x + 1), norms))

    def abc_angles_to_vec(lat: list, as_list: bool = False):
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
    
    periodic = True

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
        self.periodic = kwargs.get('periodic', True)

    def __str__(self) -> str:
        return f"""Cell
Real space lattice: \na: {self.a}, b: {self.b}, c: {self.c},
alpha: {self.alpha}, beta: {self.beta}, gamma: {self.gamma},
lattice constant: {self.lat0} Bohr/Angstrom\n
Reciprocal space lattice:
symmetry k-points: {self.sym_ks}, \npossible k-path: {self.possible_kpath},
Monkhorst-Pack mesh: {self.mpmesh_nks}\n
Structure:
coords: {self.coords}, \nvelocities: {self.vels}, \nmagnetic moments: {self.magmoms}, \noccupations: {self.mobs}
labels: {self.labels}, \nkinds: {self.kinds}, \nlabels_kinds_map: {self.labels_kinds_map}"""

    def as_dict(self) -> dict:
        return {"a": self.a, "b": self.b, "c": self.c, "alpha": self.alpha, "beta": self.beta, "gamma": self.gamma,
                "sym_ks": self.sym_ks, "possible_kpath": self.possible_kpath, "mpmesh_nks": self.mpmesh_nks,
                "coords": self.coords, "vels": self.vels, "magmoms": self.magmoms, "mobs": self.mobs,
                "labels": self.labels, "kinds": self.kinds, "labels_kinds_map": self.labels_kinds_map, "periodic": self.periodic}

class CellTransformer:
    """transform the Cell object with certain kind of perturbation.
    This collection of methods is especially useful in orbital generation
    tasks, while less useful in the general simulation tasks like 
    pseudopotential tests."""
    def _kernel(transform: str, cell: Cell, magnitude: float):
        """kernel function to call the specific transformation function
        
        Parameters
        ----------
        transform: str
            the kind of transformation, can be shear, twist, or scale with
            optional axis specification, like 'shear-a', 'twist-b', 'scale-ab',
            'scale-abc', etc.
        cell: Cell
            the cell object to be transformed
        magnitude: float
            the magnitude of the transformation
        
        Returns
        -------
        Cell
            the transformed cell object
        """
        call_map = {'shear': CellTransformer._shear, 
                    'twist': CellTransformer._twist, 
                    'scale': CellTransformer._scale}
        frags = transform.split('-')
        assert len(frags) == 2 or len(frags) == 1, \
            f'Invalid transform specified for CellTransformer: {transform}'
        transform = frags[0]
        axis = 'a' if len(frags) == 1 else frags[1]
        assert transform in call_map, f'transform should be one of shear, twist, scale: {transform}'
        print(f"CellTransformer setup\n\
Transformation: {transform}\n\
Magnitude: {magnitude}\n\
Axis: {axis}\n")
        return call_map[transform](cell, magnitude, axis)

    def _shear(cell: Cell, shear: float, axis: str = 'a'):
        """perform shear transformation keeping one axis fixed. Implemented from
        https://en.wikipedia.org/wiki/Shear_mapping

        Parameters
        ----------
        cell: Cell
            the cell object to be transformed
        shear: float
            the shear magnitude
        axis: str
            the axis that is fixed during the shear transformation
        
        Returns
        -------
        Cell
            the transformed cell object
        """
        import numpy as np
        assert axis in ['a', 'b', 'c'], f'axis should be one of a, b, c: {axis}'
        idx = ['a', 'b', 'c'].index(axis)
        shear_mat = np.eye(3)
        shear_mat[(idx+1)%3, (idx+2)%3] = shear
        cell.coords = np.dot(cell.coords, shear_mat)
        return cell
    
    def _twist(cell: Cell, twist: float, axis: str = 'a'):
        """perform twist on cell along one axis. From the bottom to the middle of the cell,
        the atomic layers are rotated in two dimension with a certain degree from 0 to deg.
        With this function, it is easy to generate the twisted bilayer graphene, or spirial
        structure. See unittest for details.

        Parameters
        ----------
        cell: Cell
            the cell object to be transformed
        twist: float
            the twist magnitude
        axis: str
            the axis that is fixed during the twist transformation
        
        Returns
        -------
        Cell
            the transformed cell object
        """
        import numpy as np
        assert axis in ['a', 'b', 'c'], f'axis should be one of a, b, c: {axis}'
        idx = ['a', 'b', 'c'].index(axis)
        # because the json cannot store numpy array, we need to convert it to list. But here
        # we need to convert it back to numpy array because the numpy array is more convenient
        # for matrix operation
        cell.coords = np.array(cell.coords)
        z = np.unique(cell.coords[:, idx]) # "z" does not mean narrowlly the z-axis, but for
                                           # easy understanding
        iatms = [np.where(cell.coords[:, idx] == i) for i in z]
        degs = np.linspace(0, twist, len(z)//2).tolist() + np.linspace(twist, 0, len(z)//2).tolist()
        for deg, idx_z in zip(degs, iatms):
            rot_mat = np.eye(3)
            rot_mat[(idx+1)%3, (idx+1)%3] = np.cos(np.deg2rad(deg))
            rot_mat[(idx+1)%3, (idx+2)%3] = -np.sin(np.deg2rad(deg))
            rot_mat[(idx+2)%3, (idx+1)%3] = np.sin(np.deg2rad(deg))
            rot_mat[(idx+2)%3, (idx+2)%3] = np.cos(np.deg2rad(deg))
            for iat in idx_z[0]:
                cell.coords[iat] = rot_mat @ cell.coords[iat]
        # convert back to list
        cell.coords = cell.coords.tolist()
        return cell

    def _scale(cell: Cell, scale: float, axis: str = 'a'):
        """perform scale transformation keeping one axis fixed.
        
        Parameters
        ----------
        cell: Cell
            the cell object to be transformed
        scale: float
            the scale magnitude
        axis: str
            the axis that is perform the scale transformation on. It can be a, b, c, ab, bc, ac, abc
            . For example if it is 'a', and the scale is given a value larger than 1, the cell will
            be stretched along the a axis. If 'ab', the cell will be stretched along the a and b axis
            with the same scale factor. If 'abc", it is the commonly seen volume scaling.
        
        Returns
        -------
        Cell
            the transformed cell object
        """
        assert axis in ['a', 'b', 'c', 'ab', 'bc', 'ac', 'abc'], \
        f'axis should be one of a, b, c, ab, bc, ac, abc: {axis}'
        if axis == 'a':
            cell.a *= scale
        elif axis == 'b':
            cell.b *= scale
        elif axis == 'c':
            cell.c *= scale
        elif axis == 'ab':
            cell.a *= scale**(1/2)
            cell.b *= scale**(1/2)
        elif axis == 'bc':
            cell.b *= scale**(1/2)
            cell.c *= scale**(1/2)
        elif axis == 'ac':
            cell.a *= scale**(1/2)
            cell.c *= scale**(1/2)
        elif axis == 'abc':
            cell.a *= scale**(1/3)
            cell.b *= scale**(1/3)
            cell.c *= scale**(1/3)
        return cell

import unittest
class TestAtomSpeciesGenerator(unittest.TestCase):

    def test_constructor(self):
        asg = AtomSpeciesGeneartor('Si', "", [])
        self.assertEqual(asg.symbol, 'Si')
        self.assertEqual(asg.pseudo_dir, "")
        self.assertEqual(asg.pptags, [])
        self.assertEqual(asg.orbital_dir, None)
        self.assertEqual(asg.naotags, None)
    
    def test_oncall(self):
        import os
        import json
        # prepare a fake database
        database = {"the_file_name_with_path": ["tag1", "tag2", "Si"]}
        fdatabase = "./database.json"
        with open(fdatabase, "w") as f:
            json.dump(database, f)
        asg = AtomSpeciesGeneartor('Si', "./", ["tag1", "tag2"])
        times = 0
        for as_ in asg():
            times += 1
            self.assertEqual(as_.pp, "the_file_name_with_path")
        self.assertEqual(times, 1)
        os.remove(fdatabase)

        asg = AtomSpeciesGeneartor('Si', "/root/abacus-develop/pseudopotentials/", ["Si", "PBE", "NC", "sg15", "sr"])
        times = 0
        ref = {'/root/abacus-develop/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.2.upf', 
               '/root/abacus-develop/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.0.upf', 
               '/root/abacus-develop/pseudopotentials/sg15_oncv_upf_2020-02-06/Si_ONCV_PBE-1.1.upf'}
        result = []
        for as_ in asg():
            result.append(as_.pp)
            times += 1
        result = set(result)
        self.assertEqual(result, ref)
        self.assertEqual(times, 3)

    def test_inflate(self):
        folded = {"__element__": ["Si", "Ge"], "__tags__": ["tag1", "tag2"]}
        plain = AtomSpeciesGeneartor.inflate(folded)
        self.assertEqual(plain, {"Si": ["tag1", "tag2"], "Ge": ["tag1", "tag2"]})

        folded = {"__element__": ["Si", "Ge"], "__tags__": ["tag1", "tag2"], "Si": ["tag3", "tag4"]}
        plain = AtomSpeciesGeneartor.inflate(folded)
        self.assertEqual(plain, {"Si": ["tag3", "tag4"], "Ge": ["tag1", "tag2"]})

        folded = {"Si": ["tag3", "tag4"], "Ge": ["tag1", "tag2"]}
        plain = AtomSpeciesGeneartor.inflate(folded)
        self.assertEqual(plain, {"Si": ["tag3", "tag4"], "Ge": ["tag1", "tag2"]})

        folded = {"Si": ["tag3", "tag4"], "Ge": ["tag1", "tag2"], "__element__": ["Si", "Ge"], "__tags__": ["tag1", "tag2"]}
        plain = AtomSpeciesGeneartor.inflate(folded)
        self.assertEqual(plain, {"Si": ["tag3", "tag4"], "Ge": ["tag1", "tag2"]})

class TestCellGeneartor(unittest.TestCase):

    def test_constructor(self):
        with self.assertRaises(AssertionError):
            CellGenerator('Si', 1.0) # raise error because scales must be list
        with self.assertRaises(AssertionError):
            CellGenerator('Si', [1.0])
        cg = CellGenerator('Si_dimer', [1.0])
        self.assertEqual(cg.identifier, 'molecule')
        self.assertEqual(cg.config, 'Si_dimer')
        self.assertEqual(cg.pertmags, [1.0])

    def test_kspacing(self):
        result = CellGenerator.kmeshgen(4.22798145, 4.22798145, 4.22798145, 60, 60, 60, 0.03*1.889725989)
        self.assertEqual(result, [33, 33, 33])

    def test_build_from_cif(self):
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
        with open(fcif, "w") as f:
            f.write(cif)
        cg = CellGenerator(fcif, [1.0])
        times = 0
        for cell in cg():
            self.assertAlmostEqual(cell.a, 4.84451701, delta=1e-6)
            self.assertAlmostEqual(cell.b, 4.84451701, delta=1e-6)
            self.assertAlmostEqual(cell.c, 4.84451701, delta=1e-6)
            self.assertAlmostEqual(cell.alpha, 60.0, delta=1e-6)
            self.assertAlmostEqual(cell.beta, 60.0, delta=1e-6)
            self.assertAlmostEqual(cell.gamma, 60.0, delta=1e-6)
            self.assertEqual(cell.labels, ["Ac", "Ac", "O"])
            self.assertEqual(cell.kinds, ["Ac", "O"])
            self.assertEqual(cell.labels_kinds_map, [0, 0, 1])
            self.assertEqual(cell.coords, [[0.75, 0.75, 0.75], [0.25, 0.25, 0.25], [0, 0, 0]])
            times += 1
        self.assertEqual(times, 1)

        cg = CellGenerator(fcif, [1.0, 1.02, 1.04, 1.06, 1.08])
        times = 0
        for cell in cg():
            self.assertAlmostEqual(cell.a, 4.84451701 * cg.pertmags[times]**(1/3), delta=1e-6)
            self.assertAlmostEqual(cell.b, 4.84451701 * cg.pertmags[times]**(1/3), delta=1e-6)
            self.assertAlmostEqual(cell.c, 4.84451701 * cg.pertmags[times]**(1/3), delta=1e-6)
            self.assertAlmostEqual(cell.alpha, 60.0, delta=1e-6)
            self.assertAlmostEqual(cell.beta, 60.0, delta=1e-6)
            self.assertAlmostEqual(cell.gamma, 60.0, delta=1e-6)
            self.assertEqual(cell.labels, ["Ac", "Ac", "O"])
            self.assertEqual(cell.kinds, ["Ac", "O"])
            self.assertEqual(cell.labels_kinds_map, [0, 0, 1])
            self.assertEqual(cell.coords, [[0.75, 0.75, 0.75], [0.25, 0.25, 0.25], [0, 0, 0]])
            times += 1
        self.assertEqual(times, 5)

        # test with magmom
        cg = CellGenerator(fcif, [1.0], magmoms = [1.0, -1.0, 0.0])
        times = 0
        for cell in cg():
            self.assertAlmostEqual(cell.a, 4.84451701, delta=1e-6)
            self.assertAlmostEqual(cell.b, 4.84451701, delta=1e-6)
            self.assertAlmostEqual(cell.c, 4.84451701, delta=1e-6)
            self.assertAlmostEqual(cell.alpha, 60.0, delta=1e-6)
            self.assertAlmostEqual(cell.beta, 60.0, delta=1e-6)
            self.assertAlmostEqual(cell.gamma, 60.0, delta=1e-6)
            self.assertEqual(cell.labels, ["Ac1", "Ac2", "O1"])
            self.assertEqual(cell.kinds, ["Ac", "O"])
            self.assertEqual(cell.labels_kinds_map, [0, 0, 1])
            self.assertEqual(cell.coords, [[0.75, 0.75, 0.75], [0.25, 0.25, 0.25], [0, 0, 0]])
            self.assertEqual(cell.magmoms, [1.0, -1.0, 0.0])
            times += 1

        os.remove(fcif)

    def test_build_bravis(self):
        cg = CellGenerator('Si_sc', [1.0]) # for bravis, the scale ought to be scaling factor but presently
        # is directly the volume of the cell
        times = 0
        for cell in cg():
            self.assertEqual(cell.a, 2.5318206866989112)
            self.assertEqual(cell.b, 2.5318206866989112)
            self.assertEqual(cell.c, 2.5318206866989112)
            self.assertEqual(cell.alpha, 90.0)
            self.assertEqual(cell.beta, 90.0)
            self.assertEqual(cell.gamma, 90.0)
            self.assertEqual(cell.labels, ["Si"])
            self.assertEqual(cell.kinds, ["Si"])
            self.assertEqual(cell.labels_kinds_map, [0])
            self.assertEqual(cell.coords[0], [0, 0, 0])
            times += 1
        self.assertEqual(times, 1)

        cg = CellGenerator('Si_sc', [1.0, 1.02, 1.04, 1.06, 1.08])
        times = 0
        for cell in cg():
            self.assertEqual(cell.a, 2.5318206866989112 * cg.pertmags[times]**(1/3))
            self.assertEqual(cell.b, 2.5318206866989112 * cg.pertmags[times]**(1/3))
            self.assertEqual(cell.c, 2.5318206866989112 * cg.pertmags[times]**(1/3))
            self.assertEqual(cell.alpha, 90.0)
            self.assertEqual(cell.beta, 90.0)
            self.assertEqual(cell.gamma, 90.0)
            self.assertEqual(cell.labels, ["Si"])
            self.assertEqual(cell.kinds, ["Si"])
            self.assertEqual(cell.labels_kinds_map, [0])
            self.assertEqual(cell.coords[0], [0, 0, 0])
            times += 1
        self.assertEqual(times, 5)

    def test_build_molecule(self):
        cg = CellGenerator('Ba_dimer', [1.0])
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
            self.assertEqual(cell.coords, [[10.0, 10.0, 10.0], [11.0, 10.0, 10.0]])
            times += 1
        self.assertEqual(times, 1)

        cg = CellGenerator('Ba_dimer', [1.0, 1.02, 1.04, 1.06, 1.08])
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
            self.assertEqual(cell.coords, [[10, 10, 10], [cg.pertmags[times]+10, 10, 10]])
            times += 1
        self.assertEqual(times, 5)

    def test_transform_cif(self):
        import uuid, os
        import numpy as np
        origin = """data_image0
_chemical_formula_structural       Pd40
_chemical_formula_sum              "Pd40"
_cell_length_a       5
_cell_length_b       5
_cell_length_c       50
_cell_angle_alpha    90
_cell_angle_beta     90
_cell_angle_gamma    90

_space_group_name_H-M_alt    "P 1"
_space_group_IT_number       1

loop_
  _space_group_symop_operation_xyz
  'x, y, z'

loop_
  _atom_site_type_symbol
  _atom_site_label
  _atom_site_symmetry_multiplicity
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  _atom_site_occupancy
  Pd  Pd1       1.0  0.00000  0.00000  0.00000  1.0000
  Pd  Pd2       1.0  0.50000  0.00000  0.00000  1.0000
  Pd  Pd3       1.0  0.00000  0.50000  0.00000  1.0000
  Pd  Pd4       1.0  0.50000  0.50000  0.00000  1.0000
  Pd  Pd5       1.0  0.00000  0.00000  0.10000  1.0000
  Pd  Pd6       1.0  0.50000  0.00000  0.10000  1.0000
  Pd  Pd7       1.0  0.00000  0.50000  0.10000  1.0000
  Pd  Pd8       1.0  0.50000  0.50000  0.10000  1.0000
  Pd  Pd9       1.0  0.00000  0.00000  0.20000  1.0000
  Pd  Pd10      1.0  0.50000  0.00000  0.20000  1.0000
  Pd  Pd11      1.0  0.00000  0.50000  0.20000  1.0000
  Pd  Pd12      1.0  0.50000  0.50000  0.20000  1.0000
  Pd  Pd13      1.0  0.00000  0.00000  0.30000  1.0000
  Pd  Pd14      1.0  0.50000  0.00000  0.30000  1.0000
  Pd  Pd15      1.0  0.00000  0.50000  0.30000  1.0000
  Pd  Pd16      1.0  0.50000  0.50000  0.30000  1.0000
  Pd  Pd17      1.0  0.00000  0.00000  0.40000  1.0000
  Pd  Pd18      1.0  0.50000  0.00000  0.40000  1.0000
  Pd  Pd19      1.0  0.00000  0.50000  0.40000  1.0000
  Pd  Pd20      1.0  0.50000  0.50000  0.40000  1.0000
  Pd  Pd21      1.0  0.00000  0.00000  0.50000  1.0000
  Pd  Pd22      1.0  0.50000  0.00000  0.50000  1.0000
  Pd  Pd23      1.0  0.00000  0.50000  0.50000  1.0000
  Pd  Pd24      1.0  0.50000  0.50000  0.50000  1.0000
  Pd  Pd25      1.0  0.00000  0.00000  0.60000  1.0000
  Pd  Pd26      1.0  0.50000  0.00000  0.60000  1.0000
  Pd  Pd27      1.0  0.00000  0.50000  0.60000  1.0000
  Pd  Pd28      1.0  0.50000  0.50000  0.60000  1.0000
  Pd  Pd29      1.0  0.00000  0.00000  0.70000  1.0000
  Pd  Pd30      1.0  0.50000  0.00000  0.70000  1.0000
  Pd  Pd31      1.0  0.00000  0.50000  0.70000  1.0000
  Pd  Pd32      1.0  0.50000  0.50000  0.70000  1.0000
  Pd  Pd33      1.0  0.00000  0.00000  0.80000  1.0000
  Pd  Pd34      1.0  0.50000  0.00000  0.80000  1.0000
  Pd  Pd35      1.0  0.00000  0.50000  0.80000  1.0000
  Pd  Pd36      1.0  0.50000  0.50000  0.80000  1.0000
  Pd  Pd37      1.0  0.00000  0.00000  0.90000  1.0000
  Pd  Pd38      1.0  0.50000  0.00000  0.90000  1.0000
  Pd  Pd39      1.0  0.00000  0.50000  0.90000  1.0000
  Pd  Pd40      1.0  0.50000  0.50000  0.90000  1.0000
"""
        fcif = uuid.uuid4().hex + '.cif'
        with open(fcif, "w") as f:
            f.write(origin)

        # twist
        cg = CellGenerator(fcif, [45.0], pertkind='twist-c')
        spiral = [
            [0.00000,  0.00000,  0.00000],
            [0.50000,  0.00000,  0.00000],
            [0.00000,  0.50000,  0.00000],
            [0.50000,  0.50000,  0.00000],
            [0.00000,  0.00000,  0.10000],
            [0.49039,  0.09755,  0.10000],
            [0.90245,  0.49039,  0.10000],
            [0.39285,  0.58794,  0.10000],
            [0.00000,  0.00000,  0.20000],
            [0.46194,  0.19134,  0.20000],
            [0.80866,  0.46194,  0.20000],
            [0.27060,  0.65328,  0.20000],
            [0.00000,  0.00000,  0.30000],
            [0.41573,  0.27779,  0.30000],
            [0.72221,  0.41573,  0.30000],
            [0.13795,  0.69352,  0.30000],
            [0.00000,  0.00000,  0.40000],
            [0.35355,  0.35355,  0.40000],
            [0.64645,  0.35355,  0.40000],
            [0.00000,  0.70711,  0.40000],
            [0.00000,  0.00000,  0.50000],
            [0.35355,  0.35355,  0.50000],
            [0.64645,  0.35355,  0.50000],
            [0.00000,  0.70711,  0.50000],
            [0.00000,  0.00000,  0.60000],
            [0.41573,  0.27779,  0.60000],
            [0.72221,  0.41573,  0.60000],
            [0.13795,  0.69352,  0.60000],
            [0.00000,  0.00000,  0.70000],
            [0.46194,  0.19134,  0.70000],
            [0.80866,  0.46194,  0.70000],
            [0.27060,  0.65328,  0.70000],
            [0.00000,  0.00000,  0.80000],
            [0.49039,  0.09755,  0.80000],
            [0.90245,  0.49039,  0.80000],
            [0.39285,  0.58794,  0.80000],
            [0.00000,  0.00000,  0.90000],
            [0.50000,  0.00000,  0.90000],
            [0.00000,  0.50000,  0.90000],
            [0.50000,  0.50000,  0.90000]
        ]
        # the above is the reference value, now we compare between
        # the reference and the transformed one
        times = 0
        for cell in cg():
            for i, j in zip(cell.coords, spiral):
                for k, l in zip(i, j):
                    self.assertAlmostEqual(k%1, l%1, delta=1e-5)
            times += 1

        self.assertEqual(times, 1)

        # shear
        cg = CellGenerator(fcif, [0.5], pertkind='shear-c')

        refxy = [[0.5, 0.25], [0.5, 0.75], [0.0, 0.0], [0.0, 0.5]]
        refxy = np.array(refxy).T.tolist()
        times = 0
        for cell in cg():
            for c in cell.coords:
                self.assertTrue(c[0]%1 in refxy[0])
                self.assertTrue(c[1]%1 in refxy[1])
            times += 1
        os.remove(fcif)

        # scale is too simple to test

if __name__ == "__main__":
    unittest.main(exit=False)
