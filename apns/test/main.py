def run(finp: str):
    import json
    import apns.test.citation as amc
    with open(finp, "r") as f:
        inp = json.load(f)
    main(inp)
    amc.citation()

def main(inp: dict):
    # then collect structures, update the descriptor in inp
    import apns.test.apns_io as apns_io
    inp = apns_io.prepare(inp)
    print("* * * Preparation done * * *".center(100, " "))
    # then loop over struset
    import apns.test.dft_param_set as dftparam
    import apns.test.atom_species_and_cell as struparam
    import itertools as it
    import os
    out_dir = inp["global"].get("out_dir", ".")
    folders = []
    for iset, struset in enumerate(inp["strusets"]):
        software = struset.get("calculator", "abacus")
        calcset = struset.get("calcset", 0)
        atomset = struset.get("atomset", 0)
        print(f"""Struset {iset} module binding information
DFT package: {software}
Calcset: {calcset}
Atomset: {atomset}
Number of structure descriptors: {len(struset["desc"])}
""")
        # the dft parameters
        dftgen = dftparam.DFTParamSetGenerator(inp[software][calcset])
        # the atom species
        asgens = [struparam.AtomSpeciesGeneartor(symbol=s, 
                  pseudo_dir=inp["global"].get("pseudo_dir", "."), pptags=tags[0], 
                  orbital_dir=inp["global"].get("orbital_dir", "."), naotags=tags[1]) \
                  for s, tags in inp["atomsets"][atomset].items()]
        # connect the asgens with converged ecutwfc (if possible), to support the arg `ecutwfc = "auto"`
        if software == "abacus":
            for asgen in asgens: 
                asgen.connect_ecut_db(os.path.join(inp["global"].get("cache_dir", "."), "ecutwfc.json"))
        # the cell
        for desc in struset["desc"]:
            # first extract the atomspecies that really needed
            cellgen = struparam.CellGenerator(**desc)
            for paramset, cell in it.product(dftgen(), cellgen()):
                asgens_subset, cell = bind_atom_species_with_cell(asgens, cell)
                for atom_species in it.product(*[asgen() for asgen in asgens_subset]):
                    standard_fname_with_contents = export(paramset, atom_species, cell, fmt=software)
                    cache = dict(zip(["AtomSpecies", "Cell", "DFTParamSet", "CellGenerator"], \
                    [[as_.as_dict() for as_ in atom_species], cell.as_dict(), paramset, 
                     dict(zip(["identifier", "config"], [cellgen.identifier, cellgen.config]))]))
                    folders.append(write_and_move(standard_fname_with_contents, cache, out_dir))
    print("* * * All structures generated * * *".center(100, " "))
    return folders

def write_and_move(standard_fname_with_contents: dict, cache: dict, out_dir: str):
    from os import makedirs
    from os.path import basename
    import uuid
    folder = f"{out_dir}/{str(uuid.uuid4())}"
    makedirs(folder, exist_ok=True)
    for fname, content in standard_fname_with_contents.items():
        with open(f"{folder}/{fname}", 'w') as f:
            f.write(content)
    with open(f"{folder}/description.json", 'w') as f:
        import json
        json.dump(cache, f, indent=4)
    # copy the pseudopotentials and numerical orbitals
    import shutil
    for as_ in cache["AtomSpecies"]:
        if as_["pp"] is not None:
            shutil.copy(f"{as_['pp']}", f"{folder}/{basename(as_['pp'])}")
        if as_["nao"] is not None:
            shutil.copy(f"{as_['nao']}", f"{folder}/{basename(as_['nao'])}")
    return folder

from apns.test.atom_species_and_cell import AtomSpeciesGeneartor, AtomSpecies, Cell
def bind_atom_species_with_cell(atom_species_generators: list[AtomSpeciesGeneartor], cell: Cell):
    assert all([isinstance(asgen, AtomSpeciesGeneartor) for asgen in atom_species_generators])
    those_wanted = [asgen for asgen in atom_species_generators if asgen.symbol in cell.kinds]
    assert len(those_wanted) == len(cell.kinds), "Not all AtomSpecies are found in the cell"
    return those_wanted, cell

def export(paramset: dict, atomset: list, cell: Cell, fmt = "abacus") -> dict:
    """with paramset, atomset and Cell, build input files for specific DFT software.
    Development:
    Once add new support, update the in-built dict"""
    call_map = {"abacus": write_abacus, "qespresso": write_qespresso}
    assert fmt in call_map, "Unsupported DFT software"
    return call_map[fmt](paramset, atomset, cell)

def write_qespresso(paramset: dict, atomset: list, cell: Cell):
    return {"pwscf.in": write_qespresso_in(paramset, atomset, cell)}

def write_qespresso_in(paramset: dict, atomset: list, cell: Cell):
    from apns.test.atom_species_and_cell import CellGenerator
    import os
    sections = ["control", "system", "electrons", "ions", "cell"] # seems must be in this order
    nat = len(cell.labels)
    ntyp = len(set(cell.labels))
    result = ""
    for section in sections:
        result += f"&{section.upper()}\n"
        source = paramset.get(section, {})
        source.update({"nat": nat, "ntyp": ntyp}) if section == "system" else None
        for k, v in source.items():
            v = v if v in [".true.", ".false."] else v if not isinstance(v, str) else f"'{v}'"
            result += f"{k:<20s} = {v}\n"
        result += "/\n\n"
    result += f"CELL_PARAMETERS (angstrom)\n"
    latvec = CellGenerator.abc_angles_to_vec([cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma], True)
    for i in range(3):
        result += f"{latvec[i][0]:>20.10f} {latvec[i][1]:>20.10f} {latvec[i][2]:>20.10f}\n"
    result += "\nATOMIC_SPECIES\n"
    for label in dict.fromkeys(cell.labels):
        ind = [i for i, aspec in enumerate(atomset) if aspec.symbol == cell.kinds[cell.labels_kinds_map[cell.labels.index(label)]]]
        result += f"{label:<4s} {atomset[ind[0]].mass:<8.4f} {os.path.basename(atomset[ind[0]].pp)}\n"
    result += "\nATOMIC_POSITIONS (crystal)\n"
    for i, l in enumerate(cell.labels):
        result += f"{l:<4s} {cell.coords[i][0]:>20.10f} {cell.coords[i][1]:>20.10f} {cell.coords[i][2]:>20.10f}\n"
    result += f"\nK_POINTS (automatic)\n{cell.mpmesh_nks[0]} {cell.mpmesh_nks[1]} {cell.mpmesh_nks[2]} 0 0 0\n"

    return result

def write_abacus(paramset: dict, atomset: list, cell: Cell):
    keys = ["INPUT", "STRU", "KPT"]
    # here it is possible to support the converged ecutwfc value auto-set for INPUT.
    # first get the fpp from atomset, then get the max, set to ecutwfc in paramset
    ecut_set = paramset.get("ecutwfc", None) # if ecut_set is None or "auto", then set it to the max of fpp
    if ecut_set is None or ecut_set == "auto":
        ecutwfc = max([as_.ecutwfc for as_ in atomset if as_.ecutwfc is not None])
        paramset["ecutwfc"] = ecutwfc
    vals = [write_abacus_input(paramset), write_abacus_stru(atomset, cell), write_abacus_kpt(cell)]
    return dict(zip(keys, vals))

def write_abacus_input(paramset: dict):
    # first precondition such that all parameters are set as string
    paramset = {k: str(v) if not isinstance(v, list) else " ".join([str(_v) for _v in v]) for k, v in paramset.items()}
    result = "INPUT_PARAMETERS\n"
    result += "\n".join([f"{k:<20s} {v:<s}" for k, v in paramset.items()])
    return result

def write_abacus_stru(atomset: list[AtomSpecies], cell: Cell):
    from os.path import basename
    result = "ATOMIC_SPECIES\n"
    # need a map from cell.labels to index of AtomSpecies in atomset list!
    temp_ = [a.symbol for a in atomset]
    uniquelabels_atomspecies_map = [temp_.index(cell.kinds[cell.labels_kinds_map[cell.labels.index(ulbl)]]) for ulbl in dict.fromkeys(cell.labels)]
    result += "\n".join([f"{label:<4s} \
{atomset[uniquelabels_atomspecies_map[il]].mass:<8.4f} {basename(atomset[uniquelabels_atomspecies_map[il]].pp)}" \
            for il, label in enumerate(dict.fromkeys(cell.labels))])
    
    if not all([as_.nao is None for as_ in atomset]):
        result += "\n\nNUMERICAL_ORBITAL\n"
        result += "\n".join([f"{basename(atomset[uniquelabels_atomspecies_map[il]].nao)}" for il in range(len(set(cell.labels)))])

    result += f"\n\nLATTICE_CONSTANT\n{cell.lat0:<20.10f}\n"
    result += "\nLATTICE_VECTORS\n"
    from apns.test.atom_species_and_cell import CellGenerator
    latvec = CellGenerator.abc_angles_to_vec([cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma], True)
    result += "\n".join([f"{latvec[i][0]:<20.10f} {latvec[i][1]:<20.10f} {latvec[i][2]:<20.10f}" for i in range(3)])
    
    coord = "Direct" if cell.periodic else "Cartesian"
    result += f"\n\nATOMIC_POSITIONS\n{coord}\n"
    for label in dict.fromkeys(cell.labels):
        ind = [i for i, l in enumerate(cell.labels) if l == label]
        result += f"{label}\n{cell.magmoms[ind[0]]:<4.2f}\n{len(ind)}\n"
        for i in ind:
            result += f"{cell.coords[i][0]:<20.10f}{cell.coords[i][1]:<20.10f}{cell.coords[i][2]:<20.10f} \
m {cell.mobs[i][0]:<2d}{cell.mobs[i][1]:<2d}{cell.mobs[i][2]:<2d}\n"
    return result

def write_abacus_kpt(cell: Cell):
    """it is not clarified how to perform a two-step calculation like scf-nscf, for example the
    band structure calculation, therefore, the band structure calculation is not supported yet."""
    result = f"K_POINTS\n0\nGamma\n\
{cell.mpmesh_nks[0]} {cell.mpmesh_nks[1]} {cell.mpmesh_nks[2]} 0 0 0"
    return result

import unittest
class TestExportFunctions(unittest.TestCase):

    def test_write_abacus_input(self):
        result = write_abacus_input({
            "calculation": "scf",
            "basis_type": "pw",
            "nspin": 1,
            "ecutwfc": 60,
            "cal_force": 1,
            "symmetry": 0
        })
        ref = """INPUT_PARAMETERS
calculation          scf
basis_type           pw
nspin                1
ecutwfc              60
cal_force            1
symmetry             0"""
        self.assertEqual(result, ref)

    def test_write_abacus_stru(self):
        # periodic cell case, bravis or cif
        cell = Cell(a=10, b=10, c=10, alpha=60, beta=60, gamma=60, 
                    coords=[[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
                    kinds=["Si", "O"], labels=["Si1", "O1"], labels_kinds_map=[0, 1],
                    magmoms=[0.0, 0.0], mobs=[[1, 1, 1], [1, 1, 1]], periodic=True)
        Si = AtomSpecies(symbol="Si", mass=28.0855, pp="Si.pz-vbc.UPF")
        O = AtomSpecies(symbol="O", mass=15.9994, pp="O.pbe-rrkjus.UPF")
        result = write_abacus_stru([Si, O], cell)
        ref = """ATOM_SPECIES
Si1  28.0855  Si.pz-vbc.UPF
O1   15.9994  O.pbe-rrkjus.UPF

LATTICE_CONSTANT
1.8897259886        

LATTICE_VECTORS
10.0000000000        0.0000000000         0.0000000000        
5.0000000000         8.6602540378         0.0000000000        
5.0000000000         2.8867513459         8.1649658093        

ATOM_POSITIONS
Direct
Si1
1
0.00
0.0000000000        0.0000000000        0.0000000000                         m 1 1 1 
O1
1
0.00
0.2500000000        0.2500000000        0.2500000000                         m 1 1 1"""
        self.assertTrue(ref in result)

        # test with nao
        Si = AtomSpecies(symbol="Si", mass=28.0855, pp="Si.pz-vbc.UPF", nao="Si_gga_6au_100Ry_2s2p1d.orb")
        O = AtomSpecies(symbol="O", mass=15.9994, pp="O.pbe-rrkjus.UPF", nao="O_gga_6au_100Ry_2s2p1d.orb")
        result = write_abacus_stru([Si, O], cell)
        ref = """ATOM_SPECIES
Si1  28.0855  Si.pz-vbc.UPF
O1   15.9994  O.pbe-rrkjus.UPF

NUMERICAL_ORBITAL
Si_gga_6au_100Ry_2s2p1d.orb
O_gga_6au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897259886        

LATTICE_VECTORS
10.0000000000        0.0000000000         0.0000000000        
5.0000000000         8.6602540378         0.0000000000        
5.0000000000         2.8867513459         8.1649658093        

ATOM_POSITIONS
Direct
Si1
1
0.00
0.0000000000        0.0000000000        0.0000000000                         m 1 1 1 
O1
1
0.00
0.2500000000        0.2500000000        0.2500000000                         m 1 1 1"""
        self.assertTrue(ref in result)
        # isolated case, molecule
        cell = Cell(a=20, b=20, c=20, alpha=90, beta=90, gamma=90, 
                    coords=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], 
                    kinds=["H"], labels=["H", "H"], labels_kinds_map=[0, 0], 
                    magmoms=[0.0, 0.0], mobs=[[1, 1, 1], [1, 1, 1]], periodic=False)
        H = AtomSpecies(symbol="H", mass=1.0079, pp="H.pz-vbc.UPF")
        result = write_abacus_stru([H], cell)
        ref = """ATOM_SPECIES
H    1.0079   H.pz-vbc.UPF

LATTICE_CONSTANT
1.8897259886        

LATTICE_VECTORS
20.0000000000        0.0000000000         0.0000000000        
0.0000000000         20.0000000000        0.0000000000        
0.0000000000         0.0000000000         20.0000000000       

ATOM_POSITIONS
Cartesian
H
2
0.00
0.0000000000        0.0000000000        0.0000000000                         m 1 1 1 
0.0000000000        0.0000000000        1.0000000000                         m 1 1 1"""
        self.assertTrue(ref in result)
        # test with nao
        H = AtomSpecies(symbol="H", mass=1.0079, pp="H.pz-vbc.UPF", nao="H_gga_6au_100Ry_2s2p1d.orb")
        result = write_abacus_stru([H], cell)
        ref = """ATOM_SPECIES
H    1.0079   H.pz-vbc.UPF

NUMERICAL_ORBITAL
H_gga_6au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897259886        

LATTICE_VECTORS
20.0000000000        0.0000000000         0.0000000000        
0.0000000000         20.0000000000        0.0000000000        
0.0000000000         0.0000000000         20.0000000000       

ATOM_POSITIONS
Cartesian
H
2
0.00
0.0000000000        0.0000000000        0.0000000000                         m 1 1 1 
0.0000000000        0.0000000000        1.0000000000                         m 1 1 1"""
        self.assertTrue(ref in result)

if __name__ == "__main__":
    unittest.main(exit=False)
    inp = {
        "global": {
            "mode": "test",
            "pseudo_dir": "./download/pseudopotentials",
            "orbital_dir": "./download/numerical_orbitals",
            "cache_dir": "./apns_cache",
            "out_dir": "./apns_out"
        },
        "credentials": {
            "materials_project": {"api_key": ""}
        },
        "abacus": [
            {
                "ecutwfc": 30,
                "calculation": "cell-relax",
                "nspin": 2
            }
        ],
        "atomsets": [
            {
                "Ba": [["sg15", "1.0"], None],
                "O": [["GBRV"], None],
                "Ti": [["rrkjus"], None],
                "H": [["PD", "04"], None],
                "Fe": [["sg15", "1.0"], None],
                "Mn": [["sg15", "1.0"], None],
                "Y": [["sg15", "1.0"], None],
            }
        ],
        "strusets": [
            {
                "calculator": "abacus",
                "calcset": 0, "atomset": 0,
                "desc": [["search", "Y2O3", [0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06], [0.05]]]
            }
        ]
    }
    main(inp)