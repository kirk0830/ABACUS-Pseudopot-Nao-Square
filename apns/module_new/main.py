def main(finp: str):
    import json
    with open(finp, 'r') as f:
        inp = json.load(f)
    
    # first check

    # then collect structures, update the descriptor in inp

    # then loop over struset
    import apns.module_new.dft_param_set as dftparam
    import apns.module_new.atom_species_and_cell as struparam
    import itertools as it
    for struset in inp["strusets"]:
        # the dft parameters
        software = struset.get("calculator", "abacus")
        calcset = struset.get("calcset", 0)
        dftgen = dftparam.DFTParamSetGenerator(inp[software][calcset])
        # the atom species
        atomset = struset.get("atomset", 0)
        asgens = [struparam.AtomSpeciesGeneartor(s, **tags) for s, tags in inp["atomsets"][atomset]]
        # the cell
        for desc in struset["desc"]:
            # first extract the atomspecies that really needed
            cellgen = struparam.CellGenerator(desc[0], desc[1], desc[2]) # kspacing? magmom?
            for paramset, cell in it.product(dftgen, cellgen):
                asgens_subset, cell = bind_atom_species_with_cell(asgens, cell)
                for atom_species in it.product(*asgens_subset):
                    export(paramset, atom_species, cell, fmt=software)

from apns.module_new.atom_species_and_cell import AtomSpeciesGeneartor, AtomSpecies, Cell
def bind_atom_species_with_cell(atom_species_generators: list[AtomSpecies], cell: Cell):
    assert all([isinstance(asgen, AtomSpeciesGeneartor) for asgen in atom_species_generators])
    those_wanted = [asgen for asgen in atom_species_generators if asgen.symbol in cell.kinds]
    return those_wanted, cell

def export(paramset: dict, atomset: list, cell: Cell, fmt = "abacus"):
    pass

def write_abacus_input(paramset: dict):
    # first precondition such that all parameters are set as string
    paramset = {k: str(v) if not isinstance(v, list) else " ".join([str(_v) for _v in v]) for k, v in paramset.items()}
    result = "INPUT_PARAMETERS\n"
    result += "\n".join([f"{k:<20s} {v:<s}" for k, v in paramset.items()])
    return result

def write_abacus_stru(atomset: list[AtomSpecies], cell: Cell):
    from os.path import basename
    result = "ATOM_SPECIES\n"
    # need a map from cell.labels to index of AtomSpecies in atomset list!
    temp_ = [a.symbol for a in atomset]
    labels_atomspecies_map = [temp_.index(cell.kinds[ik]) for ik in cell.labels_kinds_map]
    result += "\n".join([f"{label:<4s} {atomset[labels_atomspecies_map[il]].mass:<8.4f} {basename(atomset[labels_atomspecies_map[il]].pp)}" \
            for il, label in enumerate(dict.fromkeys(cell.labels))])
    result = "\nNUMERICAL_ORBITAL\n"
    result += "\n".join([f"{basename(atomset[ias].nao)}" for ias in labels_atomspecies_map])
    result += f"\nLATTICE_CONSTANT\n{cell.lat0:<20.10f}\n"
    result += "\nLATTICE_VECTORS\n"
    from apns.module_new.atom_species_and_cell import CellGenerator
    latvec = CellGenerator.abc_angles_to_vec([cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma], True)
    result += "\n".join([f"{latvec[i][0]:<20.10f} {latvec[i][1]:<20.10f} {latvec[i][2]:<20.10f}" for i in range(3)])
    
    coord = "Direct" if cell.periodic else "Cartesian"
    result += f"\nATOM_POSITIONS\n{coord}\n"
    for label in dict.fromkeys(cell.labels):
        ind = [i for i, l in enumerate(cell.labels) if l == label]
        result += f"{label}\n{len(ind)}\n{cell.magmoms[ind[0]]}\n"
        for i in ind:
            result += f"{cell.coords[i][0]:<20.10f}{cell.coords[i][1]:<20.10f}{cell.coords[i][2]:<20.10f} \
                m {cell.mobs[i][0]:<2d}{cell.mobs[i][1]:<2d}{cell.mobs[i][2]:<2d}\n"
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
        cell = Cell(a=10, b=10, c=10, alpha=60, beta=60, gamma=60, 
                    coords=[[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
                    kinds=["Si", "O"], labels=["Si1", "O1"], labels_kinds_map=[0, 1],
                    magmoms=[0.0, 0.0])
        Si = AtomSpecies()

if __name__ == "__main__":
    unittest.main()