class OrbgenParamSetGenerator:

    inp = None
    pseudo_dir, cache_dir = None, None
    def __init__(self, finp: str, pseudo_dir: str, cache_dir: str):
        import json
        with open(finp, 'r') as f:
            self.inp = json.load(f)
        self.pseudo_dir, self.cache_dir = pseudo_dir, cache_dir

    def nested_parse(nested: dict, excluded: list = None):
        """return a generator based on present nested dictionary, with consideration
        of keywords in excluded """
        import itertools as it
        excluded = [] if excluded is None else excluded
        scalar = {k: v for k, v in nested.items() if k in excluded or not isinstance(v, list)}
        scalar_joint = {k: v for k, v in nested.items() if \
        ("|" in k) and isinstance(v, list) and (len(v) == k.count("|") + 1) and not all([isinstance(i, list) for i in v])}
        scalar.update(scalar_joint)
        calculate = {k: v for k, v in nested.items() if k not in scalar.keys()}
        
        cal_keys = calculate.keys()
        cal_vals = calculate.values()
        for val in it.product(*cal_vals):
            joint = dict(zip(cal_keys, val))
            joint.update(scalar)
            yield OrbgenParamSetGenerator.dejoint(joint)
    
    def dejoint(joint: dict):
        """convert joint dictionary to normal one"""
        normal = {}
        for k, v in joint.items():
            if "|" in k:
                keys = k.split("|")
                for key in keys:
                    normal[key] = v[keys.index(key)]
            else:
                normal[k] = v
        return normal

    def characterize(pseudo_dir: str, cache_dir: str, atomspecies):
        from apns.test.atom_species_and_cell import AtomSpecies
        import json
        import os
        assert isinstance(atomspecies, AtomSpecies), "atomspecies should be an instance of AtomSpecies"
        with open(os.path.join(cache_dir, "ecutwfc.json"), "r") as f:
            ecutwfc = json.load(f)
        pp = atomspecies.pp
        if pp not in ecutwfc.keys(): 
            print(f"pseudo potential {pp} not found in ecutwfc.json")
            return None
        return {"pseudo_dir": pseudo_dir, "pseudo_name": os.path.basename(pp), "ecutwfc": ecutwfc[pp]}

    def __call__(self):
        from apns.test.atom_species_and_cell import AtomSpeciesGeneartor
        import itertools as it
        for iorb, orbset in enumerate(self.inp["orbsets"]):
            print(f"* * * Generating orbgen input for orbset {iorb} * * *".center(100))
            unpacked = OrbgenParamSetGenerator.nested_parse(orbset)
            for orbs in OrbgenParamSetGenerator.recomb_orbset(unpacked): # prototype of "orbitals" section
                print(f"Unpack nested orbset: \n{orbs}")
                ia, isp, istru = orbs[0]["atomset"], orbs[0]["spillset"], orbs[0]["struset"]
                atomset, spillset, struset = self.inp["atomsets"][ia], self.inp["spillsets"][isp], self.inp["strusets"][istru]
                for element, tags in atomset.items(): # generate for each element...
                    asgen = AtomSpeciesGeneartor(element, self.pseudo_dir, tags[0])
                    for as_, spill_, stru_ in it.product(asgen(), 
                    OrbgenParamSetGenerator.nested_parse(spillset, ["spill_coefs"]), 
                    OrbgenParamSetGenerator.nested_parse(struset)):
                        print(f"as: \n{as_}, \nspill: \n{spill_}, \nstru: \n{stru_}")

    def recomb_orbset(unpacked: list):
        import itertools as it
        groups = {}
        for u in unpacked:
            groups.setdefault(u["nzeta"], []).append(u)
        yield from it.product(*groups.values())

    def build(orb_combs: list, strusets: list):
        """because the orbital is somehow bound with the structure"""
        result = {"orbitals": [], "reference_systems": []}
        for orbs in orb_combs: # get the trivial concenation of orbital settings, no fold
            for orb in orbs: # for each level of orbitals
                folded_stru = strusets[orb["struset"]]

import unittest
class TestOrbgenParamSetGenerator(unittest.TestCase):
    def test_dejoint(self):
        joint = {"a|b|c": [1, 2, 3]}
        normal = {"a": 1, "b": 2, "c": 3}
        self.assertEqual(OrbgenParamSetGenerator.dejoint(joint), normal)
    
        joint = {"a|b|c": [1, 2, 3], "d": 4}
        normal = {"a": 1, "b": 2, "c": 3, "d": 4}
        self.assertEqual(OrbgenParamSetGenerator.dejoint(joint), normal)

        joint = {"a|b|c": [1, [2, 3], 4], "d": 5}
        normal = {"a": 1, "b": [2, 3], "c": 4, "d": 5}
        self.assertEqual(OrbgenParamSetGenerator.dejoint(joint), normal)
    
    def test_nested_parse(self):
        nested = {"a": [1, 2], "b": [3, 4], "c": [5, 6]}
        excluded = ["a"]
        normal = [{"a": [1, 2], "b": 3, "c": 5}, {"a": [1, 2], "b": 3, "c": 6}, \
        {"a": [1, 2], "b": 4, "c": 5}, {"a": [1, 2], "b": 4, "c": 6}]
        self.assertEqual(list(OrbgenParamSetGenerator.nested_parse(nested, excluded)), normal)

        nested = {"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": [7, 8]}
        excluded = ["a", "b"]
        normal = [{"a": [1, 2], "b": [3, 4], "c": 5, "d": 7}, {"a": [1, 2], "b": [3, 4], "c": 5, "d": 8}, \
        {"a": [1, 2], "b": [3, 4], "c": 6, "d": 7}, {"a": [1, 2], "b": [3, 4], "c": 6, "d": 8}]
        self.assertEqual(list(OrbgenParamSetGenerator.nested_parse(nested, excluded)), normal)

        nested = {"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": [7, 8], "e|f": [9, 10]}
        excluded = ["a", "b"]
        normal = [{"a": [1, 2], "b": [3, 4], "c": 5, "d": 7, "e": 9, "f": 10}, \
        {"a": [1, 2], "b": [3, 4], "c": 5, "d": 8, "e": 9, "f": 10}, \
        {"a": [1, 2], "b": [3, 4], "c": 6, "d": 7, "e": 9, "f": 10}, \
        {"a": [1, 2], "b": [3, 4], "c": 6, "d": 8, "e": 9, "f": 10}]
        self.assertEqual(list(OrbgenParamSetGenerator.nested_parse(nested, excluded)), normal)

        nested = {"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": [7, 8], "e|f": [9, [10, 11]]}
        excluded = ["a", "b"]
        normal = [{"a": [1, 2], "b": [3, 4], "c": 5, "d": 7, "e": 9, "f": [10, 11]}, \
        {"a": [1, 2], "b": [3, 4], "c": 5, "d": 8, "e": 9, "f": [10, 11]}, \
        {"a": [1, 2], "b": [3, 4], "c": 6, "d": 7, "e": 9, "f": [10, 11]}, \
        {"a": [1, 2], "b": [3, 4], "c": 6, "d": 8, "e": 9, "f": [10, 11]}]
        self.assertEqual(list(OrbgenParamSetGenerator.nested_parse(nested, excluded)), normal)

        nested = {"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": [7, 8], "e|f": [[9, 10], [11, 12]]}
        excluded = ["a", "b", "c"]
        normal = [{"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": 7, "e": 9, "f": 10}, \
        {"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": 7, "e": 11, "f": 12}, \
        {"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": 8, "e": 9, "f": 10}, \
        {"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": 8, "e": 11, "f": 12}]
        self.assertEqual(list(OrbgenParamSetGenerator.nested_parse(nested, excluded)), normal)

    def test_recomb_orbset(self):
        unpacked = [{"nzeta": "SZ", "struset": "dimer"}, 
                    {"nzeta": "DZP", "struset": "trimer"}, 
                    {"nzeta": "TZDP", "struset": "dimer"},
                    {"nzeta": "TZDP", "struset": "trimer"}]
        result = [({"nzeta": "SZ", "struset": "dimer"}, 
                   {"nzeta": "DZP", "struset": "trimer"}, 
                   {"nzeta": "TZDP", "struset": "dimer"}),
                  ({"nzeta": "SZ", "struset": "dimer"}, 
                   {"nzeta": "DZP", "struset": "trimer"}, 
                   {"nzeta": "TZDP", "struset": "trimer"})]
        self.assertEqual(list(OrbgenParamSetGenerator.recomb_orbset(unpacked)), result)

    def test_on_call(self):
        #return
        import json
        import uuid
        import os
        test = {"orbsets": [{"atomset": 0, "spillset": 0, "nzeta|struset|band_range|deps": [
                ["SZ", 0, 1.0, None], ["DZP", 0, [1.0, 1.5, 2.0], 0], ["TZDP", 1, [1.0, 1.5, 2.0], 1]]}],
                "atomsets": [{"Be": [["PD"], None]}],
                "spillsets": [{"spill_coefs": [0.0, 1.0], "optimizer": ["bfgs", "pytorch.SWAT"]}],
                "strusets": [{"shape": "dimer", "nbands": "auto", "nspin": 1, "bond_lengths": "auto"},
                             {"shape": "trimer", "nbands": "auto", "nspin": 1, "bond_lengths": "auto"}]}
        finp = "test_on_call_" + str(uuid.uuid4()) + ".json"
        with open(finp, 'w') as f:
            json.dump(test, f)
        opsg = OrbgenParamSetGenerator(finp, "./download/pseudopotentials/", "./apns_cache/")
        os.remove(finp)

        orbgen = opsg()
        #print(list(orbgen))

if __name__ == "__main__":
    unittest.main()