"""APNS IO module for Python."""

def check(inp: dict):
    assert isinstance(inp, dict), 'Input must be a dictionary.'
    assert 'global' in inp, 'Global settings must be provided.'
    mode = inp['global'].get('mode', 'test')
    if mode == 'test':
        assert 'abacus' in inp or 'qespresso' in inp, 'At least one of the two DFT calculator must be provided.'
        assert isinstance(inp['abacus'], list) or isinstance(inp['qespresso'], list), 'DFT settings must be a list.'
        assert all([isinstance(i, dict) for i in inp['abacus']]) or all([isinstance(i, dict) for i in inp['qespresso']]), 'DFT settings must be a dictionary.'
        assert 'atomsets' in inp, 'Atoms must be defined in atomsets section.'
        assert isinstance(inp['atomsets'], list), 'Atomsets must be a list.'
        assert all([isinstance(i, dict) for i in inp['atomsets']]), 'Organize pptags and naotags information in a dictionary for each element.'
        assert 'structures' in inp, 'Structures must be defined in structures section.'
        assert isinstance(inp['structures'], list), 'Structures must be provided as a list.'
        assert all([isinstance(i, dict) for i in inp['structures']]), 'Structure settings should be organized as a dict.'
        assert all(['calculator' in i for i in inp['structures']]), 'Calculator must be defined for each structure.'
        assert all(['calcset' in i for i in inp['structures']]), 'DFT calculator setting must be defined for each structure.'
        assert all(['atomset' in i for i in inp['structures']]), 'Atomset must be defined for each structure.'
        assert all(['desc' in i for i in inp['structures']]), 'Structure descriptor must be defined for each structure.'
        assert all([isinstance(i['desc'], list) for i in inp['structures']]), 'Structure descriptor should be organized in one list.'
        assert all([all([len(j) == 3 for j in i['desc']]) for i in inp['structures']]), 'Structure descriptor must be a list of 3 elements.'
        assert all([j[0] in ["mp", "scratch", "local"] for i in inp['structures'] for j in i['desc']]), 'Structure type can only be "mp", "scratch" or "local".'
        assert all([isinstance(j[1], str) for i in inp['structures'] for j in i['desc']]), 'Structure must be specified as a string.'
        assert all([all([isinstance(d[2], list) and all([isinstance(scale, float) for scale in d[2]]) for d in s['desc']]) for s in inp['structures']])
        
def read(finp):
    import json
    with open(finp, 'r') as f:
        data = json.load(f)
    check(data)
    return data

import unittest
class TestAPNSIO(unittest.TestCase):
    def test_check(self):
        test = {
            "global": {
                "mode": "test",
                "pseudo_dir": "./download/pseudopotentials",
                "orbital_dir": "./download/numerical_orbitals"
            },
            "abacus": [
                {
                    "ecutwfc": 30,
                    "calculation": "cell-relax",
                    "nspin": 2
                },
                {
                    "ecutwfc": [30, 40, 50],
                    "calculation": "scf",
                    "cal_force": 1,
                    "cal_stress": 1
                },
                {
                    "ecutwfc|ecutrho": [[30, 300], [40, 400], [50, 500]],
                    "calculation": "scf"
                },
                {
                    "ecutwfc": [30, 40, 50],
                    "ecutrho": [300, 400, 500],
                    "calculation": "scf"
                }
            ],
            "atomsets": [
                {
                    "H": [[],
                        []],
                    "O": [[],
                        []]
                },
                {
                    "H": [[],
                        []],
                    "O": [[],
                        []],
                    "Na": [[],
                        []]
                }
            ],
            "structures": [
                {
                    "calculator": "abacus",
                    "calcset": 3,
                    "atomset": 2,
                    "desc": [["mp", "Li3ErCl6", [1.0]],
                             ["scratch", "ErCl_xy3", [0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06]],
                             ["scratch", "Er_dimer", [2.1, 2.4, 2.7, 3.0, 3.3]]]
                },
                {
                    "calculator": "abacus",
                    "calcset": 1,
                    "atomset": 0,
                    "desc": [["mp", "Li3DyCl6", [0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06]],
                             ["local", "./apns_cache/Li3LuCl6-Pna21.cif", [1.0]]]
                }
            ]
        }
        # will not raise any exception
        check(test)
        # wrong input type
        with self.assertRaises(Exception):
            check(1)
        # wrong input type
        with self.assertRaises(Exception):
            check("test")
        # wrong input type
        with self.assertRaises(Exception):
            check(1.0)
        # wrong input content
        with self.assertRaises(Exception):
            check({
                "global": {
                    "mode": "test",
                    "pseudo_dir": "./download/pseudopotentials",
                    "orbital_dir": "./download/numerical_orbitals"
                }
            })
        # wrong input content
        with self.assertRaises(Exception):
            check({
                "global": {
                    "mode": "test",
                    "pseudo_dir": "./download/pseudopotentials",
                    "orbital_dir": "./download/numerical_orbitals"
                },
                "abacus": {
                    "ecutwfc": 30,
                    "calculation": "cell-relax",
                    "nspin": 2
                }
            })
        # wrong input content
        with self.assertRaises(Exception):
            check({
                "global": {
                    "mode": "test",
                    "pseudo_dir": "./download/pseudopotentials",
                    "orbital_dir": "./download/numerical_orbitals"
                },
                "abacus": [
                    {
                        "ecutwfc": 30,
                        "calculation": "cell-relax",
                        "nspin": 2
                    }
                ]
            })
        # wrong input content
        with self.assertRaises(Exception):
            check({
                "global": {
                    "mode": "test",
                    "pseudo_dir": "./download/pseudopotentials",
                    "orbital_dir": "./download/numerical_orbitals"
                },
                "abacus": [
                    {
                        "ecutwfc": 30,
                        "calculation": "cell-relax",
                        "nspin": 2
                    },
                    {
                        "ecutwfc": 30,
                        "calculation": "cell-relax",
                        "nspin": 2
                    }
                ]
            })
        # wrong input content
        with self.assertRaises(Exception):
            check({
                "global": {
                    "mode": "test",
                    "pseudo_dir": "./download/pseudopotentials",
                    "orbital_dir": "./download/numerical_orbitals"
                },
                "abacus": [
                    {
                        "ecutwfc": 30,
                        "calculation": "cell-relax",
                        "nspin": 2
                    },
                    {
                        "ecutwfc": 30,
                        "calculation": "cell-relax",
                        "nspin": 2
                    }
                ],
                "atomsets": {}
            })
        # wrong input content
        with self.assertRaises(Exception):
            check({
                "global": {
                    "mode": "test",
                    "pseudo_dir": "./download/pseudopotentials",
                    "orbital_dir": "./download/numerical_orbitals"
                },
                "abacus": [
                    {
                        "ecutwfc": 30,
                        "calculation": "cell-relax",
                        "nspin": 2
                    },
                    {
                        "ecutwfc": 30,
                        "calculation": "cell-relax",
                        "nspin": 2
                    }
                ],
                "atomsets": []
            })

if __name__ == '__main__':
    unittest.main()