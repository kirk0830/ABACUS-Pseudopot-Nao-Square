class Calculator:
    """Density functional theory based calculator (like ABACUS, Quantum ESPRESSO, VASP, etc.)
    input generation"""
    settings = []
    def __init__(self) -> None:
        pass

    def build(self, calcset: list):
        """build/extend to all settings: for any key has value as list, convert to list of 
        inputs like
        ```python
        {
            "ecutwfc": [100, 200, 300],
            "nspin": 2
        }
        ```
        convert to
        ```python
        [
            {"ecutwfc": 100, "nspin": 2},
            {"ecutwfc": 200, "nspin": 2},
            {"ecutwfc": 300, "nspin": 2}
        ]
        ```
        More complex example:
        ```python
        {
            "ecutwfc": [100, 200, 300],
            "ecutrho": [800, 1600],
            "nspin": 2
        }
        ```
        convert to
        ```python
        [
            {"ecutwfc": 100, "ecutrho": 800, "nspin": 2},
            {"ecutwfc": 100, "ecutrho": 1600, "nspin": 2},
            {"ecutwfc": 200, "ecutrho": 800, "nspin": 2},
            {"ecutwfc": 200, "ecutrho": 1600, "nspin": 2},
            {"ecutwfc": 300, "ecutrho": 800, "nspin": 2},
            {"ecutwfc": 300, "ecutrho": 1600, "nspin": 2}
        ]
        ```
        """
        assert isinstance(calcset, list), "calcset must be a list"
        assert all([isinstance(i, dict) for i in calcset]), "calcset must be a list of dictionaries"
        import itertools as it
        self.settings = []
        for cs in calcset:
            keys = [k for k in cs.keys() \
                    if ("|" not in k and isinstance(cs[k], list))\
                    or ("|" in k and isinstance(cs[k], list) \
                        and all([isinstance(i, list) and len(i) == k.count("|") + 1 for i in cs[k]]))]
            vals = [cs[k] for k in keys]
            res = {k: v for k, v in cs.items() if k not in keys}
            self.settings.append([Calculator.disjoint_keywords(res|dict(zip(keys, v))) for v in it.product(*vals)])

    def disjoint_keywords(calcset: dict):
        """disjoint keywords in one calcset, like if there is one/or more keywords appear like
        ```json
        {
            "ecutwfc|ecutrho": [100, 800],
        }
        ```
        then convert to
        ```json
        {
            "ecutwfc": 100,
            "ecutrho": 800
        }
        ```
        """
        assert isinstance(calcset, dict), "calcset must be a dictionary"
        disjoint = {}
        for k, v in calcset.items():
            if '|' in k:
                ks, vs = k.split('|'), v
                assert len(ks) == len(vs), "number of keywords and values must be the same"
                disjoint |= dict(zip(ks, vs))
            else:
                disjoint[k] = v
        return disjoint

class ABACUS(Calculator):
    def __init__(self) -> None:
        super().__init__()

class QuantumESPRESSO(Calculator):
    def __init__(self) -> None:
        super().__init__()

    def flatten(calcset: dict):
        """flatten one calcset, from nested dict to only one-layer dict by:
        input:
        ```bash
        &system
        calculation = 'relax'
        prefix = 'calc'
        outdir = './'
        wf_collect = .true.
        /
        &electrons
        mixing_beta = 0.7
        conv_thr = 1.0d-8
        /
        &ions
        ion_dynamics = 'bfgs'
        /
        &cell
        press = 0.0
        press_conv_thr = 0.5
        /
        ```
        input as
        ```json
        {
            "system": {
                "calculation": "relax",
                "prefix": "calc",
                "outdir": "./",
                "wf_collect": ".true."
            },
            "electrons": {
                "mixing_beta": 0.7,
                "conv_thr": 1.0e-8
            },
            "ions": {
                "ion_dynamics": "bfgs"
            },
            "cell": {
                "press": 0.0,
                "press_conv_thr": 0.5
            }
        }
        ```
        convert to
        ```json
        {
            "system.calculation": "relax",
            "system.prefix": "calc",
            "system.outdir": "./",
            "system.wf_collect": ".true.",
            "electrons.mixing_beta": 0.7,
            "electrons.conv_thr": 1.0e-8,
            "ions.ion_dynamics": "bfgs",
            "cell.press": 0.0,
            "cell.press_conv_thr": 0.5
        }
        /
        """
        assert isinstance(calcset, dict), "calcset must be a dictionary"
        flat = {}
        for k, v in calcset.items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    flat["|".join([f"{k}.{kk_}" for kk_ in kk.split("|")])] = vv
            else:
                flat[k] = v
        return flat
    def make_tree(calcset: dict):
        """do the thing in opposite way as `flatten`"""
        assert isinstance(calcset, dict), "calcset must be a dictionary"
        tree = {}
        for k, v in calcset.items():
            if '.' in k:
                k1, k2 = k.split('.')
                if k1 not in tree:
                    tree[k1] = {}
                tree[k1][k2] = v
            else:
                tree[k] = v
        return tree
    
    def build(self, calcset: list):
        """override the build method because QuantumESPRESSO need to flatten the input"""
        assert isinstance(calcset, list), "calcset must be a list"
        assert all([isinstance(i, dict) for i in calcset]), "calcset must be a list of dictionaries"
        calcset = [QuantumESPRESSO.flatten(i) for i in calcset]
        super().build(calcset)
        self.settings = [[QuantumESPRESSO.make_tree(j) for j in i] for i in self.settings]

import unittest
class TestCalculator(unittest.TestCase):
    def test_disjoint_keywords(self):
        self.assertEqual(Calculator.disjoint_keywords(
            {"ecutwfc|ecutrho": [100, 800]}), {"ecutwfc": 100, "ecutrho": 800})
    def test_flatten(self):
        self.assertEqual(QuantumESPRESSO.flatten(
            {
                "system": {
                    "calculation": "relax",
                    "prefix": "calc",
                    "outdir": "./",
                    "wf_collect": ".true."
                },
                "electrons": {
                    "mixing_beta": 0.7,
                    "conv_thr": 1.0e-8
                },
                "ions": {
                    "ion_dynamics": "bfgs"
                },
                "cell": {
                    "press": 0.0,
                    "press_conv_thr": 0.5
                }
            }
        ), {
            "system.calculation": "relax",
            "system.prefix": "calc",
            "system.outdir": "./",
            "system.wf_collect": ".true.",
            "electrons.mixing_beta": 0.7,
            "electrons.conv_thr": 1.0e-8,
            "ions.ion_dynamics": "bfgs",
            "cell.press": 0.0,
            "cell.press_conv_thr": 0.5
        })
        # test flatten with joint case
        self.assertEqual(QuantumESPRESSO.flatten(
            {
                "system": {
                    "calculation": "relax",
                    "prefix": "calc",
                    "outdir": "./",
                    "wf_collect": ".true.",
                    "ecutwfc|ecutrho": [100, 800]
                },
                "electrons": {
                    "mixing_beta": 0.7,
                    "conv_thr": 1.0e-8
                },
                "ions": {
                    "ion_dynamics": "bfgs"
                },
                "cell": {
                    "press": 0.0,
                    "press_conv_thr": 0.5
                }
            }
        ), {
            "system.calculation": "relax",
            "system.prefix": "calc",
            "system.outdir": "./",
            "system.wf_collect": ".true.",
            "system.ecutwfc|system.ecutrho": [100, 800],
            "electrons.mixing_beta": 0.7,
            "electrons.conv_thr": 1.0e-8,
            "ions.ion_dynamics": "bfgs",
            "cell.press": 0.0,
            "cell.press_conv_thr": 0.5
        })
    def test_make_tree(self):
        self.assertEqual(QuantumESPRESSO.make_tree(
            {
                "system.calculation": "relax",
                "system.prefix": "calc",
                "system.outdir": "./",
                "system.wf_collect": ".true.",
                "electrons.mixing_beta": 0.7,
                "electrons.conv_thr": 1.0e-8,
                "ions.ion_dynamics": "bfgs",
                "cell.press": 0.0,
                "cell.press_conv_thr": 0.5
            }
        ), {
            "system": {
                "calculation": "relax",
                "prefix": "calc",
                "outdir": "./",
                "wf_collect": ".true."
            },
            "electrons": {
                "mixing_beta": 0.7,
                "conv_thr": 1.0e-8
            },
            "ions": {
                "ion_dynamics": "bfgs"
            },
            "cell": {
                "press": 0.0,
                "press_conv_thr": 0.5
            }
        })
    def test_build(self):
        calc = Calculator()
        calc.build([{"ecutwfc": [100, 200, 300], "ecutrho": [800, 1600], "nspin": 2}])
        self.assertEqual(calc.settings, [
            [
                {"ecutwfc": 100, "ecutrho": 800, "nspin": 2},
                {"ecutwfc": 100, "ecutrho": 1600, "nspin": 2},
                {"ecutwfc": 200, "ecutrho": 800, "nspin": 2},
                {"ecutwfc": 200, "ecutrho": 1600, "nspin": 2},
                {"ecutwfc": 300, "ecutrho": 800, "nspin": 2},
                {"ecutwfc": 300, "ecutrho": 1600, "nspin": 2}
            ]
        ])
        # test joint case
        calc.build([{"ecutwfc|ecutrho": [100, 800], "nspin": 2}])
        self.assertEqual(calc.settings, [
            [
                {"ecutwfc": 100, "ecutrho": 800, "nspin": 2}
            ]
        ])
        # joint plus direct product
        calc.build([{"ecutwfc|ecutrho": [[100, 200], [800, 1600]], "nspin": [1, 2]}])
        self.assertEqual(calc.settings, [
            [
                {"ecutwfc": 100, "ecutrho": 200, "nspin": 1},
                {"ecutwfc": 100, "ecutrho": 200, "nspin": 2},
                {"ecutwfc": 800, "ecutrho": 1600, "nspin": 1},
                {"ecutwfc": 800, "ecutrho": 1600, "nspin": 2}
            ]
        ])
        # test quantum espresso
        qe = QuantumESPRESSO()
        qe.build([{
            "system": {
                "calculation": "relax",
                "prefix": "calc",
                "outdir": "./",
                "wf_collect": ".true."
            },
            "electrons": {
                "mixing_beta": 0.7,
                "conv_thr": 1.0e-8
            },
            "ions": {
                "ion_dynamics": "bfgs"
            },
            "cell": {
                "press": 0.0,
                "press_conv_thr": 0.5
            }
        }])
        self.assertEqual(qe.settings, [
            [
                {
                    "system": {
                        "calculation": "relax",
                        "prefix": "calc",
                        "outdir": "./",
                        "wf_collect": ".true."
                    },
                    "electrons": {
                        "mixing_beta": 0.7,
                        "conv_thr": 1.0e-8
                    },
                    "ions": {
                        "ion_dynamics": "bfgs"
                    },
                    "cell": {
                        "press": 0.0,
                        "press_conv_thr": 0.5
                    }
                }
            ]
        ])
        # test qe with both joint and direct product
        qe.build([{
            "system": {
                "calculation": "relax",
                "prefix": "calc",
                "outdir": "./",
                "wf_collect": ".true.",
                "ecutwfc|ecutrho": [[100, 800], [200, 1600]]
            },
            "electrons": {
                "mixing_beta": [0.7, 0.8],
                "conv_thr": 1.0e-8
            },
            "ions": {
                "ion_dynamics": "bfgs"
            },
            "cell": {
                "press": 0.0,
                "press_conv_thr": 0.5
            }
        }])
        self.assertEqual(qe.settings, [
            [
                {
                    "system": {
                        "calculation": "relax",
                        "prefix": "calc",
                        "outdir": "./",
                        "wf_collect": ".true.",
                        "ecutwfc": 100,
                        "ecutrho": 800
                    },
                    "electrons": {
                        "conv_thr": 1.0e-8,
                        "mixing_beta": 0.7,
                    },
                    "ions": {
                        "ion_dynamics": "bfgs"
                    },
                    "cell": {
                        "press": 0.0,
                        "press_conv_thr": 0.5
                    }
                },
                {
                    "system": {
                        "calculation": "relax",
                        "prefix": "calc",
                        "outdir": "./",
                        "wf_collect": ".true.",
                        "ecutwfc": 100,
                        "ecutrho": 800
                    },
                    "electrons": {
                        "conv_thr": 1.0e-8,
                        "mixing_beta": 0.8,
                    },
                    "ions": {
                        "ion_dynamics": "bfgs"
                    },
                    "cell": {
                        "press": 0.0,
                        "press_conv_thr": 0.5
                    }
                },
                {
                    "system": {
                        "calculation": "relax",
                        "prefix": "calc",
                        "outdir": "./",
                        "wf_collect": ".true.",
                        "ecutwfc": 200,
                        "ecutrho": 1600
                    },
                    "electrons": {
                        "conv_thr": 1.0e-8,
                        "mixing_beta": 0.7,
                    },
                    "ions": {
                        "ion_dynamics": "bfgs"
                    },
                    "cell": {
                        "press": 0.0,
                        "press_conv_thr": 0.5
                    }
                },
                {
                    "system": {
                        "calculation": "relax",
                        "prefix": "calc",
                        "outdir": "./",
                        "wf_collect": ".true.",
                        "ecutwfc": 200,
                        "ecutrho": 1600
                    },
                    "electrons": {
                        "conv_thr": 1.0e-8,
                        "mixing_beta": 0.8,
                    },
                    "ions": {
                        "ion_dynamics": "bfgs"
                    },
                    "cell": {
                        "press": 0.0,
                        "press_conv_thr": 0.5
                    }
                }
            ]
        ])
        # test ABACUS
        abacus = ABACUS()
        abacus.build([{
            "ecutwfc": [100, 200, 300],
            "ecutrho": [800, 1600],
            "nspin": 2
        }])
        self.assertEqual(abacus.settings, [
            [
                {"ecutwfc": 100, "ecutrho": 800, "nspin": 2},
                {"ecutwfc": 100, "ecutrho": 1600, "nspin": 2},
                {"ecutwfc": 200, "ecutrho": 800, "nspin": 2},
                {"ecutwfc": 200, "ecutrho": 1600, "nspin": 2},
                {"ecutwfc": 300, "ecutrho": 800, "nspin": 2},
                {"ecutwfc": 300, "ecutrho": 1600, "nspin": 2}
            ]
        ])
        # test joint case
        abacus.build([{"ecutwfc|ecutrho": [100, 800], "nspin": 2}])
        self.assertEqual(abacus.settings, [
            [
                {"ecutwfc": 100, "ecutrho": 800, "nspin": 2}
            ]
        ])
        # joint plus direct product
        abacus.build([{"ecutwfc|ecutrho": [[100, 200], [800, 1600]], "nspin": [1, 2]}])
        self.assertEqual(abacus.settings, [
            [
                {"ecutwfc": 100, "ecutrho": 200, "nspin": 1},
                {"ecutwfc": 100, "ecutrho": 200, "nspin": 2},
                {"ecutwfc": 800, "ecutrho": 1600, "nspin": 1},
                {"ecutwfc": 800, "ecutrho": 1600, "nspin": 2}
            ]
        ])

if __name__ == "__main__":
    unittest.main()