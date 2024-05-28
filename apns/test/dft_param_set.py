class DFTParamSetGenerator:

    folded = None
    keys, vals = None, None
    residual = None
    def __init__(self, folded: dict) -> None:
        self.folded = DFTParamSetGenerator.flatten(folded)
        keys = [k for k, v in self.folded.items()\
            if ("|" not in k and isinstance(v, list))\
                or ("|" in k and isinstance(v, list)\
                    and all([isinstance(_v, list) and len(_v) == k.count("|") + 1 for _v in v]))]
        self.keys = keys
        self.vals = [self.folded[k] for k in keys]
        self.residual = {k: v for k, v in self.folded.items() if k not in keys}
        print(f"""DFTParamSetGenerator setup:
Keys to be expanded: {", ".join(keys)}
""")

    def __call__(self):
        import itertools as it
        for genvals in it.product(*self.vals): # generator provided by itertools.product
            genkv = dict(zip(self.keys, genvals)) | self.residual
            # joint expansion
            dejoint = {}
            for k, v in genkv.items():
                dejoint |= dict(zip(k.split("|"), v)) if "|" in k else {k: v}
            yield DFTParamSetGenerator.make_tree(dejoint)

    def flatten(param_set: dict):
        """flatten one param_set, from nested dict to only one-layer dict by:
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
        assert isinstance(param_set, dict), "param_set must be a dictionary"
        flat = {}
        for k, v in param_set.items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    flat["|".join([f"{k}.{kk_}" for kk_ in kk.split("|")])] = vv
            else:
                flat[k] = v
        return flat

    def make_tree(param_set: dict):
        """do the thing in opposite way as `flatten`"""
        assert isinstance(param_set, dict), "param_set must be a dictionary"
        tree = {}
        for k, v in param_set.items():
            if '.' in k:
                k1, k2 = k.split('.')
                if k1 not in tree:
                    tree[k1] = {}
                tree[k1][k2] = v
            else:
                tree[k] = v
        return tree

import unittest
class TestABACUSParamSetGenerator(unittest.TestCase):
    def test_flatten(self):
        self.assertEqual(DFTParamSetGenerator.flatten(
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
        self.assertEqual(DFTParamSetGenerator.flatten(
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
        self.assertEqual(DFTParamSetGenerator.make_tree(
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

    def test_oncall(self):
        dftgen = DFTParamSetGenerator({"ecutwfc": [100, 200, 300], "ecutrho": [800, 1600], "nspin": 2})
        result = list(dftgen())
        ref = [
            {"ecutwfc": 100, "ecutrho": 800, "nspin": 2},
            {"ecutwfc": 100, "ecutrho": 1600, "nspin": 2},
            {"ecutwfc": 200, "ecutrho": 800, "nspin": 2},
            {"ecutwfc": 200, "ecutrho": 1600, "nspin": 2},
            {"ecutwfc": 300, "ecutrho": 800, "nspin": 2},
            {"ecutwfc": 300, "ecutrho": 1600, "nspin": 2}
        ]
        self.assertEqual(len(result), len(ref))
        for test in ref:
            self.assertTrue(test in result)

        dftgen = DFTParamSetGenerator({"ecutwfc|ecutrho": [100, 800], "nspin": 2})
        result = list(dftgen())
        self.assertEqual(result, [{"ecutwfc": 100, "ecutrho": 800, "nspin": 2}])

        dftgen = DFTParamSetGenerator({"ecutwfc|ecutrho": [[100, 200], [800, 1600]], "nspin": [1, 2]})
        result = list(dftgen())
        ref = [
            {"ecutwfc": 100, "ecutrho": 200, "nspin": 1},
            {"ecutwfc": 100, "ecutrho": 200, "nspin": 2},
            {"ecutwfc": 800, "ecutrho": 1600, "nspin": 1},
            {"ecutwfc": 800, "ecutrho": 1600, "nspin": 2}
        ]
        self.assertEqual(len(result), len(ref))
        for test in ref:
            self.assertTrue(test in result)

        dftgen = DFTParamSetGenerator({
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
        result = list(dftgen())
        ref = [{
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
        }]
        self.assertEqual(result, ref)

        dftgen = DFTParamSetGenerator({
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
        })
        ref = [
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
        result = list(dftgen())
        self.assertEqual(len(result), len(ref))
        for test in ref:
            self.assertTrue(test in result)

if __name__ == "__main__":
    unittest.main()