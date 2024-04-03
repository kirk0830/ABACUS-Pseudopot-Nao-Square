import unittest
import apns.module_structure.basic as amsb

class TestBasic(unittest.TestCase):
    def test_basic(self):
        """magnetization is given atom-wise"""
        self.assertEqual(amsb.expand_atomic_species(
            symbols=["H", "H", "O"],
            pseudopotentials={
                "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
                "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF"
            },
            starting_magnetization=[0.0, 0.0, 0.0],
            atomic_positions=[
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 2.0]
            ]
        ), {
            "H1": {
                "pseudopotential": "H.pbe-kjpaw_psl.1.0.0.UPF",
                "starting_magnetization": 0.0,
                "atomic_positions": [
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0]
                ]
            },
            "O2": {
                "pseudopotential": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
                "starting_magnetization": 0.0,
                "atomic_positions": [
                    [0.0, 0.0, 2.0]
                ]
            }
        })
        """magnetization is given atom-type-wise"""
        self.assertEqual(amsb.expand_atomic_species(
            symbols=["H", "H", "O"],
            pseudopotentials={
                "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
                "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF"
            },
            starting_magnetization={
                "H": 0.0,
                "O": 1.2345
            },
            atomic_positions=[
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 2.0]
            ]
        ), {
            "H": {
                "pseudopotential": "H.pbe-kjpaw_psl.1.0.0.UPF",
                "starting_magnetization": 0.0,
                "atomic_positions": [
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0]
                ]
            },
            "O": {
                "pseudopotential": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
                "starting_magnetization": 1.2345,
                "atomic_positions": [
                    [0.0, 0.0, 2.0]
                ]
            }
        })
        """magnetization is given atom-type-wise, but it is Fe2 anti-ferromagnetic"""
        self.assertEqual(amsb.expand_atomic_species(
            symbols=["Fe", "Fe", "O"],
            pseudopotentials={
                "Fe": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
                "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF"
            },
            starting_magnetization=[2.0, -2.0, 0.0],
            atomic_positions=[
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 2.0]
            ]
        ), {
            "Fe1": {
                "pseudopotential": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
                "starting_magnetization": 2.0,
                "atomic_positions": [
                    [0.0, 0.0, 0.0]
                ]
            },
            "Fe2": {
                "pseudopotential": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
                "starting_magnetization": -2.0,
                "atomic_positions": [
                    [0.0, 0.0, 1.0]
                ]
            },
            "O3": {
                "pseudopotential": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
                "starting_magnetization": 0.0,
                "atomic_positions": [
                    [0.0, 0.0, 2.0]
                ]
            }
        })

    def test_init_magmom(self):

        magnetization = amsb.init_magmom("mp-35.cif", magnetism="materials_project")
        self.assertListEqual(magnetization, [4.1, -3.6, -3.6, -3.6, -3.6, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3])
        magnetization = amsb.init_magmom("dimer", magnetism="default")
        self.assertListEqual(magnetization, [0.0, 0.0])
        magnetization = amsb.init_magmom("dimer", magnetism="ferromagnetic")
        self.assertListEqual(magnetization, [1.0, 1.0])
        magnetization = amsb.init_magmom("dimer", magnetism="antiferromagnetic")
        self.assertListEqual(magnetization, [1.0, -1.0])
        magnetization = amsb.init_magmom("trimer", magnetism="default")
        self.assertListEqual(magnetization, [0.0, 0.0, 0.0])
        magnetization = amsb.init_magmom("trimer", magnetism="ferromagnetic")
        self.assertListEqual(magnetization, [1.0, 1.0, 1.0])
        magnetization = amsb.init_magmom("trimer", magnetism="antiferromagnetic")
        self.assertListEqual(magnetization, [0.0, 0.0, 0.0])
        magnetization = amsb.init_magmom("tetramer", magnetism="default")
        self.assertListEqual(magnetization, [0.0, 0.0, 0.0, 0.0])
        magnetization = amsb.init_magmom("tetramer", magnetism="ferromagnetic")
        self.assertListEqual(magnetization, [1.0, 1.0, 1.0, 1.0])
        magnetization = amsb.init_magmom("tetramer", magnetism="antiferromagnetic")
        self.assertListEqual(magnetization, [1.0, -1.0, 1.0, -1.0])

if __name__ == "__main__":
    unittest.main()