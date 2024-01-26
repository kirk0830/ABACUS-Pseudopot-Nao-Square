import unittest
import apns.module_structure.CifParser_Pymatgen as amscp
import math

def listvec2norm(vec):
    return math.sqrt(sum([x**2 for x in vec]))

def listvecdot(vec1, vec2):
    return sum([vec1[i]*vec2[i] for i in range(3)])

def listvecangle(vec1, vec2):
    return math.acos(listvecdot(vec1, vec2)/listvec2norm(vec1)/listvec2norm(vec2))

class TestCifParserPymatgen(unittest.TestCase):

    def test_structure(self):
        fname = "./apns/module_structure/test/support/mp-8.cif"
        elements, positions = amscp.structure(fname)
        self.assertListEqual(elements, ['Re', 'Re'])
        ref = [[0.33333300, 0.66666700, 0.25000000],
                [0.66666700, 0.33333300, 0.75000000]]
        for ia in range(2):
            for i in range(3):
                self.assertAlmostEqual(positions[ia][i], ref[ia][i], places=6)
    
        positions = amscp.structure(fname, as_dict=True)
        for ia in range(2):
            for i in range(3):
                self.assertAlmostEqual(positions["Re"][ia][i], ref[ia][i], places=6)

        fname = "./apns/module_structure/test/support/mp-13.cif"
        elements, positions = amscp.structure(fname)
        ref = [[0.50000002, 0.50000008, 0.49999999],
               [0.99999998, 0.99999992, 0.00000001]]
        for ia in range(2):
            for i in range(3):
                self.assertAlmostEqual(positions[ia][i], ref[ia][i], places=6)
        
        positions = amscp.structure(fname, as_dict=True)
        for ia in range(2):
            for i in range(3):
                self.assertAlmostEqual(positions["Fe"][ia][i], ref[ia][i], places=6)

    def test_lattice(self):
        fname = "./apns/module_structure/test/support/mp-8.cif"
        a, b, c, alpha, beta, gamma, lattice_vectors = amscp.lattice(fname)
        self.assertAlmostEqual(a, 2.76533582, places=4)
        self.assertAlmostEqual(b, 2.76533582, places=4)
        self.assertAlmostEqual(c, 4.47357663, places=4)
        self.assertAlmostEqual(alpha, 90.0, places=4)
        self.assertAlmostEqual(beta, 90.0, places=4)
        self.assertAlmostEqual(gamma, 120.0, places=4)
        _a = listvec2norm(lattice_vectors[0])
        _b = listvec2norm(lattice_vectors[1])
        _c = listvec2norm(lattice_vectors[2])
        self.assertEqual(_a, a)
        self.assertEqual(_b, b)
        self.assertEqual(_c, c)

        angle_ab = listvecangle(lattice_vectors[0], lattice_vectors[1])
        angle_bc = listvecangle(lattice_vectors[1], lattice_vectors[2])
        angle_ca = listvecangle(lattice_vectors[2], lattice_vectors[0])
        self.assertEqual(angle_ab, math.radians(gamma))
        self.assertEqual(angle_bc, math.radians(alpha))
        self.assertEqual(angle_ca, math.radians(beta))
    
        fname = "./apns/module_structure/test/support/mp-13.cif"
        a, b, c, alpha, beta, gamma, lattice_vectors = amscp.lattice(fname)
        self.assertAlmostEqual(a, 2.47781303, places=4)
        self.assertAlmostEqual(b, 2.47780913, places=4)
        self.assertAlmostEqual(c, 4.05429521, places=4)
        self.assertAlmostEqual(alpha, 90.00003845, places=4)
        self.assertAlmostEqual(beta, 90.00031951, places=4)
        self.assertAlmostEqual(gamma, 109.46983668, places=4)
        _a = listvec2norm(lattice_vectors[0])
        _b = listvec2norm(lattice_vectors[1])
        _c = listvec2norm(lattice_vectors[2])
        self.assertEqual(_a, a)
        self.assertEqual(_b, b)
        self.assertEqual(_c, c)

        angle_ab = listvecangle(lattice_vectors[0], lattice_vectors[1])
        angle_bc = listvecangle(lattice_vectors[1], lattice_vectors[2])
        angle_ca = listvecangle(lattice_vectors[2], lattice_vectors[0])
        self.assertEqual(angle_ab, math.radians(gamma))
        self.assertEqual(angle_bc, math.radians(alpha))
        self.assertEqual(angle_ca, math.radians(beta))

if __name__ == "__main__":
    unittest.main()