import unittest
import apns.module_database.database as db

class TestDatabase(unittest.TestCase):

    def test_get_element_label(self):
        self.assertEqual(db.get_element_label(1), "1_H")
        self.assertEqual(db.get_element_label(2), "2_He")
        self.assertEqual(db.get_element_label(3), "3_Li")
        self.assertEqual(db.get_element_label(4), "4_Be")
        self.assertEqual(db.get_element_label(5), "5_B")

    def test_get_element_index(self):
        self.assertEqual(db.get_element_index("H"), 1)
        self.assertEqual(db.get_element_index("He"), 2)
        self.assertEqual(db.get_element_index("Li"), 3)
        self.assertEqual(db.get_element_index("Be"), 4)
        self.assertEqual(db.get_element_index("B"), 5)

    def test_sublayer_to_l(self):
        self.assertEqual(db.sublayer_to_l("s"), 0)
        self.assertEqual(db.sublayer_to_l("p"), 1)
        self.assertEqual(db.sublayer_to_l("d"), 2)
        self.assertEqual(db.sublayer_to_l("f"), 3)
        self.assertEqual(db.sublayer_to_l("g"), 4)

    def test_l_to_sublayer(self):
        self.assertEqual(db.l_to_sublayer(0), "s")
        self.assertEqual(db.l_to_sublayer(1), "p")
        self.assertEqual(db.l_to_sublayer(2), "d")
        self.assertEqual(db.l_to_sublayer(3), "f")
        self.assertEqual(db.l_to_sublayer(4), "g")

    def test_element_mass(self):
        self.assertEqual(db.element_mass("H"), 1.00794)
        self.assertEqual(db.element_mass("He"), 4.002602)
        self.assertEqual(db.element_mass("Li"), 6.941)
        self.assertEqual(db.element_mass("Be"), 9.012182)
        self.assertEqual(db.element_mass("B"), 10.811)

    def test_number_to_multiplicity(self):
        self.assertEqual(db.number_to_multiplicity(1), "s")
        self.assertEqual(db.number_to_multiplicity(2), "d")
        self.assertEqual(db.number_to_multiplicity(3), "t")
        self.assertEqual(db.number_to_multiplicity(4), "q")
    
    def test_multiplicity_to_number(self):
        self.assertEqual(db.multiplicity_to_number("s"), 1)
        self.assertEqual(db.multiplicity_to_number("d"), 2)
        self.assertEqual(db.multiplicity_to_number("t"), 3)
        self.assertEqual(db.multiplicity_to_number("q"), 4)

if __name__ == "__main__":

    unittest.main()