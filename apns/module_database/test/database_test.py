import unittest
import apns.module_database.database as db

class TestDatabase(unittest.TestCase):

    def test_element_index_tolabel(self):
        self.assertEqual(db.element_index_tolabel(1), "1_H")
        self.assertEqual(db.element_index_tolabel(2), "2_He")
        self.assertEqual(db.element_index_tolabel(3), "3_Li")
        self.assertEqual(db.element_index_tolabel(4), "4_Be")
        self.assertEqual(db.element_index_tolabel(5), "5_B")

    def test_element_label_toindex(self):
        self.assertEqual(db.element_label_toindex("H"), 1)
        self.assertEqual(db.element_label_toindex("He"), 2)
        self.assertEqual(db.element_label_toindex("Li"), 3)
        self.assertEqual(db.element_label_toindex("Be"), 4)
        self.assertEqual(db.element_label_toindex("B"), 5)

    def test_symbol_tol(self):
        self.assertEqual(db.symbol_tol("s"), 0)
        self.assertEqual(db.symbol_tol("p"), 1)
        self.assertEqual(db.symbol_tol("d"), 2)
        self.assertEqual(db.symbol_tol("f"), 3)
        self.assertEqual(db.symbol_tol("g"), 4)

    def test_l_tosymbol(self):
        self.assertEqual(db.l_tosymbol(0), "s")
        self.assertEqual(db.l_tosymbol(1), "p")
        self.assertEqual(db.l_tosymbol(2), "d")
        self.assertEqual(db.l_tosymbol(3), "f")
        self.assertEqual(db.l_tosymbol(4), "g")

    def test_element_mass(self):
        self.assertEqual(db.element_label_tomass("H"), 1.00794)
        self.assertEqual(db.element_label_tomass("He"), 4.002602)
        self.assertEqual(db.element_label_tomass("Li"), 6.941)
        self.assertEqual(db.element_label_tomass("Be"), 9.012182)
        self.assertEqual(db.element_label_tomass("B"), 10.811)

    def test_number_tomultiplicity(self):
        self.assertEqual(db.number_tomultiplicity(1), "s")
        self.assertEqual(db.number_tomultiplicity(2), "d")
        self.assertEqual(db.number_tomultiplicity(3), "t")
        self.assertEqual(db.number_tomultiplicity(4), "q")
    
    def test_multiplicity_tonumber(self):
        self.assertEqual(db.multiplicity_tonumber("s"), 1)
        self.assertEqual(db.multiplicity_tonumber("d"), 2)
        self.assertEqual(db.multiplicity_tonumber("t"), 3)
        self.assertEqual(db.multiplicity_tonumber("q"), 4)

if __name__ == "__main__":

    unittest.main()