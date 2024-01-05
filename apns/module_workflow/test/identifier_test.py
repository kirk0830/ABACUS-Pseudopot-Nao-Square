import unittest
from apns.module_workflow import identifier as id

class IdentifierTest(unittest.TestCase):
    """Test the identifier module
    """
    def test_pseudopotential(self):

        self.assertEqual(id.pseudopotential(kind="GGA", version="PBE", appendix=""), "GGA_PBE")
        self.assertEqual(id.pseudopotential(kind="GGA", version="", appendix=""), "GGA")
        self.assertEqual(id.pseudopotential(kind="GGA", version="PBE", appendix="nc"), "GGA_PBE_nc")
        self.assertEqual(id.pseudopotential(kind="GGA", version="", appendix="nc"), "GGA_nc")

    def test_numerical_orbital(self):

        self.assertEqual(id.numerical_orbital(type="DZP", rcut=0, appendix=""), "DZP_0")
        self.assertEqual(id.numerical_orbital(type="DZP", rcut=0, appendix="nc"), "DZP_0_nc")
        self.assertEqual(id.numerical_orbital(type="DZP", rcut=3, appendix="nc"), "DZP_3_nc")
        self.assertEqual(id.numerical_orbital(type="DZP", rcut=3, appendix=""), "DZP_3")
    
    def test__pseudopotential(self):

        self.assertEqual(id._pseudopotential(identifier="GGA_PBE"), ("GGA", "PBE", ""))
        self.assertEqual(id._pseudopotential(identifier="GGA"), ("GGA", "", ""))
        self.assertEqual(id._pseudopotential(identifier="GGA_PBE_nc"), ("GGA", "PBE", "nc"))
        
    def test__numerical_orbital(self):
            
        self.assertEqual(id._numerical_orbital(identifier="DZP"), ("DZP", "", ""))
        self.assertEqual(id._numerical_orbital(identifier="DZP_3_nc"), ("DZP", "3", "nc"))
        self.assertEqual(id._numerical_orbital(identifier="DZP_3"), ("DZP", "3", ""))

    def test_folder(self):

        self.assertEqual(id.folder(functional="PBE", system="mp-1", specific_test=""), "t_pbe_mp-1_")
        self.assertEqual(id.folder(functional="PBE", system="mp-1", specific_test="test"), "t_pbe_mp-1_test")
        self.assertEqual(id.folder(functional="PBE", system="mp-1", specific_test="test_1"), "t_pbe_mp-1_test_1")

if __name__ == "__main__":
    unittest.main()