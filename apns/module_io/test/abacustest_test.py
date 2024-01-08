import unittest
import apns.module_io.abacustest as abacustest

class TestAbacustest(unittest.TestCase):

    def test_image_information(self):
        abacustest.image_information(software="ABACUS")
        abacustest.image_information(software="qespresso")

if __name__ == "__main__":
    unittest.main()