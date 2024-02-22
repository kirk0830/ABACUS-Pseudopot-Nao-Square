import unittest
import apns.module_nao.orbgen as amno

class TestOrbgen(unittest.TestCase):

    def test_find_ecutwfc(self):
        result = amno.find_ecutwfc(element="Fe",
                                   ecutwfc=30,
                                   pspot_id="sg15_10")
        self.assertTupleEqual(result, ([30], ['sg15_10']))
        result = amno.find_ecutwfc(element="Zn",
                                   ecutwfc=None,
                                   pspot_id="sg15_10")
        self.assertTupleEqual(result, ([150.0], ['sg15_10']))
        result = amno.find_ecutwfc(element="Zn",
                                   ecutwfc=None,
                                   pspot_id=None)
        self.assertTupleEqual(result, 
                              ([200.0, 200.0, 200.0, 200.0, 200.0, 
                                200.0, 200.0, 150.0, 150.0, 150.0], 
                                ['dojo03', 'dojo04', 'dojo04fr', 'dojo05', 'pd03', 
                                 'pd04', 'pd04sp', 'sg1510', 'sg1510fr', 'sg1512']))
        
    def test_find_fpseudo(self):
        result = amno.find_fpseudo(element="Fe",
                                   pspot_id="sg15_10")
        self.assertTupleEqual(result, 
                              ('D:/abacus-pseudopot-nao-square/download/pseudopotentials/sg15_oncv_upf_2020-02-06/1.0/', 
                               'Fe_ONCV_PBE-1.0.upf'))
        result = amno.find_fpseudo(element="Zn",
                                   pspot_id="sg15_10")
        self.assertTupleEqual(result, 
                              ('D:/abacus-pseudopot-nao-square/download/pseudopotentials/sg15_oncv_upf_2020-02-06/1.0/', 
                               'Zn_ONCV_PBE-1.0.upf'))

    def test_siab_generator(self):
        counter_siab = 0
        ref_ecutwfc = [200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 150.0, 150.0, 150.0]
        for isiab, siab_input in enumerate(amno.siab_generator(element="Zn",
                                                               rcuts=[7, 8, 9, 10],
                                                               ecutwfc=None,
                                                               pspot_id=None)):
            counter_siab += 1
            self.assertEqual(siab_input["ecutwfc"], ref_ecutwfc[isiab])

if __name__ == "__main__":
    unittest.main()