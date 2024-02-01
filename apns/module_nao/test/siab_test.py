import unittest
import apns.module_nao.siab as siab

class SiabTest(unittest.TestCase):

    def test_merge_dicts(self):
        """test normal cases"""
        source = [
            {
                "key1": "value1",
                "key2": "value2"
            },
            {
                "key1": "value1",
                "key2": "value2"
            }
        ]
        result = siab.merge_dicts(source)
        self.assertEqual(result, {
            "key1": ["value1", "value1"],
            "key2": ["value2", "value2"]
        })
        """then test with value_tostr = True"""
        source = [
            {
                "key1": 1,
                "key2": 2
            },
            {
                "key1": 1,
                "key2": 2
            }
        ]
        result = siab.merge_dicts(source, value_tostr=True)
        self.assertEqual(result, {
            "key1": ["1", "1"],
            "key2": ["2", "2"]
        })
        source = [
            {
                "key1": [1, 2],
                "key2": [3, 4]
            },
            {
                "key1": [1, 2],
                "key2": [3, 4]
            }
        ]
        result = siab.merge_dicts(source, value_tostr=True, delimiter=",")
        self.assertEqual(result, {
            "key1": ["1,2", "1,2"],
            "key2": ["3,4", "3,4"]
        })
        """test exception cases"""
        source = [
            {
                "key1": "value1",
                "key2": "value2"
            },
            {
                "key1": "value1",
                "key3": "value3"
            }
        ]
        with self.assertRaises(ValueError):
            result = siab.merge_dicts(source)
        
    def test_dict_tostr(self):
        """test normal cases"""
        source = {
            "key1": "value1",
            "key2": "value2"
        }
        result = siab.dict_tostr(source=source,
                                 key_sequence=["key1", "key2"],
                                 direction="vertical",
                                 comment_sign="#",
                                 with_key=True,
                                 len_placeholder=10)
        reference = 'key1      value1    \nkey2      value2    \n'
        self.assertEqual(result, reference)
        result = siab.dict_tostr(source=source,
                                 key_sequence=["key1", "key2"],
                                 direction="horizontal",
                                 comment_sign="#",
                                 with_key=True,
                                 len_placeholder=10)
        reference = '#       key1      key2\n      value1    value2\n'
        self.assertEqual(result, reference)
        """test exception cases"""
        # exception1: source is empty
        source = {}
        with self.assertRaises(ValueError):
            result = siab.dict_tostr(source=source,
                                     key_sequence=["key1", "key2"],
                                     direction="vertical",
                                     comment_sign="#",
                                     with_key=True,
                                     len_placeholder=10)
        # exception2: direction must be "vertical" or "horizontal"
        source = {
            "key1": "value1",
            "key2": "value2"
        }
        with self.assertRaises(ValueError):
            result = siab.dict_tostr(source=source,
                                     key_sequence=["key1", "key2"],
                                     direction="vertical1",
                                     comment_sign="#",
                                     with_key=True,
                                     len_placeholder=10)
        # exception3: length of values in source is not consistent
        source = {
            "key1": ["value1", "value2"],
            "key2": ["value2"]
        }
        with self.assertRaises(ValueError):
            result = siab.dict_tostr(source=source,
                                     key_sequence=["key1", "key2"],
                                     direction="vertical",
                                     comment_sign="#",
                                     with_key=True,
                                     len_placeholder=10)
        # exception4: key_sequence contains key not in source
        source = {
            "key1": "value1",
            "key2": "value2"
        }
        with self.assertRaises(ValueError):
            result = siab.dict_tostr(source=source,
                                     key_sequence=["key1", "key3"],
                                     direction="vertical",
                                     comment_sign="#",
                                     with_key=True,
                                     len_placeholder=10)

    def test_dicts_tostr(self):
        """test normal cases"""
        source = [
            {
                "key1": "value1",
                "key2": "value2"
            },
            {
                "key1": "value1",
                "key2": "value2"
            }
        ]
        result = siab.dicts_tostr(source=source,
                                  key_sequence=["key1", "key2"],
                                  direction="vertical",
                                  comment_sign="#",
                                  with_key=True,
                                  len_placeholder=10)
        reference = 'key1      value1    value1    \nkey2      value2    value2    \n'
        self.assertEqual(result, reference)
        result = siab.dicts_tostr(source=source,
                                  key_sequence=["key1", "key2"],
                                  direction="horizontal",
                                  comment_sign="#",
                                  with_key=True,
                                  len_placeholder=10)
        reference = '#       key1      key2\n      value1    value2\n      value1    value2\n'
        self.assertEqual(result, reference)
        """test exception cases"""
        # exception1: source is empty
        source = []
        with self.assertRaises(ValueError):
            result = siab.dicts_tostr(source=source,
                                      key_sequence=["key1", "key2"],
                                      direction="vertical",
                                      comment_sign="#",
                                      with_key=True,
                                      len_placeholder=10)
        # exception2: direction must be "vertical" or "horizontal"
        source = [
            {
                "key1": "value1",
                "key2": "value2"
            },
            {
                "key1": "value1",
                "key2": "value2"
            }
        ]
        with self.assertRaises(ValueError):
            result = siab.dicts_tostr(source=source,
                                      key_sequence=["key1", "key2"],
                                      direction="vertical1",
                                      comment_sign="#",
                                      with_key=True,
                                      len_placeholder=10)

    def test_siab_program_section(self):
        result = siab.siab_program_section(mpi_command="mpirun -np 4",
                                           abacus_command="abacus")
        reference = """#EXE_env                                
EXE_mpi             mpirun -np 4        
EXE_pw              abacus              
"""
        self.assertEqual(result, reference)

    def test_siab_electronic_calculation(self):

        result = siab.siab_electronic_calculation(element="H",
                                                  ecutwfc=100,
                                                  rcut=5,
                                                  fpseudo="H_ONCV_PBE-1.0.upf",
                                                  pseudo_dir="/home/data/pseudo",
                                                  smearing_sigma=0.01)
        reference = """element             H                   
Ecut                100                 
Rcut                5                   
Pseudo_dir          /home/data/pseudo   
Pseudo_name         H_ONCV_PBE-1.0.upf  
smearing_sigma      0.01                
"""
        self.assertEqual(result, reference)

    def test_siab_reference_system(self):
        reference_systems = [
            {
                "shape": "dimer",
                "nbands": 8,
                "bond_lengths": [1.8, 2.0, 2.3, 2.8, 3.8],
            },
            {
                "shape": "trimer",
                "nbands": 10,
                "bond_lengths": [1.9, 2.1, 2.6],
            }
        ]
        orbital_configurations = [
            {
                "reference_structure": "dimer",
                "nbands_ref": 4,
                "from_to": [None, "SZ"]
            },
            {
                "reference_structure": "dimer",
                "nbands_ref": 4,
                "from_to": ["SZ", "DZP"]
            },
            {
                "reference_structure": "trimer",
                "nbands_ref": 6,
                "from_to": ["DZP", "TZDP"]
            }
        ]
        result = siab.siab_reference_system(reference_systems=reference_systems,
                                            nspin=1,
                                            lmax=2,
                                            orbital_configurations=orbital_configurations)
        reference = """# identifier     shape          nbands         lmax           nspin          bond_lengths   
  STRU1          dimer          8              3              1              1.8 2.0 2.3 2.8 3.8
  STRU2          trimer         10             3              1              1.9 2.1 2.6    
""" 
        self.assertEqual(result, reference)

    def test_siab_siabparams(self):
        orbital_configurations = [
            {
                "reference_structure": "dimer",
                "nbands_ref": 4,
                "from_to": [None, "SZ"]
            },
            {
                "reference_structure": "dimer",
                "nbands_ref": 4,
                "from_to": ["SZ", "DZP"]
            },
            {
                "reference_structure": "trimer",
                "nbands_ref": 6,
                "from_to": ["DZP", "TZDP"]
            }
        ]
        reference_systems = [
            {
                "shape": "dimer",
                "nbands": 8,
                "bond_lengths": [1.8, 2.0, 2.3, 2.8, 3.8],
            },
            {
                "shape": "trimer",
                "nbands": 10,
                "bond_lengths": [1.9, 2.1, 2.6],
            }
        ]
        result = siab.siab_siabparams(minimal_basis=[2, 1, 1],
                                      reference_systems=reference_systems,
                                      orbital_configurations=orbital_configurations,
                                      maxstep=9000)
        reference = """max_steps 9000
# orb_id         stru_id        nbands_ref     orb_ref        orb_config     
  Level1         STRU1          4              none           2s1p1d         
  Level2         STRU1          4              fix            4s2p2d1f       
  Level3         STRU2          6              fix            6s3p3d2f       
"""
        self.assertEqual(result, reference)
    
    def test_siab_save(self):

        orbital_configurations = [
            {
                "reference_structure": "dimer",
                "nbands_ref": 4,
                "from_to": [None, "SZ"]
            },
            {
                "reference_structure": "dimer",
                "nbands_ref": 4,
                "from_to": ["SZ", "DZP"]
            },
            {
                "reference_structure": "trimer",
                "nbands_ref": 6,
                "from_to": ["DZP", "TZDP"]
            }
        ]
        result = siab.siab_save(orbital_configurations=orbital_configurations)
        reference = """# save_id        orb_id         zeta_notation  
  Save1          Level1         SZ             
  Save2          Level2         DZP            
  Save3          Level3         TZDP           
"""
        self.assertEqual(result, reference)

    def test_orbitalconfig_tolist(self):
        """single orbital"""
        self.assertEqual(siab.orbitalconfig_tolist("1s"), [1])
        self.assertEqual(siab.orbitalconfig_tolist("2p"), [0, 2])
        self.assertEqual(siab.orbitalconfig_tolist("3d"), [0, 0, 3])
        """multiple orbitals"""
        self.assertEqual(siab.orbitalconfig_tolist("1s2p"), [1, 2])
        self.assertEqual(siab.orbitalconfig_tolist("1s3d"), [1, 0, 3])
        self.assertEqual(siab.orbitalconfig_tolist("2p3d"), [0, 2, 3])
        self.assertEqual(siab.orbitalconfig_tolist("1s2p3d"), [1, 2, 3])
        """multiple orbitals with different number of orbitals"""
        self.assertEqual(siab.orbitalconfig_tolist("1s2p3d4f"), [1, 2, 3, 4])

    def test_zeta_notation_toorbitalconfig(self):
        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZ", [2, 1, 1]), "2s1p1d")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZP", [2, 1, 1]), "2s1p1d1f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZ9P", [2, 1, 1]), "2s1p1d9f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("DZ", [2, 1, 1]), "4s2p2d")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("DZP", [2, 1, 1]), "4s2p2d1f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("TZDP", [2, 1, 1]), "6s3p3d2f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("QZDP", [2, 1, 1]), "8s4p4d2f")

        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZ", [2, 0, 1, 1]), "2s1d1f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZP", [2, 0, 1, 1]), "2s1d1f1g")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("SZ9P", [2, 0, 1, 1]), "2s1d1f9g")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("DZ", [2, 0, 1, 1]), "4s2d2f")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("DZP", [2, 0, 1, 1]), "4s2d2f1g")
        self.assertEqual(siab.zeta_notation_toorbitalconfig("TZDP", [2, 0, 1, 1]), "6s3d3f2g")                         

    def test_SIAB_INPUT(self):
        reference_systems = [
            {
                "shape": "dimer",
                "nbands": 8,
                "bond_lengths": [1.8, 2.0, 2.3, 2.8, 3.8],
            },
            {
                "shape": "trimer",
                "nbands": 10,
                "bond_lengths": [1.9, 2.1, 2.6],
            }
        ]
        orbital_configurations = [
            {
                "reference_structure": "dimer",
                "nbands_ref": 4,
                "from_to": [None, "SZ"]
            },
            {
                "reference_structure": "dimer",
                "nbands_ref": 4,
                "from_to": ["SZ", "DZP"]
            },
            {
                "reference_structure": "trimer",
                "nbands_ref": 6,
                "from_to": ["DZP", "TZDP"]
            }
        ]
        result = siab.SIAB_INPUT(element="Fe",
                                 ecutwfc=100,
                                 nspin=1,
                                 rcut=[6, 7, 8, 9, 10],
                                 fpseudo="Fe_ONCV_PBE-1.2.upf",
                                 pseudo_dir="./download/pseudopotentials/sg15_oncv_upf_2020-02-06/1.2",
                                 minimal_basis=None,
                                 smearing_sigma=0.015,
                                 reference_systems=reference_systems,
                                 orbital_configurations=orbital_configurations,
                                 maxstep=9000
                                 )
        reference = """# Refactor version of SIAB_INPUT of PTG_dpsi for generating numerical atomic orbitals of ABACUS
# from ABACUS-Planewave DFT calculations - For high throughput auto-generation and test of
# pseudopotentials and numerical atomic orbitals. This is included in project ABACUS Pseudopot-
# Nao Square (APNS). Visit related Github repos for more information:
# ABACUS (deepmodeling) Github repo: https://github.com/deepmodeling/abacus-develop
# PTG_dpsi (abacusmodeling/ABACUS-orbitals) Github repo: https://github.com/abacusmodeling/ABACUS-orbitals
# APNS Github repo: https://github.com/kirk0830/ABACUS-Pseudopot-Nao-Square
# APNS Github Pages: https://kirk0830.github.io/ABACUS-Pseudopot-Nao-Square
# APNS is mainly developed and maintained by ABACUS-AISI developer team

# PROGRAM CONFIGURATION
#EXE_env                                
EXE_mpi             mpirun -np 1        
EXE_pw              abacus              

# ELECTRONIC STRUCTURE CALCULATION
element             Fe                  
Ecut                100                 
Rcut                6 7 8 9 10          
Pseudo_dir          ./download/pseudopotentials/sg15_oncv_upf_2020-02-06/1.2
Pseudo_name         Fe_ONCV_PBE-1.2.upf 
smearing_sigma      0.015               

# REFERENCE SYSTEMS
# identifier     shape          nbands         lmax           nspin          bond_lengths   
  STRU1          dimer          8              3              1              1.8 2.0 2.3 2.8 3.8
  STRU2          trimer         10             3              1              1.9 2.1 2.6    

# SIAB PARAMETERS
max_steps 9000
# orb_id         stru_id        nbands_ref     orb_ref        orb_config     
  Level1         STRU1          4              none           2s1p1d         
  Level2         STRU1          4              fix            4s2p2d1f       
  Level3         STRU2          6              fix            6s3p3d2f       

# SAVE
# save_id        orb_id         zeta_notation  
  Save1          Level1         SZ             
  Save2          Level2         DZP            
  Save3          Level3         TZDP           
"""
        self.assertEqual(result, reference)
        
    def test_siab_input_parse(self):
        result = siab.siab_input_parse("./apns/module_nao/test/support/SIAB_INPUT")
    
    def test_set_bond_length_fromfile(self):
        result = siab.set_bond_length_fromfile("./apns/module_nao/test/support/SIAB_INPUT", [
            {"shape": "dimer", "nbands": 8},
            {"shape": "trimer", "nbands": 10}
        ])
        reference = [1.8, 2.0, 2.3, 2.8, 3.8]
        self.assertEqual(result[0]["bond_lengths"], reference)

if __name__ == "__main__":
    unittest.main()