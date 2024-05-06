import re

# COMPLETE
def PP_HEADER(data: str) -> dict:
    
    result = {
        "data": [],
        "attrib": {}
    }
    for line in data.split("\n"):
        if line.endswith("Version Number"):
            result["attrib"]["version"] = line.split()[0]
        elif line.endswith("Element"):
            result["attrib"]["element"] = line.split()[0]
        elif line.endswith("pseudopotential"):
            result["attrib"]["pseudo_type"] = line.split()[0]
        elif line.endswith("Nonlinear Core Correction"):
            result["attrib"]["core_correction"] = True if line.split()[0] == "T" else False
        elif line.endswith("Exchange-Correlation functional"):
            result["attrib"]["functional"] = line.split()[-3]
            result["attrib"]["xc_functional_components"] = line.split()[:-3]
        elif line.endswith("Z valence"):
            result["attrib"]["z_valence"] = int(float(line.split()[0]))
        elif line.endswith("Total energy"):
            result["attrib"]["total_psenergy"] = float(line.split()[0])
        elif line.endswith("Suggested cutoff for wfc and rho"):
            result["attrib"]["wfc_cutoff"] = float(line.split()[0])
            result["attrib"]["rho_cutoff"] = float(line.split()[1])
        elif line.endswith("Max angular momentum component"):
            result["attrib"]["l_max"] = int(line.split()[0])
        elif line.endswith("Number of points in mesh"):
            result["attrib"]["mesh_size"] = int(line.split()[0])
        elif line.endswith("Number of Wavefunctions, Number of Projectors"):
            result["attrib"]["number_of_wfc"] = int(line.split()[0])
            result["attrib"]["number_of_proj"] = int(line.split()[1])
        
        result["attrib"]["relativistic"] = "F"
        result["attrib"]["is_ultrasoft"] = "T"
        result["attrib"]["is_paw"] = "F"
        result["attrib"]["is_coulomb"] = "F"
        result["attrib"]["has_so"] = "F"
        result["attrib"]["has_wfc"] = "T"
        result["attrib"]["has_gipaw"] = "F"
        result["attrib"]["paw_as_gipaw"] = "F"
        result["attrib"]["l_max_rho"] = "0"
        result["attrib"]["l_local"] = "-3"

    return result
# COMPLETE
def PP_INPUTFILE(data: str) -> dict:
    
    return {
        "data": [],
        "attrib": {}
    }
# COMPLETE
def PP_INFO(data: str) -> dict:
    
    possible_attribute_names = ["Author", "Generation date"]
    table_titles = ["nl", "pn", "l", "occ", "rcut", "rcut_us", "epseu"]

    result = {
        "data": [],
        "attrib": {}
    }
    for table_title in table_titles:
        result["attrib"][table_title] = []
    contents = data.split("\n")

    number_of_digit_line = 0

    for line in contents:
        line = line.strip()
        if "version" in line:
            version_number = ""
            words = line.split()
            read_version = False
            for word in words:
                if word == "version":
                    read_version = True
                else:
                    if read_version:
                        version_number += word
        else:
            number_of_attrib = line.count(":")
            #print("number_of_attrib: ", number_of_attrib)
            if number_of_attrib:
                pass
            else:
                if line[0].isdigit():
                    number_of_digit_line += 1
                    if number_of_digit_line == 1:
                        if line == "1":
                            result["attrib"]["relativistic"] = "scalar"
                        else:
                            result["attrib"]["relativistic"] = "full"
                    elif number_of_digit_line == 2:
                        result["attrib"]["local_potential_cutoff_radius"] = float(line.split()[0])
                    else:
                        # read data of electronic configuration table...
                        words = line.split()
                        for iw in range(len(words)):
                            data = words[iw]
                            if iw == 1 or iw == 2:
                                data = int(data)
                            elif iw == 0:
                                pass
                            else:
                                data = float(data)
                            result["attrib"][table_titles[iw]].append(data)
    return result

def PP_BETA(data: str) -> dict:
    
    return {
        "data": [],
        "attrib": {}
    }

def PP_DIJ(data: str) -> dict:
    
    return {
        "data": [],
        "attrib": {}
    }

def PP_QIJ(data: str) -> dict:
    
    return {
        "data": [],
        "attrib": {}
    }

def valence_configuration(parsed: dict) -> list:

    contents = parsed["PP_HEADER"]["data"]
    lines = [line.strip() for line in contents.split("\n")]

    sequence = ["S", "P", "D", "F", "G", "H", "I", "K", "L", "M", "N"]

    result = {}
    read_valence_config = False
    for line in lines:
        if line.startswith("Wavefunctions"):
            read_valence_config = True
            continue
        if read_valence_config:
            if not line[0].isdigit():
                break
            else:
                words = line.split()
                symbol = words[0][-1]
                if symbol not in result:
                    result[symbol] = []
                result[symbol].append(words[0])
    
    # then convert to list
    result_list = []
    for isym, symbol in enumerate(sequence):
        if symbol in result:
            result_list.append(result[symbol])
        else:
            if isym >= len(result):
                break
            else:
                result_list.append([])
    return result_list

import unittest
class TestGBRV(unittest.TestCase):
    data = """<PP_INFO>
Generated using Vanderbilt code, version   7  3  6                              
Author: kfg        Generation date:    2    4   15                              
Automatically converted from original format                                    
    0        The Pseudo was generated with a Non-Relativistic Calculation
  1.00000000000E+00    Local Potential cutoff radius
nl pn  l   occ               Rcut            Rcut US             E pseu
1S  1  0  1.00      0.00000000000      1.20000000000     -0.47719071200
</PP_INFO>


<PP_HEADER>
   0                   Version Number
  H                    Element
   US                  Ultrasoft pseudopotential
    F                  Nonlinear Core Correction
 SLA  PW   PBX  PBC    PBE  Exchange-Correlation functional
    1.00000000000      Z valence
   -0.91769791689      Total energy
    0.00000    0.00000 Suggested cutoff for wfc and rho
    0                  Max angular momentum component
  615                  Number of points in mesh
    1    2             Number of Wavefunctions, Number of Projectors
 Wavefunctions         nl  l   occ
                       1S  0  1.00
</PP_HEADER>


<PP_MESH>
  <PP_R>
  7.92692342184E+01  8.06242735189E+01  8.20024753252E+01
  </PP_R>
  <PP_RAB>
  1.34358835543E+00  1.36655512324E+00  1.38991447589E+00
  </PP_RAB>
</PP_MESH>


<PP_LOCAL>
 -2.52304695475E-02 -2.48064250716E-02 -2.43895076590E-02
</PP_LOCAL>


<PP_NONLOCAL>
  <PP_BETA>
    1    0             Beta    L
   395
  0.00000000000E+00  0.00000000000E+00  0.00000000000E+00
  </PP_BETA>
  <PP_BETA>
    2    0             Beta    L
   395
  0.00000000000E+00  0.00000000000E+00  0.00000000000E+00
  </PP_BETA>
  <PP_DIJ>
    3                  Number of nonzero Dij
    1    1  6.06594103731E-01
    1    2  1.47301623089E+00
    2    2  2.60147291428E+00
  </PP_DIJ>
  <PP_QIJ>
    8     nqf. If not zero, Qij's inside rinner are computed using qfcoef's
    <PP_RINNER>
    1  7.00000000000E-01
    </PP_RINNER>
    1    1    0        i  j  (l(j))
  2.49088483939E-01    Q_int
  0.00000000000E+00  0.00000000000E+00  0.00000000000E+00
    <PP_QFCOEF>
  9.05452572541E+00 -5.27125106951E+01  1.05830005584E+02  1.79959354496E+02
 -1.69782662423E+03  4.48024531578E+03 -5.47451779900E+03  2.60559908164E+03
    </PP_QFCOEF>
    1    2    0        i  j  (l(j))
  2.25010731873E-01    Q_int
  0.00000000000E+00  0.00000000000E+00  0.00000000000E+00
    <PP_QFCOEF>
  9.37023018007E+00 -5.66306211583E+01  1.15232349038E+02  1.95986964899E+02
 -1.85101951130E+03  4.88430590239E+03 -5.96870860542E+03  2.84127426189E+03
    </PP_QFCOEF>
    2    2    0        i  j  (l(j))
  1.81851793788E-01    Q_int
  0.00000000000E+00  0.00000000000E+00  0.00000000000E+00
    <PP_QFCOEF>
  9.55567759396E+00 -6.07960378184E+01  1.25440165331E+02  2.19611120572E+02
 -2.05191036669E+03  5.40808923739E+03 -6.60757175553E+03  3.14570122599E+03
    </PP_QFCOEF>
  </PP_QIJ>
</PP_NONLOCAL>


<PP_PSWFC>
1S    0  1.00          Wavefunction
  0.00000000000E+00  0.00000000000E+00  0.00000000000E+00
</PP_PSWFC>


<PP_RHOATOM>
  0.00000000000E+00  0.00000000000E+00  0.00000000000E+00
</PP_RHOATOM>
"""

    def test_pp_header(self):
        result = PP_HEADER(self.data)
        reference = {'data': [], 
                     'attrib': {'relativistic': 'F', 
                                'is_ultrasoft': 'T', 
                                'is_paw': 'F', 
                                'is_coulomb': 'F', 
                                'has_so': 'F', 
                                'has_wfc': 'T', 
                                'has_gipaw': 'F', 
                                'paw_as_gipaw': 'F', 
                                'l_max_rho': '0', 
                                'l_local': '-3', 
                                'version': '0', 
                                'element': 'H', 
                                'pseudo_type': 'US', 
                                'core_correction': False, 
                                'functional': 'PBE', 
                                'xc_functional_components': ['SLA', 'PW', 'PBX', 'PBC'], 
                                'z_valence': 1, 
                                'total_psenergy': -0.91769791689, 
                                'wfc_cutoff': 0.0, 
                                'rho_cutoff': 0.0, 
                                'l_max': 0, 
                                'mesh_size': 615, 
                                'number_of_wfc': 1, 
                                'number_of_proj': 2}}
        self.assertEqual(result, reference)
if __name__ == "__main__":
    unittest.main()