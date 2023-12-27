import GBRV_Vanderbilt as GBRV

test_str="""0                   Version Number
As                   Element
US                  Ultrasoft pseudopotential
T                  Nonlinear Core Correction
SLA  PW   PBE  PBE     PBE  Exchange-Correlation functional
5.00000000000      Z valence
-39.93046173160      Total energy
0.00000    0.00000 Suggested cutoff for wfc and rho
2                  Max angular momentum component
875                  Number of points in mesh
2    6             Number of Wavefunctions, Number of Projectors
Wavefunctions         nl  l   occ
                       4S  0  2.00
                       4P  1  3.00"""

print(GBRV._PP_HEADER_(test_str))

test_str="""Generated using Vanderbilt code, version   7  3  6                              
Author: kfg        Generation date:   13    5 2013                              
Automatically converted from original format                                    
    1        The Pseudo was generated with a Scalar-Relativistic Calculation
  1.50000000000E+00    Local Potential cutoff radius
nl pn  l   occ               Rcut            Rcut US             E pseu
4S  4  0  2.00      0.00000000000      1.50000000000     -1.06527890600
4P  4  1  3.00      0.00000000000      1.50000000000     -0.38182955700"""

print(GBRV._PP_INFO_(test_str))