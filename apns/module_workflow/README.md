<p align="center">
    <img src="../../docs/apns.svg">
</p>

# Workflow design of ABACUS Pseudopot-Nao square
## Data carrier
There are mainly two data carrier (Python dictionary) in this package: `work_status` and `test_status`.
### `work_status`
`work_status` is used to define the workflow, or say a full-keyword input script to describe workflow to perform. `work_status` transformed from input script is straightforward, more specifically it includes following steps:
1. paths information convert to absolute path
2. expand list into element-by-element keys and values defined dictionaries of pseudopotentials and numerical atomic orbitals, like:
    ```json5
    // transform this
    {
        "systems": ["TiO2", "Pt"],
        "pseudopotentials": {
            "kinds": ["sg15", "pd", "dojo"],
            "versions": ["10", "04", "05"],
            "appendices": ["", "fr"]
        }
    }
    // into this (systems decomposed into elements is done by other function)
    {
        "systems": ["TiO2_", "Pt"],
        "pseudopotentials": {
            "Ti": {
                    "kinds": ["sg15", "pd", "dojo"],
                    "versions": ["10", "04", "05"],
                    "appendices": ["", "fr"]
                },
            "O": {
                    "kinds": ["sg15", "pd", "dojo"],
                    "versions": ["10", "04", "05"],
                    "appendices": ["", "fr"]
                },
            "Pt": {
                    "kinds": ["sg15", "pd", "dojo"],
                    "versions": ["10", "04", "05"],
                    "appendices": ["", "fr"]
                }
        }
    }
    ```
### `test_status`
`test_status` contains information of expanded parameter-determined tests, it is more about specific test information,
and organizes information test-by-test. Commonly it is organized as follows:
```json5
{
    "tests": {
        "ErO": {
            "pd04sg1510_T7D7": { // test identifier
                "elements": ["Er", "O"],
                "pseudopotentials": {
                    "files": {
                        "Er": "Er-sp.PD04.PBE.UPF",
                        "O": "O_ONCV_PBE-1.0.upf"
                    },
                    "info": {
                        "Er": {"kind": "pd", "version": "04", "appendix": ""},
                        "O": {"kind": "sg15", "version": "10", "appendix": ""}
                    }
                },
                "numerical_orbitals": {
                    "files": {
                        "Er": "Er_gga_7au_200Ry_6s3p3d3f2g.orb",
                        "O": "O_gga_7au_100Ry_2s2p1d.orb"
                    },
                    "info": {
                        "Er": {"type": "TZDP", "rcut": 7, "appendix": ""},
                        "O": {"type": "DZP", "rcut": 7, "appendix": ""}
                    }
                }
            },
            "pd04sg1510_T7D8": {
                "elements": ["Er", "O"],
            /*...*/
            }
        },
    }
}
```
All tests input scripts will generate according to `test_status`.
## Identifiers
Because we have multi-dimensional tests, the test space is spanned, or a direct Cartesian product of all parameters. It is needed to identify each test uniquely by its name, therefore identifier is designed for this task.
### Pseudopotential: `kind`\_`version`\_`appendix`
#### Domain1: `type`
Commonly pseudopotential contains information about its *type*, like *SG15*, *PD*, *ATOMPAW*, *GBRV*, and pseudopotentials collected in Quantum ESPRESSO's pslibrary, and with RRKJ scheme, the ultrasoft pseudopotential is named as *rrkjus*, and PAW one as *kjpaw*.  
#### Domain2: `version`
The second domain of identifier of pseudopotential is its version. For *SG15*, it has *1.0*, *1.1* and *1.2* versions, for *PD*, it has at least *03* and *04* versions. For norm-conserving pseudopotential collected in pseudojo.org, there are *0.3*, *0.4*, *0.5* and even one "plus" versions.
#### Domain3: `appendix`
However, even for one certain type and its specific version, there may also be different type of pseudopotentials provided for one element, like *SG15* has full-relativistic and its scalar-relativistic version, *PD* has *core* and *core-icmod1* additionally versions for rare-earth elements. Therefore, the third domain of identifier is designed to distinguish these different types of pseudopotentials.
### Numerical atomic orbital: (`pseudopotential identifier`_)`type`\_`rcut`\_`appendix`
For numerical atomic orbital, it is, tightly bind with pseudopotential. Therefore the identifier of pseudopotential it corresponds, is also part of identifier of itself. 
#### Domain1: `type`
Additionally, or say besides the specific pseudopotential, the *orbital configuration*, *cutoff radius* are two commonly referred parameters to identify a numerical atomic orbital. The first, describing zeta or say orbital configuration, is named as `type` of orbital. *DZP*, with means *double-zeta plus polarization*, is a common type of orbital, and *TZDP* is also a common type of orbital, which means *triple-zeta double polarization*, they two are denoted as *D* and *T* respectively.
#### Domain2: `rcut`
Closely behind `type`, `rcut` is another important parameter to identify a numerical atomic orbital. It is the cutoff radius of numerical atomic orbital, and is denoted as `rcut` in identifier. For example, *DZP* with *rcut* of 7 Bohr is denoted as *D7*, and *TZDP* with *rcut* of 7 Bohr is denoted as *T7*.
#### Domain3: `appendix`
It is also possible to compare orbitals before and after some special treatment like *high frequency oscillation remove*, then appendix will be, named like `osci_rm`.
### Test
With identifiers of pseudopotential and numerical atomic orbital, one test is certainly capable to uniquely be identified. Simply arrange them in the order of element declaration of system. For ErO, as shown above, it is (pseudopotential identifer short as PI, numerical atomic orbital identifier short as NAOI):
```
pd04sg1510_T7D7 -> identifier of a test
^     ^    ^ ^
PI_Er |    | |
      PI_O | |
           NAOI_Er
             NAOI_O
```
### Folder
Furtherly, the test is defined within the domain of a system, therefore system itself is also a part of identifier to identify test among all tests submit finally. The system identifier is append at the end of test identifier, and, DFT xc functional is an identifier added at the beginning of test identifier. For example, for ErO, it is:
```
t_pbe_pd04sg1510_T7D7_ErO -> identifier of a folder
  ^   ^     ^    ^ ^  ^
  xc  |     |    | |  |
      PI_Er |    | |  |
            PI_O | |  |
                 NAOI_Er
                   NAOI_O
                      system
```
## Workflow
There are multiple test workflows designed, see Feishu Documents:  
[ABACUS Pseudopot-Nao square](https://ucoyxk075n.feishu.cn/docx/DiaqdLxmyoltgLxU4iEcdzoVnwb)