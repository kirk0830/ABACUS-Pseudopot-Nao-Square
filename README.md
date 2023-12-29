<p align="center">
    <img src="docs/apns.svg">
</p>  

# ABACUS Pseudopot-Nao Square  
## About  
**ABACUS** (Atomic-orbital Based Ab-initio Computation at UStc) **Pseudopot** (pseudopotential)-**Nao** (numerical atomic orbital) **Square** (**APNS**) is a project for continuously generating pseudopotential and numerical atomic orbital test data for ABACUS user. The project is based on the [ABACUS](https://github.com/deepmodeling/abacus-develop) project, which is a high-performance ab initio simulation software for electronic structure calculation.  
## Features
- **APNS** is designed as a test data generator for ABACUS user. It can generate pseudopotential and numerical atomic orbital test data for ABACUS user.
## To run on your local machine
**WARNING:** In principle APNS is not designed to be run on your local machine. However, if you want to run it on your local machine, you can follow the instructions below.
### Prerequisites
    - ABACUS
    - Python 3.6 or higher
    - numpy
    - matplotlib
    - mp_api (Materials Project API, Python-end)
    - pymatgen
### Installation
```bash
python setup.py install
```
### Usage
#### Very first configuration
1. download pseudopotentials in directory `./download/pseudopotentials`, decompress or unzip them if necessary. Sort folders like this:
```bash
download
├── pseudopotentials
│   ├── NCPP-PBE-PD04
│   │   ├── ...upf
│   │   └── ...upf
│   ├── NCPP-PBE-PD03
│   │   ├── ...upf
│   │   └── ...upf
│   ├── sg15_oncv_upf_2020-02-06
│   │   ├── ...upf
│   │   └── ...upf
```
2. run pseudopotential archiving task like:
```python
import apns.module_pseudo.upf_archive as arch

if __name__ == '__main__':
    arch.archive(pseudo_dir = "./download/pseudopotentials",
                 only_scan = False)
```
Then in `download/pseudopotentials`, all pseudopotential files will be moved into specific folders and in each folder there will be a `description.json` created, which contains the information of the pseudopotential.
#### Run APNS
prepare `input.json` like this:
```json
{
    "global": {
        "test_mode": "pseudopotential",
        "software": "ABACUS",
        "work_dir": "./",
        "pseudo_dir": "./download/pseudopotentials",
        "orbital_dir": "./download/numerical_orbitals",
        "save_log": true
    },
    "calculation": {
        "basis_type": "pw",
        "functionals": ["PBE"],
        "ecutwfc": [100],
        "cell_scaling": [0.0]
    },
    "additional_keywords": {
    },
    "systems": ["H", "He"],
    "materials_project": {
        "api_key": "your Materials Project API key",
        "n_structures": 1,
        "theoretical": false,
        "most_stable": true
    },
    "pseudopotentials": {
        "kinds": ["all"],
        "versions": [""],
        "appendices": [""]
    },
    "numerical_orbitals": {
        "types": ["DZP"],
        "rcuts": [7, 8, 9, 10],
        "appendices": [""]
    }
}
```
For detailed description of all keywords, see [keywords.md](https://github.com/kirk0830/ABACUS-Pseudopot-Nao-Square/blob/main/docs/keywords.md). After preparing `input.json`, remember to check `main.py` if there is anything you should modify, especially the `input.json` path.  
Then run:
```bash
python main.py
```
#### Output
The output will be an `*.zip` file, contains all the test data.
#### Note
1. The `input.json` does not support comment.
2. make sure you have your own Materials Project API. For details, see `README.md` in `apns/module_structures/`
## Author information  
**APNS** is mainly developed and maintained by the ABACUS-AISI (Artificial Intelligence for Science Institute, BEIJING) team.  