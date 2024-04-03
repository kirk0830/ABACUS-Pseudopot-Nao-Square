<p align="center">
    <img src="docs/assets/images/apns.svg">
</p>  

# ABACUS Pseudopot-Nao Square  
## About  
**ABACUS** (Atomic-orbital Based Ab-initio Computation at UStc) **Pseudopot** (pseudopotential)-**Nao** (numerical atomic orbital) **Square** (**APNS**) is a project for continuously generating pseudopotential and numerical atomic orbital test data for ABACUS user. The project is based on the [ABACUS](https://github.com/deepmodeling/abacus-develop) project, which is a high-performance ab initio simulation software for electronic structure calculation.  
## Features
- **APNS** is designed as a test data generator for ABACUS user. It can generate pseudopotential and numerical atomic orbital test data for ABACUS user.
## To run on your local machine
**WARNING:** In principle APNS is not designed to be run on your local machine. However, if you want to run it on your local machine, you can follow the instructions below.
### Prerequisites
Prerequisites differ from workflows, currently there are three workflows. Mutual prerequisites are listed below:
- Python 3.6 or higher
- numpy
- scipy

Following lists prerequisites for each workflow:
#### Test
- urllib3 (< 2.1.0)
- mp_api
- pymatgen
- seekpath
#### Analysis
- matplotlib
- lbg (lebegue, developed by Deeptechnology Inc., for downloading groups of jobs from Bohrium Supercomputing Cloud platform)
### Orbgen
- pytorch
- torch_optimizer
- torch_complex
### Installation
*BEFORE INSTALLATION, WE STRONGLY RECOMMEND YOU TO CREATE A NEW VIRTUAL ENVIRONMENT.*
```bash
python3 -m virtualenv apnsvenv
source apnsvenv/bin/activate
```
in `.gitignore`, a line is already set for ignoring virtual environment folder entitled with substr `venv`. Then you can install the package.
To install, run:
```bash
pip install .
```
, add `-e` to install in editable mode.
```bash
pip install -e .
```
### Usage
#### Very first configuration
For APNS >= 1.0.0, no additional configuration is needed anymore, just run APNS, it will refresh the archive for pseudopotential and numerical atomic orbital (not implemented yet) automatically. However, you need to pay attention to the file `./download/pseudopotentials/rules.json`, it is for the configuration of pseudopotential archive. Once you add new kind of pseudopotential, you should always add a new rule to distinguish those pseudopotentials, like: 
```json
{
    "rules": [
        {
            "kind": "sg15",
            "version": "1.2",
            "appendix": "fr",
            "re.folder": ".*",
            "re.file": "^([A-Z][a-z]?)(_ONCV_PBE_FR-1\\.0\\.upf)$"
        }
    ]
}
```
, in which `re.folder` and `re.file` are regular expressions for matching the folder and file name of the pseudopotential.
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
        "cal_force": 1,
        "cal_stress": 1,
        "scf_nmax": 500
    },
    "extensive": {
        "characteristic_lengths": [0.0]
    },
    "systems": ["Cs", "Ba", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"],
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
APNS also supports specifying `pseudopotentials` section for each element:
```json
    "pseudopotentials": {
        "kinds": {"Cs": ["sg15"], "Ba": ["sg15", "dojo"]},
        "versions": {"Cs": ["1.2"], "Ba": ["all"]},
        "appendices": {"Cs": ["fr"], "Ba": [""]},
    },
```
Then run:
```bash
python main.py -i input.json
```
After jobs are done, you can download with `lbg` developed by DPTechnology Inc.:
```bash
lbg jobgroup download <group_id>
```
However if it is the first time you run `lbg`, you need to configure it first:
```bash
lbg config account
```
, your Bohrium account and password are needed. Then edit another file named `input_analysis.json`:
```json
{
    "global": {
        "test_mode": "analysis",
        "pseudo_dir": "./download/pseudopotentials",
        "orbital_dir": "./download/numerical_orbitals"
    },
    "analysis": {
        "search_domain": "./11845898",
        "items": [
            "driver_EcutwfcConv_20240319"
        ]
    }
}
```
, in which `search_domain` should be the folder (either relative or absolute path) in which APNS will search jobs to analysis, and `items` are the items you want to analyze. APNS encourages users to write their own analysis workflow, once complete, add to the `list` items, then driver will run them one-by-one. Then run:
```bash
python main.py -i input_analysis.json
```
The element in `items` should be coded with respect to some workflow regulations, they are:
1. always have an enter like
```python
import argparse
def entry():
    parser = argparse.ArgumentParser(description="APNS pseudopotential convergence test")
    # add -i
    parser.add_argument("-i", "--input", type=str, help="input json file")
    args = parser.parse_args()
    return args.input

def run():
    path = entry()
    # your code here

if __name__ == "__main__":
    run()
```
Then driver of analysis workflow will call the analyzer with
```python
def run(finp: str):
    with open(finp, "r") as f:
        inp = json.load(f)

    for item in inp["analysis"]["items"]:
        
        item = item.replace("\\", "/").split("/")[-1]
        item = item + ".py" if not item.endswith(".py") else item
        item = "apns/module_analysis/drivers/" + item
        if os.path.exists(item):
            os.system("python " + item + " -i" + inp["analysis"]["search_domain"])
        else:
            print("Warning: user-defined analysis item \"", item, "\" not found, skip.")
```
## Author information  
**APNS** is mainly developed and maintained by the ABACUS-AISI (Artificial Intelligence for Science Institute, BEIJING) team.  
