<p align="center">
    <img src="docs/assets/images/apns.svg">
</p>  

# ABACUS Pseudopot-Nao Square

## About

**ABACUS** (Atomic-orbital Based Ab-initio Computation at UStc) **Pseudopot**ential-**N**umerical **a**tomic **o**rbital **Square** (**APNS**) is a project for continuously generating pseudopotential and numerical atomic orbital test data for ABACUS user. The project is based on the [ABACUS](https://github.com/deepmodeling/abacus-develop) project, which is a high-performance ab initio simulation software for electronic structure calculation.  

### How/why this projecty is initiated?

This project is initially designed to answer the practical question that "**in practice**, how accurate can ABACUS calculate the electronic structure of a material with **currently available** given pseudopotential and numerical atomic orbital?"

That is to say, we have tested commonly-seen all pseudopotentials that supported by ABACUS.

On the other hand, we also want to provide an error-standardized reference for users who use ABACUS to train their machine-learning forcefields.

### What is obtained now?

We have tested the efficiency and precision of norm-conserving and ultra-soft pseudopotentials for ABACUS, and the results are published online in [AIS-Square](https://aissquare.com/) which is maintained by the [Artificial Intelligence for Science Institute of BEIJING](https://aisi.ac.cn/).

**Parts** of norm-conserving pseudopotential results has been plotted and shown on website: [APNS: Pseudopotentials](https://kirk0830.github.io/ABACUS-Pseudopot-Nao-Square/pseudopotential/pseudopotential.html).

### What is this repository for?

What is under development in this repository is the high-throughput and highly automated workflow for generating pseudopotential and numerical atomic orbital test suite for ABACUS development. For most of ABACUS users, this repository is not necessary.

### Who supports this project?

This project is fully supported by the [Artifical Intelligence for Science Institute of BEIJING](https://aisi.ac.cn/), as a part of [DeepModeling](https://github.com/deepmodeling) open-source development community.

### Future development plan

We have noticed the excellent works of:

- G. Prandini et al. published in 2018 ([Prandini, Gianluca, et al. "Precision and efficiency in solid-state pseudopotential calculations." npj Computational Materials 4.1 (2018): 72.](https://www.nature.com/articles/s41524-018-0127-2)), along with the [Standard Solid-State Pseudopotential library](https://www.materialscloud.org/discover/sssp/table/efficiency) online database.

- E. Bosoni et al. published in 2024 ([Bosoni, Emanuele, et al. "How to verify the precision of density-functional-theory implementations via reproducible and universal workflows." Nature Reviews Physics 6.1 (2024): 45-58.](https://www.nature.com/articles/s42254-023-00655-3)), along with the [acwf-verification](https://acwf-verification.materialscloud.org/) online database.

, where the [AiiDA](https://www.aiida.net/) provides high-throughput computing workflow and platform. We plan to substitute part of kernel of APNS with AiiDA in the future, to further standardize the workflow and improve the quality of test data.

## APNS backend workflow configuration

**WARNING:** In principle APNS is not designed to be run on local machine of common ABACUS users. However, if you want to run it on your local machine, you can follow the instructions below.

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

#### ~~Analysis~~

- ~~matplotlib~~
- ~~lbg (lebegue, developed by [DP technology Inc.](https://www.dp.tech/), for downloading groups of jobs from [Bohrium(R) platform](https://bohrium.dp.tech/))~~

#### ~~Orbgen~~

~~Please see the sub-project [ABACUS-ORBGEN](https://github.com/kirk0830/ABACUS-ORBGEN) for more details.~~
*Note: the Orbgen workflow is to be refactored now*

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

APNS will manage local pseudopotentials and numerical atomic orbitals on its own, please read the annotation in `/download/upf/update.py` and `/download/orb/update.py` carefully before running the following commands.

Once you update your local pseudopotential or numerical atomic orbital library, you should update the tag database accordingly.

#### Prepare input JSON file

prepare a JSON file like the following (but please do not include any comments in your JSON file):

```jsonc
{
    "global": {
        "mode": "test",
        "pseudo_dir": "/path/to/your/folder/stores/psp",
        "orbital_dir": "/path/to/your/folder/stores/nao",
        "cache_dir": "/path/where/apns/cache",
        "out_dir": "/path/to/your/output"
    },
    "credentials": {
        "materials_project": {"api_key": "your-materials-project-api-key"}
        // other credential settings not stably-implemented yet
    },
    "abacus": [
        // a list of ABACUS DFT calculation settings, see below for details
    ],
    "atomsets": [
        // a list of definition of atomic species, see below for details
    ],
    "strusets": [
        // a list of definition of structures, see below for details
    ],
}
```

- In `credentials` section, you can configure the structural automatic download from [Materials Project](https://next-gen.materialsproject.org/). Register, login and get your API key from the [dashboard](https://www.materialsproject.org/dashboard).

- In `abacus` section, any number of dicts are allowed, each dict contains the content exactly the same as the input of ABACUS (commonly the `INPUT` file). For a complete list of keywords supported, see the [ABACUS online manual](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html). Beyond the basic setting manner, APNS further supports the `iterator`

    ```json
    "ecutwfc": [20, 30, 40, 50, 60, 70, 80, 90, 100]
    ```

    and `joint`

    ```json
    "ecutwfc|ecutrho": [[30, 240]]
    ```

    and even their combined mode:

    ```json
    "ecutwfc|ecutrho": [[30, 240], [40, 320], [50, 400]]
    ```

    . If more than one `iterator` is given, then Cartesian product will be made. This means that if you have two `iterator`:

    ```json
    "ecutwfc": [20, 30, 40, 50, 60, 70, 80, 90, 100],
    "ecutrho": [240, 320, 400]
    ```

    , then there will be in total 27 calculation will be defined.
    Moreover, the `abacus` section is merely one of the DFT code supported by APNS. APNS also partially supports the `qespresso` (Quantum ESPRESSO), e.g.:

    ```json
    "qespresso": [
        {
            "control": {
                "outdir": "./out",
                "prefix": "test",
                "pseudo_dir": "./",
                "verbosity": "high",
                "restart_mode": "from_scratch",
                "calculation": "scf",
                "tstress": ".true.",
                "tprnfor": ".true."
            },
            "system": {
                "ibrav": 0,
                "ecutwfc|ecutrho": [[20, 160], [30, 240], [40, 320], [50, 400], [60, 480]],
                "occupations": "smearing",
                "smearing": "gaussian",
                "degauss": 0.01
            },
            "electrons": {
                "mixing_mode": "plain",
                "conv_thr": 1.0e-6,
                "diagonalization": "david"
            },
            "ions": {},
            "cell": {}
        }
    ]
    ```

    . More DFT code support can be easily added by modifying the `apns/test/main.py` file (but all relevant implementation will be moved to other files later).

    Both the `abacus` and `qespresso` is called the `calculator` in APNS.

- In `atomsets` section, any number of dicts are allowed, each dict contains the definition of atomic species. A typical definition can be:

    ```json
        {
            "Co": [["NC", "pslibrary"], null],
            "Ag": [["NC"], null],
            "Cd": [["NC"], null],
            "In": [["NC"], null]
        }
    ```

    , in which for Cobalt, all pseudopotentials satisfying the search tags `NC` and `pslibrary` will be used. For other elements, Ag, Cd and In, as long as there is any pseudopotential is tagged with `NC`, it will be used. The second element in the list is for numerical atomic orbitals, which is not needed for PW DFT calculation.

    In an ABACUS LCAO calculation, the atomic species should be defined like:

    ```json
        {
            "H": [["NC", "sg15", "1.0"], ["SimulatedAnnealing"]],
            "O": [["NC", "sg15", "1.0"], ["PTG_dpsi"]]
        }
    ```

    , which means based on searching the pseudopotential(s) that matches the tags `NC` and `sg15`, the corresponding numerical atomic orbitals will be searched with the tag `SimulatedAnnealing` for Hydrogen and `PTG_dpsi` for Oxygen.

    When it is the case that large number of elements with the same tags, you can use:

    ```json
        {
            "__element__": ["Cs", "Ba", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi"],
            "__tags__": [["US"], null]
        }
    ```

    , in which all elements in the list will be searched with the same tags. Here it is all the ultrasoft pseudopotentials available will be used.

- In `strusets` section, any number of dicts are allowed, each dict defines a series of structure that is calculated with the same `calcualtor` and `atomset`. An example is:

    ```json
        {
            "calculator": "abacus", "calcset": 0,
            "atomset": 0,
            "database": "materials_project",
            "desc": [["file", "/mnt/e/Downloads/Co3O4.cif", [1.00], [0.08]]]
        }
    ```

    APNS backend workflow will parse the cif file specified (here the file path is `/mnt/e/Downloads/Co3O4.cif`), prepare input files for `calculator` `abacus`, with the first set of DFT parameters defined in `abacus` section, the first set of atomic species defined in `atomsets` section. the third element (the list) is the characteristic length of the structure, for Co3O4 case it is the scaling factor of the crystal volume. The fourth element is the `kspacing` (see ABACUS online manual for explanation: [kspacing](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#kspacing), here APNS only supports one number, instead of three numbers for each reciprocal axis).

    Once in the `credential` section the Materials Project API key is given, APNS can try to connect with the Materials Project online database, search for the *most stable* structure of Co3O4:

    ```json
        {
            "calculator": "abacus", "calcset": 0,
            "atomset": 0,
            "database": "materials_project",
            "desc": [["search", "Co3O4", [1.00], [0.08]]]
        }
    ```

    APNS also supports building some ideal structures, e.g. for the face-centered cubic Aluminum:

    ```json
        {
            "calculator": "abacus", "calcset": 0,
            "atomset": 0,
            "desc": [["from_scratch", "Al_fcc", [1.00], [0.08]]]
        }
    ```

    A cell volume database should be correctly set in file `apns/test/bravis_and_molecule.py`, please read this file for more details. The second element also supports the oxide building, e.g. for the Al2O3:

    ```json
        {
            "calculator": "abacus", "calcset": 0,
            "atomset": 0,
            "desc": [["from_scratch", "AlO_x2y3", [1.00], [0.08]]]
        }
    ```

    , in which the Al will be at the site 'x' and O will be at the site 'y', and their mole ratio is 2:3. This also relies on the correct setting of the cell volume database. We also support the bulding of simple molecule:

    ```json
        {
            "calculator": "abacus", "calcset": 0,
            "atomset": 0,
            "desc": [["from_scratch", "H_dimer", [0.94], [0.08]]]
        }
    ```

    , in which the H2 molecule will be built with the bond length 0.94 Angstrom. The kspacing here will not be read, the calculation will be automatically set as Gamma-point only.

#### Run APNS-test workflow

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

, your Bohrium account and password are needed.

#### Run APNS-analysis workflow

Since APNS-2, the `analysis` workflow will not be executed from the `main.py`, instead, there are example scripts (for greping, analyzing and ploting data that shown in the APNS website) in `analysis/` folder, in which there are also many useful fundamental postprocessing utilities function implemented. Users are encouraged to write their own scripts based on these examples.
