<p align="center">
    <img src="apns.svg">
</p>  

This file introduces basic usage of APNS implemented presently. 
# ABACUS Pseudopot-Nao Square
## Example `input.json`
```json5
{
    // some global settings placed here
    "global": {
        "test_mode": "pseudopotential",           // "pseudopotential" or "numerical_orbital"
        "software": "ABACUS",                     // "ABACUS" or "qespresso"
        "work_dir": "./",                         // the working directory
        "pseudo_dir": "./module_pseudo/download", // always fixed, pseudopotential now is 
                                                  // designed to distribute along with the code
        "orbital_dir": "./module_nao/download",   // will be the directory to store numerical orbital files
        "save_log": true                          // will dump json
    },
    "calculation": {
        "basis_type": "pw",     // "pw" or "lcao". "pw" is supported by both ABACUS and qespresso, 
                                // "lcao" is only supported by ABACUS
        "functionals": ["PBE"], // enumerate! xc functional
        "ecutwfc": [100],       // enumerate! ecutwfc that will impose for all systems
        "cell_scaling": [0.00]  // enumerate! cell scaling factor, NOT IMPLEMENTED YET
    },
    "systems": ["Er2O3"],       // enumerate! system name, will be used to find CIF via Materials Project API.
                                // or specify not symmetrized CIF file path in structures list
    "materials_project": {
        "api_key": "xxx",       // your Materials Project API
        "n_structures": 1,      // number of structures want to download, for each system or specify a list for systems
        "theoretical": false,   // always keep false to get only experimentally observed structures
        "most_stable": true     // get the most stable structure
    },
    // psedupotential section supports two ways input, the first is list, the second is dict.
    // For a list, the specified kinds, versions and appendices will be used for all elements in
    // all systems. For a dict, must specify kinds, versions and appendices for each element.
    "pseudopotentials": {
        "kinds": ["sg15", "pd"],
        "versions": ["10", "04"],
        "appendices": [""]
    },
    // similarly for numerical orbitals.
    "numerical_orbitals": {
        "types": ["DZP"],
        "rcuts": [7, 8, 9, 10],
        "appendices": [""]
    }
}
```

## Full input keywords list
This document provides a full list of keywords for ABACUS Pseudopot-Nao Square. The keywords are divided into several sections, each section is for a specific purpose. The keywords are listed in alphabetical order within each section.

---
### `global` section

#### `test_mode`
* **Description:** mode to test, can be `pseudopotential` or `numerical_orbital`. For `pseudopotential` mode, if `software` specified as `ABACUS`, then `basis_type pw` will be automatically set, if `software` specified as `qespresso`, nothing would happen. `pseudopotential` test mode will only test the pseudopotential, and `numerical_orbital` mode need user specify information of pseudopotential and numerical orbitals both.
* **Type:** string
* **Default:** `pseudopotential`

#### `software`
* **Description:** software to test, can be `ABACUS` or `qespresso`. If `software` specified as `ABACUS`, then `basis_type pw` will be automatically set, if `software` specified as `qespresso`, nothing would happen.
* **Type:** string
* **Default:** `ABACUS`

#### `work_dir`
* **Description:** working directory, all the files will be generated in this directory.
* **Type:** string
* **Default:** `'./'`

#### `pseudo_dir`
* **Description:** directory to store pseudopotential files.
* **Type:** string
* **Default:** `'./module_pseudo/resources/'`

#### `orbital_dir`
* **Description:** directory to store numerical orbital files.
* **Type:** string
* **Default:** `'./module_nao/resources/'`

#### `save_log`
* **Description:** whether to save log (json)
* **Type:** bool
* **Default:** `True`
---
### `calculation` section

#### `basis_type`
* **Description:** basis type, can be `pw` or `lcao`. `software ABACUS` supports both `pw` and `lcao`, `software qespresso` only supports `pw`.
* **Type:** string
* **Default:** `pw`

#### `functionals`
* **Description:** DFT functionals to test, specified as list.
* **Type:** list
* **Default:** `['PBE']`

#### `ecutwfc`
* **Description:** kinetic energy cutoffs for wavefunctions, in Rydberg, specified as list, but instead this keyword can be set via ABACUSTEST application, see https://labs.dp.tech/projects/abacustest/
* **Type:** list
* **Default:** `[100]`

#### `cell_scaling`
* **Description:** scaling factor for cell, specified as list.
* **Type:** list
* **Default:** `[0.00]`
---
### `systems` section
This section should be organized in the following way:
```json
{
    "systems": ["Yb2O3", "Er", "TiO2"]
}
```
Elements will be recognized from this section, then `pseudopotentials` and `numerical_orbitals` section will search available files for these elements. For analysis task, this section will have more options, e.g.:
```json
{
    "systems": {
        "Yb2O3": {
            "band_gap": 5.0,
            "lattice_constant": 10.0,
        },
        "Er": {
            "band_gap": 5.0,
            "lattice_constant": 10.0,
        },
    }
}
```

---
### `materials_project` section

#### `api_key`
* **Description:** Materials Project Application Program Interface (API) Key. Get your API key following tutorial in https://docs.materialsproject.org/downloading-data/using-the-api
* **Type:** str
*This keyword is compulsory for running APNS test suite generation, no default value is given. If not specified, will raise error*

#### `n_structures`
* **Description:** For each system specified, how many structures are needed. If a list is given, then the length of the list should be the same as the length of `systems` section. If an integar is given, then will use the same number for all systems.
* **Type:** int or list
* **Default:** `1`

#### `theoretical`
* **Description:** Whether to get theoretical structures. If `True`, then will get theoretical structures, if `False`, then will only get experimental structures.
* **Type:** boolean
* **Default:** `False`

#### `most_stable`
* **Description:** If only get structures which are marked as `is_stable` in Materials Project. If `True`, then will only get the most stable structure, if `False`, then will get all structures. The `is_stable` here sometimes means *the most stable*, some sub-stable structures will be filtered out, like TiO2 anatase.
* **Type:** boolean
* **Default:** `False`

---
### `pseudopotentials` section

#### `kinds`
* **Description:** kinds of pseudopotentials to test (`sg15`, `pd`, `rrkjus` or something other), specified as list. If contents are not organized element-by-element, then will use the same pseudopotential for all elements.
* **Type:** list or dict
* **Default:** `[""]`

#### `versions`
* **Description:** versions of pseudopotentials to test, specified as list. If contents are not organized element-by-element, then will use the same version for all elements.
* **Type:** list or dict
* **Default:** `[""]`

#### `appendices`
* **Description:** appendices of pseudopotentials to test, specified as list. If contents are not organized element-by-element, then will use the same appendix for all elements.
* **Type:** list or dict
* **Default:** `[""]`

---
### `numerical_orbitals` section

#### `types`
* **Description:** types of numerical orbitals to test (`DZP` or `TZDP`, conventionally), specified as list. If contents are not organized element-by-element, then will use the same type for all elements.
* **Type:** list or dict
* **Default:** `[""]`

#### `rcuts`
* **Description:** cutoff radii of numerical orbitals to test, specified as list. If contents are not organized element-by-element, then will use the same cutoff radius for all elements.
* **Type:** list or dict
* **Default:** `[""]`

#### `appendices`
* **Description:** appendices of numerical orbitals to test, specified as list. If contents are not organized element-by-element, then will use the same appendix for all elements.
* **Type:** list or dict
* **Default:** `[""]`