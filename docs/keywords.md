# ABACUS Pseudopot-Nao Square
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