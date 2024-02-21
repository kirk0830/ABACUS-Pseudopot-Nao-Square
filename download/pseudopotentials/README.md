<p align="center">
    <img src="../../docs/assets/images/apns.svg">
</p>  

# ABACUS Pseudopot-Nao Square  
## In brief  
This folder stores downloaded and generated pseudopotentials.  
## Nomenclature  
Following the rule in `sources/module_workflow/README.md`, the nomenclature of pseudopotentials is:  
```text
<kind>_<version>_<appendix>
```
where:  
- `<kind>` is the kind of pseudopotential, like `sg15`, `pd`, `dojo`, etc.  
- `<version>` is the version of pseudopotential, like `10`, `04`, `05`, etc.  
- `<appendix>` is the appendix of pseudopotential, like `fr`, `alt`, etc.  
## To update
To update the pseudopotentials, please follow the steps below:
1. Download the pseudopotentials from the official website or other sources.
2. Store the downloaded pseudopotentials in this folder.
3. Extend functions in `apns/module_pseudo/archive.py`, especially function `determine_kind` and `archive`. If necessary, write a new function whose name startswith `op_` to handle the new kind of pseudopotential, for more details, refer to already implemented functions like `op_SG15_`, `op_PD03_`, `op_PD04_`, etc.
4. Delete the `description.json` in this folder and run `apns/module_pseudo/archive.py` to generate a new `description.json`.