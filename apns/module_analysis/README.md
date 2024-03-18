---
published: false
---

<p align="center">
    <img src="../../docs/assets/images/apns.svg">
</p>  

# ABACUS Pseudopot-Nao Square
Module: analysis
This module is for greping results and perform various kinds of postprocess, rendering results for presenting.

## Design notes
### 2024.03.16
For kinds of analysis are relatively unpredictable, it might be better to let user specify a code, for specifying what kind of analysis want to perform, and in folders of this module, the analysis scripts should be dynamically managed and archived. If the driver really find the analysis related script in folder, then execute, otherwise raise an error and exit.  
To dynamically manage the analyzer, a json file is needed. The code for specifying analyzer would be arbitrary, like `DOSAnalyzer_20240316` or something, only a uniformed interface is compulsory for volunteer developer to write, so that the main driver can find and execute desired analysis.
```bash
module_analysis
|- description.json # in this file, code -> analysis script mapping stores
|- drivers_collection/ # in this folder, store implemented scripts
   |- APNS_EcutwfcConv.py
   |- APNS_StressConv.py
   |- APNS_DOSIntegral.py
   |- ...
|- ... # do not change things in other folder
```
In corresponding `input.json` input script, specify target folder where results stores, and analysis drivers want to execute:
```json
{
    "global": {
        "test_mode": "analysis",
        "pseudo_dir": "./download/pseudopotentials",
        "orbital_dir": "./download/numerical_orbitals"
    },
    "analysis": {
        "search_domain": "/root/apns/BohriumJob-123456/",
        "items": [
                    "APNS_EcutwfcConv", 
                    "APNS_StressConv", 
                    "APNS_DOSIntegral"
                ]
    }
}
```