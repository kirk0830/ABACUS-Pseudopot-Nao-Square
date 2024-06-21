"""with abacustest Reuse workflow, eos_pw_vs_lcao module,
can have following properties:
1. energies
2. volumes
3. band structures

with these quantities, can the following properties be calculated:
1. equation-of-state (EOS)
    by non-linearly fitting the Birch-Murnaghan equation of state,
    can get the equilibrium volume, bulk modulus, and its pressure derivative.
    With the fit data, calculate delta value
2. basis set completeness
    by directly comparing the energy minimum of the PW-LCAO pair,
    the energy difference is a good indicator of the basis set completeness.
    The lower the energy difference, the better the basis set completeness.
3. band structure similarity
    similarity is measured by the band structure difference.

For one shot of run, abacustest will produce following files:
abacustest file structure:
```bash
- task-main
- inputs
- outputs\
    - results\
        - test_system_folder_1
            - PW
                - eos-3 <- one single job folder, can be used for band structure similarity test
                - eos-2 <- EOS point "another"
                - ...
                - eos0
                - ...
                - eos3
            - LCAO1
                - eos-3
                - ...
            - ...
            - LCAOn
        - test_system_folder_2
        - ...
        - supermetrics.json
        - post.py
        - metrics.json <- it is where the results are stored (vols, eners, etc.)
        - test_system_result_pic_1
        - test_system_result_pic_2
        - ...
```
In metrics.json, the content is like:
```json
{
    "test_system_folder_1/PW/eos0": {...},
    "test_system_folder_1/PW/eos1": {...},
    //...
}
```
"""