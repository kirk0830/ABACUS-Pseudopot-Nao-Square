---
published: false
---

# ABACUS Pseudopot-Nao Square  
This folder contains tools distributed for user.  
  
## Tools list  
### Orbital oscillation remove  
#### In brief  
Presently due to inability of opitmizer, some truncated spherical Bessel functions are failed to remove, or say their coefficients are not successfully optimized due to ill-condition of the spillage optimization problem.  
Therefore there will be unphysical oscillation in orbital and can be removed by inversely solve the coefficients of truncated spherical Bessel function and seletively remove some of them, then generate new orbitals for calculation.  
#### Usage  
1. run `nao_oscillation_remove.py` in interactive mode, find the optimal numbers of Bessel functions to remove.  
2. run again with numbers obtained in the previous step, generate orbitals.  
