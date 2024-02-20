---
published: false
---

<p align="center">
    <img src="../../docs/apns.svg">
</p>

# Module description
This module provides basic, elemental information, to use data stored, 
```python
from apns.module_database import database as db

db.element_label_toindex('H')
```

# Available data
## Elements
### Fundamental
- Get *element symbol* by *atomic number*:  
    `get_element_symbol(1)`
- Get *atomic number* by *element symbol*:  
    `element_label_toindex('H')`
### NIST data
- Get *atom mass* by *element symbol*:  
    `element_label_tomass('H')`
## Quantum chemistry
- Get *angular momentum number* by *sublayer symbol* (s, p, d, etc):  
    `l_tosymbol(1)`
- Get *sublayer symbol* by *angular momentum number* (0, 1, 2, etc):  
    `symbol_tol('s')`
- Get *number notation of multiplicity* by *multiplicity symbol* (s, d, t, etc):  
    `multiplicity_tonumber('singlet')`
- Get *multiplicity symbol* by *number notation of multiplicity* (1, 2, 3, etc):  
    `number_tomultiplicity(1)`