# Module description
This module provides basic, elemental information, to use data stored, 
```python
from module_database import database as db

db.get_element_index('H')
```

# Available data
## Elements
### Fundamental
- Get *element symbol* by *atomic number*:  
    `get_element_symbol(1)`
- Get *atomic number* by *element symbol*:  
    `get_element_index('H')`
### NIST data
- Get *atom mass* by *element symbol*:  
    `element_mass('H')`
## Quantum chemistry
- Get *angular momentum number* by *sublayer symbol* (s, p, d, etc):  
    `l_to_sublayer(1)`
- Get *sublayer symbol* by *angular momentum number* (0, 1, 2, etc):  
    `sublayer_to_l('s')`
- Get *number notation of multiplicity* by *multiplicity symbol* (s, d, t, etc):  
    `multiplicity_to_number('singlet')`
- Get *multiplicity symbol* by *number notation of multiplicity* (1, 2, 3, etc):  
    `number_to_multiplicity(1)`