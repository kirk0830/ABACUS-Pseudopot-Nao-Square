"""Cartesian direct product is crucial for generating high-dimensional cross test.

Workflow overall:

Parse user setting, line1:
                    system -> element
                    pseudopotentials + element -> valid_pseudopotential
                    (numerical_orbitals + valid_pseudopotential -> valid_numerical_orbital)
                    (valid_pseudopotential + valid_numerical_orbital -> pseudopot_nao)
                    
                    line2:
                    calculation -> calculation_settings

                    line3:
                    system -> cif

for system in systems:
    for system_pseudopot_nao_setting in system_pseudopot_nao_settings:
        for calculation_setting in calculation_settings:
            [create a folder for this test] <- system, system_pseudopot_nao_setting, calculation_setting
            [generation of INPUT file] <- calculation_setting
                    calculation_settings -> input
            [generation of STRU and KPT file] <- system_pseudopot_nao_setting, calculation_setting
                    cif + pseudopot_nao + calculation_settings -> stru
                    cif -> kpt
            [copy of pseudopotential (and numerical orbital) files] <- system_pseudopot_nao_setting

"""

import itertools as it

"""pseudopotential/pseudopot-nao combination for one specific system"""

"""This function provides pseudopot-nao combination for one single element.

once valid pseudpotentials are scanned and archived, the rcuts, ecuts and configurations are based on
apns user settings, but can also based on validity scan results, say they could be intersection between
user setting and actual validity scan results, or just the actual validity scan results.
The function pseudopot_nao provides Cartesian direct product between pseudopotential settings and 
numerical orbital setting.

pseudopotentials param here is the list of pseudopotential identifiers.

In principle the other param can also be list of numerical orbital identifiers, but will implement in
the future."""
def pseudopot_nao(pseudopotentials: list, numerical_orbitals: list = None) -> list:
    """seperate valid pseudopotentials (and numerical orbitals) into different combinations for one single element
    
    Args:
        pseudopotentials (list|dict): list of pseudopotential identifiers or dict of numerical orbital identifiers whose keys are pseudopotential identifiers
        
    Returns:
        list: list of dictionaries containing pseudopotential and numerical orbital identifiers
        
    Examples:
    ```python
    >>> pseudopot_nao(pseudopotentials=["sg15_10", "pd_04"])
    # returns:
    [
        {'pseudopotential': 'sg15_10'}, 
        {'pseudopotential': 'pd_04'}
    ]
    >>> pseudopot_nao(pseudopotentials=["sg15_10", "pd_04"],
                      numerical_orbitals=["DZP_6@sg15_10", "TZDP_6@sg15_10", "TZDP_10@pd_04"])
    # returns:
    [
        {'pseudopotential': 'sg15_10', 'numerical_orbital': 'DZP_6'}, 
        {'pseudopotential': 'sg15_10', 'numerical_orbital': 'TZDP_6'}, 
        {'pseudopotential': 'pd_04', 'numerical_orbital': 'TZDP_10'}
    ]
    ```
    """
    if numerical_orbitals is None:
        return [{"pseudopotential": pseudopotential} for pseudopotential in pseudopotentials]
    else:
        return [{"pseudopotential": numerical_orbital.split("@")[1],
                 "numerical_orbital": numerical_orbital.split("@")[0]} for numerical_orbital in numerical_orbitals]

"""This function provides element-wise combination of pseudopotential-nao combination yielded by
function pseudopot_nao. Elements' pseudopotential-nao combination should be stored in a list.

According to the return value of this function, combined with the dictionary "valid_pseudopotential",
It is easy to quickly gather all pseudopotential and numerical orbital files because the dictionary
"valid_pseudopotential"'s first level key is element and second level key is pseudopotential identifier.
With these two key, the returned dictionary has a key "file", which is the absolute path of the file.
So it is for numerical orbitals.

Methodologically, a loop over the returned list will quickly adjust the template file, especially for
ABACUS the "X_pseudopot" and "X_numerical_orbital" keywords in STRU file.
"""
def system(elements: list) -> list:
    """generate all possible combination of pseudopotential and numerical orbital for each element

    Args:
        elements (list): elements' packed pseudopotential and numerical orbital combination, returned by pseudopot_nao

    Returns:
        list: a list of dictionaries containing pseudopotential and numerical orbital identifiers
    
    Examples:
    ```python
    >>> system(elements=[[{"pseudopotential": "sg15_10"}, {"pseudopotential": "pd_04_spd"}],
                         [{"pseudopotential": "dojo_03"}, {"pseudopotential": "pd_03"}]])
    # returns:              element1,   element2
    [
        {'pseudopotential': ['sg15_10', 'dojo_03']},
        {'pseudopotential': ['sg15_10', 'pd_03']},
        {'pseudopotential': ['pd_04_spd', 'dojo_03']},
        {'pseudopotential': ['pd_04_spd', 'pd_03']}
    ]
    
    >>> system(elements=[[{"pseudopotential": "sg15_10", "numerical_orbital": "DZP_6"}, {"pseudopotential": "pd_04_spd", "numerical_orbital": "TZDP_6"}],
                         [{"pseudopotential": "dojo_03", "numerical_orbital": "QZTP_10"}, {"pseudopotential": "pd_03", "numerical_orbital": "SZ_10"}]])
    # returns:              element1,   element2                          element1,   element2
    [
        {'pseudopotential': ['sg15_10', 'dojo_03'], 'numerical_orbital': ['DZP_6', 'QZTP_10']},
        {'pseudopotential': ['sg15_10', 'pd_03'], 'numerical_orbital': ['DZP_6', 'SZ_10']},
        {'pseudopotential': ['pd_04_spd', 'dojo_03'], 'numerical_orbital': ['TZDP_6', 'QZTP_10']},
        {'pseudopotential': ['pd_04_spd', 'pd_03'], 'numerical_orbital': ['TZDP_6', 'SZ_10']}
    ]
    """
    indices = [list(range(len(element))) for element in elements]
    combinations = list(it.product(*indices))
    return [{key: [elements[i][combination[i]][key]\
                    for i in range(len(elements))] \
                    for key in elements[0][0].keys()}\
                    for combination in combinations]


import apns.module_structure.basic as amsb
def systems(systems: dict, vupfs: dict, vorbs: dict = None) -> list:
    """generate all possible combinations for all systems, elemental information such as valid pseudopotentials and
    numerical orbitals should be provided
    
    Args:
        systems (list): a list of system_with_mpid
        valid_pseudopotentials (dict): a dictionary of valid pseudopotentials, which is generated from function in apns.module_pseudo
        valid_numerical_orbitals (dict, optional): a dictionary of valid numerical orbitals, which is generated from function in apns.module_pseudo
        
    Returns:
        list: a list of lists of dictionaries containing pseudopotential and numerical orbital identifiers for each pseudopot-nao combination for each
        system, in sequence of system_list
    """
    #vorbs = {} if vorbs is None else vorbs
    upforb_settings = []
    for formula in systems.keys():
        elements = amsb.scan_elements(formula)
        # then make Cartesian direct product of pseudopotential and numerical orbital
        # for each element, return a list by function pseudopot_nao
        # the .keys() method is used to get all pseudopotentials and numerical orbitals'
        # identifiers, which are used to make Cartesian direct product
        vupf_ids = [list(vupfs[element].keys()) for element in elements]
        vorb_ids = [list(vorbs[element].keys()) if len(vorbs[element].keys()) > 0 \
                    else None for element in elements] if vorbs is not None else None
        element_wise_combinations = [pseudopot_nao(pseudopotentials=vupf_ids[i], 
                                                   numerical_orbitals=None) 
                                     for i in range(len(elements))]
        upforb_settings.append(system(element_wise_combinations))
    return upforb_settings

"""calculation parameters combination for global """

"""This function provides Cartesian direct product of calculation parameters.

All parameters specified in calculation section of work_status or input.json will be expanded to list
of param suites. If a parameter is specified as list, then will be a dimension for direct product,
otherwise will be a scalar and its value will be copied to all param suites.

Methodologically, a loop over the returned list will quickly adjust the template file, especially for
ABACUS the INPUT file. There is a parameter not in INPUT, it is, the characteristic_lengths, is used to scale
LATTICE_CONSTANT in STRU file."""
def calculation(calculation_settings: dict) -> list:

    """Gnerate all possible combination of parameters whose values specified as list to enumerate in 
    work_status/input.json calculation section
    
    Args:
        work_status_calculation_section (dict): calculation section of work_status
        
    Returns:
        list: a list of calculation parameters dictionaries
        
    Examples:
    >>> calculation(
        ...     work_status_calculation_section={
            ...         "basis_type": "pw",
            ...         "kpoints": [1, 2],
            ...         "kpoints_type": ["automatic", "gamma"],
            ...         }
        ...     )
    [
        {'kpoints': 1, 'kpoints_type': 'automatic', 'basis_type': 'pw'}, 
        {'kpoints': 1, 'kpoints_type': 'gamma', 'basis_type': 'pw'},
        {'kpoints': 2, 'kpoints_type': 'automatic', 'basis_type': 'pw'},
        {'kpoints': 2, 'kpoints_type': 'gamma', 'basis_type': 'pw'}
    ]
    """
    # separate list and scalar parameters, for any variable which is given as list, it will be
    # expanded to a list of dictionaries, each dictionary contains all scalar parameters and one
    # list parameter
    lists_names = [key for key in calculation_settings.keys() if isinstance(calculation_settings[key], list)]
    lists_tocombine = [calculation_settings[key] for key in lists_names]
    scalar_names = [key for key in calculation_settings.keys() if key not in lists_names]
    # make Cartesian direct product
    result = [dict(zip(lists_names, combination)) for combination in it.product(*lists_tocombine)]
    # add scalar parameters to each dictionary
    for scalar_name in scalar_names:
        for i in range(len(result)):
            result[i][scalar_name] = calculation_settings[scalar_name]
    return result

def extensive(extensive_settings: dict) -> list:
    """generate all possible combinations of extensive settings set in
    `extensive` section of input.json. Extensive setting may actually
    affect STRU (for example the cell scaling), KPT (for example the
    kpoints grid)
    
    Args:
        extensive_settings (dict): extensive settings in input.json
        
    Returns:
        list: a list of extensive settings dictionaries
    """
    result = []
    
    iterable_keys = ["characteristic_lengths"]
    other_keys = [key for key in extensive_settings.keys() if key not in iterable_keys]
    iterables = [extensive_settings[key] for key in iterable_keys]
    for combination in list(it.product(*iterables)):
        result.append(dict(zip(iterable_keys, combination)))
    for key in other_keys:
        for i in range(len(result)):
            result[i][key] = extensive_settings[key]
    return result

def setup_iterables(**kwargs) -> tuple[list, list, list]:
    """setup iterables, return three iterables representing three dimension to iterate:
    calculation settings, extensive settings, and systems"""
    return systems(kwargs.get("systems", {}), 
                   kwargs.get("pseudopotentials", {}), 
                   kwargs.get("numerical_orbitals", {})), \
           calculation(kwargs.get("calculation_settings", {})), \
           extensive(kwargs.get("extensive_settings", {}))