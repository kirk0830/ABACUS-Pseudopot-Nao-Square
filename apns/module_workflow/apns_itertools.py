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
def pseudopot_nao(pseudopotentials: list, numerical_orbitals: list = []) -> list:
    """seperate valid pseudopotentials (and numerical orbitals) into different combinations for one single element
    
    Args:
        pseudopotentials (list|dict): list of pseudopotential identifiers or dict of numerical orbital identifiers whose keys are pseudopotential identifiers
        
    Returns:
        list: list of dictionaries containing pseudopotential and numerical orbital identifiers
        
    Examples:
    >>> pseudopot_nao(
        ...     pseudopotentials=["sg15_10", "pd_04"]
        ... )
    [
        {'pseudopotential': 'sg15_10'}, 
        {'pseudopotential': 'pd_04'}
    ]
    >>> pseudopot_nao(
        ...     pseudopotentials=["sg15_10", "pd_04"],
        ...     numerical_orbitals=["DZP_6@sg15_10", "TZDP_6@sg15_10", "TZDP_10@pd_04"]
        ... )
    [
        {'pseudopotential': 'sg15_10', 'numerical_orbital': 'DZP_6'}, 
        {'pseudopotential': 'sg15_10', 'numerical_orbital': 'TZDP_6'}, 
        {'pseudopotential': 'pd_04', 'numerical_orbital': 'TZDP_10'}
    ]
    """
    if len(numerical_orbitals) == 0:
        return [
            { "pseudopotential": pseudopotential } for pseudopotential in pseudopotentials
            ]
    else:
        return [
            {
                "pseudopotential": numerical_orbital.split("@")[1],
                "numerical_orbital": numerical_orbital.split("@")[0]
            } for numerical_orbital in numerical_orbitals
            ]

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
    >>> system(
        ...     elements=[
        ...                 [{"pseudopotential": "sg15_10"}, {"pseudopotential": "pd_04_spd"}],
        ...                 [{"pseudopotential": "dojo_03"}, {"pseudopotential": "pd_03"}]
        ...             ]
        ... )
    [
        {'pseudopotential': ['sg15_10', 'dojo_03']},
        {'pseudopotential': ['sg15_10', 'pd_03']},
        {'pseudopotential': ['pd_04_spd', 'dojo_03']},
        {'pseudopotential': ['pd_04_spd', 'pd_03']}
    ]
    
    >>> system(
        ...     elements=[
        ...                 [{"pseudopotential": "sg15_10", "numerical_orbital": "DZP_6"}, {"pseudopotential": "pd_04_spd", "numerical_orbital": "TZDP_6"}],
        ...                 [{"pseudopotential": "dojo_03", "numerical_orbital": "QZTP_10"}, {"pseudopotential": "pd_03", "numerical_orbital": "SZ_10"}]
        ...             ]
        ... )
    [
        {'pseudopotential': ['sg15_10', 'dojo_03'], 'numerical_orbital': ['DZP_6', 'QZTP_10']},
        {'pseudopotential': ['sg15_10', 'pd_03'], 'numerical_orbital': ['DZP_6', 'SZ_10']},
        {'pseudopotential': ['pd_04_spd', 'dojo_03'], 'numerical_orbital': ['TZDP_6', 'QZTP_10']},
        {'pseudopotential': ['pd_04_spd', 'pd_03'], 'numerical_orbital': ['TZDP_6', 'SZ_10']}
    ]
    """
    indices = [
        list(range(len(element))) for element in elements
    ]
    combinations = list(it.product(*indices))
    result = []
    for combination in combinations:
        result.append(
            { key: [elements[i][combination[i]][key] for i in range(len(elements))] for key in elements[0][0].keys() }
        )
    return result

import apns.module_structure.basic as amsb
def systems(system_list: list, valid_pseudopotentials: dict = {}, valid_numerical_orbitals: dict = {}) -> list:
    print("valid_pseudopotentials: ", valid_pseudopotentials)
    pseudopot_nao_settings = []
    for _system in system_list:
        elements = amsb.scan_elements(_system)
        elements = [
            pseudopot_nao(
                pseudopotentials=list(valid_pseudopotentials[element].keys()),
                numerical_orbitals=list(valid_numerical_orbitals[element].keys())
                ) for element in elements
                ]
        pseudopot_nao_settings.append(system(elements=elements))
    
    print("pseudopot_nao_settings: ", pseudopot_nao_settings)
    return pseudopot_nao_settings
"""calculation parameters combination for global """

"""This function provides Cartesian direct product of calculation parameters.

All parameters specified in calculation section of work_status or input.json will be expanded to list
of param suites. If a parameter is specified as list, then will be a dimension for direct product,
otherwise will be a scalar and its value will be copied to all param suites.

Methodologically, a loop over the returned list will quickly adjust the template file, especially for
ABACUS the INPUT file. There is a parameter not in INPUT, it is, the characteristic_lengths, is used to scale
LATTICE_CONSTANT in STRU file."""
def calculation(work_status_calculation_section: dict) -> list:

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

    wscs = work_status_calculation_section
    lists_to_combine = []
    lists_names = []
    scalar_names = []
    for key in wscs.keys():
        if isinstance(wscs[key], list):
            lists_to_combine.append(wscs[key])
            lists_names.append(key)
        else:
            scalar_names.append(key)

    result = []
    for combination in list(it.product(*lists_to_combine)):
        result.append(dict(zip(lists_names, combination)))
    for scalar_name in scalar_names:
        for i in range(len(result)):
            result[i][scalar_name] = wscs[scalar_name]

    return result


