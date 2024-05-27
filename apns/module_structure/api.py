from apns.module_workflow.identifier import TEMPORARY_FOLDER as cache_dir
def download_all(formula: dict,
                 api_keys: dict,
                 n_structures: dict = None,
                 stru_cache_dir = cache_dir):
    """# API of online crystal structure database
    General interface to download structures from different databases.
    
    ### Usage
    a complete form of function call would be:
    ```python
    from apns.module_structure.api import download
    log = download(formula, n_structures, api_key, database, stru_cache_dir)
    ```
    The `formula` is a dict and should be like:
    ```python
    formula = {"mp": ["Er2O3", "TiO2"],
               "cod": ["V2O5", "Bk"],
               "optimade": ["Ti4O9", "TiO2"]}
    ```
    The `n_structures` is also a dict and should be like:
    ```python
    n_structures = {"mp": [1, 2], "cod": [1, 1], "optimade": [1, 1]}
    ```
    , which means download 1 structure for "Er2O3" and 2 structures for "TiO2" from Materials Project,
    1 structure for "V2O5" and 1 structure for "Bk" from Crystallography Open Database, and 1 structure
    for "Ti4O9" and 1 structure for "TiO2" from OPTIMADE API.
    or the simplest form:
    ```python
    from apns.module_structure.api import download
    log = download(formula, n_structures)
    ```
    Function will return a dict whose keys are formula and values are tuple, the first
    is the file path of downloaded structure, the second is the magnetism information.
    ### Database documentation
    - mp: Materials Project
    - mv: Materials Virtual Lab
    - cod: Crystallography Open Database
    - optimade: OPTIMADE API
    ### Online database credentials
    See their official website for more information. For example the Materials Project, you should
    have an API key to download structures.
    """
    from apns.module_structure.materials_project import download as materials_project
    from apns.module_structure.cod import download as cod
    from apns.module_structure.optimade import download as optimade
    n_structures = n_structures if n_structures is not None else {db: [1] * len(formula[db]) for db in formula.keys()}
    dbs = formula.keys()
    assert all([db in ["materials_project", "matterverse", "cod", "optimade"] for db in dbs]), f'database not supported: {dbs}'
    log = {}
    for db in dbs:
        if db == "materials_project":
            log[db] = materials_project(api_keys[db], formula[db], n_structures[db], stru_cache_dir)
        elif db == "matterverse":
            raise NotImplementedError("Materials Virtual Lab is not implemented yet.")
        elif db == "cod":
            log[db] = cod(formula[db])
        elif db == "optimade":
            log[db] = optimade(formula[db])
    return log

def download(formula: list|str, n_structures: list|int, api_key: str, database: str = "materials_project", cache_dir: str = cache_dir):
    """after all structures needed are collected, in-one-shot download structures as cif, return
    the formula as dict key, value as list of tuples (fname, magmoms)"""
    return download_all({database: formula}, {database: api_key}, {database: n_structures}, cache_dir)[database]