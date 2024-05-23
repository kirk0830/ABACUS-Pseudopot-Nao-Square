from apns.module_structure.materials_project import download as mpdwld
from apns.module_structure.cod import download as coddwld
from apns.module_structure.optimade import download as optimadedwld
from apns.module_workflow.identifier import TEMPORARY_FOLDER as cache_dir
def download(formula: dict,
             n_structures: dict = None,
             api_keys: str = "",
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
    n_structures = n_structures if n_structures is not None else {db: [1] * len(formula[db]) for db in formula.keys()}
    api_keys = api_keys if api_keys != "" else {db: "" for db in formula.keys()}
    dbs = formula.keys()
    assert all([db in ["mp", "mv", "cod", "optimade"] for db in dbs]), f'database not supported: {dbs}'
    log = {}
    for db in dbs:
        if db == "mp":
            log[db] = mpdwld(api_keys[db], formula[db], n_structures[db], stru_cache_dir)
        elif db == "mv":
            raise NotImplementedError("Materials Virtual Lab is not implemented yet.")
        elif db == "cod":
            log[db] = coddwld(formula[db])
        elif db == "optimade":
            log[db] = optimadedwld(formula[db])
    return log