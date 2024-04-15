import apns.module_workflow.identifier as amwi
import json
import os

def lookup(formula: str, n_structures: int, with_magmom: bool,
           flocal_db: str = amwi.TEMPORARY_FOLDER + "/mpapi_search_cache.json"):
    """look up the structure search log cached locally, the file
    temporary_dir/mpapi_search_cache.json, return queried results.

    example:
    ```python
    result = lookup("Fe2O3", 5, True)
    ```
    return:
    ```python
    [("./apns_cache/mp-1234.cif", [5.0, 5.0, 5.0]), 
     ("./apns_cache/mp-5678.cif", [4.0, 4.0, 4.0]),
     ("./apns_cache/mp-9101.cif", [3.0, 3.0, 3.0]), 
     ("./apns_cache/mp-1121.cif", [2.0, 2.0, 2.0]),
     ("./apns_cache/mp-3141.cif", [1.0, 1.0, 1.0])]
    # the sequence of elements in the list following file path would be in accord with
    # the sequence of atoms in cif file, this is guaranteed by Materials Project API.
    # therefore, in example result the lengths of magmom are uniformly to be 3 indicates
    # there are uniformly 3 atoms in each cif file.
    ```
    """
    if not os.path.exists(flocal_db):
        with open(flocal_db, "w") as f:
            json.dump({}, f)
    with open(flocal_db, "r") as f:
        local_db = json.load(f)
    if formula in local_db:
        mpids = list(local_db[formula].keys()) # would be mpid
        mpids = [mpid for mpid in mpids if local_db[formula][mpid].get("magmom", None) is not None]\
            if with_magmom else mpids
        mpids = mpids[:n_structures] if len(mpids) > n_structures else mpids
        result = [(local_db[formula][mpid]["file"], local_db[formula][mpid]["magmom"]) for mpid in mpids]
        return result
    else:
        return []

def update(formula: str, mpid: str|list[str], fcif_magmom: tuple[str, list]|list[tuple[str, list]],
           flocal_db: str = amwi.TEMPORARY_FOLDER + "/mpapi_search_cache.json"):
    """update the structure search log cached locally, the file
    temporary_dir/mpapi_search_cache.json, with new search results.

    example:
    ```python
    update("Fe2O3", ("./apns_cache/mp-1234.cif", [5.0, 5.0, 5.0]))
    ```
    """
    if not os.path.exists(flocal_db):
        with open(flocal_db, "w") as f:
            json.dump({}, f)
    with open(flocal_db, "r") as f:
        local_db = json.load(f)
    
    if isinstance(mpid, str) and isinstance(fcif_magmom, tuple):
        mpid, fcif_magmom = [mpid], [fcif_magmom]
    elif isinstance(mpid, list) and isinstance(fcif_magmom, list):
        pass
    else:
        raise ValueError("mpid and fcif_magmom should be both list or both str and tuple")
    for m, fcif_m in zip(mpid, fcif_magmom):
        local_db.setdefault(formula, {}).update({m: dict(zip(["file", "magmom"], fcif_m))})
    
    with open(flocal_db, "w") as f:
        json.dump(local_db, f, indent=4)

def delete(formula: str, mpid: str,
           flocal_db: str = amwi.TEMPORARY_FOLDER + "/mpapi_search_cache.json"):
    """delete the structure search log cached locally, the file
    temporary_dir/mpapi_search_cache.json, with new search results.

    example:
    ```python
    delete("Fe2O3", "mp-1234")
    ```
    """
    if not os.path.exists(flocal_db):
        with open(flocal_db, "w") as f:
            json.dump({}, f)
    with open(flocal_db, "r") as f:
        local_db = json.load(f)
    
    if formula in local_db:
        if mpid in local_db[formula]:
            local_db[formula].pop(mpid)
    
    with open(flocal_db, "w") as f:
        json.dump(local_db, f, indent=4)

import apns.module_structure.materials_project as amsmp
def search(api_key: str,
           formula: str|list[str],
           n_structures: int|list[int],
           **kwargs):
    """search structures from Materials Project API, and cache the results locally
    
    Compulsory domains:
    formula: str|list[str], the formula of the compound
    n_structures: int|list[int], the number of structures to search

    return:
    ```python
    {
        "formula": [
            ("./apns_cache/mp-1234.cif", [5.0, 5.0, 5.0]),
            ("./apns_cache/mp-5678.cif", [4.0, 4.0, 4.0]),
            ("./apns_cache/mp-9101.cif", [3.0, 3.0, 3.0]),
            ("./apns_cache/mp-1121.cif", [2.0, 2.0, 2.0]),
            ("./apns_cache/mp-3141.cif", [1.0, 1.0, 1.0])
        ],
        "formula": [
            ("./apns_cache/mp-1234.cif", [5.0, 5.0, 5.0]),
            ("./apns_cache/mp-5678.cif", [4.0, 4.0, 4.0]),
            ("./apns_cache/mp-9101.cif", [3.0, 3.0, 3.0]),
            ("./apns_cache/mp-1121.cif", [2.0, 2.0, 2.0]),
            ("./apns_cache/mp-3141.cif", [1.0, 1.0, 1.0])
        ],
        #...
    }
    ```
    """
    assert api_key is not None, "api_key is required, for more information, see Materials Project API guideline"

    if isinstance(formula, str):
        formula = [formula]
    if isinstance(n_structures, int):
        n_structures = [n_structures]
    with_magmom = kwargs.get("with_magmom", True)
    if isinstance(with_magmom, bool):
        with_magmom = [with_magmom]
    is_stable = kwargs.get("is_stable", True)
    if isinstance(is_stable, bool):
        is_stable = [is_stable]
    theoretical = kwargs.get("theoretical", False)
    if isinstance(theoretical, bool):
        theoretical = [theoretical]

    assert len(formula) == len(n_structures) == len(with_magmom) == len(is_stable) == len(theoretical)
    # because use Materials Project API, must decrease the frequency of query as much as possible
    # therefore, the search results should be cached locally, and print warning if the search
    # is just 1.
    if len(formula) == 1:
        # the input() function used here is risky when APNS workflow is not run
        # interactively. 
        answer = input("""Materials Project API Usage Warning:
Due to avoid unexpected block of your computer IP, please make sure presently you only have one search.
Otherwise APNS would suggest you search all needed structures at once. To do this,
feed list instead of single string of formula into this function.
Continue with ONE SINGLE search? (y/n)""")
        exit(0) if answer != "y" else None
    result = amsmp.download(api_key=api_key, 
                            formula=formula, 
                            n_structures=n_structures, 
                            stru_cache_dir=amwi.TEMPORARY_FOLDER,
                            is_stable=is_stable, 
                            theoretical=theoretical)
    for formula in result.keys():
        for impid in range(len(result[formula])):
            fcif, magmom = result[formula][impid]
            mpid = result[formula][impid][0].split("/")[-1].split(".")[0]
            update(formula, mpid, (fcif, magmom))
    return result
    
import unittest
class TestStructureManage(unittest.TestCase):

    def test_lookup_update_delete(self):
        result = lookup("Fe2O3", 5, True, "./test_cache.json")
        self.assertEqual(result, [])
        update("Fe2O3", "mp-1234", ("./apns_cache/mp-1234.cif", [5.0, 5.0, 5.0]), "./test_cache.json")
        result = lookup("Fe2O3", 5, True, "./test_cache.json")
        self.assertEqual(result, [("./apns_cache/mp-1234.cif", [5.0, 5.0, 5.0])])
        update("Fe2O3", "mp-5678", ("./apns_cache/mp-5678.cif", [4.0, 4.0, 4.0]), "./test_cache.json")
        result = lookup("Fe2O3", 5, True, "./test_cache.json")
        self.assertEqual(result, [("./apns_cache/mp-1234.cif", [5.0, 5.0, 5.0]), 
                                  ("./apns_cache/mp-5678.cif", [4.0, 4.0, 4.0])])
        update("Fe2O3", "mp-9101", ("./apns_cache/mp-9101.cif", [3.0, 3.0, 3.0]), "./test_cache.json")
        result = lookup("Fe2O3", 5, True, "./test_cache.json")
        self.assertEqual(result, [("./apns_cache/mp-1234.cif", [5.0, 5.0, 5.0]), 
                                  ("./apns_cache/mp-5678.cif", [4.0, 4.0, 4.0]), 
                                  ("./apns_cache/mp-9101.cif", [3.0, 3.0, 3.0])])
        delete("Fe2O3", "mp-5678", "./test_cache.json")
        result = lookup("Fe2O3", 5, True, "./test_cache.json")
        self.assertEqual(result, [("./apns_cache/mp-1234.cif", [5.0, 5.0, 5.0]), 
                                  ("./apns_cache/mp-9101.cif", [3.0, 3.0, 3.0])])
        delete("Fe2O3", "mp-1234", "./test_cache.json")
        result = lookup("Fe2O3", 5, True, "./test_cache.json")
        self.assertEqual(result, [("./apns_cache/mp-9101.cif", [3.0, 3.0, 3.0])])
        delete("Fe2O3", "mp-9101", "./test_cache.json")
        result = lookup("Fe2O3", 5, True, "./test_cache.json")
        self.assertEqual(result, [])
        os.remove("./test_cache.json")

if __name__ == "__main__":
    unittest.main()