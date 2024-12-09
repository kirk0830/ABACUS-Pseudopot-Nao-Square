def lookup_cache(formula: str, n_structures: int = 1, cache_dir: str = "./apns_cache"):
    """# Lookup structures in the cache.
    Args:
        cache_dir (str): Path to the cache directory.
        formula (str): Chemical formula.
        n_structures (int): Number of structures to return.
        database (str): Database to search in.
    Returns:
        list: List of structures.
    """
    import os
    import json
    
    if not os.path.exists(cache_dir):
        return []
    fcache = os.path.join(cache_dir, "structures.json")
    if not os.path.exists(fcache):
        return []
    with open(fcache, "r") as f:
        structures = json.load(f)
    if formula not in structures:
        return []
    if n_structures > len(structures[formula]):
        print(f"Warning: Requested {n_structures} structures, but only {len(structures[formula])} structures are available.")
        return structures[formula]
    
    return [(data["file"], data["magmom"]) for data in structures[formula][:n_structures]]

def update_cache(search_result, cache_dir: str = "./apns_cache"):
    """# Update local cache by search result
    write back the search result returned by for instance the `search()` function. 
    Search result is organized as:
    ```python
    {
        "formula1": [
            (fname11, magmom11),
            (fname12, magmom12), 
            ...
        ],
        "formula2": [
            (fname21, magmom21),
            (fname22, magmom22), 
            ...
        ],
        ...
    }
    ```
    the record in cache is organized as:
    ```python
    {
        "formula1": [
            [fname11, magmom11],
            [fname12, magmom12],
            ...
        ],
        "formula2": [
            [fname21, magmom21],
            [fname22, magmom22],
            ...
        ],
        ...
    }
    ```
    """
    import os
    import json
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)
    fcache = os.path.join(cache_dir, "structures.json")
    if os.path.exists(fcache):
        with open(fcache, "r") as f:
            structures = json.load(f)
    else:
        structures = {}
    for formula in search_result.keys():
        if formula not in structures:
            structures[formula] = []
        for value in search_result[formula]:
            structures[formula].append({"file": value[0], "magmom": value[1]})
    with open(fcache, "w") as f:
        json.dump(structures, f, indent=4)

def search(formula: list|str, n_structures: list|int, api_key: str, cache_dir: str = "./apns_cache", database: str = "materials_project"):
    """# Search
    First look up the local cache to check if there are already enough number of structures of 
    the formula in query, if so, give the first `n_structures`, otherwise will deposite
    from online database, such as `materials_project`, `matterverse`, `cod` or `optimade`.
    Return the research result as a dict, whose keys are formula and values are list of exact
    structures: tuple (fname, magmom), where `fname` is the absolute path to the structure file
    and `magmom` is the magnetic moment of each site.
    
    ### Example
    ```python
    formula = ["Li2O", "Li2O2"]
    n_structures = [1, 2]
    api_key = "your_api_key"
    cache_dir = "./apns_cache"
    database = "materials_project"
    result = search(formula, n_structures, api_key, cache_dir, database)
    ```
    """
    from apns.structure.api import download

    if isinstance(formula, str):
        print("Warning: only one formula as query, high-frequency query is not recommended.")
        assert isinstance(n_structures, int), "n_structures should be an integer."
        formula = [formula]
        n_structures = [n_structures]
    # then handle with the general case
    assert len(formula) == len(n_structures), "Length of formula and n_structures should be the same."
    result = {}
    ifo_to_download = []
    for i in range(len(formula)):
        structures = lookup_cache(formula[i], n_structures[i], cache_dir)
        if len(structures) < n_structures[i]:
            print(f"""Warning: Only {len(structures)} structures of {formula[i]} are available in the cache.
          will search {n_structures[i]} structures from {database} database.""")
            ifo_to_download.append(i)
        else:
            print(f"{n_structures[i]} structures of {formula[i]} are available in local cache.")
            result[formula[i]] = [(s[0], s[1]) for s in structures]
    if len(ifo_to_download) > 0:
        formula_to_download = [formula[i] for i in ifo_to_download]
        n_structures_to_download = [n_structures[i] for i in ifo_to_download]
        log = download(formula_to_download, n_structures_to_download, api_key, database, cache_dir)
        update_cache(log, cache_dir)
        print("""Search result has been updated to the cache, next time you can directly use it without download again.
Privacy: please delete your api_key information in input files in case you leave it in your pull request.""")
        result.update(log)
    return result

def initialize(inp: dict):
    """# Initialize structure download query
    collect the demands of structure search, return two dicts, one for demands, one for coordinates of demands in the input dict.
    
    Args:
        inp (dict): input dict.
    
    Returns:
        tuple: two dicts, one for demands, one for coordinates of demands in the input dict.
    
    Example:
    ```python
    # given the inp
    inp = {
        "strusets": [
            {"database": "mp",
             "desc": [["search", "Li2O", [1.0], [0.08]],
                      ["search", "Li2O2", [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.10]]]},
            {"database": "mp",
             "desc": [["search", "Li2O3", [1.0], [0.08]],
                      ["local", "/path/to/structure", [1.0], [0.08]]]},
            {"database": "mp",
             "desc": [["search", "Li2O", [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.05]]]}
        ]
    }
    # the function will return
    demands = {"mp": ["Li2O", "Li2O2", "Li2O3"]}
    coords = {("mp", "Li2O"):  [(0, [1.0], [0.08]), (2, [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.05])], 
              ("mp", "Li2O2"): [(0, [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.10])], 
              ("mp", "Li2O3"): [(1, [1.0], [0.08])]}
    ```
    """
    demands = {}
    db_formu_traits = {} # (database, formula) -> (istruset, scales, kspacing) mapping
    for i, struset in enumerate(inp["strusets"]): # will raise if there is no this section
        database = struset.get("database", "materials_project")
        for desc in struset["desc"]: # will raise if there is no this section
            if desc[0] == "search":
                demands.setdefault(database, []).append(desc[1])
                db_formu_traits.setdefault((database, desc[1]), []).append((i, desc[2], desc[3])) # istruset and kspacing
    return {k: list(dict.fromkeys(v)) for k, v in demands.items()}, db_formu_traits

def prepare(inp: dict):
    demands, traits = initialize(inp)
    cache_dir = inp.get("cache_dir", "./apns_cache")
    print("* * * Structure search and download * * *".center(100))
    print(f"Cache directory where structures are download, stored and local searched: {cache_dir}")
    search_result = {}
    for db, formula in demands.items():
        print(f"""Structures queried from `{db}` database:
{", ".join(formula)}""")
        n_structures = [1] * len(formula)
        api_key = inp.get("credentials", {}).get(db, {}).get("api_key", "")
        search_result[db] = search(formula, n_structures, api_key, cache_dir, db)
    # remove all desc starting with "search" by overwrite
    for struset in inp["strusets"]:
        struset["desc"] = [{k: v for k, v in dict(zip(["config", "pertmags", "kspacing", "magmoms"], desc[1:])).items() if v is not None}\
             for desc in struset["desc"] if desc[0] != "search"]
    # append the search result to the input dict
    for k, v in traits.items(): # (database, formula) -> [(istruset, scaling, kspacing), ...] mapping
        db, formula = k
        result = search_result[db][formula] # [(fname, magmom), ...]
        for i, scales, kspacing in v:
            inp["strusets"][i]["desc"].extend(\
                [dict(zip(["config", "pertmags", "kspacing", "magmoms"], [r[0], scales, kspacing, r[1]])) for r in result])
    return inp

import unittest
class TestAPNSIO(unittest.TestCase):
    def test_initialize(self):
        inp = {
            "strusets": [
                {"database": "mp",
                 "desc": [["search", "Li2O", [1.0], [0.08, 0.05]],
                          ["search", "Li2O2", [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.10, 0.08, 0.06]]]},
                {"database": "mp",
                 "desc": [["search", "Li2O3", [1.0], [0.1]],
                          ["local", "/path/to/structure", [1.0], [0.1]]]},
                {"database": "mp",
                 "desc": [["search", "Li2O", [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.05, 0.08, 0.1]]]},
                {"database": "pymatgen.cod",
                 "desc": [["search", "Li2O3", [1.0], [0.08]]]},
            ]
        }
        demands, traits = initialize(inp)
        self.assertEqual(demands, {"mp": ["Li2O", "Li2O2", "Li2O3"], "pymatgen.cod": ["Li2O3"]})
        self.assertEqual(traits, {("mp", "Li2O"): [(0, [1.0], [0.08, 0.05]), (2, [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.05, 0.08, 0.1])], 
                                 ("mp", "Li2O2"): [(0, [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.1, 0.08, 0.06])], 
                                 ("mp", "Li2O3"): [(1, [1.0], [0.1])],
                                 ("pymatgen.cod", "Li2O3"): [(3, [1.0], [0.08])]})
    def test_prepare(self):
        """this testcase illustrates the way how to convert user setting to backend function used parameters.
        user can specify descriptor of structure in order like, 
        - identifier, can be `search`, `from_scratch` and `local`. `search` means search structure from online database, in this case,
        user should explicitly define which database to search with, or default to `materials_project`. `from_scratch` means draw structure from scratch,
        can be a molecule or a crystal. `local` means user has a structure file and want to use it.
        - formula/fcif/config: for `search`, type chemical formula, for `from_scratch`, type like `Si_dimer`, `Fe_bcc`, `FeO_x2y3`. for `local`, type the path to the structure file.
        - scales: for periodic system cases, this list will be used to rescale the vol of cell size, useful when perform EOS test. for isolated system, this list will be taken as
        bond lengths, useful when scanning PES of molecules.
        - kspacing: for periodic system cases, this list will be used to set the spacing between kpoints. For isolated system, this keyword should not be set."""
        inp = {
            "credentials": {"materials_project": {"api_key": ""}}, # unless you have an API key and fill in for the first time, you cannot pass this test
            "strusets": [
                {"database": "materials_project",
                 "desc": [["search", "Li2O", [1.0], [0.08, 0.05]],  # the case user only wants to search one structure, scale vol with 1, kspacing can be 0.08 or 0.05
                          ["from_scratch", "Na_dimer", [2.1, 2.3, 2.5, 2.7, 3.0, 3.5, 4.0]], # the case user draw dimer of Na, bond length can be 2.1, 2.3, 2.5, 2.7, 3.0, 3.5, 4.0
                          ["search", "Na2O2", [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.10, 0.08, 0.06]], # the case user wants to search Na2O2, scale vol with 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, kspacing can be 0.10, 0.08, 0.06
                          ["from_scratch", "NaO_xy2", [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.10, 0.08, 0.06]]]}, # the case user draw NaO2, scale vol with 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, kspacing can be 0.10, 0.08, 0.06
                {"database": "materials_project", 
                 "desc": [["from_scratch", "VO_x2y5", [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.12]], # the case user draw VO2.5, scale vol with 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, kspacing can be 0.12
                          ["search", "Al2O3", [1.0], [0.1]], # the case user wants to search Al2O3, scale vol with 1.0, kspacing = 0.1
                          ["local", "/path/to/structure", [1.0], [0.1]]]}, # the case user has a structure file, scale vol with 1.0, kspacing = 0.1
                {"database": "materials_project",
                 "desc": [["search", "Li2O", [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], [0.05, 0.08, 0.1]]]}, # the case user wants to search Li2O, scale vol with 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, kspacing can be 0.05, 0.08, 0.1
                {"database": "materials_project",
                 "desc": [["search", "Li2O", [1.0], [0.08]]]}, # the case user wants to search Li2O, scale vol with 1.0, kspacing = 0.08
            ]
        }
        result = prepare(inp)
        ref = {'credentials': {'materials_project': {'api_key': ''}}, 
               'strusets': [{'database': 'materials_project', 
                             'desc': [{'config': 'Na_dimer', 'scales': [2.1, 2.3, 2.5, 2.7, 3.0, 3.5, 4.0]}, 
                                      {'config': 'NaO_xy2', 'scales': [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], 'kspacing': [0.1, 0.08, 0.06]}, 
                                      {'config': './apns_cache/mp-1960.cif', 'scales': [1.0], 'kspacing': [0.08, 0.05], 'magmoms': [0.0, 0.0, 0.0]}, 
                                      {'config': './apns_cache/mp-2340.cif', 'scales': [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], 'kspacing': [0.1, 0.08, 0.06], 'magmoms': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, 
                            {'database': 'materials_project', 
                             'desc': [{'config': 'VO_x2y5', 'scales': [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], 'kspacing': [0.12]}, 
                                      {'config': '/path/to/structure', 'scales': [1.0], 'kspacing': [0.1]}, 
                                      {'config': './apns_cache/mp-1143.cif', 'scales': [1.0], 'kspacing': [0.1], 'magmoms': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}]}, 
                            {'database': 'materials_project', 
                             'desc': [{'config': './apns_cache/mp-1960.cif', 'scales': [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06], 'kspacing': [0.05, 0.08, 0.1], 'magmoms': [0.0, 0.0, 0.0]}]}, 
                            {'database': 'materials_project', 
                             'desc': [{'config': './apns_cache/mp-1960.cif', 'scales': [1.0], 'kspacing': [0.08], 'magmoms': [0.0, 0.0, 0.0]}]}]}
        self.assertEqual(result, ref)

if __name__ == "__main__":
    unittest.main()