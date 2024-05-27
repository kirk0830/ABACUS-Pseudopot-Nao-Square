"""this script is for updating the pseudopotential database
After 2024 04 30, each pseudopotential file is described by
a series of tags. Then accessing pseudopotential file would
be the process of roll and roll tag-filtering task"""
PSEUDO_DIR = "./download/pseudopotentials"

TAGRULES = PSEUDO_DIR + "/rules.json"
"""this file saves rules to add tags onto various pseudopotential
files. The rules are in the form of a dictionary, where the value
of key `rules` is a list of dictionaries, each dictionary contains
`re.folder`, `re.file` and `tags` keys. The `re.folder` and `re.file`
are regular expressions to match the folder and file name, respectively.
If two regular expressions are matched, the tags will be added to the
pseudopotential file. The `tags` key is a list of tags to be added."""
FDATABASE = PSEUDO_DIR + "/database.json"
"""this file to save all descriptions of pseudopotential files
is hard-coded as database.json"""

import os
import json
import apns.pspot.parse as ampp
import apns.pspot.parse_special.GBRV_Vanderbilt as amppsg
def initialize(refresh: bool = False) -> list[str]:
    """initialize will create a database file if not exists"""
    if not os.path.exists(FDATABASE):
        with open(FDATABASE, "w") as f:
            json.dump({}, f)
    # refresh the database with TAGRULES
    if not refresh:
        return []
    upfs_unclassified = []
    with open(FDATABASE) as f:
        database = json.load(f)
    if os.path.exists(TAGRULES):
        with open(TAGRULES) as f:
            rules = json.load(f)
        for root, dirs, files in os.walk(PSEUDO_DIR): # big triangle code is a bad practice...
            for file in files:
                if file.lower().endswith(".upf"): # if there is any pseudopotential file, then add tags
                    key = os.path.abspath(os.path.join(root, file)) # key is directly the full path of file
                    for rule in rules["rules"]: # add tags according to rules, iterate over all rules
                        re_folder, re_file, tags = rule["re.folder"], rule["re.file"], rule["tags"]
                        if re.search(re_folder, root) and re.match(re_file, file):
                            if not key in database: # it is the first time to add tags, so add element tag
                                # the parse of pseudopotential will be slow, so need to reduce this operation
                                # as much as possible
                                pp = ampp.as_dict(key)
                                ppgencode = ampp.determine_code(pp)
                                if ppgencode == "GBRV":
                                    parsed = amppsg.PP_HEADER(pp["PP_HEADER"]["data"])
                                    element = parsed["attrib"]["element"].strip().lower().capitalize()
                                else:
                                    element = pp["PP_HEADER"]["attrib"]["element"].strip().lower().capitalize()
                                database.setdefault(key, []).append(element)
                            # check if tags are already in the database by comparing sets
                            if not set(tags) <= set(database[key]):
                                print(f"Appending tags {tags} to {key}...")
                                database[key] = list(set(database[key] + tags))
                    upfs_unclassified.append(key) if key not in database else None
        # save the database
        with open(FDATABASE, "w") as f:
            json.dump(database, f, indent=4)
        complete_tag(["fr", "full-relativistic"], ["sr", "scalar-relativistic"], FDATABASE)
    else:
        raise FileNotFoundError("Rules file not found")
    print(f"""there are {len(upfs_unclassified)} unclassified pseudopotential files, see returned value for details.
you can update the rules file to include these files, or manually add tags to them.""")
    return upfs_unclassified

def complete_tag(if_without: list, add: list, fdb: str):
    """for all pseudopotential files, if without tags in `if_without`, add tags in `add`
    useful for some pseudopotentials marked as "fr" and "full-relativistic", while those
     not fr, add "sr" and "scalar-relativistic" """
    assert os.path.exists(fdb), "database file not found"
    with open(fdb) as f:
        database = json.load(f)
    
    record = []
    for fpp in database.keys():
        if not set(if_without) <= set(database[fpp]):
            print(f"Appending tags {add} to {fpp}...")
            record.append(fpp)
            database[fpp] = list(set(database[fpp] + add))

    with open(fdb, "w") as f:
        json.dump(database, f, indent=4)

    return record 

import unittest
import re
class TestRegularExpression(unittest.TestCase):
    def test_refolder(self):
        """test re.folder key in rules"""
        PseudoDojov10 = "/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/pseudopotentials/pseudos_ac_she/116_Lv/Lv-6spd_r.upf"
        refolder = r"pseudos_ac_she/\d+_[A-Z][a-z]?"
        refile = r"[A-Z][a-z]?-\d+[a-z]*_r\.upf" 
        # will have tags ["Lv", 
        #                 "PseudoDojo", "DOJO", "abinit", 
        #                 "v1.0", "1.0", 
        #                 "PBE", "Perdew-Burke-Ernzerhof",
        #                 "NC", "norm-conserving",
        #                 "full-relativistic", "rel",
        #                 ...]
        self.assertTrue(re.search(refolder, PseudoDojov10))
        fupf = PseudoDojov10.split("/")[-1]
        print(fupf)
        self.assertTrue(re.match(refile, PseudoDojov10.split("/")[-1]))

if __name__ == "__main__":

    # unittest.main()
    # exit()
    # fail_upfs = initialize(True)
    # print(fail_upfs)
    # exit()
    import apns.new.tag_search as ts
    searcher = ts.TagSearcher(FDATABASE)
    print("\n".join(searcher(False, False, "Sr", "sr", "DOJO")))