"""NOTE: this script is for updating the local (user-specific) library of numerical orbitals.
Therefore the following path is not general. specify with the one fit your own file structure."""
def initialize(orbital_dir: str = "/root/abacus-develop/numerical_orbitals", 
               pseudo_dir: str = "/root/abacus-develop/pseudopotentials", 
               refresh: bool = False):
    """(re-)initialize the database.json file in one folder, for describing the numerical orbitals
    files within. Based on tag searching system implemented in APNS. For rules how to set tags
    on the numerical orbitals, see the rules.json file in the same folder.
    
    Args:
        orbital_dir (str): the folder containing the numerical orbitals
        pseudo_dir (str): the folder containing the pseudopotentials
        refresh (bool): whether to refresh the database file
    
    Returns:
        list[str]: a list of unclassified numerical orbital files
    
    Note:
        always take care on the returned value, because it may contain some unclassified files.
        If so, you need to manually tag them or update the rules file.
    """
    import os, json, re
    from apns.test.tag_search import TagSearcher
    
    # check if the database file exists, if not, create one
    fdb_orb = os.path.join(orbital_dir, "database.json")
    if not os.path.exists(fdb_orb):
        with open(fdb_orb, "w") as f:
            json.dump({}, f)
    
    # no matter what workflow is (refresh or not), always create a database file.
    # while if not refresh, return an empty list
    if not refresh:
        return []
    
    # else, if refresh, update the database file with the rules. But due to there
    # would be some unclassified files, return them for further manual tagging.
    orbs_unclassified = []
    with open(fdb_orb) as f:
        orb_db = json.load(f)
    
    # because orb is based on the pseudopotential, so first initialize the pseudopotential database
    fdb_pp = os.path.join(pseudo_dir, "database.json")
    pp_searcher = TagSearcher(fdb_pp) # with this, can search local available pseudopotential

    # orbpat, for matching the numerical orbital file name
    orbpat = r"([A-Z][a-z]?)_gga_(\d+(\.\d+)?)au_(\d+(\.\d+)?)Ry_(.*)\.orb"

    # read the rules file, according to user preset rules, tag the numerical orbital files
    ftag = os.path.join(orbital_dir, "rules.json")
    if os.path.exists(ftag):
        with open(ftag) as f:
            rules = json.load(f)
        
        # loop over all the numerical orbital files, and tag them according to the rules
        for root, _, files in os.walk(orbital_dir):
            for file in files:
                # match the file name with the pattern
                m = re.match(orbpat, file)
                if m:
                    element, rcut, _, ecut, _, conf = m.groups()
                    key = os.path.abspath(os.path.join(root, file))

                    # scan the rules, if the file matches the rule, add the tags to the database
                    for rule in rules["rules"]:
                        re_folder, re_file, tags = rule["re.folder"], rule["re.file"], rule["tags"]
                        if len(tags) == 2 and all([isinstance(taggrp, list) for taggrp in tags]):

                            # sometimes it is needed to add its corresponding upf file as a tag of the orbital file
                            # then there will be two groups of tags, the first is for pseudopotential, the second is for orbital
                            pptags, orbtags = tags[0], tags[1]
                            pptags = [element] + pptags
                            fupf = pp_searcher(False, False, *pptags)
                            assert len(fupf) > 0, f"No pseudopotential found for {element} with tags {pptags}"
                            assert len(fupf) < 2, f"Multiple pseudopotentials found for {element} with tags {pptags}:\n{fupf}"
                            orbtags = orbtags + list(fupf)
                        else:
                            # otherwise the tags are for the orbital file only
                            orbtags = tags

                        tags_add = [rcut + "au", ecut + "Ry", conf] + orbtags
                        if re.search(re_folder, root) and re.match(re_file, file):
                            if key not in orb_db:
                                orb_db.setdefault(key, []).extend(tags_add)
                            if not set(tags_add) <= set(orb_db[key]):
                                print(f"Appending tags {tags_add} to {key}...")
                                orb_db[key] = list(set(orb_db[key] + tags_add))
                        # else:
                        #     print(f"WARNING: Regular expression does not match:\
                        #           \n{re_folder} <= {root},\n{re_file} <= {file}")
                    orbs_unclassified.append(key) if key not in orb_db else None
        with open(fdb_orb, "w") as f:
            json.dump(orb_db, f, indent=4)
    else:
        raise FileNotFoundError("Rules file not found")
    print(f"""there are {len(orbs_unclassified)} unclassified numerical orbital files, see returned value for details.
you can update the rules file to include these files, or manually add tags to them.""")
    return orbs_unclassified

if __name__ == "__main__":
    orbs = initialize(refresh=True)
    print(orbs)