ORBITAL_DIR = "download/orb"
TAGRULES = ORBITAL_DIR + "/rules.json"
FORBDATABASE = ORBITAL_DIR + "/database.json"
PSEUDO_DIR = "download/upf"
FPPDATABASE = PSEUDO_DIR + "/database.json"

def initialize(refresh: bool = False):
    import os, json, re
    from apns.test.tag_search import TagSearcher
    if not os.path.exists(FORBDATABASE):
        with open(FORBDATABASE, "w") as f:
            json.dump({}, f)
    if not refresh:
        return []
    orbs_unclassified = []
    with open(FORBDATABASE) as f:
        orb_db = json.load(f)
    pp_searcher = TagSearcher(FPPDATABASE)
    orbpat = r"([A-Z][a-z]?)_gga_(\d+(\.\d+)?)au_(\d+(\.\d+)?)Ry_(.*)\.orb"
    if os.path.exists(TAGRULES):
        with open(TAGRULES) as f:
            rules = json.load(f)
        for root, _, files in os.walk(ORBITAL_DIR):
            for file in files:
                if file.lower().endswith(".orb"):
                    element, rcut, _, ecut, _, conf = re.match(orbpat, file).groups()
                    key = os.path.abspath(os.path.join(root, file))
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
                            if not key in orb_db:
                                orb_db.setdefault(key, []).extend(tags_add)
                            if not set(tags_add) <= set(orb_db[key]):
                                print(f"Appending tags {tags_add} to {key}...")
                                orb_db[key] = list(set(orb_db[key] + tags_add))
                        # else:
                        #     print(f"WARNING: Regular expression does not match:\
                        #           \n{re_folder} <= {root},\n{re_file} <= {file}")
                    orbs_unclassified.append(key) if key not in orb_db else None
        with open(FORBDATABASE, "w") as f:
            json.dump(orb_db, f, indent=4)
    else:
        raise FileNotFoundError("Rules file not found")
    print(f"""there are {len(orbs_unclassified)} unclassified numerical orbital files, see returned value for details.
you can update the rules file to include these files, or manually add tags to them.""")
    return orbs_unclassified

if __name__ == "__main__":
    orbs = initialize(True)
    print(orbs)