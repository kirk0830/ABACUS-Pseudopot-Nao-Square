"""this script is for updating the pseudopotential database
After 2024 04 30, each pseudopotential file is described by
a series of tags. Then accessing pseudopotential file would
be the process of roll and roll tag-filtering task"""
FDATABASE = "./download/pseudopotentials/database.json"
"""this file to save all descriptions of pseudopotential files
is hard-coded as database.json"""
TAGRULES = "./download/pseudopotentials/tag_rules.json"
"""this file saves rules to add tags onto various pseudopotential
files. The rules are in the form of a dictionary, where the value
of key `rules` is a list of dictionaries, each dictionary contains
`re.folder`, `re.file` and `tags` keys. The `re.folder` and `re.file`
are regular expressions to match the folder and file name, respectively.
If two regular expressions are matched, the tags will be added to the
pseudopotential file. The `tags` key is a list of tags to be added."""
import os
import json
def initialize():
    """initialize will create a database file if not exists"""
    if not os.path.exists(FDATABASE):
        with open(FDATABASE, "w") as f:
            json.dump({}, f)

import re
def mark_mutual_tags(folder: str, tags: list, regex: str = r".*", unique_tag: str = "from_name"):
    """mark mutual tags for all files in the folder,
    and add unique tag to the tags list. The unique tag
    can be set to `from_name` or `from_file`, for the former,
    will grep the first one or two elements from the file name (split
    the file name by `_`, `-`, or `.`), for the latter, will grep element
    information from file."""
    assert os.path.exists(folder), "Folder not found"
    assert os.path.isdir(folder), "Not a folder"
    assert os.path.exists(FDATABASE), "Database file not found"
    assert unique_tag in ["from_name", "from_file"], "Invalid unique tag"

    with open(FDATABASE) as f:
        database = json.load(f)
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.lower().endswith(".upf") and re.match(regex, file):
                # use the full path of file as the key
                key = os.path.join(root, file)
                original_tags = database.get(key, [])
                database[key] = list(set(original_tags + tags))
                if unique_tag == "from_name":
                    element = file.split("_")[0].split("-")[0].split(".")[0]
                    database[key] = list(set(database[key] + [element.lower().capitalize()]))
                elif unique_tag == "from_file":
                    raise NotImplementedError("Not implemented yet")

    with open(FDATABASE, "w") as f:
        json.dump(database, f, indent=4)

def overwrite_tags(folder: str, tag_to_overwrite: str, tags: list, regex: str = r".*"):
    """overwrite tags for all files in the folder"""
    assert os.path.exists(folder), "Folder not found"
    assert os.path.isdir(folder), "Not a folder"
    assert os.path.exists(FDATABASE), "Database file not found"

    with open(FDATABASE) as f:
        database = json.load(f)
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.lower().endswith(".upf") and re.match(regex, file):
                key = os.path.join(root, file)
                original_tags = database.get(key, [])
                database[key] = list(set(original_tags + tags))
                database[key] = [tag for tag in database[key] if tag != tag_to_overwrite]

    with open(FDATABASE, "w") as f:
        json.dump(database, f, indent=4)

def remove_tags(folder: str, tags: list, regex: str = r".*"):
    """remove tags for all files in the folder"""
    assert os.path.exists(folder), "Folder not found"
    assert os.path.isdir(folder), "Not a folder"
    assert os.path.exists(FDATABASE), "Database file not found"

    with open(FDATABASE) as f:
        database = json.load(f)
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.lower().endswith(".upf") and re.match(regex, file):
                key = os.path.join(root, file)
                original_tags = database.get(key, [])
                database[key] = [tag for tag in original_tags if tag not in tags]

    with open(FDATABASE, "w") as f:
        json.dump(database, f, indent=4)

if __name__ == "__main__":
    # initialize()
    # mark_mutual_tags("./download/pseudopotentials/pbe_s_sr", 
    #                  ["ONCVPSP", "ONCV"], 
    #                  regex=r".*",
    #                  unique_tag="from_name")
    # #overwrite_tags("./download/pseudopotentials/pslibrary-pbe.0.3.1", "NC", ["rrkjus", "US", "ultrasoft", "rrkj"], regex=r".*rrkjus.*")
    # print("Database updated")
    import apns.module_pseudo.tag_search as ts
    searcher = ts.TagSearcher(FDATABASE)
    print(searcher(True, False, "Mn", "psl", "US"))