"""because the search task of pseudopotential and/or 
numerical atomic orbitals is becoming more and more
complicated, it should be based on tag search instead
of simply filtering on names anymore.
Note that this searcher is both needed by pseudopotential
and numerical atomic orbitals management, it should be
placed elsewhere other than module_pseudo"""

def tag_search(tag: str, database: dict, ignore_case: bool = True):
    """with one single tag, return set of pseudopotential"""
    result = []
    for key in database:
        if ignore_case:
            if tag.lower() in [x.lower() for x in database[key]]:
                result.append(key)
        else:
            if tag in database[key]:
                result.append(key)
    return set(result)

import os
import json
class TagSearcher:
    """A class to search pseudopotential and/or numerical atomic orbitals
    based on specified tags. The database is a dictionary with keys as
    file full-path and values as tags. The search is based on set intersection
    of tags. The database is loaded from a json file. The search is performed
    by calling the class with tags as arguments. The result is a set of
    file full-path."""
    result = set()
    database = dict()
    def __init__(self, fdatabase: str):
        assert os.path.exists(fdatabase), "Database file not found"
        with open(fdatabase) as f:
            self.database = json.load(f)

    def __call__(self, ignore_case: bool = True, return_with_tags: bool = False, *tags):
        self.result = set(self.database.keys())
        for tag in tags:
            self.result = set.intersection(self.result, tag_search(tag, self.database, ignore_case))
        if return_with_tags:
            return {key: self.database[key] for key in self.result}
        return self.result
    
import unittest
class TestTagSearcher(unittest.TestCase):

    def test_tag_search(self):
        database = {
            "file1": ["tag1", "tag2"],
            "file2": ["tag2", "tag3"],
            "file3": ["tag3", "tag4"]
        }
        self.assertEqual(tag_search("tag2", database), {"file1", "file2"})
        self.assertEqual(tag_search("tag2", database, ignore_case=False), {"file1", "file2"})
        self.assertEqual(tag_search("TAG2", database, ignore_case=False), set())
        self.assertEqual(tag_search("tag5", database), set())

    def test_tag_searcher(self):
        database = {
            "file1": ["tag1", "tag2"],
            "file2": ["tag2", "tag3"],
            "file3": ["tag3", "tag4"]
        }
        fdatabase = "./test_database.json"
        with open(fdatabase, "w") as f:
            json.dump(database, f)
        searcher = TagSearcher(fdatabase)
        self.assertEqual(searcher(True, False, "tag2"), {"file1", "file2"})
        self.assertEqual(searcher(False, False, "tag2"), {"file1", "file2"})
        self.assertEqual(searcher(False, False, "TAG2"), set())
        self.assertEqual(searcher(False, False, "tag5"), set())
        self.assertEqual(searcher(True, True, "tag2"), {"file1": ["tag1", "tag2"], "file2": ["tag2", "tag3"]})
        os.remove(fdatabase)

if __name__ == "__main__":
    unittest.main()
