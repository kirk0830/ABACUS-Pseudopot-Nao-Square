'''to face with the demand of file searching and managing in the project, we
develop a class to search files based on tags. The class is called TagSearcher.

The class is designed to search pseudopotential and/or numerical atomic orbitals
based on specified tags. The database is a dictionary with keys as file full-path
and values as tags. The search is based on set intersection of tags. The database
is loaded from a json file. The search is performed by calling the class with tags
as arguments. The result is a set of file full-path.
'''

import os
import json
class TagSearcher:
    '''
    A class to search pseudopotential and/or numerical atomic orbitals
    based on specified tags. The database is a dictionary with keys as
    file full-path and values as tags. The search is based on set intersection
    of tags. The database is loaded from a json file. The search is performed
    by calling the class with tags as arguments. The result is a set of
    file full-path.
    '''
    database = {} # the database to store the tags
    
    @staticmethod
    def impl(tag: str, database: dict, ignorecase: bool = True):
        '''search all keys in database whose values contain the tag
        
        Parameters
        ----------
        tag : str
            the tag to search
        database : dict
            the database to search, in which keys have one or more values
            as tags
        ignorecase : bool
            whether to ignore case when searching
        
        Returns
        -------
        set
            the set of keys whose values contain the tag
        '''
        result = []
        
        # we convert to lowercase for ignorecase mode
        tag = tag.lower() if ignorecase else tag
        
        for key in database:
            character = database[key]
            if ignorecase:
                character = [c.lower() for c in character]
            if tag in character:
                result.append(key)
                    
        return set(result)
        
    def __init__(self, fdatabase: str):
        '''instantiate with the database that stores as a JSON file,
        in this way, the file I/O is minimized
        
        Parameters
        ----------
        fdatabase : str
            the file path of the database
        '''
        with open(fdatabase) as f:
            self.database = json.load(f)
        # this function is expected to raise:
        # FileNotFoundError if the file does not exist
        # JSONDecodeError if the file is not a valid JSON file

    def search(self, 
               ignorecase: bool = True, 
               return_with_tags: bool = False, 
               *tags):
        '''
        search with tags, also support the excluded tags by using '^' as prefix
        '''
        result = set(self.database.keys()) # each time we search, we start with all keys
        excluded = set() # we start with empty set for excluded tags
        # do iteratively search for each tag
        for tag in [t for t in tags if not t.startswith('^')]:
            # iteratively find the intersection between all tags...
            result = set.intersection(result, 
                                      TagSearcher.impl(tag, self.database, ignorecase))
        # those tags with '^' as prefix are excluded
        for tag in [t[1:] for t in tags if t.startswith('^')]:
            excluded = set.union(excluded,
                                 TagSearcher.impl(tag, self.database, ignorecase))
        result = set.difference(result, excluded)
    
        return result if not return_with_tags else {key: self.database[key] for key in result}
    
import unittest
class TestTagSearcher(unittest.TestCase):

    def test_kernel_impl(self):
        # normal mode
        database = {
            'file1': ['tag1', 'tag2'],
            'file2': ['tag2', 'tag3'],
            'file3': ['tag3', 'tag4']
        }
        self.assertEqual(TagSearcher.impl('tag2', database), 
                         {'file1', 'file2'})
        self.assertEqual(TagSearcher.impl('tag2', database, ignorecase=False), 
                         {'file1', 'file2'})
        self.assertEqual(TagSearcher.impl('TAG2', database, ignorecase=False), 
                         set())
        self.assertEqual(TagSearcher.impl('tag5', database), 
                         set())

    def test_tagsearcher(self):
        database = {
            'file1': ['tag1', 'tag2'],
            'file2': ['tag2', 'tag3'],
            'file3': ['tag3', 'tag4']
        }
        fdatabase = 'test_database.json'
        with open(fdatabase, 'w') as f:
            json.dump(database, f)
        mysearcher = TagSearcher(fdatabase)
        os.remove(fdatabase) # we can remove this file at once
        
        # then we test
        # normal mode
        self.assertEqual(mysearcher.search(True, False, 'tag2'), 
                         {'file1', 'file2'})
        self.assertEqual(mysearcher.search(False, False, 'tag2'), 
                         {'file1', 'file2'})
        self.assertEqual(mysearcher.search(False, False, 'TAG2'), 
                         set())
        self.assertEqual(mysearcher.search(False, False, 'tag5'), 
                         set())
        self.assertEqual(mysearcher.search(True, True, 'tag2'), 
                         {'file1': ['tag1', 'tag2'], 'file2': ['tag2', 'tag3']})
        # exclusive
        self.assertEqual(mysearcher.search(True, False, '^tag2'), 
                         {'file3'})

if __name__ == '__main__':
    unittest.main()
