'''
update the database indexing the numerical atomic orbitals locally

Details
-------
The ABACUS-Pseudopot-Nao-Square backend high-throughput workflows
manage all pseudopotentials and numerical atomic orbitals downloaded
locally in a tag-search manner. Each numerical atomic orbital will be
tagged with a list of strings, including the corresponding 
pseudopotential, the element symbol, the cutoff radius, orbital
configuration, etc. The tags are used to search for the numerical
atomic orbitals in the database.

Usage
-----
NOTE: you should read the `usage` of the `update.py` in the `download/orb` 
directory first.
The run of tagging on numerical atomic orbitals requires the rules for
tagging. The rules should be stored in a JSON file, besides the format
of the rules is the same as the pseudopotential rules, there is an
additional pseudopotential binding rule:
```json
{
    "description": "***",
    "rules": [
        {
            "re.folder": "SG15\\-Version1p0__AllOrbitals\\-Version2p1",
            "re.file": "^([A-Z][a-z]?)(_gga_\\d+(\\.\\d+)?au_\\d+(\\.\\d+)?Ry.*\\.orb)$",
            "tags": [["sg15", "1.0", "sr"], 
                     ["2.1", "bfgs"]]
        },
        {
            "re.folder": "20240618",
            "re.file": "^([A-Z][a-z]?)(_gga_\\d+(\\.\\d+)?au_\\d+(\\.\\d+)?Ry.*\\.orb)$",
            "tags": ["20240618"]
        }
    ]
}
```
, in which the first rule is for binding the pseudopotential to the
orbitals, and add two tags `2.1`, `bfgs` to the orbital file. The tags
in the first list are tags can be used to search a *unique* pseudopotential.
Likewise, once there is a file satisfies the regex in `re.folder` and
`re.file`, the tags in `tags` will be added to the file. The tags are
appended to the file, instead of replacing the original tags.
Then you can safely run the following command to update the database
(if the database file is not found, it will be created. But this cannot
be applied to the non-existing case of rules file because it is
considered as a fatal error):
```bash
python update.py
```

Version
-------
20241230
'''

# in-built modules
import os
import re
import json
import logging
import time

# APNS special modules
from apns.test.tag_search import TagSearcher

def _impl(fn_psp_db, fn_orb_db, fn_rules, orbital_dir):
    '''
    initialize the database of pseudopotentials and numerical orbitals
    
    Parameters
    ----------
    fn_psp_db: str
        the path to the database of pseudopotentials
    
    fn_orb_db: str
        the path to the database of numerical orbitals
    
    fn_rules: str
        the path to the rules for tagging the numerical orbitals
    
    orbital_dir: str
        the root directory of the numerical orbital files.
    
    Returns
    -------
    list: a list containing the paths of numerical orbital files that not recognized by rules
    '''
    unlabelled = []

    logging.info(f'Initialize the orbital tag-search system database...')

    # load the databases and rules
    with open(fn_orb_db) as f:
        orb_db = json.load(f)
    logging.info(f'Orbital database loaded/opened from {fn_orb_db}.')

    with open(fn_rules) as f:
        rules = json.load(f)['rules']
    logging.info(f'Rules loaded from {fn_rules}.')

    # load the pseudopotential database
    psplib = TagSearcher(fn_psp_db) # read database by opening with TagSearcher
    logging.info(f'Pseudopotential database loaded from {fn_psp_db} to APNS-TagSearcher.')

    for root, _, files in os.walk(orbital_dir):
        logging.info(f'Walking through {root}...')
        for file in files:
            m = re.match(r'^([A-Z][a-z]?)_gga_(\d+(\.\d+)?)au_(\d+(\.\d+)?)Ry_(.*)\.orb$', 
                         os.path.basename(file))
            if m:
                elem, rcut, _, ecut, _, component = m.groups()
                logging.info(f'Found orbital file {file} for element {elem} with rcut {rcut} au and ecut {ecut} Ry.')
                fn = os.path.abspath(os.path.join(root, file))

                logging.info(f'Transversing rules to add tags to {fn}...')
                for rule in rules:
                    redir, refn, add = rule['re.folder'], rule['re.file'], rule['tags']
                    mm = re.search(redir, root) and re.search(refn, file)
                    if mm:
                        logging.info(f'Valid rule found for {fn}.')
                    else:
                        continue # no need to waste time on invalid rules

                    # there are two possiblities: add pseudopotential tags or not. For the former, the tags will be
                    # in two groups, the first is for pseudopotential, the second is for orbital. For the latter, the
                    # tags are for the orbital file only.
                    temp = None
                    if all([isinstance(t, str) for t in add]): # the direct case: only tags for the orbital file
                        temp = add.copy()
                    else: # instead, assume the `add` is a list of two lists of tags, the first is the corresponding psp
                        pspt, orbt = add
                        fupf = psplib.search(False, False, elem, *pspt)
                        if len(fupf) != 1:
                            errmsg = f'Each orbital file should correspond to ONLY ONE pseudopotential file. Found {len(fupf)} for {elem} with tags {pspt}.'
                            logging.error(errmsg)
                            raise ValueError(errmsg)
                        temp = orbt + list(fupf)
                    temp += [rcut + 'au', ecut + 'Ry', component]

                    # the database is a dictionary, the key is the absolute path of the orbital file, and the value is
                    # a list of tags
                    if set(orb_db.get(fn, [])) != set(temp):
                        logging.info(f'Refresh. Appending tags {temp} to {fn} and updating the database.')
                        orb_db[fn] = list(set(orb_db.get(fn, []) + temp))
                    else:
                        logging.info(f'No need to refresh, have a good day: tags {temp} already in {fn}.')
                    
                unlabelled.append(fn) if fn not in orb_db else None
                    
    with open(fn_orb_db, "w") as f:
        json.dump(orb_db, f, indent=4)
    logging.info(f'Orbital database updated and saved to {fn_orb_db}.')

    logging.info(f'Orbital tag-search system database initialized, {len(unlabelled)} files not recognized by rules.')
    return unlabelled

def main(prefix, fn_psp_db, fn_orb_db, fn_rules, orbital_dir = None):
    '''
    update the database indexing the numerical atomic orbitals locally
    
    Parameters
    ----------
    fn_psp_db: str
        the path to the database of pseudopotentials
    
    fn_orb_db: str
        the path to the database of numerical orbitals
    
    fn_rules: str
        the path to the rules for tagging the numerical orbitals
    
    orbital_dir: str
        the root directory of the numerical orbital files. Numerical orbital files
        are allowed to be placed in any level of subdirectories. If this is
        None, the directory of the orbital database will be used.
    
    Returns
    -------
    list: a list containing the paths of numerical orbital files that not recognized by rules
    '''
    # init log
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    fnlog = f'{prefix}@{timestamp}.log'
    print(f'Logging to {fnlog}...')
    logging.basicConfig(filename=fnlog, level=logging.INFO)
    logging.info(f'Update started at {timestamp}...')

    # sanity
    if not os.path.exists(fn_psp_db):
        raise FileNotFoundError(f'Pseudopotential database {fn_psp_db} not found.')
    if not os.path.exists(fn_orb_db):
        # create the database if not exists
        with open(fn_orb_db, 'w') as f:
            json.dump({}, f)
        logging.info(f'Orbital database created at {fn_orb_db}.')

    # kernel
    orbital_dir = os.path.dirname(fn_orb_db) if orbital_dir is None else orbital_dir
    files = _impl(fn_psp_db, fn_orb_db, fn_rules, orbital_dir)

    # close log
    logging.info(f'Update completed at {time.strftime("%Y%m%d-%H%M%S")}.')
    logging.shutdown()
    return files

if __name__ == "__main__":
    
    fdb_psp = '/root/abacus-develop/pseudopotentials/database.json'
    fdb_orb = '/root/abacus-develop/numerical_orbitals/database.json'
    frules  = '/root/abacus-develop/numerical_orbitals/rules.json'

    # run
    unlabelled = main('local_naolib_update', fdb_psp, fdb_orb, frules)
    n = len(unlabelled) # number of unlabelled files
    if n > 0:
        print('The unlabelled files are:')
        for fn in unlabelled:
            print(f'  {fn}')
        print('Make sure it is what you expected. Otherwise please check the rules in the rules.json file.')

