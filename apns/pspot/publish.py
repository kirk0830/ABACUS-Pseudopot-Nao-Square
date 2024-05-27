"""
Collect eos test result and pack the pseudopotential file.
The time to publish recommended pseudopotentials is coming.
"""
import apns.module_workflow.identifier as amwi
import os
def initialize():
    """initialize by asserting the existence of all files"""
    APNS_CACHE_DIR = amwi.TEMPORARY_FOLDER
    EOS_DB = "apns_eos_result.json"
    ECUTWFC_DB = "apns_ecutwfc_db.json"
    PSEUDO_DIR = "./download/pseudopotentials"
    PSEUDO_DB = "pseudo_db.json"
    feos = os.path.join(APNS_CACHE_DIR, EOS_DB)
    assert os.path.exists(feos), f"{feos} does not exist"
    fecut = os.path.join(APNS_CACHE_DIR, ECUTWFC_DB)
    assert os.path.exists(fecut), f"{fecut} does not exist"
    fpspot = os.path.join(PSEUDO_DIR, PSEUDO_DB)
    assert os.path.exists(fpspot), f"{fpspot} does not exist"
    return feos, fecut, fpspot

import json
def calculate_delta_peratom(feos, excluded: list = None):
    """index the feos by element, bravais lattice, pspotid,
    then will get natoms and delta from each entry"""
    with open(feos, "r") as f:
        data = json.load(f)
    
    result = {}
    for element in data.keys():
        if element in excluded:
            continue
        for brav in data[element].keys():
            pspotids = list(data[element][brav].keys())
            pspotids = [pspotid for pspotid in pspotids if pspotid != "AEref"]
            for pspotid in pspotids:
                delta = data[element][brav][pspotid]["delta"]
                natoms = data[element][brav][pspotid]["natoms"]
                result.setdefault(element, {}).setdefault(pspotid, []).append(delta/natoms)
        # then average over all pspotids/bravis lattices/structures/phases for each element
        for pspotid in result[element].keys():
            vals = result[element][pspotid]
            print(f"Element: {element}, Pseudopotential: {pspotid}, delta_peratom: {vals}")
            avg = sum(vals)/len(vals)
            result[element][pspotid] = avg

    return result

def sort_delta_peratom(data):
    """sort the delta per atom by average"""
    result = {}
    for element in data.keys():
        for pspotid in data[element].keys():
            result.setdefault(element, {})[pspotid] = data[element][pspotid]
        result[element] = dict(sorted(result[element].items(), key=lambda x: x[1]))
    return result

def load_ecutwfc(fecut):
    """load the ecutwfc database"""
    with open(fecut, "r") as f:
        data = json.load(f)
    return data

def load_pseudo(fpspot):
    """load the pseudopotential database"""
    with open(fpspot, "r") as f:
        data = json.load(f)
    return data

def pspot_eq(val1: str, val2: str):
    """due to test name encoding, all underlines and dots are removed in the test name,
    however in pseudopotenital database, all underlines and dots are unchanged. It is
    needed to identify if the val1 and val2 indicate the same pseudopotential"""

    abbr = val1 if len(val1) < len(val2) else val2
    full = val1 if len(val1) > len(val2) else val2

    return full.replace("_", "").replace(".", "") == abbr

def attributes(val):
    """there are three attributes, kind, version and appendix"""
    words = val.split("_")
    words = [words[0], words[1], "_".join(words[2:])] if len(words) > 2 else words
    # map words to results with length 3, absent will be filled with empty string
    results = words + [""] * (3 - len(words))
    return results

import apns.pspot.parse as ampp
def pack(eos_db, ecutwfc_db, pseudo_db, n_lowest: int = 5):
    """for elements in eos_db, get n_lowest delta per atom-corresponding pseudopotentials,
    then collect ecutwfc and its path"""
    result = {}
    for element in eos_db.keys():
        ecutwfcs = ecutwfc_db[element]
        deltas_peratom = eos_db[element]
        pseudos = pseudo_db[element]
        for pspotid in list(deltas_peratom.keys())[:min(n_lowest, len(deltas_peratom))]:
            for _pspotid in pseudos.keys():
                if pspot_eq(pspotid, _pspotid):
                    kind, version, appendix = attributes(_pspotid)
                    ecutwfc = ecutwfcs[pspotid]
                    fpseudo = pseudos[_pspotid]
                    zval = ampp.z_valence(fpseudo)
                    delta_peratom = deltas_peratom[pspotid]
                    print(f"Element: {element}, Pseudopotential: {_pspotid}, ecutwfc: {ecutwfc}, delta_peratom: {delta_peratom}")
                    result.setdefault(element, []).append(
                        dict(zip(["kind", "version", "appendix", "ecutwfc", "fpseudo", "delta_peratom", "zval"],
                                [kind, version, appendix, ecutwfc, fpseudo, delta_peratom, zval])))
    
    return result

def file_pack(pack):
    """according to dict returned by function pack, copy the pseudopotential files to a new directory"""
    database = {}
    for element in pack.keys():
        os.makedirs(element, exist_ok=True)
        # create a dict storing the pseudopotential files, ecutwfc and delta per atom
        
        for item in pack[element]: # loop over all pseudopotentials
            fpseudo = item["fpseudo"]
            if not os.path.exists(fpseudo):
                print(f"{fpseudo} does not exist")
                continue
            kind = item["kind"]
            version = item["version"]
            appendix = item["appendix"]
            fname = f"{element}_{kind}_{version}_{appendix}.UPF"
            fnew = os.path.join(element, fname)
            os.system(f"cp {fpseudo} {fnew}")
            print(f"{fpseudo} copied to {fnew}")
            
            database.setdefault(element, []).append(
                dict(zip(["kind", "version", "appendix", "ecutwfc", "file", "delta_peratom", "zval"],
                         [kind, version, appendix, item["ecutwfc"], fnew, item["delta_peratom"], item["zval"]]))
            )
            
    with open("database.json", "w") as f:
        json.dump(database, f, indent=4)
    
    return list(database.keys())

import apns.new.compress as compress
import time
def main(excluded: list = None, n_lowest: int = 5):
    feos, fecut, fpspot = initialize()
    eos_db = calculate_delta_peratom(feos, excluded=excluded)
    eos_db = sort_delta_peratom(eos_db)
    ecutwfc_db = load_ecutwfc(fecut)
    pseudo_db = load_pseudo(fpspot)
    pack_db = pack(eos_db, ecutwfc_db, pseudo_db, n_lowest)
    folders_and_files = file_pack(pack_db)
    compress.pack(folders_and_files, f"apns_{time.strftime('%Y%m%d%H%M%S')}.zip")
    for f in folders_and_files:
        os.system(f"rm -rf {f}")

if __name__ == "__main__":
    main(excluded=["Nb", "W", "Al", "Cu", "Hf", "Mg", "Mo", "Ta", "Ti", "V", "Zr"],
         n_lowest=5)
