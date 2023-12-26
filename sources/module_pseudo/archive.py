import os
import re
import json

"""This is for archiving pseudopotentials element-wise.
_scan_ returns a dictionary of pseudopotentials, element-wise, recursively including subdirectories.
"""

"""Extension for this module
Once new kind of pseudopotential is added into folder where all pseudopotentials are stored, 
update _folder_match_ and archive.
"""

"""Variable"""
only_scan = True

"""Functions"""
def _folder_match_(folder: str):

    if folder.startswith("nc-"):
        return "dojo"
    elif folder.startswith("pbe_s_sr"):
        return "dojo"
    elif folder.startswith("NCPP-PD04"):
        return "pd04"
    elif folder.startswith("NCPP-PD03"):
        return "pd03"
    elif folder.startswith("sg15"):
        return "sg15"
    else:
        print("Current folder name is: ", folder)
        raise ValueError("Folder name not recognized, add a logic branch in this function"
                      +"\nand create a function to create description.json to describe the folder.")

def _scan_(pseudo_dir: str) -> dict:
    """Returns a dictionary of pseudopotentials, element-wise.

    Args:
        pseudo_dir (str): path to the directory containing pseudopotentials

    Returns:
        dict: a dictionary of pseudopotentials, element-wise. Values are lists of pseudopotentials including path.
        list: a list contains folder in which pseudopotentials are stored.
    """
    result = {}
    folders = []
    for dir in list(os.walk(pseudo_dir)):
        # print(dir)
        # print(dir[0]) # present directory, str
        # print(dir[1]) # all subfolders names in list
        # print(dir[2]) # all files in list
        for pseudopotential in dir[2]:
            match = re.match(r"^([A-Za-z]{1,2})([-._]?.*)(.upf)$", pseudopotential, re.IGNORECASE)
            if match:
                # find dir contains pseudopotential
                element = match.group(1)[0].upper() + match.group(1)[1:].lower()
                if element not in result:
                    result[element] = []
                _dir = os.path.abspath(dir[0])
                result[element].append(
                    _dir+('\\' if _dir.count('\\') > 0 else _dir+'/')+pseudopotential
                )
                if dir[0] not in folders:
                    folder = dir[0].split("/")[-1] if dir[0].count("/") > 0 else dir[0].split("\\")[-1]
                    folders.append(folder)
    return result, folders

def _PD04_(path: str):
    """archive PD04 pseudopotentials. Becuase PD04 pseudopotential has more than 3 subsets.
    With more details, will create different subfolders and move upfs into them.
    Args:
        path (str): path of PD04 pseudopotentials
    """
    pseudopotential_kinds = {"default": []}
    path_backup = os.path.abspath(os.getcwd())
    os.chdir(path)
    files = os.listdir()
    for file in files:
        _match = re.match(r"^([A-Za-z]{1,2})(.*\.PD04)(.*)$", file)
        if _match:
            kind = _match.group(2).replace(".PD04", "")
            if kind == "":
                pseudopotential_kinds["default"].append(file)
            else:
                kind = kind[1:] if kind.startswith("-") else kind
                if kind not in pseudopotential_kinds.keys():
                    pseudopotential_kinds[kind] = []
                pseudopotential_kinds[kind].append(file)
    for kind in pseudopotential_kinds.keys():
        os.mkdir(kind)
        for file in pseudopotential_kinds[kind]:
            os.rename(file, kind+'/'+file)
        description = {"kind": "pd", "version": "04"}
        if kind != "default":
            description["appendix"] = kind
        else:
            description["appendix"] = ""
        with open(kind+'/'+"description.json", "w") as json_f:
            json.dump(description, json_f, indent=4)

    os.chdir(path_backup)

    return pseudopotential_kinds

def _PD03_(path: str):
    """add description.json in PD03 pseudopotential folder 
    """
    path_backup = os.path.abspath(os.getcwd())
    os.chdir(path)
    description = {"kind": "pd", "version": "03", "appendix": ""}
    with open("description.json", "w") as json_f:
        json.dump(description, json_f, indent=4)
    os.chdir(path_backup)
    return description

def _SG15_(path: str):
    """archive SG15 pseudopotentials. Becuase SG15 pseudopotential has more than 3 subsets.
    With more details, will create different subfolders and move upfs into them.
    Args:
        path (str): path of SG15 pseudopotentials
    """
    pseudopotential_kinds = {"1.0": []} # at least one version
    path_backup = os.path.abspath(os.getcwd())
    os.chdir(path)
    files = os.listdir()
    for file in files:
        _match = re.match(r"^(.*)(-)([0-9]\.[0-9])(.upf)$", file, re.IGNORECASE)
        if _match:
            version = _match.group(3)
            version += "" if "FR" not in _match.group(1) else "_fr"
            if version not in pseudopotential_kinds.keys():
                pseudopotential_kinds[version] = []
            pseudopotential_kinds[version].append(file)
    for kind in pseudopotential_kinds.keys():
        os.mkdir(kind)
        for file in pseudopotential_kinds[kind]:
            os.rename(file, kind+'/'+file)
        description = {"kind": "sg15"}
        _match = re.match(r"^([0-9]\.[0-9])(_)?(FR)?", kind)
        description["version"] = _match.group(1)
        description["appendix"] = "fr" if _match.group(3) else ""
        with open(kind+'/'+"description.json", "w") as json_f:
            json.dump(description, json_f, indent=4)

    os.chdir(path_backup)

    return pseudopotential_kinds

def _DOJO_(path: str):
    """add description.json in DOJO pseudopotential folder  
    special case: DOJO v0.3: pbe_s_sr  
    general pattern: r"^(nc-)([s|f]r)(-)([0-9]{1,2})([-_])([.*])?(_pbe_standard_upf)(.*)$"  
    """
    description = {"kind": "dojo"}
    path_backup = os.path.abspath(os.getcwd())
    folder = path.split("\\")[-1] if path.count("\\") > 0 else path.split("/")[-1]
    os.chdir(path)
    if folder == "pbe_s_sr":
        description["version"] = "03"
        description["appendix"] = ""
    else:
        _match = re.match(r"^(nc-)([sf]r)(-)([0-9]{1,2})([-_]*)(.*)(_pbe_standard)(.*)$", folder)
        if not _match:
            print("Cannot recognize arbitrary named pseudopotential folder name: ", folder)
            raise ValueError("Not standard folder name of DOJO pseudopotential.")
        if not _match.group(1).startswith("nc"):
            raise ValueError("Non norm-conserving pseudopotential is not supported by ABACUS yet.")
        description["version"] = _match.group(4)
        if _match.group(2) == "sr":
            description["appendix"] = _match.group(6)
        else:
            description["appendix"] = _match.group(2) + "_" + _match.group(6) if _match.group(6) != "" else _match.group(2)
    with open("description.json", "w") as json_f:
        json.dump(description, json_f, indent=4)

    os.chdir(path_backup)
    return description

def archive(pseudo_dir: str = "./download/pseudopotentials/", only_scan: bool = only_scan):
    """archive pseudopotential files, will also create description.json in each folder created

    Args:
        pseudo_dir (str, optional): folder where all kinds of pseudopotentials are stored folder-by-folder. Defaults to "./download/pseudopotentials/".
        only_scan (bool, optional): if only scan without further moving upf files. Defaults to True.

    Returns:
        dict: available pseudopotentials for each element is stored in list and act as value, whose corresponding
        key is the element symbol.
    """
    if not only_scan:
        _, folders = _scan_(pseudo_dir=pseudo_dir)
        pseudo_dir = pseudo_dir[:-1] if (pseudo_dir.endswith("/") or pseudo_dir.endswith("\\")) else pseudo_dir
        
        for folder in folders:
            folder_path = pseudo_dir + "/" + folder
            if _folder_match_(folder) == "pd04":
                _PD04_(folder_path)
            elif _folder_match_(folder) == "pd03":
                _PD03_(folder_path)
            elif _folder_match_(folder) == "dojo":
                _DOJO_(folder_path)
            elif _folder_match_(folder) == "sg15":
                _SG15_(folder_path)
    
    results, _ = _scan_(pseudo_dir=pseudo_dir)
    return results

def description(upf_path: str):
    """Get the description.json contents that generated by archive function.
    via the full path of upf file
    Args:
        upf_path (str): upf file path

    Raises:
        FileNotFoundError: description.json not generated

    Returns:
        dict: contents, dict saved in description.json
    """
    path_backup = os.path.abspath(os.getcwd())
    path = ".\\"
    if upf_path.count("/") > 0:
        path = "\\".join(upf_path.split("/")[:-1])
    elif upf_path.count("\\") > 0:
        path = "\\".join(upf_path.split("\\")[:-1])
    else:
        path = path_backup
    os.chdir(path)
    description = {}
    if os.path.isfile("description.json"):
        with open("description.json", "r") as json_f:
            description = json.load(json_f)
    else:
        print("Current directory: ", path)
        raise FileNotFoundError("description.json is not found.")

    os.chdir(path_backup)
    return description

if __name__ == "__main__":

    print(archive())