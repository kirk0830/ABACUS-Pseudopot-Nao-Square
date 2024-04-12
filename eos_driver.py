elements = ["Li", "Be", "B", "C", "Na", "Mg", "Al", "Si", "S"]
pseudo_db_path = "./download/pseudopotentials/pseudo_db.json"
ecutwfc_db_path = "./apns_cache/ecutwfc_conv.json"
template_json = "./ignorethis_input.json"

import json

with open(pseudo_db_path, "r") as f:
    pseudo_db = json.load(f)
with open(ecutwfc_db_path, "r") as f:
    ecutwfc_db = json.load(f)

import os
import time
for element in elements:
    fjson = "temp.json"
    with open(template_json, "r") as f:
        data = json.load(f)
    data["systems"] = [element]
    for pseudo in pseudo_db[element]:
        result = pseudo.split("_")
        if len(result) == 3:
            kind, version, appendix = result
        elif len(result) == 2:
            kind, version = result
            appendix = ""
        elif len(result) == 1:
            kind = result[0]
            version = ""
            appendix = ""
        else:
            print("Error: ", pseudo)
            continue
        if kind == "hgh":
            print("Skip: ", pseudo)
            continue
        data["pseudopotentials"]["kinds"] = [kind]
        data["pseudopotentials"]["versions"] = [version]
        data["pseudopotentials"]["appendices"] = [appendix]

        for pseudo_ in ecutwfc_db[element]:
            if pseudo.replace("_", "").replace(".", "") == pseudo_:
                data["calculation"]["ecutwfc"] = ecutwfc_db[element][pseudo_]
                break
        with open(fjson, "w") as f:
            json.dump(data, f, indent=4)
        os.system("python3 main.py -i {}".format(fjson))
        time.sleep(10)