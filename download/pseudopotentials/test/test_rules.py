import os
import re
import json

path = "./download/pseudopotentials/"
with open(path + "rules.json") as f:
    rules = json.load(f)

nrules = len(rules["rules"])
print("Number of rules: ", nrules)

for root, dirs, files in os.walk(path):
    for file in files:
        if file.lower().endswith(".upf"):
            folder = root.replace("\\", "/").split("/")[-1]
            print("Folder: ", folder)
            for i in range(nrules):
                match_file = re.match(rules["rules"][i]["re.file"], file)
                match_folder = re.match(rules["rules"][i]["re.folder"], folder)
                if match_file and match_folder:
                    print("File: ", file, " matches rule: ", rules["rules"][i]["re.file"], " in folder: ", folder)
                    print("kind: ", rules["rules"][i]["kind"], " version: ", rules["rules"][i]["version"], " appendix: ", rules["rules"][i]["appendix"])
                    element = match_file.group(1)
                    break
            if match_file is None or match_folder is None:
                raise ValueError("No rule found for file: ", file, " in folder: ", folder)
