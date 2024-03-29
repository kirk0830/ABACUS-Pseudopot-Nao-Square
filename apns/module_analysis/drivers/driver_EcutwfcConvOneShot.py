import os
import json

path = "../apns_conv_testresult/"
folders = os.listdir(path)
folders = [path + f for f in folders if os.path.isdir(path + f)]

for folder in folders:
    with open("input_analysis.json", "r") as f:
        data = json.load(f)
        data["analysis"]["search_domain"] = folder
    with open("input_analysis.json", "w") as f:
        json.dump(data, f, indent=4)
    os.system("python.exe main.py -i input_analysis.json")

    