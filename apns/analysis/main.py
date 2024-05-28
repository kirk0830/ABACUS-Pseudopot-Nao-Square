import json
import os

def run(finp: str):
    with open(finp, "r") as f:
        inp = json.load(f)

    for item in inp["analysis"]["items"]:
        
        item = item.replace("\\", "/").split("/")[-1]
        item = item + ".py" if not item.endswith(".py") else item
        item = "apns/module_analysis/drivers/" + item
        if os.path.exists(item):
            os.system("python " + item + " -i" + inp["analysis"]["search_domain"])
        else:
            print("Warning: user-defined analysis item \"", item, "\" not found, skip.")

import unittest
import sys
class TestRun(unittest.TestCase):
    def test_run(self):
        path_backup = os.getcwd()
        os.chdir("apns/module_analysis/drivers")
        contents = """# only for unittest of apns.module_analysis.workflow_analysis.driver
if __name__ == "__main__":
    print("Hello, world!")
"""
        with open("test_run.py", "w") as f:
            f.write(contents)
        os.chdir(path_backup)
        inp = {"analysis": {"items": ["test_run"]}}
        with open("test_run.json", "w") as f:
            json.dump(inp, f)
        # catch stdout of test_run.py
        f = open("apns/module_analysis/drivers/stdout", "w")
        sys.stdout = f
        run("test_run.json")
        f.close()
        sys.stdout = sys.__stdout__
        with open("apns/module_analysis/drivers/stdout", "r") as f:
            self.assertEqual(f.read(), "Hello, world!\n")
        os.remove("apns/module_analysis/drivers/stdout")
        os.remove("test_run.json")
        os.remove("apns/module_analysis/drivers/test_run.py")

if __name__ == "__main__":
    unittest.main()