
# ABACUS-Pseudopot-Nao-Sqaure
## Orbital generation
### Description
Files in this folder is from ABACUS official Numerical atomic orbital for SG15-1.0 pseudopotential. The file manage script is written like:
```python
import os
import re

path_backup = os.path.abspath(os.getcwd())
origin_path = "D:/abacus-pseudopot-nao-square/download/numerical_orbitals/sg15_10/SG15-Version1p0__AllOrbitals-Version2p0"
os.chdir(origin_path)
os.chdir("..")

visit_history = []
for root, dirs, files in os.walk(os.getcwd()):
    
    if "SIAB_INPUT" in files:
        _match = re.match(r"(.*)([0-9]+\_[A-Z][a-z]?\_(TZDP|DZP|SZ))(.*)", root)
        if _match:
            print(_match.groups())
            element = _match.group(2).split("_")[1]
            if element in visit_history:
                continue
            else:
                visit_history.append(element)
                print(f"element: {element}, find reference file for generating orbitals")
                os.system(f"cp {root}/SIAB_INPUT {path_backup}/{element}.SIAB_INPUT")
    continue # comment this line to enable the following code, moving numerical orbitals to folders
    for file in files:
        _match = re.match(r"^([A-Z][a-z]?)(\_gga\_)([0-9]+)(au\_)([0-9]+)(Ry\_)([\w]+)(.orb)$", file)
        if _match:
            element = _match.group(1)
            rcut = float(_match.group(3))
            ecutwfc = float(_match.group(5))
            zeta = root.replace("\\", "/").split("/")[-1].split("_")[-1]
            print(f"element: {element}, rcut: {rcut}, ecutwfc: {ecutwfc}, zeta: {zeta}")
            folder = f"{rcut}au_{ecutwfc}Ry_{zeta}"
            os.makedirs(folder, exist_ok=True)
            os.system(f"mv {root}/{file} {folder}/{file}")

os.chdir(path_backup)
```