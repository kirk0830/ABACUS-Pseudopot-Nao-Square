"""generate series of SIAB_INPUT in-one-shot

Only following sections in input are allowed to be generator to iterate on:
pwsets"""

def convert_orbset(orbset: list):
    """convert the given orbset to the one in SIAB_INPUT"""
    result = []
    for orb in orbset:
        keys = ["zeta_notation", "shape", "orb_ref", "nbands_ref"]
        vals = [orb["conf"], orb["shape"], orb["dep"], orb["states"]]
        result.append(dict(zip(keys, vals)))
    return result

def convert_calcset(calcset: dict, pp: str, rcuts: list):
    from apns.test.dft_param_set import DFTParamSetGenerator
    result = {"pseudo_name": pp, "bessel_nao_rcut": rcuts}
    calcgen = DFTParamSetGenerator(calcset)
    for calc in calcgen():
        yield {**result, **calc}

def build(calcset: dict,
          siabset: dict,
          struset: list,
          orbset: list,
          environment: str = "", 
          mpi_command: str = "mpirun -np 16", 
          abacus_command: str = "abacus"):
    result = {}
    result.update(calcset)
    result.update({k: v for k, v in siabset.items() if k != "rcuts"})
    result["reference_systems"] = struset
    result["orbitals"] = convert_orbset(orbset)
    result.update({"environment": environment, "mpi_command": mpi_command, "abacus_command": abacus_command})
    return result

def gen(finp: str):
    import json, uuid, os
    from apns.test.atom_species_and_cell import AtomSpeciesGeneartor
    with open(finp, 'r') as f:
        inp = json.load(f)
    with open(os.path.join(inp["global"]["cache_dir"], "ecutwfc.json"), "r") as f:
        ecuts = json.load(f)
    out_dir = inp["global"].get("out_dir", "./")
    for task in inp["tasks"]:
        atomset = inp["ppsets"][task["pp"]]
        calcset = inp["pwsets"][task["pw"]]
        siab = inp["siabsets"][task["siab"]]
        tags = atomset["tags"]
        orb = inp["orbsets"][task["orb"]]
        struset = inp["strusets"][task["stru"]]
        for element in atomset["elements"]:
            print(f"* * * Generate SIAB_INPUT for {element} * * *".center(100))
            asgen = AtomSpeciesGeneartor(element, inp["global"]["pseudo_dir"], tags)
            for atom in asgen():
                ecut = max(100, ecuts.get(atom.pp, 100)) # let it no larger than 100 Ry
                for calc in convert_calcset(calcset, atom.pp, siab["rcuts"]):
                    folder = str(uuid.uuid4())
                    out = {"pseudo_dir": "./", "ecutwfc": ecut}
                    out.update({k: v for k, v in build(calc, siab, struset, orb).items()\
                                if k not in ["pseudo_dir", "ecutwfc"]})
                    yield os.path.join(out_dir, folder), out

def write_autorun(finp: str):
    import json
    import os
    with open(finp, 'r') as f:
        inp = json.load(f)
    out_dir = inp["global"].get("out_dir", "./")
    siab_dir = inp["global"].get("siab_dir", "./")
    siab = os.path.join(siab_dir, "SIAB_nouvelle.py")
    assert os.path.exists(siab), \
        f"SIAB_nouvelle.py not found in {siab_dir}"
    
    script = f"""import os
cwd = os.path.abspath(os.getcwd())
folders = [f for f in os.listdir(cwd) if os.path.isdir(f)]
for folder in folders:
    if "SIAB_INPUT.json" in os.listdir(folder):
        os.chdir(folder)
        print(f"Running SIAB in {folder}", flush=True)
        os.system("nohup python3 /root/deepmodeling/abacus_orbital_generation/SIAB/SIAB_nouvelle.py -i SIAB_INPUT.json > log")
        os.chdir(cwd)
        
"""
    with open(os.path.join(out_dir, "autorun.py"), 'w') as f:
        f.write(script)

def run(finp: str):
    import json, shutil, os
    from apns.test.citation import citation
    for folder, task in gen(finp):
        os.makedirs(folder, exist_ok=True)
        # copy the pseudopotential file into the directory
        shutil.copy(task["pseudo_name"], folder)
        task["pseudo_name"] = os.path.basename(task["pseudo_name"])
        with open(os.path.join(folder, "SIAB_INPUT.json"), 'w') as f:
            json.dump(task, f, indent=4)
    write_autorun(finp)
    citation()

if __name__ == "__main__":
    run("./orbgen.json")
