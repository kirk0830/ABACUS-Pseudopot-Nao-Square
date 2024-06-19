"""this driver is for gathering Cohesive energy test data from abacus"""
def collect_jobs(folder: str):
    """Collect APNS jobs for cohesive energy calculation"""
    import os, re, json
    from apns.analysis.postprocess.read_abacus_out import read_e_fromlog
    print("* * * Collect ABACUS result * * *".center(100))
    result = []
    for root, _, files in os.walk(folder):
        for file in files:
            if re.match(r"(running_)(\w+)(\.log)", file):
                # this seems to be a outdir of abacus, but the present folder is the work folder
                # also need something in the work folder to identify the job
                parent = os.path.dirname(root)
                # all calculation setting can be found in the description.json
                with open(os.path.join(parent, "description.json"), "r") as f:
                    desc = json.load(f)
                eks = read_e_fromlog(os.path.join(root, file))
                if eks is None:
                    print(f"WARNING: Present APNS job is broken: {parent}")
                    continue
                else:
                    result.append((desc, eks))
    return result

def main(folder: str):
    import os
    from apns.analysis.drivers.apns2_ecoh_utils import pair, cal_e_cohesive
    data = collect_jobs(folder)
    paired = pair(data)
    for desc, e_b, e_a in paired:
        natom = len(desc["Cell"]["coords"])
        e_coh = cal_e_cohesive(e_b, e_a, natom)
        system = os.path.basename(desc["CellGenerator"]["config"])
        pps = [os.path.basename(a["pp"]) for a in desc["AtomSpecies"]]
        pps = "|".join(pps)
        orbs = [os.path.basename(a["nao"]) for a in desc["AtomSpecies"]]
        orbs = "|".join(orbs)
        print(f"""For system {system}, 
Number of atoms: {natom},
Pseudopotentials:\n{pps},
Numerical atomic orbitals:\n{orbs},
Bulk energy: {e_b} eV,
Atom energy: {e_a} eV,
Cohesive energy: {e_coh} eV/atom.""")
        