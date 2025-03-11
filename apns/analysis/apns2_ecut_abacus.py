"""this file is for greping and rendering svg plots of ecutwfc convergence test for new version of APNS"""
import os
import re
import json
import unittest

from apns.analysis.apns2_ecut_utils \
    import update_ecutwfc, build_sptc_from_nested, plot_log, plot_stack
from apns.analysis.apns2_utils import stru_rev_map
from apns.analysis.postprocess.read_abacus_out \
    import read_natom_fromlog, read_e_fromlog, read_pressure_fromlog, \
        read_istate, read_keyvals_frominput
from apns.analysis.apns2_utils import \
    read_apnsjob_desc, convert_fpp_to_ppid
from apns.pspot.parse import z_valence

def collect(folder: str):
    print("* * * Collect ABACUS result * * *".center(100))
    
    result = {}
    for root, _, files in os.walk(folder):
        for file in files:
            if re.match(r"running_(\w+)(\.log)", file): # reach the inner most folder like OUT.ABACUS
                natom = read_natom_fromlog(os.path.join(root, file))
                eks = read_e_fromlog(os.path.join(root, file))
                pressure = read_pressure_fromlog(os.path.join(root, file))
                bs = read_istate(os.path.join(root, "istate.info"))
                # continue if there is None among eks, pressure and bs
                parent = os.path.dirname(root)
                if None in [eks, pressure, bs]:
                    print(f"""WARNING: Present APNS job is broken: {parent}""")
                    continue
                bs = bs[0]
                atom_species, cellgen = read_apnsjob_desc(os.path.join(parent, "description.json"))
                system = os.path.basename(cellgen["config"])
                abacus_input = read_keyvals_frominput(os.path.join(parent, "INPUT"))
                ecutwfc = abacus_input["ecutwfc"]
                pps = [a["pp"] for a in atom_species]
                ppids = [": ".join(convert_fpp_to_ppid(pp)) for pp in pps]
                zvals = [float(z_valence(os.path.join(parent, pp))) for pp in pps]
                s = "\n".join(ppids)
                print(f"""In folder {parent}
Structure tested: {system}
Number of atoms: {natom}
ecutwfc: {ecutwfc}
Final Kohn-Sham energy: {eks}
Pressure: {pressure}
Pseudopotentials are used:\n{s}
""")
                data = {"ecutwfc": ecutwfc, "eks": eks, "pressure": pressure, "istate": bs, "natom": natom, "z_valence": zvals}
                # band structure is not easy to print, therefore omitted
                idx = -1 if result.get(system) is None \
                    or result[system].get("ppcases") is None \
                        or result[system]["ppcases"].count(pps) == 0 \
                    else result[system]["ppcases"].index(pps)
                if idx == -1:
                    result.setdefault(system, {"ppcases": [], "pptests": []}).setdefault("ppcases", []).append(pps)
                    result[system]["pptests"].append([data])
                else:
                    result[system]["pptests"][idx].append(data)
                #result[(system, "|".join(pps), ecutwfc)] = (natom, zvals, eks, pressure, bs)
    return result

class TestAPNS2EcutABACUS(unittest.TestCase):

    sysrevmap_ = stru_rev_map("./apns_cache/structures.json", True)
    jobpath = "/path/to/your/abacus/job/folder"

    @unittest.skip('we do not read')
    def test_read(self):
        collected = collect(self.jobpath)
        system_and_stpcs = build_sptc_from_nested(collected)
        result = []
        for s, stpcs in system_and_stpcs.items():
            for stpc in stpcs:
                pp, data = stpc()
                temp = {"name": self.sysrevmap_.get(s), "fcif": s, "pp": pp}
                temp['name'] = temp['fcif'] if temp['name'] is None else temp['name'] # use fcif if name is None
                temp.update(data)
                result.append(temp)
                ecut_conv = stpc.ecuts[stpc.iconv]
                pp = stpc.pp(as_list=True)
                
                # NOTE: for some hard-to-converge cases, other auxillary elements may be
                #       added to the pseudopotential list. In this case, the number of 
                #       psp is no longer 1. For this case, comment out the following
                #       assertion
                assert len(pp) == 1, "The pseudopotential should be unique for each test case"

                # update the local database if needed.
                update_ecutwfc(pp[0], ecut_conv)

        with open(os.path.basename(self.jobpath)+".json", "w") as f:
            json.dump(result, f, indent=4)

    @unittest.skip('we do not plot')
    def test_plot(self):
        with open('/the/json/dumped/by/test_read/func.json') as f:
            result = json.load(f)

        # result = sorted(result, key=lambda x: x["pp"])[:8]
        flogs = plot_log(result, 'png')
        fstacks = plot_stack(result, 'png')

if __name__ == "__main__":
    unittest.main()
