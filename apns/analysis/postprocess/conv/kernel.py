import apns.analysis.postprocess.read_abacus_out as amarao
import numpy as np
import os

def cal_diff(val, val_ref):
    """a default calculator, to pass as argument to amack.calculate,
    if no calculator is provided"""
    if isinstance(val, list):
        if not isinstance(val_ref, list):
            return (np.array(val) - val_ref).tolist()
        else:
            return (np.array(val) - np.array(val_ref)).tolist()
    else:
        if not isinstance(val_ref, list):
            return val - val_ref
        else:
            raise TypeError("Tend to scalarize one value respect to one list.")

def Yb_special_case(path):
    """this is a special treatment for the ecutwfc convergence test on Yb, which is, not
    recored in Materials Project but available in Crystallography Open Database (COD). The
    convergence test was happen to be carried out during the refactor of APNS, therefore
    there is no source information where the structure from, and lack of the pseudopotential
    information.
    Luckily, for Yb, there are only finite cases so we can hard-code.
    
    The work's Bohrium jobgroup id is: 12310698
    Download via Lesbegue's API:

    ```bash
    lbg jobgroup download 12310698
    ```
    
    each job folder is like:
    "/root/documents/simulation/abacus/12310698/12902679/tmp
    /outputs/artifacts/outputs/8b044239-815e-416d-ac14-7a23030e8cad/00003/OUT.ABACUS"
    """
    import os
    print(f"Parsing {path}")
    frags = path.replace("\\", "/").split("/")
    assert frags[-1] == "OUT.ABACUS", f"Not expected path: {path}"
    basename = "/".join(frags[:-1])
    print(f"Searching {basename}")
    files = os.listdir(basename) # list all files in the parent dir
    pnid = {"Yb3+_f--core-icmod1.PD04.PBE.UPF": "pd043+f--core-icmod1",
            "Yb.pbe-n-nc.UPF": "pslnc031",
            "Yb.upf": "dojo043plus",
            "Yb_GTH_NC_LnPP1.upf": "gth",
            "Yb-sp.PD04.PBE.UPF": "pd04sp",
            "Yb3+_f--core.PD04.PBE.UPF": "pd043+f--core"} # there are all possible pseudopotentials
    fupf = [f for f in files if f.upper().endswith(".UPF")][0]
    return "Yb", "cod-9010997.cif", pnid[fupf]

def search(search_domain: str, searcher: callable, scalarizer: callable = None):

    paths, vals = searcher(search_domain=search_domain)
    ecutwfc = [float(amarao.read_keyvals_frominput(os.path.join(path, "INPUT"), "ecutwfc")) for path in paths]
    systems, mpids, pnids = [], [], []
    for path in paths:
        _, sys, mpid, pnid, _ = amarao.read_testconfig_fromBohriumpath(path)
        if all([not sys, not mpid, not pnid]): # Yb special case
            sys, mpid, pnid = Yb_special_case(path)
        systems.append(sys)
        mpids.append(mpid)
        pnids.append(pnid)
    
    first = []
    second = []
    for s, m, p, e, v in zip(systems, mpids, pnids, ecutwfc, vals):
        # ecutwfc and vals are distinct for each (s, m, p), while
        # the tuple (s, m, p) might be repeated, because it is for
        # one system, one mpid and one pnid

        # combine all (e, ener) pair with identical (s, m, p)
        if (s, m, p) not in first:
            first.append((s, m, p))
            second.append([])
            index = -1
        else:
            index = first.index((s, m, p))
        second[index].append((e, v))
        # therefore "second" is a list indiced by unique (s, m, p)
        # for each (s, m, p), it is a list of (e, val) pair, no matter
        # what exactly the val is, it is the value to be compared

    # sort "second" by e
    for i in range(len(second)): # loop over all (s, m, p)...
        second[i] = sorted(second[i], key=lambda x: x[0])
        # tranverse to list of two tuples
        second[i] = list(zip(*second[i]))
        second[i][0] = list(second[i][0])
        second[i][1] = list(second[i][1])

    # normalize data
    scalarizer = cal_diff if scalarizer is None else scalarizer
    for i in range(len(second)): # loop over all (s, m, p)...
        # normalize the energy to the last value
        second[i][1] = scalarizer(second[i][1], second[i][1][-1])

    return first, second

def calculate(first: list, second: list, thr: float = None):

    thr = 1e-3 if thr is None else thr

    nsuite = len(first)
    assert nsuite == len(first)
    assert nsuite == len(second)

    conv = []
    for i in range(nsuite):
        ecutwfc = second[i][0]
        vals = second[i][1]
        assert len(ecutwfc) == len(vals)    
        for j in range(len(ecutwfc)):
            if abs(vals[j]) < thr:
                #            system       mpid         pnid         ecutwfc
                conv.append((first[i][0], first[i][1], first[i][2], ecutwfc[j]))
                break
        else:
            conv.append((first[i][0], first[i][1], first[i][2], None))
    return conv
