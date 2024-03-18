import apns.module_analysis.read_abacus_out as amarao
import numpy as np
import os

def search(search_domain: str, searcher: callable):

    paths, vals = searcher(search_domain=search_domain)
    ecutwfc = [float(amarao.read_keyvals_frominput(os.path.join(path, "INPUT"), "ecutwfc")) for path in paths]
    systems, mpids, pnids = [], [], []
    for path in paths:
        _, sys, mpid, pnid, _ = amarao.read_testconfig_fromBohriumpath(path)
        systems.append(sys)
        mpids.append(mpid)
        pnids.append(pnid)
    
    first = []
    second = []
    for s, m, p, e, v in zip(systems, mpids, pnids, ecutwfc, vals):
        # combine all (e, ener) pair with identical (s, m, p)
        if (s, m, p) not in first:
            first.append((s, m, p))
            second.append([])
            index = -1
        else:
            index = first.index((s, m, p))
        second[index].append((e, v))
    # sort second by e
    for i in range(len(second)):
        second[i] = sorted(second[i], key=lambda x: x[0])
        # tranverse to list of two tuples
        second[i] = list(zip(*second[i]))
        second[i][1] = np.array(second[i][1])
        second[i][0] = list(second[i][0])
        second[i][1] = (second[i][1] - second[i][1][-1]).tolist()

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
                conv.append((first[i][0], first[i][1], first[i][2], ecutwfc[j]))
                break
        else:
            conv.append((first[i][0], first[i][1], first[i][2], None))
    return conv
