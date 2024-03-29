import apns.module_analysis.postprocess.read_abacus_out as amarao
import numpy as np
import os

def default_calculator(val, val_ref):
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

def search(search_domain: str, searcher: callable, scalarizer: callable = None):

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
    scalarizer = default_calculator if scalarizer is None else scalarizer
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
