"""Driver for get converged ecutwfc with respect to pressure"""
# tool collection for reading ABACUS output files,
# make sure the out dir is complete
import apns.module_analysis.postprocess.read_abacus_out as amarao
import apns.module_analysis.postprocess.conv.kernel as amack
import numpy as np
import os

def search_pressure(search_domain: str):
    path, pressures = [], []
    for root, dirs, files in os.walk(search_domain):
        for file in files:
            if file.startswith("running_") and file.endswith(".log"):
                # get the pressure trajectory
                press = amarao.read_pressure_fromlog(os.path.join(root, file))
                if press is not None:
                    path.append(root)
                    pressures.append(press)
    return path, np.array(pressures)

def run(search_domain: str, thr: float = 1e-1):
    first, second = amack.search(search_domain=search_domain,
                                 searcher=search_pressure)
    conv = amack.calculate(first, second, thr=thr)
    return conv, second

if __name__ == "__main__":
    search_domain = "../11549321"
    conv, _ = run(search_domain)
    print(conv)