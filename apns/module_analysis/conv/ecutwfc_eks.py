"""Driver for get converged ecutwfc with respect to energy per atom"""
# tool collection for reading ABACUS output files,
# make sure the out dir is complete
import apns.module_analysis.read_abacus_out as amarao
import apns.module_analysis.conv.kernel as amack
import numpy as np
import os

def search_eks(search_domain: str):
    paths, eks = [], []
    for root, dirs, files in os.walk(search_domain):
        for file in files:
            if file.startswith("running_") and file.endswith(".log"):
                # get the energy trajectory
                eks_traj = amarao.read_etraj_fromlog(os.path.join(root, file), unit="eV", term="EKS")
                if eks_traj is not None and len(eks_traj) > 0:
                    paths.append(root)
                    eks.append(eks_traj[-1])
    return paths, np.array(eks)

def run(search_domain: str):
    first, second = amack.search(search_domain=search_domain,
                                 searcher=search_eks)
    conv = amack.calculate(first, second, thr=1e-3)
    return conv, second

if __name__ == "__main__":
    search_domain = "../11549321"
    conv, _ = run(search_domain)
    print(conv)