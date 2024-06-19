####
# utilties for calculating the band structure similarity between two band structures
####

def cal_nelec(occ: list, kwt: list) -> float:
    """integrate the occupation number over kpoints.
    Args:
        occ (list): the occupation number, in the form of [ispin][ik][iband]
        kwt (list): the weight of kpoints, in the form of [ik]
    """
    nelec = 0.0
    for i in range(len(occ)): # loop over spin
        for j in range(len(occ[i])): # loop over kpoints
            nelec += sum([occ[i][j][k] for k in range(len(occ[i][j]))]) * kwt[j]
    return nelec

def cal_bs_dist(band_struct_1, band_struct_2, smear: str, sigma: float, v: float) -> float:
    """
    Args:
        band_struct_1 (list): a nested list storing band structure and k-weighted occupation
        information: 
        [ispin][iks][ibands][1] -> energy (ekb),
        [ispin][iks][ibands][2] -> wk*occ (weight of kpoint multiplied with occupation)
        band_struct_2 (list): the same as band_struct_1 
        smear (str): smearing function type, can be "gaussian" or "fermi-dirac"
        sigma(float): smearing_sigma, in unit eV
        v (float): energy shift of fermi level
    Return:
        float: the "eta" defined above
    """
    return 0

