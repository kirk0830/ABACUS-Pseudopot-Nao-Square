####
# utilties for calculating the band structure similarity between two band structures
####

def cal_nelec(occ: list, kwt: list = None) -> float:
    """integrate the occupation number over kpoints.
    Args:
        occ (list): the occupation number, in the form of [ispin][ik][iband]
        kwt (list, optional): the weight of kpoints, in the form of [ik]. Defaults to a list of ones if not provided.
    """
    nelec = 0.0
    kwt = kwt or [1.0] * len(occ[0])
    for i in range(len(occ)): # loop over spin
        for j in range(len(occ[i])): # loop over kpoints
            nelec += sum([occ[i][j][k] for k in range(len(occ[i][j]))]) * kwt[j]
    return nelec

def delta_band(band_energy1, band_energy2, n_elec, wk, smearing, smearing_sigma, efermi_shift = 0, return_all = False):
    '''
    Calculate the "distance" between two band structures.

    Parameters
    ----------
        band_energy1, band_energy2 : list
            Nested lists that contain the band energies.
            band_energy[ispin][ik][iband] specifies the band energy of a state.
        wk : list
            Weight of k-points.
        n_elec : float or tuple
            Total number of electrons used to determine the Fermi level.
            If it is a tuple, the first element is for band_energy1 and
            the second is for band_energy2.
        smearing : str
            Smearing method, can be 'gaussian' or 'fermi-dirac'
        smearing_sigma : float
            Smearing parameter.
        efermi_shift : float
            Energy shift of the Fermi level.
        return_all : bool
            If True, return a tuple (eta, efermi1, efermi2, omega) where
            eta is the distance, efermi1 and efermi2 are the Fermi levels,
            and omega is the optimized energy shift.
            If False, return only eta.

    '''
    import numpy as np
    from scipy.special import erf
    from scipy.optimize import brentq, minimize_scalar
    # occupation function
    def f_occ(x, x0):
        if smearing == 'gaussian':
            return 0.5 * (1.0 - erf((x - x0) / smearing_sigma)) \
                if smearing_sigma > 0 else 0.5 * (1 - np.sign(x - x0))
        elif smearing == 'fermi-dirac':
            return 1.0 / (1.0 + np.exp((x - x0) / smearing_sigma)) \
                if smearing_sigma > 0 else 0.5 * (1 - np.sign(x - x0))
        else:
            raise ValueError('Unknown smearing method: %s'%smearing)

    def efermi(wk, be, n_elec):
        _nmax = np.sum(wk * f_occ(be, np.max(be)))
        _delta = (_nmax - n_elec) / n_elec
        if np.abs(_delta) < 1e-4 or np.abs(_delta * n_elec) <= 0.01: # 0.1% error: the case where all bands are occupied
            print(f"WARNING: all bands are occupied in band_energy1, error of this estimation: {_delta:.4%}")
            return np.max(be)
        else: # if error is too large, let it directly fail
            if _delta < 0:
                raise ValueError(f"""WARNING: maximum possible number of electrons in band structure not-enough:
{n_elec:.4f} vs. {_nmax:.4f} (nelec vs. nmax). This is always because of too small basis size and all
bands are occupied, otherwise please check your data.""")
            return brentq(lambda x: np.sum(wk * f_occ(be, x)) - n_elec, np.min(be), np.max(be))
        
    # convert to arrays for convenience
    be1 = np.array(band_energy1)
    be2 = np.array(band_energy2)
    
    # convert spinless weight to the one with spin
    nspin = len(be1)
    wk = np.array(wk).reshape(1, len(wk), 1) * (2 / nspin)
    wk = 2 * wk/np.sum(wk) # normalize the weight

    n_elec1, n_elec2 = n_elec if isinstance(n_elec, tuple) else (n_elec, n_elec)

    if be1.shape != be2.shape:
        raise TypeError(f'Error: Inconsistent shape between two band structures: {be1.shape} vs {be2.shape}.')
    assert be1.shape[1] == wk.shape[1]
    assert smearing_sigma >= 0 and n_elec1 > 0 and n_elec2 > 0

    # determine the Fermi levels for two band structures by root finding
    efermi1 = efermi(wk, be1, n_elec1)
    efermi2 = efermi(wk, be2, n_elec2)

    # geometrically averaged occupation (under shifted Fermi level)
    f_avg = np.sqrt(f_occ(be1, efermi1 + efermi_shift) * f_occ(be2, efermi2 + efermi_shift))

    res = minimize_scalar(lambda omega: np.sum(wk * f_avg * (be1 - be2 + omega)**2), \
            (-10, 10), method='brent')
    
    omega = res.x
    eta = np.sqrt(res.fun / np.sum(wk * f_avg))
    eta_max = np.max(np.abs(be1 - be2 + omega))

    # for ispin in range(nspin):
    #     for ik in range(len(be1[ispin])):
    #         delta = np.array(be1[ispin][ik]) - np.array(be2[ispin][ik]) + omega
    #         # zval = 19, natom = 4, nocc = 19*2 = 38
    #         if np.linalg.norm(delta, np.inf)*1e3 > 100:
    #             print(f_occ(be1, efermi1 + efermi_shift)[ispin][ik], delta)
    #             print(ispin, ik, np.linalg.norm(delta, np.inf))

    return (eta, eta_max) if not return_all else (eta, eta_max, efermi1, efermi2, omega)

import unittest
class TestAPNS2EtaUtilities(unittest.TestCase):
    def test_cal_nelec(self):
        pass

    def test_delta_band(self):
        pass

if __name__ == "__main__":
    unittest.main()