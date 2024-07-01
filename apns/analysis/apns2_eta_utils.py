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

    # convert to arrays for convenience
    be1 = np.array(band_energy1)
    be2 = np.array(band_energy2)
    
    # convert spinless weight to the one with spin
    nspin = len(be1)
    wk = np.array(wk).reshape(1, len(wk), 1) * (2 / nspin)

    n_elec1, n_elec2 = n_elec if isinstance(n_elec, tuple) else (n_elec, n_elec)

    assert be1.shape == be2.shape and be1.shape[1] == wk.shape[1]
    assert smearing_sigma >= 0 and n_elec1 > 0 and n_elec2 > 0

    # determine the Fermi level
    efermi1 = brentq(lambda x: np.sum(wk * f_occ(be1, x)) - n_elec1, np.min(be1), np.max(be1))
    efermi2 = brentq(lambda x: np.sum(wk * f_occ(be2, x)) - n_elec2, np.min(be2), np.max(be2))

    # geometrically averaged occupation (under shifted Fermi level)
    f_avg = np.sqrt(f_occ(be1, efermi1 + efermi_shift) * f_occ(be2, efermi2 + efermi_shift))

    res = minimize_scalar(lambda omega: np.sum(wk * f_avg * (be1 - be2 + omega)**2), \
            (-10, 10), method='brent')

    omega = res.x
    eta = np.sqrt(res.fun / np.sum(wk * f_avg))
    eta_max = np.max(np.abs(be1 - be2 + omega))

    return (eta, eta_max) if not return_all else (eta, eta_max, efermi1, efermi2, omega)