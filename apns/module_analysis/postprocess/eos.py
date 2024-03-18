import scipy.optimize as opt

def birch_murnaghan_eos(volumes, energies):
    """Fit the Birch-Murnaghan equation of state to the given volumes and energies."""
    def birch_murnaghan(v, e0, b0, b0p, v0):
        return e0 + 9 * v0 * b0 / 16 * (b0p * ((v0 / v)**(2/3) - 1)**3 + ((v0 / v)**(2/3) - 1)**2 * (6 - 4 * (v0 / v)**(2/3)))
    popt, _ = opt.curve_fit(birch_murnaghan, volumes, energies, p0=(energies[0], 1, 1, volumes[0]))
    return popt[3], popt[0], popt[1], popt[2]