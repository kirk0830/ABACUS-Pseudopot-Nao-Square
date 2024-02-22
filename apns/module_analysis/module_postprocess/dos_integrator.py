"""old version of driver for analysis of band structure of qespresso non-band structure calculation"""

import numpy as np
import matplotlib.pyplot as plt

def cal_bandgap(band_energies, nelec: int, efermi: float, metal = True):
    """find band gap, only for nspin = 1"""

    if len(band_energies) < nelec:
        raise ValueError("number of electrons is larger than number of bands")
    
    dtype = [('energy', float), ('wk', float)]
    to_sort = np.array(band_energies, dtype=dtype)
    sorted_ = np.sort(to_sort, order='energy')
    _nelec = 0
    idx_level = 0
    while (_nelec < nelec):
        _nelec += sorted_[idx_level][1]
        idx_level += 1
    idx_level -= 1

    homo = sorted_[idx_level][0]
    print("foldband>> find homo with energy: {0} ev, its (k,b)-weight is: {1}".format(homo, sorted_[idx_level][1]))
    if(homo > efermi):
        if (nelec > (_nelec - sorted_[idx_level][1])) and (nelec < _nelec):
            Warning("foldband>> It maybe a gapless system!")
            homo = efermi
            lumo = efermi
            band_gap = 0
            print("foldband>> ***WARNING*** you may meet a gapless system, set homo = lumo = {0} ev and band gap = 0 ev".format(efermi))
            return 0, sorted_, efermi, efermi
        else:
            raise ValueError("foldband>> homo energy is larger than efermi")
    
    lumo = sorted_[idx_level+1][0]
    band_gap = lumo - homo

    #while (band_gap < 1E-3 and metal):
    #    idx_level += 1
    #    lumo = sorted_[idx_level][0]
    #    band_gap = lumo - homo
    print("foldband>> find lumo with energy: {0} ev, its (k,b)-weight is: {1}".format(lumo, sorted_[idx_level][1]))
    return band_gap, sorted_, homo, lumo

def scan_degeneracies(band_energies):

    energies = []
    degeneracies = []
    end_indices = []

    for index, band_energy in enumerate(band_energies):
        if index != 0:
            if energies[-1] != band_energy[0]:
                energies.append(band_energy[0])
                degeneracies.append(1)
                end_indices.append(index)
            else:
                degeneracies[-1] += 1
                end_indices[-1] += 1
        else:
            energies.append(band_energy[0])
            degeneracies.append(1)
            end_indices.append(index)
    dtype = [('energy', float), ('degeneracies', int), ('end_index', int)]
    result = []
    for index, energy in enumerate(energies):
        result.append((energy, degeneracies[index], end_indices[index]))
    result = np.array(result, dtype = dtype)

    return energies, degeneracies, end_indices, result

def smooth_1d(source, lb=-5.0, ub=15.0, dx=0.01, smear=0.05):

    x = np.zeros(shape=(int((ub - lb) / dx), 1))
    result = np.zeros_like(x)
    for idx in range(result.shape[0]):
        for xypair in source:
            de = idx * dx + lb - xypair[0]
            if smear != 0:
                result[idx] += (1 / (np.sqrt(2 * np.pi) * smear)) * np.exp(-(de ** 2) / (2 * smear ** 2)) * xypair[1]
            else:
                if de < dx:
                    result[idx] = xypair[1]
            x[idx] = idx * dx + lb
    result /= np.sum(result) * dx
    return x, result

def plot_dos(dos: np.ndarray, homo: float, lumo: float, efermi: float, emin: float, emax: float, title: str = "") -> None:
    """
    Plot the Density of States (DOS) with energy levels and Fermi level.

    Args:
        dos (np.ndarray): Array of DOS values.
        homo (float): Highest Occupied Molecular Orbital (HOMO) energy level.
        lumo (float): Lowest Unoccupied Molecular Orbital (LUMO) energy level.
        efermi (float): Fermi energy level.
        emin (float): Minimum energy value for the plot.
        emax (float): Maximum energy value for the plot.
        title (str, optional): Title for the plot. Defaults to "".
    """
    x = np.linspace(emin, emax, len(dos))

    fig, ax = plt.subplots()
    ax.plot(x, dos)
    ax.axvline(x=efermi, color='b', linestyle='-.')
    ax.axvline(x=homo, color='r')
    ax.axvline(x=lumo, color='r')
    plt.xlabel("Energy (eV)")
    plt.ylabel("DOS")
    if title:
        plt.title(title)
    plt.xlim(emin, emax)
    plt.ylim(0, 1.2 * max(dos))
    ax.fill_between(x, dos.flatten(), y2=-1, where=((x >= emin) & (x <= efermi)), alpha=0.5)
    plt.show()

if __name__ == "__main__":
    #system = "t_pd_03_inas/scf.log"
    #band_energies, nelec, nband, efermi = parse_qe_output(system)
    import apns.module_analysis.module_grep.abacus as amamga
    system = "./apns/module_analysis/module_grep/test/support/Lu.log"
    band_energies, nelec, nband, efermi = amamga.grep_band(system)
    print("foldband>> SYSTEM: ", system)
    print("foldband>> nelec = ", nelec, " nband = ", nband, " efermi = ", efermi)
    window_width = 20
    emin = efermi - window_width/2
    emax = emin + window_width
    band_gap, band_energies, homo, lumo = cal_bandgap(band_energies, nelec, efermi)
    print("foldband>> band gap = ", round(band_gap, 4), " ev")
    e, dos = smooth_1d(source=band_energies, 
                       lb=emin, 
                       ub=emax, 
                       dx=0.01, 
                       smear=0.05)
    _, _, end_indices, degen_data = scan_degeneracies(band_energies)
    print("foldband>> undegenerated band energies with (k,b)-weights:")
    print(band_energies)
    print("foldband>> degenerated band energies, degeneracies and end indices:")
    print(degen_data)

    plot_dos(dos, homo, lumo, efermi, emin, emax, system)
