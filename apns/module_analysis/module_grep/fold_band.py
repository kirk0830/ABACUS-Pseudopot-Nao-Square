"""old version of driver for analysis of band structure of qespresso non-band structure calculation"""

import numpy as np
import matplotlib.pyplot as plt

def split_band(line: str) -> list:
    # first split by space:
    result = []
    words = line.split()
    for word in words:
        if(word.count("-")) > 1:
            subwords =  [f"-{x}" for x in word.split("-") if x]
            for subword in subwords:
                result.append(float(subword))
        else:
            result.append(float(word))
    return result

def calculate_band_gap(band_energies, nelec: int, efermi: float, metal = True):
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
        print("foldband>> number of electrons before adding the last: ", _nelec - sorted_[idx_level][1])
        print("foldband>> number of electrons after adding the last: ", _nelec)
        print("foldband>> total number of electrons of system: ", nelec)
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

    _energy = 0
    _degen = 1
    _end_index = 0
    
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

def parse_qe_output(filename, nband = 0):
    """
    parse QE stdout
    INPUT: stdout redirected file name
    OUTPUT: band energies, weight of kpoints
    """
    read_band_section = False
    read_kpoint_info = False
    band_energies = []
    band_index = 0
    line_index = 0
    kpoint_index = 0
    nkpoints = 0
    wk = []
    nelec = 0
    efermi = 0
    with open(filename, 'r') as f:
        line = f.readline()
        line_index += 1
        # if line != EOF
        while(True):
            clean_line = line.strip()
            """
            Read mode switch
            """
            # reach the end, jump out
            if clean_line.startswith("JOB DONE."):
                break
            # the beginning report
            if clean_line.startswith("number of electrons"):
                words = clean_line.split()
                nelec = int(words[-1].split(".")[0])
            # the beginning report
            if clean_line.startswith("number of Kohn-Sham states"):
                words = clean_line.split()
                nband = int(words[-1])
            # band structure result begin
            if clean_line.startswith("End of self-consistent calculation"):
                read_band_section = True
                read_kpoint_info = False
                band_index = 0
                kpoint_index = -1
                band_energies.clear()
            # band structure result end
            if clean_line.startswith("the Fermi energy is"):
                read_band_section = False
                read_kpoint_info = False
                words = clean_line.split()
                efermi = float(words[-2])
                #print("efermi = ", efermi)
            # the beginning report
            if clean_line.startswith("number of k points="):
                words = clean_line.split()
                for word in words:
                    if word.isdigit():
                        nkpoints = int(word)
                        read_band_section = False
                        read_kpoint_info = True
                        break
            """
            Contents parsing of each mode
            """
            # BAND STRUCTURE RESULT PARSE
            if read_band_section and len(clean_line) > 1:
                if (clean_line.find("k =") != -1) and (clean_line.find("PWs") != -1):
                    # if line starts with "k =", it is the beginning of a new kpoint
                    kpoint_index += 1
                    pass
                else:
                    # if line starts with number, it is still of present kpoint
                    if clean_line[0].isdigit() or clean_line[0] == '-':
                        # read band energy and store in band_energies
                        energies_one_line = [float(e) for e in split_band(clean_line)]
                        for e in energies_one_line:
                            band_energies.append((e, wk[kpoint_index]))
                            band_index += 1

            # KPOINT INFO PARSE
            if read_kpoint_info and len(wk) < nkpoints:
                if clean_line.startswith("k("):
                    kpoint_index += 1
                    words = clean_line.split()
                    for index, word in enumerate(words):
                        if word.startswith("wk"):
                            wk.append(float(words[index+2]))
                        
            line = f.readline() # read next line
            line_index += 1
            if (line_index % 100) == 0:
                #print("line_index = ", line_index)
                pass

    return band_energies, nelec, nband, efermi

def calculate_dos(band_energies, emin = -5.0, emax = 15.0, delta_e = 0.01, smear = 0.05):
    # draw probability distribution
    dos = np.zeros(shape = (int((emax - emin)/delta_e), 1))
    if smear != 0:
        for idx in range(dos.shape[0]):
            for band_energy in band_energies:
                dos[idx] += 1/(np.sqrt(2*np.pi)*smear)*np.exp(-(idx*delta_e+emin - band_energy[0])**2/(2*smear**2)) * band_energy[1]
    else:
        for idx in range(dos.shape[0]):
            for band_energy in band_energies:
                if (idx*delta_e+emin - band_energy[0]) < delta_e:
                    dos[idx] = band_energy[1]
    # normalize
    integral = 0
    for idx in range(dos.shape[0]):
        integral += dos[idx]*delta_e
    dos /= integral
    return dos

def plot_dos(dos: np.ndarray, homo: float, lumo: float, efermi: float, emin: float, emax: float, title = ""):

    x = np.linspace(emin, emax, dos.shape[0])

    fig, ax = plt.subplots()
    ax.plot(x, dos)
    ax.axvline(x = efermi, color = 'b', linestyle='-.')
    ax.axvline(x = homo, color = 'r')
    ax.axvline(x = lumo, color = 'r')
    plt.xlabel("energy (ev)")
    plt.ylabel("DOS")
    if title != "":
        plt.title(title)
    plt.xlim(emin, emax)
    plt.ylim(0, 2)
    ax.fill_between(x, dos.flatten(), y2 = -1, where=((x >= emin) & (x <= efermi)), alpha=0.5)
    plt.show()

if __name__ == "__main__":
    system = "t_pd_03_inas/scf.log"
    band_energies, nelec, nband, efermi = parse_qe_output(system)
    print("foldband>> SYSTEM: ", system)
    print("foldband>> nelec = ", nelec, " nband = ", nband, " efermi = ", efermi)
    window_width = 5
    emin = efermi - window_width/2
    emax = emin + window_width
    band_gap, band_energies, homo, lumo = calculate_band_gap(band_energies, nelec, efermi, True)
    print("foldband>> band gap = ", round(band_gap, 4), " ev")
    dos = calculate_dos(band_energies, emin, emax, 0.01, 0.02)
    _, _, end_indices, degen_data = scan_degeneracies(band_energies)
    print("foldband>> undegenerated band energies with (k,b)-weights:")
    print(band_energies)
    print("foldband>> degenerated band energies, degeneracies and end indices:")
    print(degen_data)

    plot_dos(dos, homo, lumo, efermi, emin, emax, system)
