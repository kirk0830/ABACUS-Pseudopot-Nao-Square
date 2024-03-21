import apns.module_analysis.postprocess.read_abacus_out as amarao
import apns.module_analysis.postprocess.utils as amapu
import numpy as np

def read(fistate, flog):
    return amarao.read_istate(fistate), amarao.read_etraj_fromlog(flog, term="efermi")[-1]
    
def integral(istate, virtual=True):
    """from one spin channel"""
    to_integral = [[], []]
    for i in range(len(istate)):
        # loop over kpoints
        e = istate[i][:, 1]
        occ = istate[i][:, 2]
        to_integral[0] += e.tolist()
        to_integral[1] += occ.tolist()

    # sort by energy
    to_integral = [list(x) for x in zip(*sorted(zip(to_integral[0], to_integral[1])))]
    x = np.array(to_integral[0])
    y = np.ones_like(x) if virtual else np.array(to_integral[1])
    # get the unique energies
    unique_x = np.unique(x)
    # accumulate the occ with the same energy
    unique_y = np.array([np.sum(y[x == e]) for e in unique_x])

    return unique_x, unique_y

def calculate_bandgap(x: np.ndarray, y: np.ndarray):
    """calculate on shifted data"""
    ind = np.argmin(np.abs(x))
    if y[ind] > 1e-6:
        return 0.0
    else:
        # -x search
        negind = ind
        while negind > 0 and abs(y[negind]) < 1e-6: # -x search for the first non-zero value
            negind -= 1
        # +x search
        posind = ind
        while posind < len(x) and abs(y[posind]) < 1e-6: # +x search for the first non-zero value
            posind += 1
        if negind == 0 or posind == len(x):
            return 0.0
        else:
            return x[posind] - x[negind]

def run(fistate, flog, de = 0.02):
    xs, ys, bandgaps = [], [], []
    istates, efermi = read(fistate, flog)
    for i in range(len(istates)):
        x, y = integral(istates[i])
        x = x - efermi
        y = y if i == 0 else -y
        x, y = amapu.zero_padding(min(x) // de * de, max(x) // de * de, de, x, y)
        xs.append(x)
        ys.append(y)
        bandgaps.append(calculate_bandgap(x, y))

    return xs, ys, bandgaps

def draw(fistate, flog, de = 0.02):
    colors = ["blue", "red"]
    xs, ys, _ = run(fistate, flog, de)
    xmin = xs[0][0]
    import matplotlib.pyplot as plt
    for i in range(len(xs)):
        y = amapu.Gaussian_filter(xs[i], ys[i], sigma=de*10)
        #y = ys[i]
        plt.plot(xs[i], y, color=colors[i], label="spin channel " + str(i))
        # fill x < 0
        plt.fill_between(xs[i], y, where=xs[i] < 0, color=colors[i], alpha=0.2)
        if xmin > xs[i][0]:
            xmin = xs[i][0]
    plt.axhline(0, color="black", linestyle="-")
    plt.axvline(0, linestyle="--")
    plt.xlim(xmin, 10)
    plt.xlabel("Energy (eV)")
    plt.ylabel("DOS")
    plt.legend()
    
    return plt
    
if __name__ == "__main__":

    root = "apns/module_analysis/postprocess/test/support"
    #xs, ys, bandgaps = run(root + "/cobalt_istate.info", root + "/cobalt_running_scf.log")
    #print(bandgaps)
    draw(root + "/cobalt_istate.info", root + "/cobalt_running_scf.log").show()