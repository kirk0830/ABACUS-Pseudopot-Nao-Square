import apns.module_analysis.read_abacus_out as amarao
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
        while negind > 0 and y[negind] < 1e-6:
            negind -= 1
        # +x search
        posind = ind
        while posind < len(x) and y[posind] < 1e-6:
            posind += 1
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
    import matplotlib.pyplot as plt
    for i in range(len(xs)):
        y = amapu.Gauss_smoothing(xs[i], ys[i], sigma=de*20)
        plt.plot(xs[i], y, color=colors[i], label="spin channel " + str(i))
        # fill x < 0
        plt.fill_between(xs[i], y, where=xs[i] < 0, color=colors[i], alpha=0.2)
    plt.axhline(0, color="black", linestyle="-")
    plt.axvline(0, linestyle="--")
    plt.show()
    

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    root = "apns/module_analysis/postprocess/test/support"
    #xs, ys, bandgaps = run(root + "/istate.info", root + "/running_scf.log")
    
    draw(root + "/istate.info", root + "/running_scf.log")