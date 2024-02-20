from tools.orbital_oscillation_remove.orb_parser import data_import, filename_parse, plot_orb_json
from tools.orbital_oscillation_remove.ovlp_jjtilde_calculator import generate_truncated_spherical_bessel
from tools.orbital_oscillation_remove.database import symbol_tol, l_tosymbol
from orb_io import dict_to_orb, dict_to_c4
import numpy as np
from scipy.integrate import simps
from scipy.special import spherical_jn as jn
from scipy.linalg import solve, cho_factor, cho_solve
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

def oscillation_remove(path: str, orbital_file: str, n_jtilde_remove: list, linalg_solver = 'numpy') -> tuple[dict, dict, dict, dict]:
    """ remove high frequency oscillation from numerical atomic orbitals """
    
    print("[nao_oscillation_remove: oscillation_remove] linear algebra solver: %s" % linalg_solver)
    numerical_orbitals = data_import(path, orbital_file) # the input numerical atomic orbitals
    _out = {}                                            # the output numerical atomic orbitals
    _c_in = {}                                           # the input coefficients of the truncated spherical bessel functions
    _c_out = {}                                          # the output coefficients of the truncated spherical bessel functions
    orbital_basic_info = filename_parse(orbital_file)
    bessel_nao_ecut = orbital_basic_info['bessel_nao_ecut']
    bessel_nao_rcut = orbital_basic_info['bessel_nao_rcut']

    r = np.array(numerical_orbitals['r'])

    # solve Ax = y equation, 
    # A is overlap between truncated spherical bessel functions, 
    # x is coefficient (of tsbf)-vector from numerical atomic orbital
    # y is overlap between numerical atomic orbital and tsbf
    # calculate A
    for key in numerical_orbitals.keys():
        if key == 'r':
            _out['r'] = numerical_orbitals['r']
        else:
            l = symbol_tol(key)
            # initialize _c_in and _c_out
            _c_in[key] = []
            _c_out[key] = []
            # generate q of spherical bessel functions for current l
            qs = generate_truncated_spherical_bessel(bessel_nao_ecut, bessel_nao_rcut, [l])
            print("|-calculating overlap matrix between truncated spherical bessel functions for l = %d" % l)
            # allocate memory for overlap matrix between truncated spherical bessel functions
            ovlp_jtildejtilde = np.zeros(shape=(len(qs[0]), len(qs[0]))) # A
            # loop over q...
            for iq in range(len(qs[0])):
                # generate one of the truncated spherical bessel functions
                jlq_i = jn(l, qs[0][iq]*r)
                for jq in range(len(qs[0])):
                    # generate another one of the truncated spherical bessel functions, but only if iq >= jq, because the overlap matrix is symmetric
                    if iq >= jq:
                        jlq_j = jn(l, qs[0][jq]*r)
                        ovlp_jtildejtilde[iq][jq] = simps(jlq_i*jlq_j*r**2, r)
                        # fill the other half of the overlap matrix
                        if iq != jq:
                            ovlp_jtildejtilde[jq][iq] = ovlp_jtildejtilde[iq][jq]
            # loop over numerical atomic orbitals for current l
            for inao in range(len(numerical_orbitals[key])):
                print("|--calculating overlap matrix between numerical atomic orbital and truncated spherical bessel functions for l = %d, inao = %d" % (l, inao))
                # for current numerical orbital...1d array
                orb = np.array(numerical_orbitals[key][inao])
                # calculate y
                # allocate memory for overlap between numerical atomic orbital and truncated spherical bessel functions
                ovlp_fjtilde = np.zeros(shape=(len(qs[0]), 1)) # y
                for iq in range(len(qs[0])):
                    # generate one of the truncated spherical bessel functions
                    jlq_i = jn(l, qs[0][iq]*r)
                    # calculate overlap between numerical atomic orbital and truncated spherical bessel functions
                    ovlp_fjtilde[iq] = simps(orb*jlq_i*r**2, r)
                # solve x
                if linalg_solver == 'numpy':
                    c = np.linalg.solve(ovlp_jtildejtilde, ovlp_fjtilde)
                elif linalg_solver == 'scipy':
                    c = solve(ovlp_jtildejtilde, ovlp_fjtilde)
                elif linalg_solver == 'scipy_cholesky':
                    cho, low = cho_factor(ovlp_jtildejtilde)
                    c = cho_solve((cho, low), ovlp_fjtilde)
                _c_in[key].append(c)
                _c_out[key].append(np.zeros_like(c))
                print("the obtained c has shape: ", c.shape)
                # re-construct orbitals
                print("|--re-constructing numerical atomic orbital for l = %d, inao = %d. number of truncated spherical bessel functions to remove = %d" % (l, inao, n_jtilde_remove[l]))
                # allocate memory for reconstructed orbital
                new_orb = np.zeros_like(orb)
                # for each r, calculate new combination of truncated spherical bessel functions
                for ir, _r in enumerate(r):
                    # loop over q, but only up to the number of truncated spherical bessel functions to remove
                    for iq in range(len(qs[0]) - n_jtilde_remove[l]):
                        new_orb[ir] += jn(l, qs[0][iq]*_r)*c[iq][0]
                # add new C information to _c_out
                for iq in range(len(qs[0]) - n_jtilde_remove[l]):
                    _c_out[key][inao][iq][0] = c[iq][0]
                # make the newly added coefficients in _c_in and _c_out to 1d array
                _c_in[key][inao] = _c_in[key][inao].flatten()
                _c_out[key][inao] = _c_out[key][inao].flatten()
                if key not in _out.keys():
                    _out[key] = []
                _out[key].append(new_orb)
    return _out, numerical_orbitals, _c_out, _c_in

def compare_plot(old_orb: dict, new_orb: dict, ncol = 3):

    ncol = 3
    nplots = 0
    for key in new_orb.keys():
        if key == 'r':
            continue
        else:
            nplots += len(new_orb[key])
    nrow = int(nplots/ncol) + 1
    # plot
    fig, axs = plt.subplots(nrow, ncol, sharex=True, sharey=True)
    # set size
    fig.set_size_inches(20, 10)
    fig.suptitle('oscillation remove of orbital file %s' % orbital_file, fontsize=20)
    irow = 0
    icol = 0
    for key in new_orb.keys():
        if key == 'r':
            continue
        else:
            for i in range(len(new_orb[key])):
                # set in as semi-transparent
                axs[irow][icol].plot(old_orb['r'], old_orb[key][i], label='in', linewidth=3, alpha=0.5)
                axs[irow][icol].plot(new_orb['r'], new_orb[key][i], label='out', linewidth=3, alpha=0.5)
                # add y = 0 line
                axs[irow][icol].axhline(y=0, color='k', linestyle='--')
                axs[irow][icol].set_title('%s-%d' % (key, i))
                axs[irow][icol].legend()
                icol += 1
                if icol == ncol:
                    icol = 0
                    irow += 1
        # set ylim according to the maximum value
        max_value = max(max(old_orb[key][i]), max(new_orb[key][i]))
        min_value = min(min(old_orb[key][i]), min(new_orb[key][i]))
        axs[irow][icol].set_ylim([min_value, max_value])

    plt.show()

def orbital_inverse(orbs_to_inverse: dict, inverse_flags = []) -> dict:
    """ inverse all orbitals if the first few values are negative """
    if len(inverse_flags) == 0:
        inverse_flags = [True for i in range(len(orbs_to_inverse.keys()))]
    for ik, key in enumerate(orbs_to_inverse.keys()):
        if key == 'r':
            continue
        else:
            for i in range(len(orbs_to_inverse[key])):
                sign = True # whether to reverse
                for ir in range(10):
                    if orbs_to_inverse[key][i][ir] > 0:
                        sign = False
                        break
                if sign and inverse_flags[ik]:
                    orbs_to_inverse[key][i] = -orbs_to_inverse[key][i]

    return orbs_to_inverse

def write_log(filename: str, n_jtilde: int, n_jtilde_remove: list, inverse: bool, inverse_flags: list):
    """ write log file """
    with open(filename, 'w') as f:
        for i in range(len(n_jtilde_remove)):
            f.write("number of j'l(qr) removed [%i] = %i\n" % (i, n_jtilde_remove[i]))
        for i in range(len(n_jtilde_remove)):
            f.write("number of j'l(qr) rest in the truncated basis [%i] = %i\n" % (i, n_jtilde - n_jtilde_remove[i]))
        if inverse:
            for i in range(len(inverse_flags)):
                f.write("orbital inverse on %i-th sublayer = %s\n" % (i, str(inverse_flags[i])))
    print("log file %s written" % filename)

def main(path: str, orbital_file: str, n_jtilde_remove: list, plot_mode: int, inverse: bool, inverse_flags: list, output: bool) -> None:
    """ main function """
    basic_information = filename_parse(orbital_file)
    print("basic information of the orbital file")
    print("filename: %s" % basic_information["filename"])
    print("element: %s" % basic_information["element"])
    print("functional: %s" % basic_information["functional"])
    print("bessel_nao_rcut: %s" % basic_information["bessel_nao_rcut"])
    print("bessel_nao_ecut: %s" % basic_information["bessel_nao_ecut"])
    print("valence split: %s" % basic_information["zeta"])
    n_jtilde = int(np.sqrt(basic_information['bessel_nao_ecut']) * basic_information['bessel_nao_rcut'] / np.pi)
    orbitals_osc_rm, orbitals, coefficients_osc_rm, _ = oscillation_remove(path, orbital_file, n_jtilde_remove)
    if inverse:
        orbitals_osc_rm = orbital_inverse(orbitals_osc_rm, inverse_flags)
    
    # ------------------- plot -------------------
    if plot_mode == 0:
        plot_orb_json(orbitals_osc_rm, 'x', -1)
        plot_orb_json(orbitals, 'x', -1)
    elif plot_mode == 1:
        compare_plot(orbitals, orbitals_osc_rm)

    # ------------------- write -------------------
    file_appendix = '_osci_rm'
    if inverse:
        file_appendix += '_inv'
    
    if output:
        new_orbital_filename = orbital_file.split('.')[0] + file_appendix + '.orb'
        coefficient_filename = orbital_file.split('.')[0] + file_appendix + '.c4'
        log_filename = orbital_file.split('.')[0] + file_appendix + '.log'
        dict_to_c4(filename = coefficient_filename, element = orbital_file[0:2], coefficients = coefficients_osc_rm)
        dict_to_orb(new_orbital_filename, orbital_file[0:2], 100, float(orbital_file.split("_")[2].replace("au", "")), orbitals_osc_rm)
        write_log(log_filename, n_jtilde, n_jtilde_remove, inverse, inverse_flags)

def interactive_solve_coefficients(path: str, orbital_file: str, linalg_solver = 'numpy') -> dict:
    """ remove high frequency oscillation from numerical atomic orbitals """
    
    numerical_orbitals = data_import(path, orbital_file) # the input numerical atomic orbitals    
    # structure:
    # "coefficients": {
    #   "s": {
    #        [...] # list of coefficients for each numerical atomic orbital of the 1st sublayer
    #        [...] # list of coefficients for each numerical atomic orbital of the 2nd sublayer
    #        ...
    #   },
    #   "p": {
    #        [...] # list of coefficients for each numerical atomic orbital of the 1st sublayer
    #        ...
    #   ...
    # }
    # "qs": {
    #   "s": [...],
    #   "p": [...],
    #   ...
    # }
    interactive_information = {
        "coefficients": {},
        "qs": {}
    }
    orbital_basic_info = filename_parse(orbital_file)
    bessel_nao_ecut = orbital_basic_info['bessel_nao_ecut']
    bessel_nao_rcut = orbital_basic_info['bessel_nao_rcut']

    r = np.array(numerical_orbitals['r'])

    # solve Ax = y equation, 
    # A is overlap between truncated spherical bessel functions, 
    # x is coefficient (of tsbf)-vector from numerical atomic orbital
    # y is overlap between numerical atomic orbital and tsbf
    # calculate A
    for key in numerical_orbitals.keys():
        if key == 'r':
            interactive_information["r"] = np.array(numerical_orbitals['r'])
        else:
            l = symbol_tol(key)
            # initialize _c_in and _c_out
            interactive_information["coefficients"][key] = []
            # generate q of spherical bessel functions for current l
            qs = generate_truncated_spherical_bessel(bessel_nao_ecut, bessel_nao_rcut, [l])
            interactive_information["qs"][key] = qs[0]
            # allocate memory for overlap matrix between truncated spherical bessel functions
            ovlp_jtildejtilde = np.zeros(shape=(len(qs[0]), len(qs[0]))) # A
            # loop over q...
            for iq in range(len(qs[0])):
                # generate one of the truncated spherical bessel functions
                jlq_i = jn(l, qs[0][iq]*r)
                for jq in range(len(qs[0])):
                    # generate another one of the truncated spherical bessel functions, but only if iq >= jq, because the overlap matrix is symmetric
                    if iq >= jq:
                        jlq_j = jn(l, qs[0][jq]*r)
                        ovlp_jtildejtilde[iq][jq] = simps(jlq_i*jlq_j*r**2, r)
                        # fill the other half of the overlap matrix
                        if iq != jq:
                            ovlp_jtildejtilde[jq][iq] = ovlp_jtildejtilde[iq][jq]
            # loop over numerical atomic orbitals for current l
            for inao in range(len(numerical_orbitals[key])):
                # for current numerical orbital...1d array
                orb = np.array(numerical_orbitals[key][inao])
                # calculate y
                # allocate memory for overlap between numerical atomic orbital and truncated spherical bessel functions
                ovlp_fjtilde = np.zeros(shape=(len(qs[0]), 1)) # y
                for iq in range(len(qs[0])):
                    # generate one of the truncated spherical bessel functions
                    jlq_i = jn(l, qs[0][iq]*r)
                    # calculate overlap between numerical atomic orbital and truncated spherical bessel functions
                    ovlp_fjtilde[iq] = simps(orb*jlq_i*r**2, r)
                # solve x
                if linalg_solver == 'numpy':
                    c = np.linalg.solve(ovlp_jtildejtilde, ovlp_fjtilde)
                elif linalg_solver == 'scipy':
                    c = solve(ovlp_jtildejtilde, ovlp_fjtilde)
                elif linalg_solver == 'scipy_cholesky':
                    cho, low = cho_factor(ovlp_jtildejtilde)
                    c = cho_solve((cho, low), ovlp_fjtilde)
                interactive_information["coefficients"][key].append(c.flatten())

    return interactive_information

def interactive_main(path: str, orbital_file: str) -> None:

    interactive_information = interactive_solve_coefficients(path, orbital_file)
    n_jtilde_remove = [0 for i in range(len(interactive_information["coefficients"].keys()))]

    r = interactive_information["r"]

    # create the sliders
    ncol = len(interactive_information["coefficients"]["s"])
    nplots = 0
    for key in interactive_information["coefficients"].keys():
        nplots += len(interactive_information["coefficients"][key])
    nrow = int(nplots/ncol) + 1
    # plot
    fig, axs = plt.subplots(nrow, ncol, sharex=True, sharey=False)
    plt.subplots_adjust(bottom=0.15)
    # set size
    fig.set_size_inches(20, 10)
    fig.suptitle('oscillation remove of orbital file %s' % orbital_file, fontsize=20)

    sliders = []
    for i in range(len(n_jtilde_remove)):
        # equally spaced
        ax = plt.axes([0.1 + 0.8*i/len(n_jtilde_remove), 0.05, 0.6/len(n_jtilde_remove), 0.03])
        slider = Slider(ax, l_tosymbol(i), 0, len(interactive_information["qs"]['s'])-1, valinit=n_jtilde_remove[i], valstep=1)

        sliders.append(slider)

    # define the update function
    def update(val):
        for i in range(len(n_jtilde_remove)):
            n_jtilde_remove[i] = int(sliders[i].val)
        # update the plot
        old_orbs = {}
        new_orbs = {}
        for key in interactive_information["coefficients"].keys():
            # keys are, 's', 'p', ...
            old_orbs[key] = []
            new_orbs[key] = []
            qs_l = interactive_information["qs"][key]
            for inao in range(len(interactive_information["coefficients"][key])):
                # inside interactive_information["coefficients"][key][inao], there are coefficients
                old_orb = []
                new_orb = []
                for ir in range(len(r)):
                    old_orb.append(0)
                    new_orb.append(0)
                    for iq in range(len(qs_l)):
                        old_orb[-1] += interactive_information["coefficients"][key][inao][iq] * jn(symbol_tol(key), qs_l[iq]*r[ir])
                    for iq in range(len(qs_l) - n_jtilde_remove[symbol_tol(key)]):
                        new_orb[-1] += interactive_information["coefficients"][key][inao][iq] * jn(symbol_tol(key), qs_l[iq]*r[ir])
                old_orbs[key].append(old_orb)
                new_orbs[key].append(new_orb)
                
        irow = 0
        icol = 0
        for key in new_orbs.keys():
            if key == 'r':
                continue
            else:
                for i in range(len(new_orbs[key])):
                    # set in as semi-transparent
                    axs[irow][icol].clear()
                    axs[irow][icol].plot(r, old_orbs[key][i], label='in', linewidth=3, alpha=0.5)
                    axs[irow][icol].plot(r, new_orbs[key][i], label='out', linewidth=3, alpha=0.5)
                    # add y = 0 line
                    axs[irow][icol].axhline(y=0, color='k', linestyle='--')
                    axs[irow][icol].set_title('%s-%d' % (key, i))
                    axs[irow][icol].legend()
                    # set ylim according to the maximum value
                    max_value = max(max(old_orbs[key][i]), max(new_orbs[key][i]))
                    min_value = min(min(old_orbs[key][i]), min(new_orbs[key][i]))
                    axs[irow][icol].set_ylim([min_value, max_value])
                    icol += 1
                    if icol == ncol:
                        icol = 0
                        irow += 1
        
        plt.draw()

    # connect the sliders to the update function
    for slider in sliders:
        slider.on_changed(update)

    # plot
    plt.show()

if __name__ == '__main__':

    path = './Er_sigma_0.1/68_Er_DZP'
    orbital_file = 'Er_gga_10au_100Ry_4s2p2d2f1g.orb'
    #interactive_main(path, orbital_file)

    n_jtilde_remove = [0, 0, 1, 1, 2]
    plot_mode = -1
    inverse = False
    inverse_flags = [True, False, True, True, False]
    output = True
    main(path, orbital_file, n_jtilde_remove, plot_mode, inverse, inverse_flags, output)
    