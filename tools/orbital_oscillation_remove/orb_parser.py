"""
Coding powered by Github.copilot
Author: Kirk0830
"""
import matplotlib.pyplot as plt
#from numerical_derivatives import numerical_first_derivative, numerical_second_derivative
from orb_io import orb_to_dict

def data_import(path, orbital_file):
    """Import data from .orb file."""
    if not path.endswith('/'):
        path += '/'
    orbital_file_path = path + orbital_file
    orbs = orb_to_dict(orbital_file_path)
    return orbs

def filename_parse(filename: str):
    words = filename.split('_')
    element = words[0]
    functional = words[1]
    bessel_nao_rcut = float(words[2][:-2])
    bessel_nao_ecut = float(words[3][:-2])
    zeta = words[4].split('.')[0]
    ns_zeta = []
    ls_zeta = []
    for letter in zeta:
        if letter.isdigit():
            ns_zeta.append(int(letter))
        else:
            if letter == 's':
                ls_zeta.append(0)
            elif letter == 'p':
                ls_zeta.append(1)
            elif letter == 'd':
                ls_zeta.append(2)
            elif letter == 'f':
                ls_zeta.append(3)
            elif letter == 'g':
                ls_zeta.append(4)
            else:
                RuntimeError("numerical atomic orbital file's name is not in standard format: [element]_[functional]_[rcut]au_[ecut]Ry_[zeta].orb")
    return {
        "filename": filename,
        "element": element,
        "functional": functional,
        "bessel_nao_rcut": bessel_nao_rcut,
        "bessel_nao_ecut": bessel_nao_ecut,
        "zeta": zeta,
        "ls": ls_zeta,
        "ns": ns_zeta
    }

def plot_orb_json(orb_json: dict, sublayer: str, index: int, save_fig: bool = False, fig_name: str = 'orbital.png'):
    
    if sublayer == 'x':
        for key in orb_json:
            if key == 'r':
                continue
            for i in range(len(orb_json[key])):
                plt.plot(orb_json['r'], orb_json[key][i], label='{}-th {} orbital'.format(i, key))
    else:
        if index == -1:
            for i in range(len(orb_json[sublayer])):
                plt.plot(orb_json['r'], orb_json[sublayer][i], label='{}-th {} orbital'.format(i, sublayer))
        else:
            plt.plot(orb_json['r'], orb_json[sublayer][index], label='{}-th {} orbital'.format(index, sublayer))

    # show y=0
    plt.plot([0, orb_json['r'][-1]], [0, 0], 'k--')
    plt.xlabel('r (a.u.)')
    plt.ylabel('radial numerical atomic orbital')
    plt.legend()
    if save_fig:
        plt.savefig(fig_name)
    plt.show()

def main(sublayer = 'x', index = -1, path = './', orbital_file = 'In_gga_6au_100Ry_3s3p3d2f.orb'):
    
    #path = './PD04-PBE-Numerical_Orbitals_v20231027-1/49_In_TZDP'
    orbs = data_import(
        path = path, 
        orbital_file = orbital_file
        )
    
    #print(filename_parse(orbital_file))
    if sublayer == 'x' or index == -1:
        for key in orbs:
            if key == 'r':
                continue
            for i in range(len(orbs[key])):
                plt.plot(orbs['r'], orbs[key][i], label='{}-th {} orbital'.format(i, key))
    else:
        plt.plot(orbs['r'], orbs[sublayer][index], label='{}-th {} orbital'.format(index, sublayer))

    # show y=0
    plt.plot([0, orbs['r'][-1]], [0, 0], 'k--')
    plt.xlabel('r (a.u.)')
    plt.ylabel('radial numerical atomic orbital')
    plt.title(orbital_file)
    plt.legend()
    plt.savefig('{}.png'.format(orbital_file[:-4]))
    plt.show()

if __name__ == '__main__':

    # to plot all orbitals, set sublayer = 'x', index = -1
    main(sublayer = 'x', 
         index = -1,
         path = 'D:/pseudopotential_test_workflow/abacus_pseudopotential_square/numerical_orbitals/resources/pslnc_031/49_In_DZP',
         orbital_file = 'In_gga_6au_100Ry_2s2p2d1f.orb'
         )
    # to plot a specific orbital, set sublayer = 's', 'p', 'd', 'f', 'g', index = 0, 1, 2, ...
    # main(sublayer = 's', index = 0)
    