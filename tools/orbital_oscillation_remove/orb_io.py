from tools.orbital_oscillation_remove.database import symbol_tol

def orb_to_dict(filename):
    """Reads an ABACUS .orb file and returns a dictionary with the data."""
    mesh = 0
    line = ''
    end_sign = False
    f = open(filename, 'r')

    element = ''
    energy_cutoff = 0
    cutoff_radius = 0
    lmax = -1
    number_of_orbitals = {}

    # realspace mapping
    mesh = 0
    dr = 0

    # result dict
    result = {}

    while not end_sign:

        line = f.readline()
        words = line.strip().split()
        if len(words) == 0:
            continue

        if line.startswith('Element'):
            element = words[1]
        if line.startswith('Energy Cutoff'):
            energy_cutoff = float(words[-1])
        if line.startswith('Radius Cutoff'):
            cutoff_radius = float(words[-1])
        if line.startswith('Lmax'):
            lmax = int(words[-1])
        if line.startswith('Number of Sorbital'):
            number_of_orbitals['s'] = int(words[-1])
            result['s'] = []
        if line.startswith('Number of Porbital'):
            number_of_orbitals['p'] = int(words[-1])
            result['p'] = []
        if line.startswith('Number of Dorbital'):
            number_of_orbitals['d'] = int(words[-1])
            result['d'] = []
        if line.startswith('Number of Forbital'):
            number_of_orbitals['f'] = int(words[-1])
            result['f'] = []
        if line.startswith('Number of Gorbital'):
            number_of_orbitals['g'] = int(words[-1])
            result['g'] = []
        if line.startswith('Mesh'):
            mesh = int(words[-1])
        if line.startswith('dr'):
            end_sign = True
            dr = float(words[-1])

    total_number_of_orbitals = 0
    for key in number_of_orbitals:
        total_number_of_orbitals += number_of_orbitals[key]

    orbital_sign = ['s', 'p', 'd', 'f', 'g']

    # assert section
    if element == '':
        assert False, 'element not found in .orb file'
    if energy_cutoff == 0:
        assert False, 'energy_cutoff not found in .orb file'
    if cutoff_radius == 0:
        assert False, 'cutoff_radius not found in .orb file'
    if lmax == -1:
        assert False, 'lmax not found in .orb file'
    if len(number_of_orbitals) == 0:
        assert False, 'number of orbitals data not found in .orb file'
    if mesh == 0:
        assert False, 'mesh not found in .orb file'
    if dr == 0:
        assert False, 'dr not found in .orb file'
    
    # print basic information
    print('[orb_io: orb_to_dict] Print basic information')
    print('|-element: {}'.format(element))
    print('|-energy_cutoff: {} Ry'.format(energy_cutoff))
    print('|-cutoff_radius: {} a.u.'.format(cutoff_radius))
    print('|-lmax: {}'.format(lmax))
    for key in number_of_orbitals:
        print('|--number of {} orbitals: {}'.format(key, number_of_orbitals[key]))
    print('|-real space mapping information:')
    print('|-mesh: {}'.format(mesh))
    print('|-dr: {} a.u.'.format(dr))

    # generate radius mesh
    result['r'] = [i * dr for i in range(mesh)]

    # read the orbital data
    end_sign = False
    read_sign = False
    number_of_mesh_read = 0
    mesh_data = []
    l = -1
    number_of_orbitals_read = 0
    print('[orb_io: orb_to_dict] Continue reading the orbital data...')

    while not end_sign:
        line = f.readline()
        words = line.strip().split()

        if words[0].startswith('Type'):
            line = f.readline()
            words = line.strip().split()
            l = int(words[1])
            index_zeta = int(words[2])

            read_sign = True
            print('|-Now reading the {}-th orbital of orbital type {}'.format(index_zeta, orbital_sign[l]))
            number_of_mesh_read = 0
            mesh_data = []
            continue

        if read_sign:
            for word in words:
                mesh_data.append(float(word))
                number_of_mesh_read += 1
                if number_of_mesh_read == mesh:
                    read_sign = False
                    result[orbital_sign[l]].append(mesh_data)
                    number_of_orbitals_read += 1
                    break

        if number_of_orbitals_read == total_number_of_orbitals:
            end_sign = True

    f.close()

    return result

def dict_to_orb(filename: str, 
                element: str, 
                bessel_nao_ecut: float, 
                bessel_nao_rcut: float, 
                numerical_orbitals: dict):
    lmax = 0
    for key in numerical_orbitals.keys():
        if key == 'r':
            continue
        else:
            lmax = max(lmax, symbol_tol(key))
    with open(filename, 'w') as f:
        f.writelines('-'*75+'\n')
        f.writelines('Element'+ ' '*21 +'{}\n'.format(element))
        f.writelines('Energy Cutoff(Ry)'+ ' '*11 +'{}\n'.format(bessel_nao_ecut))
        f.writelines('Radius Cutoff(a.u.)'+ ' '*9 +'{}\n'.format(bessel_nao_rcut))
        f.writelines('Lmax'+ ' '*24 +'{}\n'.format(lmax))
        try:
            f.writelines('Number of Sorbital-->'+' '*7+str(len(numerical_orbitals['s']))+'\n')
        except KeyError:
            f.writelines('Number of Sorbital-->'+' '*7+str(0)+'\n')
        try:
            f.writelines('Number of Porbital-->'+' '*7+str(len(numerical_orbitals['p']))+'\n')
        except KeyError:
            f.writelines('Number of Porbital-->'+' '*7+str(0)+'\n')
        try:
            f.writelines('Number of Dorbital-->'+' '*7+str(len(numerical_orbitals['d']))+'\n')
        except KeyError:
            f.writelines('Number of Dorbital-->'+' '*7+str(0)+'\n')
        try:
            f.writelines('Number of Forbital-->'+' '*7+str(len(numerical_orbitals['f']))+'\n')
        except KeyError:
            f.writelines('Number of Forbital-->'+' '*7+str(0)+'\n')
        try:
            f.writelines('Number of Gorbital-->'+' '*7+str(len(numerical_orbitals['g']))+'\n')
        except KeyError:
            f.writelines('Number of Gorbital-->'+' '*7+str(0)+'\n')
        f.writelines('-'*75+'\n')
        f.writelines('SUMMARY  END\n\n')
        f.writelines('Mesh'+' '*24+str(len(numerical_orbitals['r']))+'\n')
        f.writelines('dr'+' '*26+'{}\n'.format(numerical_orbitals['r'][1]-numerical_orbitals['r'][0]))
        for key in numerical_orbitals.keys():
            if key == 'r':
                continue
            else:
                for izeta in range(len(numerical_orbitals[key])):
                    f.writelines('%20s%20s%20s\n'%('Type', 'L', 'N'))
                    f.writelines('%20i%20i%20i\n'%(0, symbol_tol(key), izeta))
                    for ir in range(len(numerical_orbitals[key][izeta])):
                        f.writelines('%15.14e  '%numerical_orbitals[key][izeta][ir])
                        if (ir+1)%4 == 0:
                            f.writelines('\n')
                    f.writelines('\n')

def dict_to_c4(filename: str,
               element: str,
               coefficients: dict,
               spillage = 0.0):
    """Write the coefficients of the truncated spherical bessel functions to a file."""
    """file format:
<Coefficient>
	 2 Total number of radial orbitals.
	Type	L	Zeta-Orbital
	  As 	0	    1
	  -0.71046261596632
	  -0.87150454963187
	  -0.23646387299581
	   0.28776919206825
	   0.47540091427994
	   0.43982444988140
	   0.30578481766673
	   0.17343492774223
	   0.06930964764613
	   0.01996799336757
	  -0.00631726367699
	  -0.00018450887673
	  -0.00411007710240
	   0.00814400724524
	  -0.00347818500548
	   0.00651426857901
	  -0.01001816829606
	   0.00770750565032
	  -0.02051897776865
	Type	L	Zeta-Orbital
	  As 	1	    1
	   0.80052977287852
	   0.82508500983994
	   0.50929769725993
	   0.15759602831618
	  -0.00486605041634
	  -0.10009368686188
	  -0.07929118332438
	  -0.07919771572262
	  -0.02126185028817
	  -0.02902640359012
	   0.01541573671325
	  -0.01709928430423
	   0.02232237296791
	  -0.02319736484789
	   0.02520343536762
	  -0.03016309406565
	   0.03356008982830
	  -0.04074883882235
	   0.16284943964046
</Coefficient>
<Mkb>
Left spillage = 2.6553225276e-02
</Mkb>
    """
    with open(filename, 'w') as f:
        f.writelines('<Coefficient>\n')
        # total number of radial orbitals
        total_number_of_radial_orbitals = 0
        for key in coefficients.keys():
            total_number_of_radial_orbitals += len(coefficients[key])
        f.writelines('    %2i Total number of radial orbitals.\n'%total_number_of_radial_orbitals)
        for key in coefficients.keys(): # keys are s, p, d, ...
            l = symbol_tol(key)
            for inao in range(len(coefficients[key])):
                f.writelines('%8s%8s%16s\n'%('Type', 'L', 'Zeta-Orbital'))
                f.writelines('%8s%8i%9i\n'%(str(element), l, inao+1))
                for iq in range(len(coefficients[key][inao])):
                    f.writelines('%22.14f\n'%coefficients[key][inao][iq])
        f.writelines('</Coefficient>\n')
        f.writelines('<Mkb>\n')
        f.writelines('Left spillage = %14.10e\n'%spillage)
        f.writelines('</Mkb>\n')
if __name__ == "__main__":

    read_orbs = orb_to_dict('In_gga_7au_100Ry_2s2p2d1f.orb')
    dict_to_orb('In_gga_7au_100Ry_2s2p2d1f_new.orb', 'In', 100, 7, read_orbs)