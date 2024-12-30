      
'''
an individual interface to abacus-test Python package
'''
import os
import json
import time
import unittest
import argparse
import logging

# this script requires both the abacustest package and lbg, so try and catch
try:
    import abacustest
    import lbgcore
except ImportError:
    err = 'Please install abacustest and lbg packages first.\n'
    err += 'abacustest can be installed with following steps:\n'
    err += '1. git clone from official Github repo of @pxlxingliang:\n'
    err += '   git clone https://github.com/pxlxingliang/abacus-test.git\n'
    err += '2. cd into the cloned repo and install the package:\n'
    err += '   cd abacus-test\n'
    err += '   python3 setup.py install\n'
    err += '3. install lbg package:\n'
    err += '   pip install lbg -U\n'
    err += '4. configure your Bohrium account and password via lbg command.\n'
    err += '   For more information on this, refer to `lbg --help`\n'
    raise ImportError(err)

IMAGES = {'abacus': 'registry.dp.tech/dptech/abacus:3.8.1',
          'quantum-espresso': 'registry.dp.tech/dptech/prod-471/abacus-vasp-qe:20230116',
          'vasp': 'registry.dp.tech/dptech/prod-471/abacus-vasp-qe:20230116',
          'python': 'python:3.8'}

def _run_on_bohrium(program):
    '''
    due to the Bohrium enviroment, always need to add some commands
    to run the command expected
    '''
    return 'ulimit -c 0; ' + program + '| tee out.log'

def _parallel_command(ncore, maxmem, 
                      nmpi, nthreads = None):
    '''
    determine the way to run in parallel on Bohrium. The machine of
    Bohrium has super-threading (Ali), so the number of physical cores
    is the half of cores marked.
    '''
    if nmpi == 1:
        nthreads = ncore if nthreads is None else nthreads
    
    return f'OMP_NUM_THREADS={nthreads} mpirun -np {nmpi}'

def _redirect(stdout = None, stderr = None):
    '''
    redirect the output of the command to the specified files
    '''
    return f'> {stdout} 2> {stderr}' if stdout and stderr else ''

def _run_program(program, 
                 ncore, maxmem,
                 nmpi, nthreads = None,
                 stdout = None, stderr = None):
    '''
    '''
    return _run_on_bohrium(
        _parallel_command(ncore, maxmem, nmpi, nthreads) 
        + ' ' + program + _redirect(stdout, stderr))

def _bohrium_credential(**kwargs):
    '''
    Configure Bohrium account information

    Parameters
    ----------
    username : str
        example: 'username@aisi.ac.cn'

    password : str
        your Bohrium password

    project_id : int
        your Bohrium project id, for more information, see https://bohrium.dp.tech/projects -> 'Project ID'

    Returns
    -------
    bohrium settings : dict
        involved in abacustest configuration

    Note
    ----
    However, this way is highly not recommended, because the password is stored in plain text. For abacustest, 
    there is an alternative way in which the account and password can be read from two environment variables:
    ```
    export BOHRIUM_USERNAME=<bohrium-email>
    export BOHRIUM_PASSWORD=<bohrium-password>
    export BOHRIUM_PROJECT_ID=<bohrium-project-id>
    ```
    '''
    # if kwargs is empty, return an empty dictionary
    if not kwargs:
        return {}
    # if there are already set environment variables, return an empty dictionary
    if {'BOHRIUM_USERNAME', 'BOHRIUM_PASSWORD', 'BOHRIUM_PROJECT_ID'}.issubset(os.environ.keys()):
        return {} # abacustest will read from environment variables, dont worry
    
    print('Warning: the password is stored in plain text, please use environment variables instead!')
    username = kwargs.get('bohrium.account', None)
    password = kwargs.get('bohrium.password', None)
    project_id = kwargs.get('project_id', None)
    assert username and password and project_id, 'Please provide Bohrium username, password, and project id.'
    return {
        'bohrium_username': username,
        'bohrium_password': password,
        'bohrium_project_id': project_id
    }

def _bohrium_machine(ncores: int, memory: float, device: str, supplier: str):
    '''Configure Bohrium machine information
    
    Parameters
    ----------
    ncores : int
        number of cores

    memory : float
        memory size, in GB

    device : str
        device type, example: 'cpu', 'gpu'

    supplier : str
        supplier name, example: 'ali', 'para'

    Returns
    -------
    dict
        involved in abacustest configuration
    '''
    return {
        'scass_type': '_'.join(['c'+str(ncores), 'm'+str(memory), device]),
        'job_type': 'container',
        'platform': supplier
    }

def _prepare_dft(**kwargs):
    '''
    Generate 'prepare' section of abacustest configuration

    Returns
    -------
    dict
        prepare configuration
    '''
    # there are two ways to generate batch of input:
    # 1. with already-existed input files to overwrite
    # 2. generate input files from scratch
    #
    # for mode 1, a list of folders will be the value
    # of list 'example_template'
    # for mode 2, three keys 'input_template', 'stru_template',
    # and 'kpt_template' will/should be defined.
    result = {}
    for key in ABACUS_KEYS:
        if key in kwargs.keys() and not key in ['pseudo_dir', 'orbital_dir']:
            result.setdefault('mix_input', {})[key] = kwargs[key] if isinstance(kwargs[key], list) else [kwargs[key]]
    result.update({'example_template': kwargs.get('folders', [])})
    result.update(dict(zip(['input_template', 'stru_template', 'kpt_template'], 
                           [kwargs.get(key, key.upper()) for key in ['input', 'stru', 'kpt']]
                           ))) if result['example_template'] == [] else None
    result = {} if result.get('mix_input', {}) == {} and result.get('example_template', []) == [] else result

    result.update(dict(zip(['mix_kpt', 'mix_stru'], [[], []])))
    result['pp_dict'] = kwargs.get('pp_dict', {})
    result['orb_dict'] = kwargs.get('orb_dict', {})
    result['pp_path'] = kwargs.get('pseudo_dir', '')
    result['orb_path'] = kwargs.get('orbital_dir', '')
    result['dpks_descriptor'] = kwargs.get('dpks_descriptor', '')
    result['extra_files'] = kwargs.get('shared_files', [])
    result['abacus2qe'] = kwargs.get('abacus2qe', False)
    result['qe_setting'] = kwargs.get('qe_setting', {})
    result['abacus2vasp'] = kwargs.get('abacus2vasp', False)
    result['vasp_setting'] = kwargs.get('vasp_setting', {})
    result['potcar'] = kwargs.get('potcar', [])
    # seperate key-value pairs, for value who has len() attribute, save to _container, otherwise save to _rest
    _container = {key: value for key, value in result.items() if hasattr(value, '__len__')}
    _rest = {key: value for key, value in result.items() if not key in _container.keys()}
    _container = {key: value for key, value in _container.items() if len(value) > 0}
    _rest = {key: value for key, value in _rest.items() if value}
    return {**_container, **_rest}

def _setup_dft(**kwargs):
    '''
    Generate predft/rundft/postdft configuration according to given parameters

    Returns
    -------
    dict
        predft/rundft/postdft configuration
    '''
    # if kwargs is empty, return an empty dictionary
    if not kwargs:
        return {}
    def _fimg(command):
        cand = ['abacus', 'pw.x', 'vasp', 'python']
        for c in cand:
            if c in command:
                return IMAGES[c]
        return IMAGES['abacus']
    # to distinguish between predft/rundft and postdft
    metrics = kwargs.get('metrics', None)
    postdft = True if metrics else False
    # rundft specific
    sub_save_path = kwargs.get('sub_save_path', None)
    rundft = True if sub_save_path else False
    # mutual keys
    switch = kwargs.get('ifrun', False)
    command = kwargs.get('command', 'abacus --version')
    shared_files = kwargs.get('shared_files', [])
    image = _fimg(command)
    if 'image' in kwargs.keys():
        print(f'User setting image {kwargs["image"]}, default image {image} will be overwritten.')
        image = kwargs['image']
    
    example = kwargs.get('job_folders', [])
    print('Presently the key `outputs` is deprecated due to imcompleted implementation of the feature of abacustest.')
    machine = _bohrium_machine(kwargs.get('ncores', 16), 
                              kwargs.get('memory', 32), 
                              'cpu', 
                              'ali') if 'ncores' in kwargs.keys() or 'memory' in kwargs.keys() else None
    # rundft specific
    group_size = kwargs.get('njobs_node', 1) # rundft specific

    result = {
        'ifrun': switch,
        'command': command,
        'extra_files': shared_files,
        'image': image,
        'example': example
    }
    if machine is not None:
        result['bohrium'] = machine
    # on demand is always on, avoiding to be killed
    result['bohrium']['on_demand'] = 1

    if rundft:
        result['group_size'] = group_size
        result['sub_save_path'] = sub_save_path
    if postdft:
        result['metrics'] = metrics
    return result

def _write_abacustest_param(jobgroup_name, 
                            bohrium_login, 
                            save_dir = None, 
                            prepare = None,
                            predft = None, 
                            rundft = None, 
                            postdft = None, 
                            export: bool = False):
    '''
    Generate abacustest param.json contents

    Parameters
    ----------
    jobgroup_name : str
        Identifier for the jobgroup, not the lbg_jobgroup_id

    bohrium_login : dict
        A dictionary should have `username`, `password`, `project_id` keys

    save_dir : str
        Path to save the job results

    prepare : dict
        Preparation procedure, usually for setting up some parameters in batch

    predft : dict
        Preprocessing procedure, usually for setting up some parameters in batch

    rundft : list
        A list of dictionaries, each dictionary represents a manner of running DFT calculation

    postdft : dict
        Postprocessing procedure, usually for analyzing the results

    export : bool
        Whether to export the configuration to a JSON file

    Returns
    -------
    dict
        abacustest param.json contents
    '''
    save_dir = save_dir if save_dir \
        else f'abacustest-autosubmit-{time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())}'
    prepare = prepare if prepare else {}
    predft = predft if predft else {}
    rundft = rundft if rundft else []
    postdft = postdft if postdft else {}
    result = {
        'bohrium_group_name': jobgroup_name,
        'config': _bohrium_credential(**bohrium_login),
        'save_path': save_dir,
        'prepare': _prepare_dft(**prepare),
        'pre_dft': _setup_dft(**predft),
        'run_dft': [_setup_dft(**dft) for dft in rundft],
        'post_dft': _setup_dft(**postdft)
    }
    result = {k: v for k, v in result.items() if len(v) > 0}

    # add compress: True to compress job files, may accelerate the upload process...
    result.update({'compress': True}) 

    if export:
        with open('param.json', 'w') as f:
            json.dump(result, f, indent=4)
        return os.path.abspath('/'.join([os.getcwd(), 'param.json']))
    return result

def _abacustest_submit_kernel_impl(abacustest_param: dict) -> str:
    fparam = f'param-{time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())}.json'
    with open(fparam, 'w') as f:
        json.dump(abacustest_param, f, indent=4)
    
    folder = abacustest_param.get('save_path', 'result')
    flog = fparam.rsplit('.', 1)[0] + '.log'
    os.system(f'nohup abacustest submit -p {fparam} > {flog}&')
    print(f'Job submitted, log file is {flog}, results will be downloaded into {folder}')
    return folder

def abacus(usr, 
           pwd, 
           projid, 
           ncores, 
           mem, 
           jobs,
           nmpi=None,
           nomp=1):
    '''
    Manually _abacustest_submit_kernel_impl abacus jobs prepared in one folder to Bohrium platform

    Parameters
    ----------
    usr : str
        Bohrium username

    pwd : str
        Bohrium password. This is dangerous, please use environment variables instead.

    projid : str
        Bohrium project id

    ncores : int
        Number of cores of machine to run the job

    mem : float
        Memory size of the machine

    jobs : list
        List of folders, each folder contains the input files for one job

    nomp : int, optional
        Number of OpenMP threads, by default 1

    Returns
    -------
    str
        Result folder
    '''
    nmpi = nmpi if nmpi else int(ncores/2)    
    run_dft = [{'ifrun': True, 'job_folders': jobs, 
                'command': _run_program('abacus', ncores, mem, nmpi, nomp, 'out.log', 'err.log'),
                'ncores': ncores, 'memory': mem}]
    
    jobgroup = f'apns_{time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())}'
    fparam = _write_abacustest_param(jobgroup_name=jobgroup, 
                            bohrium_login={'bohrium.account': usr, 
                                           'bohrium.password': pwd, 
                                           'project_id': projid}, 
                            rundft=run_dft)
    result_folder = _abacustest_submit_kernel_impl(fparam)
    return result_folder

ABACUS_KEYS = {'suffix', 'latname', 'stru_file', 'kpoint_file', 'pseudo_dir',
               'orbital_dir', 'pseudo_rcut', 'pseudo_mesh', 'lmaxmax', 
               'dft_functional', 'xc_temperature', 'calculation', 
               'esolver_type', 'ntype', 'nspin', 'kspacing', 'min_dist_coef',
               'nbands', 'nbands_sto', 'nbands_istate', 'symmetry', 
               'init_vel', 'symmetry_prec', 'symmetry_autoclose', 'nelec', 
               'nelec_delta', 'out_mul', 'noncolin', 'lspinorb', 'kpar', 
               'bndpar', 'out_freq_elec', 'dft_plus_dmft', 'rpa', 'printe', 
               'mem_saver', 'diago_proc', 'nbspline', 'wannier_card', 
               'soc_lambda', 'cal_force', 'out_freq_ion', 'device', 
               'precision', 'ecutwfc', 'ecutrho', 'erf_ecut', 'erf_height', 
               'erf_sigma', 'fft_mode', 'pw_diag_thr', 'scf_thr', 
               'scf_thr_type', 'init_wfc', 'init_chg', 'chg_extrap', 
               'out_chg', 'out_pot', 'out_wfc_pw', 'out_wfc_r', 'out_dos', 
               'out_band', 'out_proj_band', 'restart_save', 'restart_load', 
               'read_file_dir', 'nx', 'ny', 'nz', 'ndx', 'ndy', 'ndz', 
               'cell_factor', 'pw_seed', 'method_sto', 'npart_sto', 
               'nche_sto', 'emin_sto', 'emax_sto', 'seed_sto', 'initsto_ecut', 
               'initsto_freq', 'cal_cond', 'cond_che_thr', 'cond_dw', 
               'cond_wcut', 'cond_dt', 'cond_dtbatch', 'cond_smear', 'cond_fwhm', 
               'cond_nonlocal', 'ks_solver', 'scf_nmax', 'relax_nmax', 
               'out_stru', 'force_thr', 'force_thr_ev', 'force_thr_ev2', 
               'relax_cg_thr', 'stress_thr', 'press1', 'press2', 'press3', 
               'relax_bfgs_w1', 'relax_bfgs_w2', 'relax_bfgs_rmax', 
               'relax_bfgs_rmin', 'relax_bfgs_init', 'cal_stress', 'fixed_axes', 
               'fixed_ibrav', 'fixed_atoms', 'relax_method', 'relax_new', 
               'relax_scale_force', 'out_level', 'out_dm', 'out_bandgap', 
               'use_paw', 'deepks_out_labels', 'deepks_scf', 'deepks_bandgap', 
               'deepks_out_unittest', 'deepks_model', 'basis_type', 'nb2d', 
               'gamma_only', 'search_radius', 'search_pbc', 'lcao_ecut', 
               'lcao_dk', 'lcao_dr', 'lcao_rmax', 'out_mat_hs', 'out_mat_hs2', 
               'out_mat_dh', 'out_mat_xc', 'out_interval', 'out_app_flag', 
               'out_mat_t', 'out_element_info', 'out_mat_r', 'out_wfc_lcao', 
               'bx', 'by', 'bz', 'smearing_method', 'smearing_sigma', 
               'mixing_type', 'mixing_beta', 'mixing_ndim', 'mixing_restart', 
               'mixing_gg0', 'mixing_beta_mag', 'mixing_gg0_mag', 
               'mixing_gg0_min', 'mixing_angle', 'mixing_tau', 'mixing_dftu', 
               'mixing_dmr', 'dos_emin_ev', 'dos_emax_ev', 'dos_edelta_ev', 
               'dos_scale', 'dos_sigma', 'dos_nche', 'md_type', 'md_thermostat', 
               'md_nstep', 'md_dt', 'md_tchain', 'md_tfirst', 'md_tlast', 
               'md_dumpfreq', 'md_restartfreq', 'md_seed', 'md_prec_level', 
               'ref_cell_factor', 'md_restart', 'lj_rcut', 'lj_epsilon', 
               'lj_sigma', 'pot_file', 'msst_direction', 'msst_vel', 'msst_vis', 
               'msst_tscale', 'msst_qmass', 'md_tfreq', 'md_damp', 'md_nraise', 
               'cal_syns', 'dmax', 'md_tolerance', 'md_pmode', 'md_pcouple', 
               'md_pchain', 'md_pfirst', 'md_plast', 'md_pfreq', 'dump_force', 
               'dump_vel', 'dump_virial', 'efield_flag', 'dip_cor_flag', 
               'efield_dir', 'efield_pos_max', 'efield_pos_dec', 'efield_amp', 
               'gate_flag', 'zgate', 'relax', 'block', 'block_down', 'block_up', 
               'block_height', 'out_alllog', 'nurse', 'colour', 't_in_h', 
               'vl_in_h', 'vnl_in_h', 'vh_in_h', 'vion_in_h', 'test_force', 
               'test_stress', 'test_skip_ewald', 'vdw_method', 'vdw_s6', 'vdw_s8', 
               'vdw_a1', 'vdw_a2', 'vdw_d', 'vdw_abc', 'vdw_C6_file', 
               'vdw_C6_unit', 'vdw_R0_file', 'vdw_R0_unit', 'vdw_cutoff_type', 
               'vdw_cutoff_radius', 'vdw_radius_unit', 'vdw_cn_thr', 
               'vdw_cn_thr_unit', 'vdw_cutoff_period', 'exx_hybrid_alpha', 
               'exx_hse_omega', 'exx_separate_loop', 'exx_hybrid_step', 
               'exx_mixing_beta', 'exx_lambda', 'exx_real_number', 
               'exx_pca_threshold', 'exx_c_threshold', 'exx_v_threshold', 
               'exx_dm_threshold', 'exx_cauchy_threshold', 'exx_c_grad_threshold', 
               'exx_v_grad_threshold', 'exx_cauchy_force_threshold', 
               'exx_cauchy_stress_threshold', 'exx_ccp_rmesh_times', 
               'exx_opt_orb_lmax', 'exx_opt_orb_ecut', 'exx_opt_orb_tolerence', 
               'td_force_dt', 'td_vext', 'td_vext_dire', 'out_dipole', 
               'out_efield', 'out_current', 'ocp', 'ocp_set', 'berry_phase', 
               'gdir', 'towannier90', 'nnkpfile', 'wannier_spin', 
               'wannier_method', 'out_wannier_mmn', 'out_wannier_amn', 
               'out_wannier_unk', 'out_wannier_eig', 'out_wannier_wvfn_formatted', 
               'imp_sol', 'eb_k', 'tau', 'sigma_k', 'nc_k', 'of_kinetic', 
               'of_method', 'of_conv', 'of_tole', 'of_tolp', 'of_tf_weight', 
               'of_vw_weight', 'of_wt_alpha', 'of_wt_beta', 'of_wt_rho0', 
               'of_hold_rho0', 'of_lkt_a', 'of_full_pw', 'of_full_pw_dim', 
               'of_read_kernel', 'of_kernel_file', 'dft_plus_u', 'yukawa_lambda', 
               'yukawa_potential', 'omc', 'hubbard_u', 'orbital_corr', 
               'bessel_nao_ecut', 'bessel_nao_tolerence', 'bessel_nao_rcut', 
               'bessel_nao_smooth', 'bessel_nao_sigma', 'bessel_descriptor_lmax', 
               'bessel_descriptor_ecut', 'bessel_descriptor_tolerence', 
               'bessel_descriptor_rcut', 'bessel_descriptor_smooth', 
               'bessel_descriptor_sigma', 'sc_mag_switch', 'decay_grad_switch', 
               'sc_thr', 'nsc', 'nsc_min', 'sc_scf_nmin', 'alpha_trial', 'sccut', 
               'sc_file', 'qo_switch', 'qo_basis', 'qo_thr'}

def init():
    '''
    read the input from command line
    '''
    parser = argparse.ArgumentParser(description='Submit jobs to Bohrium platform')
    parser.add_argument('--folder', '-f', default=None, help='Folder containing input files')
    parser.add_argument('--image', '-i', default=None, help='Image to run the job')
    parser.add_argument('--command', '-c', default=None, help='Command to run the job')
    
    args = parser.parse_args()
    image = IMAGES[args.image] if args.image is not None else None

    return args.folder, image, args.command

def general(image, command, machine, usr, pwd, projid, jobroot):
    '''
    '''
    cwd = os.getcwd()
    os.chdir(jobroot)
    jobs = [f for f in os.listdir() if not os.path.isfile(f)]
    print(f'Check in {jobroot}, {len(jobs)} folders found. Will attempt to submit all.')

    run_dft = [{'ifrun': True, 'job_folders': jobs, 
                'command': command,
                'ncores': machine['ncores'], 'memory': machine['memory'],
                'image': image}]
    
    jobgroup = f'apns_{time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())}'
    fparam = _write_abacustest_param(jobgroup_name=jobgroup,
                                     bohrium_login={'bohrium.account': usr, 
                                                    'bohrium.password': pwd, 
                                                    'project_id': projid},
                                     rundft=run_dft)
    result_folder = _abacustest_submit_kernel_impl(fparam)
    os.chdir(cwd)

    return result_folder
    

class TestIndividualAbacustest(unittest.TestCase):
    def test_nothing(self):
        self.assertTrue(True)

if __name__ == '__main__':
    '''
    execute the script to _abacustest_submit_kernel_impl jobs to Bohrium platform

    Usage
    -----
    it is HIGH RECOMMENDED to use environment variables to store Bohrium account information
    , but please do not export your envars ...

    specify the value of folder in which there are folders containing input files
    as the value of variable `jobdir` in the following.
    '''
    unittest.main(exit=False)

    elem = 'Se'
    item = 'JYLmaxRcutJointConv' # JYEkinConv or JYLmaxRcutJointConv

    folder_specified_here = f'{item}Test-{elem}'

    image_specified_here = 'registry.dp.tech/dptech/abacus:3.8.4' # 3.8.3 is the version that FFT is fixed

    command_specified_here  = 'ulimit -c 0; '
    command_specified_here += f'python3 {item}TestDriver.py -i driver.json | tee out.log'

    jobdir, imag, cmd = init()

    jobdir = folder_specified_here if jobdir is None else jobdir
    imag = image_specified_here if imag is None else imag
    cmd = command_specified_here if cmd is None else cmd

    # now this script can run not only abacus, so the entry here becomes much more
    # general than ever before

    _ = general(imag, cmd, {'ncores': 32, 'memory': 128}, 
                None, 
                None, 
                28682, 
                jobdir) # the root folder of jobs
    

