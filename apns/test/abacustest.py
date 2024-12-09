'''
INTERFACE TO ABACUSTEST
Author: @kirk0830
Github repo: https://github.com/pxlxingliang/abacus-test (by @pxlxingliang)

Usage:
    call function `write_abacustest_param` to generate abacustest param.json contents

Example:
    write_abacustest_param(jobgroup_name='example_run',
                           bohrium_login={'username': '***@aisi.ac.cn', 'password': '****', 'project_id': '***'},
                           save_dir='~/abacus_test',
                           prepare={'abacus2qe': True, 'folders': ['example1', 'example2']},
                           predft={'ifrun': True, 'command': 'python3 prepare_something.py', 'shared_files': ['file1', 'file2']},
                           rundft=[{'ifrun': True, 'command': 'abacus --version', 'shared_files': ['file1', 'file2']}],
                           postdft={'ifrun': True, 'command': 'python3 post_something.py', 'shared_files': ['file1', 'file2']})
    can specify `export=True` to save the param.json file to the current directory, then will return absolute path of the file
    instead of contents of the file

'''

#ABACUS_IMAGE = 'registry.dp.tech/deepmodeling/abacus-intel:latest'
ABACUS_IMAGE = 'registry.dp.tech/dptech/abacus:3.8.1' # this is for submit large batch of jobs, but need to always update
QE_IMAGE = 'registry.dp.tech/dptech/prod-471/abacus-vasp-qe:20230116'
VASP_IMAGE = 'registry.dp.tech/dptech/prod-471/abacus-vasp-qe:20230116'
PYTHON_IMAGE = 'python:3.8'
ABACUS_COMMAND = 'ulimit -c 0; OMP_NUM_THREADS=1 mpirun -n 16 abacus | tee out.log' # update to avoid core dump

import os
import json
import time

def read_apns_inp(fname: str) -> dict:
    assert os.path.exists(fname), f'File not found: {fname}'
    with open(fname, 'r') as f:
        inp = json.load(f)
    return inp.get('abacustest', {})

def abacus_command(run_abacus: str = 'abacus', 
                   stdout: str = 'out.log', 
                   omp_num_threads: int = 1,
                   mpi_np: int = 16) -> str:
    # update to avoid massive core dump
    return f'ulimit -c 0; OMP_NUM_THREADS={omp_num_threads} mpirun -n {mpi_np} {run_abacus} | tee {stdout}'

def manual_submit(usr, 
                  pwd, 
                  projid, 
                  ncores, 
                  mem, 
                  jobs,
                  nmpi=None,
                  nomp=1):
    '''
    Manually submit jobs to Bohrium platform

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
    # because hyperthreading... the c32 just means 16 cores
    options = {'run_abacus': 'abacus', 'stdout': 'out.log', 'omp_num_threads': nomp, 'mpi_np': nmpi}
    
    run_dft = [{'ifrun': True, 'job_folders': jobs, 
                'command': abacus_command(**options),
                'ncores': ncores, 'memory': mem}]
    
    jobgroup = f'apns_{time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())}'
    param = write_abacustest_param(jobgroup_name=jobgroup, 
                                   bohrium_login={'bohrium.account': usr, 'bohrium.password': pwd, 'project_id': projid}, 
                                   rundft=run_dft)
    result_folder = submit(param)
    return result_folder

def auto_api(test_setting: dict, folders: list):
    
    auto_submit = test_setting['global'].get('auto_submit', False)
    
    ncores = test_setting.get('ncores', 32)
    jobgroup = f'apns_{time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())}'
    run_dft = [{'ifrun': True, 'job_folders': folders, 
                'command': abacus_command('abacus', 'out.log', 1, int(ncores/2)), # because hyperthreading... the c32 just means 16 cores
                'ncores': ncores, 'memory': test_setting.get('memory', 64)}]
    if auto_submit:
        param = write_abacustest_param(jobgroup_name=jobgroup, 
                                       bohrium_login=test_setting['credentials']['abacustest'], 
                                       rundft=run_dft)
        result_folder = submit(param)
        return result_folder
    # else:
    #     print('Job is not submitted, compress to one zip file instead')
    #     fjob = f'{jobgroup}.zip'
    #     # use shutil to do the compress
    #     shutil.make_archive(fjob, 'zip', folders)
    #     os.system('rm -rf {}'.format(' '.join(folders)))
    #     return None

def bohrium_config(**kwargs):
    '''Configure Bohrium account information

    Args:
        username (str): example: 'username@aisi.ac.cn'
        password (str): your Bohrium password
        project_id (int): your Bohrium project id, for more information, see https://bohrium.dp.tech/projects -> 'Project ID'

    Returns:
        bohrium settings: dict, involved in abacustest configuration
    
    Note:
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
    username = kwargs.get('bohrium.account', None)
    password = kwargs.get('bohrium.password', None)
    project_id = kwargs.get('project_id', None)
    assert username and password and project_id, 'Please provide Bohrium username, password, and project id.'
    return {
        'bohrium_username': username,
        'bohrium_password': password,
        'bohrium_project_id': project_id
    }

def bohrium_machine(ncores: int, memory: float, device: str, supplier: str):
    '''Configure Bohrium machine information
    
    Args:
        ncores (int): number of cores
        memory (float): memory size
        device (str): device type, example: 'cpu', 'gpu'
        supplier (str): supplier name, example: 'ali', 'para'

    Returns:
        bohrium machine settings: dict, involved in abacustest configuration
    '''
    return {
        'scass_type': '_'.join(['c'+str(ncores), 'm'+str(memory), device]),
        'job_type': 'container',
        'platform': supplier
    }

def prepare_dft(**kwargs):
    '''Generate 'prepare' section of abacustest configuration

    Returns:
        dict: prepare configuration
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
    keys, vals = abacus_default()
    for key in keys:
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

def setup_dft(**kwargs):
    '''Generate predft/rundft/postdft configuration according to given parameters

    Returns:
        dict: predft/rundft/postdft configuration
    '''
    # if kwargs is empty, return an empty dictionary
    if not kwargs:
        return {}
    def find_image(command: str):
        if "abacus" in command:
            return ABACUS_IMAGE
        elif "pw.x" in command:
            return QE_IMAGE
        elif "vasp" in command:
            return VASP_IMAGE
        elif "python" in command:
            return PYTHON_IMAGE
        else:
            return kwargs.get("image", PYTHON_IMAGE)
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
    image = find_image(command)
    example = kwargs.get('job_folders', [])
    print('Presently the key `outputs` is deprecated due to imcompleted implementation of the feature of abacustest.')
    machine = bohrium_machine(kwargs.get('ncores', 16), 
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

def write_abacustest_param(jobgroup_name: str, bohrium_login: dict, save_dir: str = '', prepare: dict = {},
                           predft: dict = {}, rundft: list = [], postdft: dict = {}, export: bool = False):
    '''Generate abacustest param.json contents

    Args:
        jobgroup_name (str): identifier for the jobgroup, not the lbg_jobgroup_id
        bohrium_login (dict): a dictionary should have `username`, `password`, `project_id` keys
        save_dir (str): path to save the job results
        predft (dict): preprocessing procedure, usually for setting up some parameters in batch
        rundft (list): a list of dictionaries, each dictionary represents a manner of running DFT calculation
        postdft (dict): postprocessing procedure, usually for analyzing the results

    Returns:
        dict: abacustest param.json contents
    '''
    save_dir = save_dir if len(save_dir) > 0 else f'abacustest-autosubmit-{time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())}'
    result = {
        'bohrium_group_name': jobgroup_name,
        'config': bohrium_config(**bohrium_login),
        'save_path': save_dir,
        'prepare': prepare_dft(**prepare),
        'pre_dft': setup_dft(**predft),
        'run_dft': [setup_dft(**dft) for dft in rundft],
        'post_dft': setup_dft(**postdft)
    }
    result = {k: v for k, v in result.items() if len(v) > 0}
    result.update({'compress': True}) # add compress: True to compress job files, may accelerate the upload process...

    if export:
        with open('param.json', 'w') as f:
            json.dump(result, f, indent=4)
        return os.path.abspath('/'.join([os.getcwd(), 'param.json']))
    return result

def submit(abacustest_param: dict) -> str:
    fparam = f'param-{time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())}.json'
    with open(fparam, 'w') as f:
        json.dump(abacustest_param, f, indent=4)
    
    folder = abacustest_param.get('save_path', 'result')
    flog = fparam.rsplit('.', 1)[0] + '.log'
    os.system(f'nohup abacustest submit -p {fparam} > {flog}&')
    print(f'Job submitted, log file is {flog}, results will be downloaded into {folder}')
    return folder

def abacus_default():
    from apns.analysis.postprocess.read_abacus_out import read_keyvals_frominput
    fthis = os.path.abspath(__file__)
    fabacus = fthis.replace('/abacustest.py', '/abacus_input_example')
    result = read_keyvals_frominput(fabacus).items()
    keys, vals = zip(*result)
    return keys, vals

####
# Convert APNS output collection to abacustest reuse model
####

####
# Reuse: 011_LcaoVSPw_eos
####
def _group_eos_pw_vs_lcao(out_dir: str):
    '''For reuse mode, pack up the target folder to the format that abacustest support.
    For eos test between ps and lcao, should do things following:
    For each system, a folder created, INPUT, KPT and STRU files are inside.
    A pporb.json should be created and specified for each element in the following way:
    ```json
    {
        'H': {'pp': 'H_ONCV_PBE-1.0.upf', 
              'orb': [//...
                     ]},
        'He': {'pp': 'He_ONCV_PBE-1.0.upf', 
               'orb': [//...
                      ]},
        // ... other elements
    }
    ```
    and all the pseudopotentials and orbitals should be put into the same folder.
    The APNS generated out_dir has flat structure. There are folders in out_dir, each folder contains
    INPUT, KPT, STRU and pseudopotentials, if needed, also with orbital files. An additional json
    file named `description.json` is used to quicky check the identity of test case.
    '''
    import os
    import json
    from apns.analysis.apns2_utils import cal_desc_diff
    folders = os.listdir(out_dir)
    print(f'ABACUSTEST: Found {len(folders)} folders in {out_dir}, group them to abacustest reuse mode accepted structure.')
    # test caeses with the identical desc['CellGenerator'], desc['DFTParamSet'], desc['Cell'] sections, and 
    # desc['AtomSpecies'][i]['pp'], are considered to be the same test case
    groups = [[folders[0]]]
    for folder in folders[1:]:
        grouped = False
        with open(f'{out_dir}/{folder}/description.json', 'r') as f:
            desc_i = json.load(f)
        for g in groups:
            with open(f'{out_dir}/{g[0]}/description.json', 'r') as f:
                desc_j = json.load(f)
            diff = cal_desc_diff(desc_i, desc_j)
            if list(diff.keys()) == ['AtomSpecies'] and all([list(v.keys()) == ['nao'] for v in diff['AtomSpecies']]):
                g.append(folder)
                grouped = True
                break
        if not grouped:
            groups.append([folder])
    print(f'ABACUSTEST: {len(groups)} groups are generated.')
    return groups

def _pack_eos_pw_vs_lcao(grouped, src, dst):
    '''take the grouped folders and pack them into the format that abacustest support.'''
    import os
    import shutil
    import json
    pp_orb = {}
    for i, folders in enumerate(grouped):
        os.makedirs(f'{dst}/group_{i}', exist_ok=True)
        # copy all *.upf/*.UPF and *.orb files to dst folder
        for folder in folders:
            folder = os.path.join(src, folder)
            files = [file for file in os.listdir(folder) if file not in ['INPUT', 'STRU', 'KPT', 'description.json']]
            # get the pp_orb information from description.json
            with open(f'{folder}/description.json', 'r') as f:
                desc = json.load(f)
            for label in dict.fromkeys(desc['Cell']['labels']):
                il = desc['Cell']['labels'].index(label)
                aspec = desc['AtomSpecies'][desc['Cell']['labels_kinds_map'][il]]
                pp_orb.setdefault(label, {'pp': os.path.basename(aspec['pp'])}).setdefault('orb', []).append(
                    os.path.basename(aspec['nao']))
                pp_orb[label]['orb'] = list(set(pp_orb[label]['orb'])) # remove duplicates
            for file in files:
                shutil.copy(f'{folder}/{file}', f'{dst}/')
        # copy INPUT, STRU and KPT to dst/group_i
        for file in ['INPUT', 'STRU', 'KPT']:
            shutil.copy(f'{folder}/{file}', f'{dst}/group_{i}/')
    with open(f'{dst}/pporb.json', 'w') as f:
        json.dump(pp_orb, f, indent=4)
    # remove the src
    shutil.rmtree(src)

    return pp_orb

def reuse_eos_pw_vs_lcao(out_dir: str, dst: str = 'abacustest_reuse_eos_pw_vs_lcao'):
    '''abacustest reuse mode: transform the APNS generated out_dir to the format that abacustest support.
    Args:
        out_dir (str): the folder contains the APNS generated results
        dst (str): the destination folder to save the transformed results
    Returns:
        pporb (dict): the pseudopotential and orbital information, in the format of {'element': {'pp': 'pp_file', 'orb': ['orb_file1', 'orb_file2', ...]}}
    '''
    grouped = _group_eos_pw_vs_lcao(out_dir)
    pporb = _pack_eos_pw_vs_lcao(grouped, out_dir, dst)
    return pporb

import unittest
class TestABACUSTest(unittest.TestCase):
    def test_bohrium_config(self):
        # empty case
        result = bohrium_config()
        self.assertIsInstance(result, dict)
        self.assertDictEqual(result, {})
        # routine case
        result = bohrium_config(username='__unittest__@aisi.ac.cn',
                                password='UnItTeSt',
                                project_id='123456')
        self.assertIsInstance(result, dict)
        self.assertDictEqual(result, {'bohrium_username': '__unittest__@aisi.ac.cn',
                                      'bohrium_password': 'UnItTeSt',
                                      'bohrium_project_id': '123456'})

    def test_bohrium_machine(self):
        result = bohrium_machine(ncores=16, memory=32, device='cpu', supplier='ali')
        self.assertIsInstance(result, dict)
        self.assertDictEqual(result, {'scass_type': 'c16_m32_cpu',
                                      'job_type': 'container',
                                      'platform': 'ali'})
    
    def test_prepare_dft(self):
        result = prepare_dft()
        self.assertIsInstance(result, dict)
        self.assertDictEqual(result, {})

        result = prepare_dft(abacus2qe=True, abacus2vasp=True, folders=['example1', 'example2'])
        self.assertIsInstance(result, dict)
        self.assertDictEqual(result, {'example_template': ['example1', 'example2'],
                                      'abacus2qe': True,
                                      'abacus2vasp': True})
        
        result = prepare_dft(ecutwfc=[30, 40, 50], pseudo_dir='__unittest__', orbital_dir='__unittest__',
                             pp_dict={'Li': 'Li.pbe-rrkjus.UPF'}, orb_dict={'Li': 'Li.pbe-rrkjus.UPF'})
        self.assertIsInstance(result, dict)
        self.assertDictEqual(result, {'mix_input': {'ecutwfc': [30, 40, 50]},
                                      'input_template': 'INPUT',
                                      'stru_template': 'STRU',
                                      'kpt_template': 'KPT',
                                      'pp_dict': {'Li': 'Li.pbe-rrkjus.UPF'},
                                      'orb_dict': {'Li': 'Li.pbe-rrkjus.UPF'},
                                      'pp_path': '__unittest__',
                                      'orb_path': '__unittest__'})

        result = prepare_dft(shared_files=['file1', 'file2'], dpks_descriptor='descriptor', potcar=['Li.pbe-rrkjus.UPF'])
        self.assertIsInstance(result, dict)
        self.assertDictEqual(result, {'extra_files': ['file1', 'file2'],
                                      'dpks_descriptor': 'descriptor',
                                      'potcar': ['Li.pbe-rrkjus.UPF']})
        
    def test_setup_dft(self):
        # empty case
        result = setup_dft()
        self.assertIsInstance(result, dict)
        self.assertDictEqual(result, {})
        # arbitrary case1
        result = setup_dft(ifrun=True, command='abacus --version', shared_files=['file1', 'file2'])
        self.assertIsInstance(result, dict)
        self.assertDictEqual(result, {'ifrun': True,
                                      'image': 'registry.dp.tech/deepmodeling/abacus-intel:latest',
                                      'example': [],
                                      'command': 'abacus --version',
                                      'extra_files': ['file1', 'file2']})
        # rundft case
        result = setup_dft(ncores=16, memory=32, sub_save_path='sub_save', njobs_node=2)
        self.assertIsInstance(result, dict)
        self.assertDictEqual(result, {'ifrun': False,
                                      'command': 'abacus --version',
                                      'extra_files': [],
                                      'image': 'registry.dp.tech/deepmodeling/abacus-intel:latest',
                                      'example': [],
                                      'group_size': 2,
                                      'sub_save_path': 'sub_save',
                                      'bohrium': {'scass_type': 'c16_m32_cpu',
                                                  'job_type': 'container',
                                                  'platform': 'ali'}})
        # postdft case
        result = setup_dft(metrics={'something inside': 'and its key'}, command='python3 __postdft__.py')
        self.assertIsInstance(result, dict)
        self.assertDictEqual(result, {'ifrun': False,
                                      'command': 'python3 __postdft__.py',
                                      'extra_files': [],
                                      'image': 'python:3.8',
                                      'example': [],
                                      'metrics': {'something inside': 'and its key'}})

    def test_write_abacustest_param(self):
        result = write_abacustest_param(jobgroup_name='unittest', 
                                        bohrium_login={'username': '__unittest__@aisi.ac.cn',
                                                        'password': 'UnItTeSt',
                                                        'project_id': '123456'},
                                        save_dir='__unittest__',
                                        prepare={},
                                        predft={},
                                        rundft=[],
                                        postdft={})
        self.assertIsInstance(result, dict)
        ref = {
            'bohrium_group_name': 'unittest',
            'config': {'bohrium_username': '__unittest__@aisi.ac.cn',
                          'bohrium_password': 'UnItTeSt',
                          'bohrium_project_id': '123456'},
            'save_path': '__unittest__'
        }
        self.assertDictEqual(result, ref)
        # prepare task, for example, converting jobs in example to qe input
        result = write_abacustest_param(jobgroup_name='', 
                                        bohrium_login={},
                                        save_dir='',
                                        prepare={'abacus2qe': True, 'folders': ['example1', 'example2']},
                                        predft={},
                                        rundft=[],
                                        postdft={})
        ref = {
            'prepare': {'example_template': ['example1', 'example2'],
                        'abacus2qe': True}
        }
        self.assertDictEqual(result['prepare'], ref['prepare'])
        # rundft task, for example, run abacus calculation

if __name__ == '__main__':
    #unittest.main(exit=False)
    import os, time
    os.chdir('/root/documents/simulation/abacus/EkinConvTest-JYlmax2-Fe')
    folders = [f for f in os.listdir() if not os.path.isfile(f)]
    manual_submit('_', '_', '28682', 32, 256, folders)
