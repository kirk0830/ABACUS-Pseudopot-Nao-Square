import apns.module_workflow.identifier as amwi
import apns.module_software.abacus.generation as amsag
import apns.module_software.qespresso.generation as amsqg
import apns.module_structure.basic as amsb
import os

def iterate(software: str = "abacus", # which software? abacus or qespresso
            systems: list = [],
            pseudopot_nao_settings: list = [],
            calculation_settings: list = [],
            characteristic_lengths: list = [],
            valid_pseudopotentials: dict = {},
            valid_numerical_orbitals: dict = {},
            pspot_archive: dict = {},
            nao_archive: dict = {},
            test_mode: bool = True):
    """iterate over all possible combinations of input parameters and generate folders
    
    Args:
        `systems (list, optional)`: systems to iterate = [system_mpid, ...]. Defaults to empty.
        
        `pseudopot_nao_settings (list, optional)`: pseudopotential and numerical orbital settings to iterate = 
        >>> [{
            "pseudopotential": [pseudopotential_identifier, ...],
            "numerical_orbital": [numerical_orbital_identifier, ...]
        }, ...]. 

        Defaults to empty.
        
        `calculation_settings (list, optional)`: calculation settings to iterate = 
        >>> [
            {"basis_type": "pw", "dft_functional": "pbe", ...}, ...
        ]. 
        
        Defaults to empty.
        
        `characteristic_lengths (list, optional)`: characteristic_lengths settings to iterate = 
        >>> [0, 0.01, ...]. # for crystal calculating EOS
        >>> [1.5, 2.0, ...]. # for molecule specifying bond length
        
        Defaults to empty.
        
        `valid_pseudopotentials (dict, optional)`: 
        >>> valid pseudopotentials = {
            [element]: {
                [pseudopotential_identifier]: {
                    "file": [pseudopotential_file_name],
                    "kind" ..., "version": ..., "appendix": ...
                }, ...
            }, ...
        }. 
        
        Defaults to empty.

        `valid_numerical_orbitals (dict, optional)`: 
        >>> valid numerical orbitals = {
            [element]: {
                [[numerical_orbital_identifier]@[pseudopotential_identifier]]: {
                    "file": [numerical_orbital_file_name],
                    "kind" ..., "version": ..., "appendix": ...
                }, ...
            }, ...
        }. 
        
        Defaults to empty.

        `pspot_archive (dict, optional)`: pseudopotential archive = [identifier]->[folder address] map. Defaults to empty.
        `nao_archive (dict, optional)`: numerical orbital archive = [identifier]->[folder address] map. Defaults to empty.
        
    Raises:
        
    Returns:
        `list`: folders generated
    """
    folders = []
    for _system in systems: # iterate over structures
        for system_pseudopot_nao_setting in pseudopot_nao_settings: # iterate over pseudopotential and numerical orbital settings
            for calculation_setting in calculation_settings: # iterate over calculation settings
                for characteristic_length in characteristic_lengths: # iterate over cell scaling
                    # make folder
                    if "numerical_orbital" not in system_pseudopot_nao_setting.keys():
                        system_pseudopot_nao_setting["numerical_orbital"] = []
                    folder = amwi._folder_(_system, 
                                           amwi.pseudopot_nao(system_pseudopot_nao_setting["pseudopotential"],
                                                              system_pseudopot_nao_setting["numerical_orbital"]), 
                                           amwi.calculation(calculation_setting))
                    folder = amwi.folder_reduce(folder)
                    os.makedirs(folder, exist_ok=True) if not test_mode else print("mkdir {}".format(folder))
                    folders.append(folder)

                    
                    # write pseudopotential and numerical orbital over all combinations
                    _elements = amsb.scan_elements(_system)
                    # copy pseudopotential
                    pseudopotentials = {}
                    _pspotids = system_pseudopot_nao_setting["pseudopotential"] # get pseudopotential identifiers of one combination
                    for i in range(len(_pspotids)): # iterate over all elements
                        _pspotid = _pspotids[i] # for one element, get its pseudopotential identifier, however, which is this element?
                        _element = _elements[i]
                        fpseudo = pspot_archive[_pspotid] + "/" + valid_pseudopotentials[_element][_pspotid]["file"]
                        os.system("cp {} {}/".format(fpseudo, folder)) if not test_mode else print("cp {} {}/".format(fpseudo, folder))
                        pseudopotentials[_element] = valid_pseudopotentials[_element][_pspotid]["file"]
                    # copy numerical orbital
                    numerical_orbitals = {}
                    _naoids = system_pseudopot_nao_setting["numerical_orbital"] # get numerical orbital identifiers of one combination
                    for i in range(len(_naoids)):
                        _naoid = _naoids[i] + "@" + _pspotids[i]
                        _element = _elements[i]
                        fnao = nao_archive[_naoid] + "/"
                        fnao += valid_numerical_orbitals[_element][_naoid]["file"]
                        os.system("cp {} {}/".format(fnao, folder)) if not test_mode else print("cp {} {}/".format(fnao, folder))
                        numerical_orbitals[_element] = valid_numerical_orbitals[_element][_naoid]["file"]
                    isolated = False
                    for word in _system.split("_")[1]:
                        if word.isalpha():
                            isolated = True # mpid only contains numbers, so if it contains letters, it will be dimer, trimer or tetramer
                            break
                    if software == "qespresso":
                        iterate_qespresso(system_with_mpid=_system, target_folder=folder, calculation_setting=calculation_setting,
                                          pseudopotentials=pseudopotentials, 
                                          characteristic_length=characteristic_length, isolated=isolated,
                                          test_mode=test_mode)
                    elif software == "abacus":
                        iterate_abacus(system_with_mpid=_system, target_folder=folder, calculation_setting=calculation_setting,
                                       pseudopotentials=pseudopotentials, numerical_orbitals=numerical_orbitals,
                                       characteristic_length=characteristic_length, isolated=isolated,
                                       test_mode=test_mode)
    return folders

def iterate_abacus(system_with_mpid: str = "", # on which structure? system with mp-id
                   target_folder: str = "./", # in which folder? target folder
                   calculation_setting: dict ={}, # use what calculation setting?
                   pseudopotentials: dict = {}, # which set of pseudopotentials?
                   numerical_orbitals: dict = {}, # which set of numerical orbitals?
                   characteristic_length: float = 0.0, # for crystal, it is cell scaling; for molecule, it is bond length
                   isolated: bool = False,  # crystal or molecule?
                   test_mode: bool = True # for test
                   ) -> None:
    """iterate over all possible combinations of input parameters and generate folders
    especially for ABACUS

    Args:
        `system_with_mpid (str, optional)`: on which structure? system with mp-id. Defaults to "".
        `target_folder (str, optional)`: in which folder? target folder. Defaults to "./".
        `calculation_setting (dict, optional)`: use what calculation setting? Defaults to {}.
        `pseudopotentials (dict, optional)`: which set of pseudopotentials? Defaults to {}.
        `numerical_orbitals (dict, optional)`: which set of numerical orbitals? Defaults to {}.
        `characteristic_length (float, optional)`: how much distortion? Defaults to 0.0.
        `if_crystal (bool, optional)`: crystal or molecule? for crystal, distortion is cell scaling; for molecule, distortion is bond length. Defaults to True.
        `test_mode (bool, optional)`: for test. Defaults to True.
    """
    if system_with_mpid.count("_") != 1:
        raise ValueError("system_with_mpid should be like 'Er2O3_12345'")
    if calculation_setting == {}:
        print("Warning: calculation_setting is empty, use default")
    if pseudopotentials == {}:
        raise ValueError("pseudopotentials is empty")
    
    # write INPUT
    _input = amsag.INPUT(calculation_setting)
    if not test_mode:
        with open(target_folder + "/INPUT", "w") as f: f.write(_input)
    else:
        print("write %s\ncontents:\n%s"%(target_folder + "/INPUT", _input))
    # write STRU
    if not isolated:
        _stru, _cell = amsag._STRU_(
            fname=amwi.TEMPORARY_FOLDER + "/" + amwi.cif(system_with_mpid), 
            pseudopotentials=pseudopotentials,
            numerical_orbitals=numerical_orbitals,
            cell_scaling=characteristic_length)
    else:
        _stru, _cell = amsag._STRU_ISOLATED_(shape=system_with_mpid.split("_")[1],
                                             pseudopotentials=pseudopotentials,
                                             numerical_orbitals=numerical_orbitals,
                                             bond_length=characteristic_length)
    if not test_mode:
        with open(target_folder + "/STRU", "w") as f: f.write(_stru)
    else:
        print("write %s\ncontents:\n%s"%(target_folder + "/STRU", _stru))
    # write KPT
    _kpt = amsag._KPT_(isolated=isolated, cell_parameters=_cell)
    if not test_mode:
        with open(target_folder + "/KPT", "w") as f: f.write(_kpt)
    else:
        print("write %s\ncontents:\n%s"%(target_folder + "/KPT", _kpt))

def iterate_qespresso(system_with_mpid: str = "", # on which structure? system with mp-id, or element with shape (e.g. Er_dimer)
                      target_folder: str = "./", # in which folder? target folder
                      calculation_setting: dict = {}, # use what calculation setting?
                      pseudopotentials: dict = {}, # which set of pseudopotentials?
                      characteristic_length: float = 0.0, # for non isolated, it is cell scaling; for isolated, it is bond length
                      isolated: bool = False,  # crystal or molecule?
                      test_mode: bool = True # for test
                      ) -> None:
    """iterate over all possible combinations of input parameters and generate folders
    especially for Quantum Espresso

    Args:
        `system_with_mpid (str, optional)`: on which structure? system with mp-id, or element with shape (e.g. Er_dimer). Defaults to "".
        `target_folder (str, optional)`: in which folder? target folder. Defaults to "./".
        `calculation_setting (dict, optional)`: use what calculation setting? Defaults to {}.
        `pseudopotentials (dict, optional)`: which set of pseudopotentials? Defaults to {}.
        `characteristic_length (float, optional)`: for non isolated, it is cell scaling; for isolated, it is bond length. Defaults to 0.0.
        `isolated (bool, optional)`: crystal or molecule? Defaults to False.
        `test_mode (bool, optional)`: for test. Defaults to True.
    
    Raises:
        `ValueError`: system_with_mpid should be like 'Er2O3_12345'
        `ValueError`: calculation_setting is empty
        `ValueError`: pseudopotentials is empty
    
    Returns:
        `None`
    """
    _in = ""
    ntype, natom = len(pseudopotentials), 0
    if isolated:
        _in, natom = amsqg._ISOLATED(element=system_with_mpid.split("_")[0], 
                                     shape=system_with_mpid.split("_")[1],
                                     bond_length=characteristic_length)
    else:
        _in, natom = amsqg._CIF(fname=amwi.TEMPORARY_FOLDER + "/" + amwi.cif(system_with_mpid),
                                cell_scaling=characteristic_length,
                                constraints=[])
        
    _in = amsqg._ATOMIC_SEPCIES(pseudopotential=pseudopotentials) + _in

    iterate_sections = ["control", "system", "electrons", "ions", "cell"]
    for section in iterate_sections:
        _in = amsqg._calculation(section=section, ntype=ntype, natom=natom, **calculation_setting) + _in

    if not test_mode:
        with open(target_folder + "/scf.in", "w") as f: f.write(_in)
    else:
        print("write %s\ncontents:\n%s"%(target_folder + "/scf.in", _in))
        print("Contents:\n%s"%_in)
    return None

