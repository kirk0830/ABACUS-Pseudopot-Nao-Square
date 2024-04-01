import apns.module_workflow.identifier as amwi
import apns.module_software.abacus.generation as amsag
import apns.module_software.qespresso.generation as amsqg
import apns.module_structure.basic as amsb
import os

def iterate(**kwargs):
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
        
        `extensive` (dict, optional): extensive settings in integrated in apns. Defaults to empty.
        >>> {
            "characteristic_lengths": [characteristic_length, ...],
            "nkpoints_in_line": [nkpoints_in_line, ...]
            ...
        }

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
    software = kwargs.get("software", "abacus")
    systems = kwargs.get("systems", [])
    extensive_settings = kwargs.get("extensive_settings", [])
    upforb_bundles = kwargs.get("upforb_bundles", [])
    calculation_settings = kwargs.get("calculation_settings", [])
    upfs = kwargs.get("upfs", {})
    orbs = kwargs.get("orbs", {})
    test_mode = kwargs.get("test_mode", True)

    folders = []
    for _is, _system in enumerate(systems): # iterate over structures
        upforb_bundles_system = upforb_bundles[_is]
        for upforb_bundle in upforb_bundles_system: # iterate over pseudopotential and numerical orbital settings
            # spns: system_pseudopot_nao_setting
            for calculation_setting in calculation_settings: # iterate over calculation settings
                for extensive_setting in extensive_settings: # iterate over extensive settings
                    # make folder
                    if "numerical_orbital" not in upforb_bundle.keys():
                        upforb_bundle["numerical_orbital"] = []
                    folder = amwi._folder_(system=_system,
                                           pseudo_nao_identifier=amwi.pseudopot_nao(**upforb_bundle).replace(".", "").replace("_", ""),
                                           calculation_identifier=amwi.calculation(calculation_setting),
                                           extensive_identifier=amwi.extensive(extensive_setting))
                    # if name is too long, find ways to reduce the length
                    folder = amwi.folder_reduce(folder)
                    if len(folder) > 100:
                        folder = "_".join(folder.split("_")[1:])
                    os.makedirs(folder, exist_ok=True) if not test_mode else print("mkdir {}".format(folder))
                    folders.append(folder)

                    # write pseudopotential and numerical orbital over all combinations
                    _elements = amsb.scan_elements(_system)
                    # copy pseudopotential
                    pseudopotentials = {}
                    # upforb_bundle["pseudopotential"] arranges like [pseudopotential_identifier, ...]
                    # indiced by index of element
                    _pspotids = upforb_bundle["pseudopotential"] # get pseudopotential identifiers of one combination
                    for i in range(len(_pspotids)): # iterate over all elements
                        _pspotid = _pspotids[i] # for one element, get its pseudopotential identifier, however, which is this element?
                        _element = _elements[i]
                        fpseudo_withpath = upfs[_element][_pspotid]
                        fpseudo = fpseudo_withpath.replace("\\", "/").split("/")[-1]
                        os.system(f"cp {fpseudo_withpath} {folder}/{fpseudo}") if not test_mode else print(f"cp {fpseudo_withpath} {folder}/{fpseudo}")
                        pseudopotentials[_element] = fpseudo
                    # copy numerical orbital
                    numerical_orbitals = {}
                    _naoids = upforb_bundle["numerical_orbital"] # get numerical orbital identifiers of one combination
                    for i in range(len(_naoids)):
                        _naoid = _naoids[i] + "@" + _pspotids[i]
                        _element = _elements[i]
                        fnao_withpath = orbs[_element][_naoid]
                        fnao = fnao_withpath.replace("\\", "/").split("/")[-1]
                        os.system("cp {} {}/{}".format(fnao_withpath, folder, fnao)) if not test_mode else print("cp {} {}/{}".format(fnao_withpath, folder, fnao))
                        numerical_orbitals[_element] = fnao
                    for word in _system.split("_")[1]:
                        if word.isalpha():
                            extensive_setting["nkpoints_in_line"] = -1 # mpid only contains numbers, so if it contains letters, it will be dimer, trimer or tetramer
                            break
                    if software == "qespresso":
                        iterate_qespresso(system_with_mpid=_system, 
                                          target_folder=folder, 
                                          calculation_setting=calculation_setting,
                                          extensive_setting=extensive_setting,
                                          pseudopotentials=pseudopotentials,
                                          test_mode=test_mode)
                    elif software == "abacus":
                        iterate_abacus(system_with_mpid=_system, 
                                       target_folder=folder, 
                                       calculation_setting=calculation_setting,
                                       extensive_setting=extensive_setting,
                                       pseudopotentials=pseudopotentials, numerical_orbitals=numerical_orbitals,
                                       test_mode=test_mode)
    return folders

def iterate_abacus(system_with_mpid: str = "", # on which structure? system with mp-id
                   target_folder: str = "./", # in which folder? target folder
                   calculation_setting: dict = None, # use what calculation setting?
                   pseudopotentials: dict = None, # which set of pseudopotentials?
                   numerical_orbitals: dict = None, # which set of numerical orbitals?
                   extensive_setting: dict = None,
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
    if calculation_setting is None:
        print("Warning: calculation_setting is empty, use default")
    if pseudopotentials is None:
        raise ValueError("pseudopotentials is empty")
    
    # write INPUT
    _input = amsag.INPUT(calculation_setting)
    if not test_mode:
        with open(target_folder + "/INPUT", "w") as f: f.write(_input)
    else:
        print("write %s\ncontents:\n%s"%(target_folder + "/INPUT", _input))
    # write STRU
    numerical_orbitals = None if numerical_orbitals is not None and len(numerical_orbitals) == 0 else numerical_orbitals

    if extensive_setting["nkpoints_in_line"] >= 0:
        starting_magnetization = amsb.starting_magnetization(amwi.cif(system_with_mpid), magnetism=extensive_setting["magnetism"])
        _stru, _cell = amsag.STRU_Pymatgen(fname=amwi.TEMPORARY_FOLDER + "/" + amwi.cif(system_with_mpid),
                                           pseudopotentials=pseudopotentials,
                                           numerical_orbitals=numerical_orbitals,
                                           cell_scaling=extensive_setting["characteristic_lengths"],
                                           starting_magnetization=starting_magnetization)
    else:
        starting_magnetization = amsb.starting_magnetization(system_with_mpid.split("_")[1], magnetism=extensive_setting["magnetism"])
        _stru, _cell = amsag.STRU_Molecule(shape=system_with_mpid.split("_")[1],
                                           pseudopotentials=pseudopotentials,
                                           numerical_orbitals=numerical_orbitals,
                                           bond_length=extensive_setting["characteristic_lengths"])
    if not test_mode:
        with open(target_folder + "/STRU", "w") as f: f.write(_stru)
    else:
        print("write %s\ncontents:\n%s"%(target_folder + "/STRU", _stru))
    # write KPT
    if extensive_setting["nkpoints_in_line"] < 0:
        _kpt = amsag._KPT_(isolated=True, cell_parameters=_cell)
    elif extensive_setting["nkpoints_in_line"] == 0:
        _kpt = amsag._KPT_(isolated=False, cell_parameters=_cell)
    else:
        _kpt = amsag._KLINE_(fname=amwi.TEMPORARY_FOLDER + "/" + amwi.cif(system_with_mpid),
                             nkpts_in_line=extensive_setting["nkpoints_in_line"])
    if not test_mode:
        with open(target_folder + "/KPT", "w") as f: f.write(_kpt)
    else:
        print("write %s\ncontents:\n%s"%(target_folder + "/KPT", _kpt))

def iterate_qespresso(system_with_mpid: str = "", # on which structure? system with mp-id, or element with shape (e.g. Er_dimer)
                      target_folder: str = "./", # in which folder? target folder
                      calculation_setting: dict = None, # use what calculation setting?
                      pseudopotentials: dict = None, # which set of pseudopotentials?
                      extensive_setting: dict = None,
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
    if calculation_setting is None:
        print("Warning: calculation_setting is empty, use default")
    if pseudopotentials is None:
        raise ValueError("pseudopotentials is empty")
    
    _in = ""
    ntype, natom = len(pseudopotentials), 0
    if extensive_setting["nkpoints_in_line"] < 0:
        _in, natom = amsqg._ISOLATED(element=system_with_mpid.split("_")[0], 
                                     shape=system_with_mpid.split("_")[1],
                                     bond_length=extensive_setting["characteristic_lengths"])
    else:
        _in, natom = amsqg._CIF(fname=amwi.TEMPORARY_FOLDER + "/" + amwi.cif(system_with_mpid),
                                cell_scaling=extensive_setting["characteristic_lengths"],
                                constraints=[])
        
    _in = amsqg._K_POINTS(fname=amwi.TEMPORARY_FOLDER + "/" + amwi.cif(system_with_mpid), 
                          nkpoints_in_line=extensive_setting["nkpoints_in_line"]) + _in
    
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

