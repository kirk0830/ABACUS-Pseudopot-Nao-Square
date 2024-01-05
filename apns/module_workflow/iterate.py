import apns.module_workflow.identifier as amwi
import apns.module_software.abacus.generation as amsag
import apns.module_structure.basic as amsb
import os

def iterate(systems: list = [],
            pseudopot_nao_settings: list = [],
            calculation_settings: list = [],
            cell_scalings: list = [],
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
        
        `cell_scalings (list, optional)`: cell scaling settings to iterate = 
        >>> [0, 0.01, ...]. 
        
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
                for cell_scaling in cell_scalings: # iterate over cell scaling
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
                    # write INPUT
                    _input = amsag.INPUT(calculation_setting)
                    if not test_mode:
                        with open(folder + "/INPUT", "w") as f: f.write(_input)
                    else:
                        print("write %s\ncontents:\n%s"%(folder + "/INPUT", _input))
                    # write pseudopotential and numerical orbital over all combinations
                    _elements = amsb.scan_elements(_system)
                    # pseudopotential
                    pseudopotentials = {}
                    _pspotids = system_pseudopot_nao_setting["pseudopotential"] # get pseudopotential identifiers of one combination
                    for i in range(len(_pspotids)): # iterate over all elements
                        _pspotid = _pspotids[i] # for one element, get its pseudopotential identifier, however, which is this element?
                        _element = _elements[i]
                        fpseudo = pspot_archive[_pspotid] + "/"
                        fpseudo += valid_pseudopotentials[_element][_pspotid]["file"]
                        os.system("cp {} {}/".format(fpseudo, folder)) if not test_mode else print("cp {} {}/".format(fpseudo, folder))
                        pseudopotentials[_element] = valid_pseudopotentials[_element][_pspotid]["file"]
                    # numerical orbital
                    numerical_orbitals = {}
                    _naoids = system_pseudopot_nao_setting["numerical_orbital"] # get numerical orbital identifiers of one combination
                    for i in range(len(_naoids)):
                        _naoid = _naoids[i] + "@" + _pspotids[i]
                        _element = _elements[i]
                        fnao = nao_archive[_naoid] + "/"
                        fnao += valid_numerical_orbitals[_element][_naoid]["file"]
                        os.system("cp {} {}/".format(fnao, folder)) if not test_mode else print("cp {} {}/".format(fnao, folder))
                        numerical_orbitals[_element] = valid_numerical_orbitals[_element][_naoid]["file"]
                    # write STRU
                    _stru, _cell = amsag._STRU_(
                        fname=amwi.TEMPORARY_FOLDER + "/" + amwi.cif(_system), 
                        pseudopotentials=pseudopotentials,
                        numerical_orbitals=numerical_orbitals,
                        cell_scaling=cell_scaling)
                    if not test_mode:
                        with open(folder + "/STRU", "w") as f: f.write(_stru)
                    else:
                        print("write %s\ncontents:\n%s"%(folder + "/STRU", _stru))
                    # write KPT
                    _kpt = amsag.KPT_generation("crystal", cell_parameters=_cell)
                    if not test_mode:
                        with open(folder + "/KPT", "w") as f: f.write(_kpt)
                    else:
                        print("write %s\ncontents:\n%s"%(folder + "/KPT", _kpt))
    return folders

if __name__ == "__main__":

    print(iterate(
        systems=["Er2O3_679"],
        pseudopot_nao_settings=[
            {
                "pseudopotential": ["sg15_10", "sg15_10"],
                "numerical_orbital": ["TZDP_6", "TZDP_10"]
            },
            {
                "pseudopotential": ["sg15_11", "sg15_10"],
                "numerical_orbital": ["TZDP_6", "TZDP_10"]
            },
            {
                "pseudopotential": ["pd_04", "sg15_10"],
                "numerical_orbital": ["TZDP_10", "TZDP_10"]
            }
        ],
        calculation_settings=[
            {
                "basis_type": "lcao",
                "dft_functional": "pbe",
                "ecutwfc": 100
            },
            {
                "basis_type": "pw",
                "dft_functional": "pbe",
                "ecutwfc": 200
            },
            {
                "basis_type": "pw",
                "dft_functional": "pbe",
                "ecutwfc": 300
            }
        ],
        cell_scalings=[0.00],
        valid_pseudopotentials={
            "Er": {
                "sg15_10": {
                    "file": "Er_ONCV_PBE-1.0.upf",
                    "kind": "sg15",
                    "version": "10",
                    "appendix": ""
                },
                "sg15_11": {
                    "file": "Er_ONCV_PBE-1.1.upf",
                    "kind": "sg15",
                    "version": "11",
                    "appendix": ""
                },
                "pd_04": {
                    "file": "Er.PD04.PBE.UPF",
                    "kind": "pd",
                    "version": "04",
                    "appendix": ""
                }
            },
            "O": {
                "sg15_10": {
                    "file": "O_ONCV_PBE-1.0.upf",
                    "kind": "sg15",
                    "version": "10",
                    "appendix": ""
                }
            }
        },
        valid_numerical_orbitals={
            "Er": {
                "TZDP_6@sg15_10": {
                    "file": "Cr_gga_6au_100Ry_3s3p3d2f.orb",
                    "type": "TZDP",
                    "rcut": "6",
                    "appendix": ""
                },
                "TZDP_6@sg15_11": {
                    "file": "Cr_gga_6au_100Ry_3s3p3d2f.orb",
                    "type": "TZDP",
                    "rcut": "6",
                    "appendix": ""
                },
                "TZDP_10@pd_04": {
                    "file": "Cr_gga_10au_100Ry_3s3p3d2f.orb",
                    "type": "TZDP",
                    "rcut": "10",
                    "appendix": ""
                }
            },
            "O": {
                "TZDP_10@sg15_10": {
                    "file": "O_gga_10au_100Ry_3s3p3d2f.orb",
                    "type": "TZDP",
                    "rcut": "10",
                    "appendix": ""
                }
            }
        },
        pspot_archive={
            "sg15_10": "/home/zhuzhen/Work/psp/sg15/1.0/",
            "sg15_11": "/home/zhuzhen/Work/psp/sg15/1.1/",
            "pd_04": "/home/zhuzhen/Work/psp/pd/04/"
        },
        nao_archive={
            "TZDP_6@sg15_10": "/home/zhuzhen/Work/nao/sg15/1.0/",
            "TZDP_10@sg15_10": "/home/zhuzhen/Work/nao/sg15/1.0/",
            "TZDP_6@sg15_11": "/home/zhuzhen/Work/nao/sg15/1.1/",
            "TZDP_10@pd_04": "/home/zhuzhen/Work/nao/pd/04/"
        },
        test_mode=True
    ))