"""Copy, configure and write files for specific software to run tests
abbreviation used globally: 
gppt: generate_pw_pseudopotential_test
"""
import os
import apns.module_workflow.identifier as id
import apns.module_io.output_basic as ob
import json
import apns.module_io.compress as cprs
import time
import itertools as it
# NOT REFACTORED YET
def qespresso(test_status: dict, functionals: list = ["pbe"], cell_scaling: list = [0.0]) -> list:
    """
    Generate Quantum ESPRESSO tests according to test_status, based on preset templates.

    Args:
        test_status (dict): test_status dictionary

    Returns:
        None

    Notes:
        template follows name convention defined in module_workflow.identifier.qespresso
    """
    os.chdir(test_status["paths"]["work_folder"])
    # ob._mkdir_(id.TEMPORARY_FOLDER)
    test_folders = []
    for functional in functionals:
        for system in test_status["tests"].keys():
            # so only focus on information of present system
            system_test_information = test_status["tests"][system]
            for test_of_system in test_status["tests"][system].keys():
                # extract information of present test of present system
                test_information = system_test_information[test_of_system]
                # name folder with preset name convention
                folder = id.folder(functional=functional,
                                   system=system,
                                   specific_test=test_of_system)
                test_folders.append(folder)
                ob._mkdir_(folder=folder)
                # copy input file to test folder, use unified name scf.in
                os.system("cp %s/%s %s/scf.in"%(id.TEMPORARY_FOLDER, 
                                                id.qespresso(system=system, 
                                                             template=True), 
                                                folder))
                ob._sed_(folder=folder.replace("\\", "/"),
                      file_to_sed="scf.in",
                      functional=functional,
                      ecutwfc=test_information["ecutwfc"])
                # copy pseudopotential to test folder
                p_information = test_information["pseudopotentials"]["info"]
                p_file = test_information["pseudopotentials"]["files"]
                for element in test_information["elements"]:
                    pinfo_element = p_information[element]
                    pseudopot_path = test_status["paths"]["pseudo_dir"] + "/"
                    pseudopot_path += id.pseudopotential(kind=pinfo_element["kind"], 
                                                         version=pinfo_element["version"], 
                                                         appendix=pinfo_element["appendix"])
                    pseudopot_file = p_file[element]
                    print("copy " + pseudopot_path + "/" + pseudopot_file + " to " + folder)
                    os.system("cp %s/%s %s"%(pseudopot_path.replace("\\", "/"), pseudopot_file, folder.replace("\\", "/")))
                    # swap pseudopotential file name in input file
                    os.system("sed -i 's/%s_pseudopot/%s/g' %s/scf.in"%(element, pseudopot_file, folder.replace("\\", "/")))
    
    os.chdir(test_status["paths"]["work_dir"])
    print_str = """------------------------------------------------------------------------------
Generation Done.
To run Quantum ESPRESSO tests on ABACUS Test, use the following parameters:
rundft: mpirun -np 16 pw.x -i scf.in | tee scf.log; rm -rf pwscf.save
Bohrium image: registry.dp.tech/dptech/prod-471/abacus-vasp-qe:20230116
Bohrium_machine_type: c32_m128_cpu
Bohrium_platform: ali
------------------------------------------------------------------------------
"""
    print(print_str)
    return test_folders

def _qespresso_(work_status: dict, test_status: dict):

    pseudopot_folders_arch = {}
    with open(test_status["paths"]["pseudo_dir"] + "/" + "description.json") as json_f:
        pseudopot_folders_arch = json.load(json_f)
    os.chdir(test_status["paths"]["work_dir"])

    lists_to_iterate = [
        work_status["calculations"]["functionals"]
        ]
    keys_to_literate = ["input_dft", "cell_scalings"] # what is the counterpart here???
    for other in work_status["calculation"].keys():
        if isinstance(work_status["calculation"][other], list):
            lists_to_iterate.append(work_status["calculation"][other])
            keys_to_literate.append(other)
    shortkeys_to_literate = [key[0] for key in keys_to_literate]
    # generate all combinations of keywords
    combinations = list(it.product(*lists_to_iterate))

    test_folders = []
    for system in test_status["tests"].keys():
        system_test_information = test_status["tests"][system]
        for test in test_status["tests"][system].keys():
            test_information = system_test_information[test]
            # then generate all combinations of keywords for one test
            for c in combinations:
                appendix = "".join([str(item[0])+str(item[1]) for item in zip(shortkeys_to_literate, c)])
                folder = "_".join(["t", system, test, appendix])
                test_folders.append(folder)
                ob._mkdir_(folder=folder)
                # copy template files
                fin = id.qespresso(system=system, template=True)
                os.system("cp %s/%s %s/scf.in"%(id.TEMPORARY_FOLDER, fin, folder))
                for i, item in enumerate(c):
                    key = keys_to_literate[i]
                    value = c[i]
                    os.system("sed -i 's/%s_to_test/%s/g' %s/scf.in"%(key, value, folder))
                # copy pseudopotential to test folder
                p_information = test_information["pseudopotentials"]["info"]
                p_file = test_information["pseudopotentials"]["files"]
                for element in test_information["elements"]:
                    pinfo_element = p_information[element]
                    _identifier =  id.pseudopotential(kind=pinfo_element["kind"], 
                                                        version=pinfo_element["version"], 
                                                        appendix=pinfo_element["appendix"])
                    pseudopot_path = pseudopot_folders_arch[_identifier]
                    pseudopot_file = p_file[element]
                    print("copy " + pseudopot_path + "/" + pseudopot_file + " to " + folder)
                    os.system("cp %s/%s %s"%(pseudopot_path.replace("\\", "/"), pseudopot_file, folder.replace("\\", "/")))
                    # swap pseudopotential file name in input file
                    os.system("sed -i 's/%s_pseudopot/%s/g' %s/scf.in"%(element, pseudopot_file, folder.replace("\\", "/")))

    os.chdir(test_status["paths"]["work_dir"])
    print_str = """------------------------------------------------------------------------------
Generation Done.
To run Quantum ESPRESSO tests on ABACUS Test, use the following parameters:
rundft: mpirun -np 16 pw.x -i scf.in | tee scf.log; rm -rf pwscf.save
Bohrium image: registry.dp.tech/dptech/prod-471/abacus-vasp-qe:20230116
Bohrium_machine_type: c32_m128_cpu
Bohrium_platform: ali
------------------------------------------------------------------------------
"""
    print(print_str)
    return test_folders

def abacus(test_status: dict, functionals: list = ["pbe"], cell_scalings: list = [0.0]) -> list:

    pseudopot_folders_arch = {}
    with open(test_status["paths"]["pseudo_dir"] + "/" + "description.json") as json_f:
        pseudopot_folders_arch = json.load(json_f)
    os.chdir(test_status["paths"]["work_dir"])

    with_cell_appendix = True
    if len(cell_scalings) == 1 and cell_scalings[0] == 0.0:
        with_cell_appendix = False
    # ob._mkdir_(id.TEMPORARY_FOLDER)
    test_folders = []
    for functional in functionals:
        for cell_scaling in cell_scalings:
            for system in test_status["tests"].keys():
                # so only focus on information of present system
                system_test_information = test_status["tests"][system]
                for test in test_status["tests"][system].keys():
                    test_information = system_test_information[test]
                    folder = id.folder(functional=functional,
                                    system=system,
                                    specific_test=test)
                    if with_cell_appendix:
                        folder += "_cell_" + str(cell_scaling)
                    test_folders.append(folder)
                    ob._mkdir_(folder=folder) # mkdir for present test
                    # copy template files
                    fINPUT, fSTRU, fKPT = id.abacus(system=system,
                                                    template=True)
                    os.system("cp %s/%s %s/INPUT"%(id.TEMPORARY_FOLDER,
                                                    fINPUT,
                                                    folder))
                    os.system("cp %s/%s %s/STRU"%(id.TEMPORARY_FOLDER,
                                                    fSTRU,
                                                    folder))
                    os.system("cp %s/%s %s/KPT"%(id.TEMPORARY_FOLDER,
                                                fKPT,
                                                folder))
                    # whether should be support multiple ecutwfc inside APNS? ABACUSTest has implemented that.
                    ob._sed_(folder=folder.replace("\\", "/"),
                        file_to_sed="INPUT",
                        functional=functional,
                        ecutwfc=test_information["ecutwfc"])
                    lattice_constant = 1.889725989 * (1 + cell_scaling)
                    ob._sed_(folder=folder.replace("\\", "/"),
                            file_to_sed="STRU",
                            lattice_constant=lattice_constant)
                    # copy pseudopotential to test folder
                    p_information = test_information["pseudopotentials"]["info"]
                    p_file = test_information["pseudopotentials"]["files"]
                    for element in test_information["elements"]:
                        pinfo_element = p_information[element]
                        _identifier =  id.pseudopotential(kind=pinfo_element["kind"], 
                                                        version=pinfo_element["version"], 
                                                        appendix=pinfo_element["appendix"])
                        pseudopot_path = pseudopot_folders_arch[_identifier]
                        pseudopot_file = p_file[element]
                        print("copy " + pseudopot_path + "/" + pseudopot_file + " to " + folder)
                        os.system("cp %s/%s %s"%(pseudopot_path.replace("\\", "/"), pseudopot_file, folder.replace("\\", "/")))
                        # swap pseudopotential file name in input file
                        os.system("sed -i 's/%s_pseudopot/%s/g' %s/STRU"%(element, pseudopot_file, folder.replace("\\", "/")))

    os.chdir(test_status["paths"]["work_dir"])
    print_str = """------------------------------------------------------------------------------
Generation Done.
To run ABACUS tests on ABACUS Test, use the following parameters:
rundft: OMP_NUM_THREADS=16 mpirun -np 1 abacus | tee out.log 
(especially when dft_functional = HSE, arbitrary if PBE)
Bohrium image: registry.dp.tech/deepmodeling/abacus-intel:latest (this one is default)
Bohrium_machine_type: c32_m128_cpu (c32_m64_cpu has bad performance)
Bohrium_platform: ali
------------------------------------------------------------------------------
"""
    print(print_str)
    return test_folders

def _abacus_(work_status: dict, test_status: dict):
    """new version of function abacus, only take additional and enumerative keywords from work_status.
    New code has following differences with the old one:
    
    1. Folder nomenclature
    old: t_pbe_test_system_mpid_cell_0.0
    new: t_test_system_mpid_fpbec0.0

    2. Requiring changes in the way writing template files
    now the logic changes to, except functionals and cell_scalings, ecutwfc, all other additionaly keywords
    if are enumerative, then template will write as xxx_to_test, then it is sed to exact value here.
    """
    pseudopot_folders_arch = {}
    with open(test_status["paths"]["pseudo_dir"] + "/" + "description.json") as json_f:
        pseudopot_folders_arch = json.load(json_f)
    os.chdir(test_status["paths"]["work_dir"])

    lists_to_iterate = [
        work_status["calculations"]["functionals"]
        ]
    keys_to_literate = ["dft_functional", "cell_scalings"]
    for other in work_status["calculation"].keys():
        if isinstance(work_status["calculation"][other], list):
            lists_to_iterate.append(work_status["calculation"][other])
            keys_to_literate.append(other)
    shortkeys_to_literate = [key[0] for key in keys_to_literate]
    # generate all combinations of keywords
    combinations = list(it.product(*lists_to_iterate))

    test_folders = []
    for system in test_status["tests"].keys():
        system_test_information = test_status["tests"][system]
        for test in test_status["tests"][system].keys():
            test_information = system_test_information[test]
            # then generate all combinations of keywords for one test
            for c in combinations:
                appendix = "".join([str(item[0])+str(item[1]) for item in zip(shortkeys_to_literate, c)])
                folder = "_".join(["t", system, test, appendix])
                test_folders.append(folder)
                ob._mkdir_(folder=folder)
                # copy template files
                fINPUT, fSTRU, fKPT = id.abacus(system=system, template=True) # the names of template files
                for fn, ft in zip([fINPUT, fSTRU, fKPT], ["INPUT", "STRU", "KPT"]):
                    os.system("cp %s/%s %s/%s"%(id.TEMPORARY_FOLDER, fn, folder, ft))
                for i, item in enumerate(c):
                    key = keys_to_literate[i]
                    value = c[i]
                    for fn in ["INPUT", "STRU", "KPT"]:
                        os.system("sed -i 's/%s_to_test/%s/g' %s/%s"%(key, value, folder, fn))
                # copy pseudopotential to test folder
                p_information = test_information["pseudopotentials"]["info"]
                p_file = test_information["pseudopotentials"]["files"]
                for element in test_information["elements"]:
                    pinfo_element = p_information[element]
                    _identifier =  id.pseudopotential(kind=pinfo_element["kind"], 
                                                        version=pinfo_element["version"], 
                                                        appendix=pinfo_element["appendix"])
                    pseudopot_path = pseudopot_folders_arch[_identifier]
                    pseudopot_file = p_file[element]
                    print("copy " + pseudopot_path + "/" + pseudopot_file + " to " + folder)
                    os.system("cp %s/%s %s"%(pseudopot_path.replace("\\", "/"), pseudopot_file, folder.replace("\\", "/")))
                    # swap pseudopotential file name in input file
                    os.system("sed -i 's/%s_pseudopot/%s/g' %s/STRU"%(element, pseudopot_file, folder.replace("\\", "/")))
                    
    os.chdir(test_status["paths"]["work_dir"])
    print_str = """------------------------------------------------------------------------------
Generation Done.
To run ABACUS tests on ABACUS Test, use the following parameters:
rundft: OMP_NUM_THREADS=16 mpirun -np 1 abacus | tee out.log 
(especially when dft_functional = HSE, arbitrary if PBE)
Bohrium image: registry.dp.tech/deepmodeling/abacus-intel:latest (this one is default)
Bohrium_machine_type: c32_m128_cpu (c32_m64_cpu has bad performance)
Bohrium_platform: ali
------------------------------------------------------------------------------
"""
    print(print_str)
    return test_folders

def generate(test_status: dict, software: str, functionals: list = ["pbe"], cell_scalings: list = [0.0]):
    """Generate all input scripts for specific software in folders

    Args:
        test_status (dict): test_status
    """
    #software = test_status["global"]["software"]
    test_folders = []
    if software.lower().startswith("q") and software.lower().endswith("espresso"):
        test_folders = qespresso(test_status,
                                 functionals=functionals)
    elif software == "ABACUS":
        test_folders = abacus(test_status,
                              functionals=functionals,
                              cell_scalings=cell_scalings)
    else:
        print("Error: software not recognized.")
        exit(1)
    # compress test folders
    timestamp = time.strftime("%Y%m%d%H%M%S", time.localtime())
    zip_fname = "tests_" + timestamp + ".zip"
    
    cprs.pack(test_folders, zip_fname)

    print("Have compressed test suites into file ", zip_fname, ", cleaning...")
    for folder in test_folders:
        os.system("rm -rf {}".format(folder))