"""Copy, configure and write files for specific software to run tests
abbreviation used globally: 
gppt: generate_pw_pseudopotential_test
"""
import os
import module_workflow.identifier as id
import module_io.output_basic as ob

def _qespresso_(test_status: dict):
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
    ob._mkdir_(id.TEMPORARY_FOLDER)
    for functional in test_status["calculation"]["functionals"]:
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
                ob._mkdir_(folder=folder)
                # copy input file to test folder, use unified name scf.in
                os.system("cp %s\\%s %s\\scf.in"%(id.TEMPORARY_FOLDER, 
                                                  id.qespresso(system=system, 
                                                               template=True), 
                                                  folder))
                ob._sed_(folder=folder,
                      file_to_sed="scf.in",
                      functional=functional,
                      ecutwfc=test_information["ecutwfc"])
                # copy pseudopotential to test folder
                p_information = test_information["pseudopotentials"]["info"]
                p_file = test_information["pseudopotentials"]["files"]
                for element in test_information["elements"]:
                    pinfo_element = p_information[element]
                    pseudopot_path = test_status["paths"]["pseudo_dir"] + "\\"
                    pseudopot_path += id.pseudopotential(kind=pinfo_element["kind"], 
                                                         version=pinfo_element["version"], 
                                                         appendix=pinfo_element["appendix"])
                    pseudopot_path += "\\"
                    pseudopot_file = p_file[element]
                    print("copy " + pseudopot_path + pseudopot_file + " to " + folder)
                    os.system("cp %s\\%s %s"%(pseudopot_path, pseudopot_file, folder))
                    # swap pseudopotential file name in input file
                    os.system("sed -i 's/%s_pseudopot/%s/g' %s\\scf.in"%(element, pseudopot_file, folder))
    
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

def _abacus_(test_status: dict):

    os.chdir(test_status["paths"]["work_folder"])
    ob._mkdir_(id.TEMPORARY_FOLDER)
    for functional in test_status["calculation"]["functionals"]:
        for system in test_status["tests"].keys():
            # so only focus on information of present system
            system_test_information = test_status["tests"][system]
            for test in test_status["tests"][system].keys():
                test_information = system_test_information[test]
                folder = id.folder(functional=functional,
                                   system=system,
                                   specific_test=test)
                ob._mkdir_(folder=folder) # mkdir for present test
                # copy template files
                fINPUT, fSTRU, fKPT = id.abacus(system=system,
                                                template=True)
                os.system("cp %s\\%s %s\\INPUT"%(id.TEMPORARY_FOLDER,
                                                 fINPUT,
                                                 folder))
                os.system("cp %s\\%s %s\\STRU"%(id.TEMPORARY_FOLDER,
                                                fSTRU,
                                                folder))
                os.system("cp %s\\%s %s\\KPT"%(id.TEMPORARY_FOLDER,
                                               fKPT,
                                               folder))
                ob._sed_(folder=folder,
                      file_to_sed="INPUT",
                      functional=functional,
                      ecutwfc=test_information["ecutwfc"])
                # copy pseudopotential to test folder
                p_information = test_information["pseudopotentials"]["info"]
                p_file = test_information["pseudopotentials"]["files"]
                for element in test_information["elements"]:
                    pinfo_element = p_information[element]
                    pseudopot_path = test_status["paths"]["pseudo_dir"] + "\\"
                    pseudopot_path += id.pseudopotential(kind=pinfo_element["kind"], 
                                                         version=pinfo_element["version"], 
                                                         appendix=pinfo_element["appendix"])
                    pseudopot_file = p_file[element]
                    print("copy " + pseudopot_path + pseudopot_file + " to " + folder)
                    os.system("cp %s\\%s %s"%(pseudopot_path, pseudopot_file, folder))
                    # swap pseudopotential file name in input file
                    os.system("sed -i 's/%s_pseudopot/%s/g' %s\\STRU"%(element, pseudopot_file, folder))

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

def generate(test_status: dict):
    """Generate all input scripts for specific software in folders

    Args:
        test_status (dict): test_status
    """
    software = test_status["global"]["software"]
    if software.lower().startswith("q") and software.lower().endswith("espresso"):
        _qespresso_(test_status)
    elif software == "ABACUS":
        _abacus_(test_status)
    else:
        print("Error: software not recognized.")
        exit(1)
