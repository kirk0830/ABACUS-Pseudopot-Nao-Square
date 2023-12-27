import os
import apns.module_database.database as ai
import apns.module_workflow.identifier as id
import apns.module_io.output_basic as ob

def generate(test_status: dict):

    if test_status["global"]["software"] == "ABACUS":
        _abacus_(test_status)
    else:
        print("Error: software not recognized.")
        exit(1)

def _abacus_(test_status: dict):

    os.chdir(test_status["paths"]["work_folder"])
    # ob._mkdir_(id.TEMPORARY_FOLDER)
    for functional in test_status["calculation"]["functionals"]:
        for system in test_status["tests"].keys():
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
                os.system("cp %s/%s %s/INPUT"%(id.TEMPORARY_FOLDER,
                                                 fINPUT,
                                                 folder))
                os.system("cp %s/%s %s/STRU"%(id.TEMPORARY_FOLDER,
                                                fSTRU,
                                                folder))
                os.system("cp %s/%s %s/KPT"%(id.TEMPORARY_FOLDER,
                                               fKPT,
                                               folder))
                ob._sed_(folder=folder,
                      file_to_sed="INPUT",
                      functional=functional,
                      ecutwfc=test_information["ecutwfc"])
                # copy pseudopotential and numerical orbital to test folder
                p_information = test_information["pseudopotentials"]["info"]
                n_information = test_information["numerical_orbitals"]["info"]
                p_file = test_information["pseudopotentials"]["files"]
                n_file = test_information["numerical_orbitals"]["files"]

                for element in test_status["tests"][system][test]["elements"]:
                    pinfo_element = p_information[element]
                    pseudopot_path = test_status["paths"]["pseudo_dir"] + "/"
                    pseudopot_path += id.pseudopotential(kind=pinfo_element["kind"], 
                                                         version=pinfo_element["version"], 
                                                         appendix=pinfo_element["appendix"])
                    pseudopot_path += "/"
                    pseudopot_file = test_status["tests"][system][test]["pseudopotentials"]["files"][element]
                    print("copy " + pseudopot_path + pseudopot_file + " to " + folder)
                    os.system("cp %s/%s %s"%(pseudopot_path, pseudopot_file, folder))
                    # swap pseudopotential file name in input file
                    os.system("sed -i 's/{}_pseudopot/".format(element) + pseudopot_file + "/g' " + folder + "/STRU")

                    ninfo_element = n_information[element]
                    nao_path = test_status["paths"]["orbital_dir"] + "/"
                    nao_path += id.numerical_orbital(type = ninfo_element["type"],
                                                     rcut = ninfo_element["rcut"],
                                                     appendix = ninfo_element["appendix"])
                    nao_et = str(ai.get_element_index(element)) + "_" + element + "_"
                    nao_et += ninfo_element["type"]
                    nao_path += nao_et + "/"
                    nao_file = n_file[element]
                    print("copy " + nao_path + nao_file + " to " + folder)
                    os.system("cp " + nao_path + nao_file + " " + folder)
                    # swap numerical orbital file name in input file
                    os.system("sed -i 's/{}_numerical_orbital/".format(element) + nao_file + "/g' " + folder + "/STRU")

    os.chdir(test_status["paths"]["root"])
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

def packup_tests(test_status: dict):

    pass