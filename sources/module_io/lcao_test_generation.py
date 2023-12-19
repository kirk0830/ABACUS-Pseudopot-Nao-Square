import os
import module_database.database as ai

def generate_lcao_pseudopotential_test(test_status: dict):

    if test_status["global"]["software"] == "ABACUS":
        glpt_abacus(test_status)
    else:
        print("Error: software not recognized.")
        exit(1)

def glpt_abacus(test_status: dict):

    os.chdir(test_status["paths"]["work_folder"])
    for functional in test_status["calculation"]["functionals"]:
        for system in test_status["tests"].keys():
            for test in test_status["tests"][system].keys():
                folder = "t_" + functional + "_" + test + "_" + system
                if os.path.exists(folder):
                    print("Warning: " + folder + " already exists, clean.")
                    os.system("rm -rf " + folder)
                os.mkdir(folder)
                input_file_path = test_status["paths"]["numerical_orbital"]["work"] + "\\templates\\ABACUS\\INPUT"
                stru_file_path = test_status["paths"]["numerical_orbital"]["work"] + "\\templates\\ABACUS\\STRU_" + system
                kpt_file_path = test_status["paths"]["numerical_orbital"]["work"] + "\\templates\\ABACUS\\KPT"
                # copy INPUT, STRU, KPT to test folder, use unified name STRU on STRU_*
                os.system("cp " + input_file_path + " " + folder + "\\INPUT")
                # functional is configured here
                os.system("sed -i 's/functional_to_test/" + functional + "/g' " + folder + "\\INPUT")
                os.system("cp " + stru_file_path + " " + folder + "\\STRU")
                os.system("cp " + kpt_file_path + " " + folder + "\\KPT")
                # copy pseudopotential and numerical orbital to test folder
                for element in test_status["tests"][system][test]["elements"]:
                    pseudopot_path = test_status["paths"]["pseudopotential"]["resources"] + "\\"
                    pseudopot_kva = test_status["tests"][system][test]["pseudopotentials"]["info"][element]["kind"]
                    if test_status["tests"][system][test]["pseudopotentials"]["info"][element]["version"] != "":
                        pseudopot_kva += "_" + test_status["tests"][system][test]["pseudopotentials"]["info"][element]["version"]
                    if test_status["tests"][system][test]["pseudopotentials"]["info"][element]["appendix"] != "":
                        pseudopot_kva += "_" + test_status["tests"][system][test]["pseudopotentials"]["info"][element]["appendix"]
                    pseudopot_path += pseudopot_kva + "\\"
                    pseudopot_file = test_status["tests"][system][test]["pseudopotentials"]["files"][element]
                    print("copy " + pseudopot_path + pseudopot_file + " to " + folder)
                    os.system("cp " + pseudopot_path + pseudopot_file + " " + folder)
                    # swap pseudopotential file name in input file
                    os.system("sed -i 's/{}_pseudopot/".format(element) + pseudopot_file + "/g' " + folder + "\\STRU")

                    nao_path = test_status["paths"]["numerical_orbital"]["resources"] + "\\"
                    nao_path += pseudopot_kva + "\\"
                    nao_et = str(ai.get_element_index(element)) + "_" + element + "_"
                    nao_et += test_status["tests"][system][test]["numerical_orbitals"]["info"][element]["type"]
                    nao_path += nao_et + "\\"
                    nao_file = test_status["tests"][system][test]["numerical_orbitals"]["files"][element]
                    print("copy " + nao_path + nao_file + " to " + folder)
                    os.system("cp " + nao_path + nao_file + " " + folder)
                    # swap numerical orbital file name in input file
                    os.system("sed -i 's/{}_numerical_orbital/".format(element) + nao_file + "/g' " + folder + "\\STRU")

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