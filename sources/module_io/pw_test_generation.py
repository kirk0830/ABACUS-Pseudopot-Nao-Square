import os

def generate_pw_pseudopotential_test(test_status: dict):

    if test_status["global"]["software"] == "Quantum ESPRESSO":
        gppt_quantum_espresso(test_status)
    elif test_status["global"]["software"] == "ABACUS":
        gppt_abacus(test_status)
    else:
        print("Error: software not recognized.")
        exit(1)

def gppt_quantum_espresso(test_status: dict):

    os.chdir(test_status["paths"]["work_folder"])
    for functional in test_status["calculation"]["functionals"]:
        for system in test_status["tests"].keys():
            for test in test_status["tests"][system].keys():
                folder = "t_" + functional + "_" + test + "_" + system
                if os.path.exists(folder):
                    print("Warning: " + folder + " already exists, clean.")
                    os.system("rm -rf " + folder)
                os.mkdir(folder)
                input_file_path = test_status["paths"]["pseudopotential"]["work"] + "\\templates\\Quantum_ESPRESSO\\template_" + system + ".in"
                # copy input file to test folder, use unified name scf.in
                os.system("cp " + input_file_path + " " + folder + "\\scf.in")
                # adjust functional in input file
                os.system("sed -i 's/functional_to_test/" + functional + "/g' " + folder + "\\scf.in")
                # adjust ecutwfc in input file
                os.system("sed -i 's/ecut_to_test/" + str(test_status["tests"][system][test]["ecutwfc"]) + "/g' " + folder + "\\scf.in")
                # copy pseudopotential to test folder
                for element in test_status["tests"][system][test]["elements"]:
                    pseudopot_path = test_status["paths"]["pseudopotential"]["resources"] + "\\"
                    pseudopot_path += test_status["tests"][system][test]["pseudopotentials"]["info"][element]["kind"]
                    if test_status["tests"][system][test]["pseudopotentials"]["info"][element]["version"] != "":
                        pseudopot_path += "_" + test_status["tests"][system][test]["pseudopotentials"]["info"][element]["version"]
                    if test_status["tests"][system][test]["pseudopotentials"]["info"][element]["appendix"] != "":
                        pseudopot_path += "_" + test_status["tests"][system][test]["pseudopotentials"]["info"][element]["appendix"]
                    pseudopot_path += "\\"
                    pseudopot_file = test_status["tests"][system][test]["pseudopotentials"]["files"][element]
                    print("copy " + pseudopot_path + pseudopot_file + " to " + folder)
                    os.system("cp " + pseudopot_path + pseudopot_file + " " + folder)
                    # swap pseudopotential file name in input file
                    os.system("sed -i 's/{}_pseudopot/".format(element) + pseudopot_file + "/g' " + folder + "\\scf.in")
    os.chdir(test_status["paths"]["root"])
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

def gppt_abacus(test_status: dict):

    os.chdir(test_status["paths"]["work_folder"])
    for functional in test_status["calculation"]["functionals"]:
        for system in test_status["tests"].keys():
            for test in test_status["tests"][system].keys():
                folder = "t_" + functional + "_" + test + "_" + system
                if os.path.exists(folder):
                    print("Warning: " + folder + " already exists, clean.")
                    os.system("rm -rf " + folder)
                os.mkdir(folder)
                template_files_path = test_status["paths"]["pseudopotential"]["work"] + "\\templates\\ABACUS\\"
                stru_file = template_files_path + "STRU_" + system
                input_file = template_files_path + "INPUT"
                kpt_file = template_files_path + "KPT"
                # copy STRU_* to test folder, use unified name STRU
                os.system("cp " + stru_file + " " + folder + "\\STRU")
                # copy INPUT to test folder, use unified name INPUT
                os.system("cp " + input_file + " " + folder + "\\INPUT")
                # copy KPT to test folder, use unified name KPT
                os.system("cp " + kpt_file + " " + folder + "\\KPT")
                # adjust functional in input file
                os.system("sed -i 's/functional_to_test/" + functional + "/g' " + folder + "\\INPUT")
                # adjust ecutwfc in input file
                os.system("sed -i 's/ecut_to_test/" + str(test_status["tests"][system][test]["ecutwfc"]) + "/g' " + folder + "\\INPUT")
                # copy pseudopotential to test folder
                for element in test_status["tests"][system][test]["elements"]:
                    pseudopot_path = test_status["paths"]["pseudopotential"]["resources"] + "\\"
                    pseudopot_path += test_status["tests"][system][test]["pseudopotentials"]["info"][element]["kind"]
                    if test_status["tests"][system][test]["pseudopotentials"]["info"][element]["version"] != "":
                        pseudopot_path += "_" + test_status["tests"][system][test]["pseudopotentials"]["info"][element]["version"]
                    if test_status["tests"][system][test]["pseudopotentials"]["info"][element]["appendix"] != "":
                        pseudopot_path += "_" + test_status["tests"][system][test]["pseudopotentials"]["info"][element]["appendix"]
                    pseudopot_path += "\\"
                    pseudopot_file = test_status["tests"][system][test]["pseudopotentials"]["files"][element]
                    print("copy " + pseudopot_path + pseudopot_file + " to " + folder)
                    os.system("cp " + pseudopot_path + pseudopot_file + " " + folder)
                    # swap pseudopotential file name in input file
                    os.system("sed -i 's/{}_pseudopot/".format(element) + pseudopot_file + "/g' " + folder + "\\STRU")

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
