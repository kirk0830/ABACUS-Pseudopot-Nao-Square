def image_information(software: str = "abacus"):
    # software: ABACUS, qespresso
    if software == "abacus":
        abacus()
    elif software == "qespresso":
        qespresso()

def abacus():

    template(software="ABACUS",
             command="OMP_NUM_THREADS=16 mpirun -np 1 abacus | tee out.log ",
             additional_information="(especially when dft_functional = HSE, arbitrary if PBE)",
             image="registry.dp.tech/deepmodeling/abacus-intel:latest",
             machine_type="c32_m128_cpu")

def qespresso():

    template(software="qespresso",
             command="mpirun -np 16 pw.x -i scf.in | tee scf.log; rm -rf pwscf.save",
             additional_information="",
             image="registry.dp.tech/dptech/prod-471/abacus-vasp-qe:20230116",
             machine_type="c32_m128_cpu")

def template(software: str = "ABACUS",
             command: str = "OMP_NUM_THREADS=16 mpirun -np 1 abacus | tee out.log ",
             additional_information: str = "",
             image: str = "registry.dp.tech/deepmodeling/abacus-intel:latest",
             machine_type: str = "c32_m128_cpu"):
    result = TEMPLATE.replace("_software_", software)
    result = result.replace("_commands_", command)
    if additional_information == "":
        result = result.replace("_additional_information_\n", "")
    else:
        result = result.replace("_additional_information_", additional_information)
    result = result.replace("_image_", image)
    result = result.replace("_machine_type_", machine_type)
    print(result)

TEMPLATE = """
GENERATION DONE
To run _software_ tests on abacustest, use the following parameters:
rundft: _commands_ 
_additional_information_
Bohrium image: _image_
Bohrium_machine_type: _machine_type_
Bohrium_platform: ali
You can download your job results by: lbg jobgroup download <jobgroup_id>"""
