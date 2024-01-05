
def driver_v1(input_file: str):
    """new version of driver"""

    """initialize"""
    import apns.module_workflow.initialize as amwinit
    inp, vpspot, vnao, pspot_arch, nao_arch = amwinit.initialize(input_file)
    """iteratively generation"""
    import apns.module_workflow.iterate as amwi
    import apns.module_workflow.apns_itertools as amwai
    import apns.module_structure.basic as amsb
    import apns.module_io.compress as amic
    folders = amwi.iterate(systems=inp["systems"],
                           pseudopot_nao_settings=amwai.system(
                               elements=[
                                   amwai.pseudopot_nao(
                                   list(vpspot[element].keys())) for element in amsb.scan_elements(inp["systems"])
                                   ]
                           ),
                           calculation_settings=amwai.calculation(inp["calculation"]),
                           cell_scalings=inp["calculation"]["cell_scaling"],
                           valid_pseudopotentials=vpspot,
                           valid_numerical_orbitals=vnao,
                           pspot_archive=pspot_arch,
                           nao_archive=nao_arch,
                           test_mode=False)
    """compress"""
    import time
    import os
    amic.pack(folders, "apns_{}.zip".format(time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())))
    for folder in folders:
        os.system("rm -rf {}".format(folder))
    """submit?"""
    # sorry I don't connect with abacustest yet

def driver_v0(input_file: str):
    """old version of driver"""
    # Step 0: initialize
    import apns.module_workflow.initialize as init
    init.initialize_cache()
    # Step 1: input -> work_status
    import apns.module_workflow.to_work_status as tws
    work_status = tws.to(fname=input_file)
    # Step 2: work_status -> test_status
    import apns.module_workflow.to_test_status as tts
    test_status = tts.to(work_status=work_status)
    # Step 3: test_status -> test
    import apns.module_workflow.to_test as tt
    tt.to(test_status=test_status,
          software=work_status["global"]["software"],
          basis_type=work_status["calculation"]["basis_type"],
          functionals=work_status["calculation"]["functionals"],
          cell_scalings=work_status["calculation"]["cell_scaling"])

if __name__ == "__main__":
    driver_v1("input.json")
