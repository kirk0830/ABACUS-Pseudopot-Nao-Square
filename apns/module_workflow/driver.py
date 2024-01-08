"""version 1 driver is more well-designed, it seperate the whole workflow in a different way from version 0

version 1 considers more about "ITERATE" rather than "FLOW" in version 0.
Therefore version 1 will first do Cartesian direct product on calculation_settings and pseudopotential-
numerical atomic orbital settings. The first iteration layer is caculation_settings, the second iteration-
layer is pseudopotential-numerical atomic orbital settings. The third iteration layer is systems.

Therefore 
"""
def driver_v1(input_file: str):
    """new version of driver"""

    """initialize"""
    import apns.module_workflow.initialize as amwinit
    inp, vpspot, vnao, pspot_arch, nao_arch = amwinit.initialize(input_file)
    """iteratively generation"""
    import apns.module_workflow.apns_itertools as amwai
    import apns.module_io.compress as amic
    import apns.module_workflow.iterate as amwi
    folders = amwi.iterate(software=inp["global"]["software"].lower(),
                           systems=inp["systems"],
                           pseudopot_nao_settings=amwai.systems(system_list=inp["systems"], 
                                                                valid_pseudopotentials=vpspot, 
                                                                valid_numerical_orbitals=vnao),
                           calculation_settings=amwai.calculation(inp["calculation"]),
                           extensive=inp["extensive"],
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
        pass
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
          cell_scalings=work_status["calculation"]["characteristic_lengths"])

def configure(input_file: str):
    """configure the apns storing files, only run this at the first time"""
    import json
    with open(input_file, "r") as f:
        inp = json.load(f)
    import apns.module_pseudo.upf_archive as ampua
    ampua.archive(pseudo_dir=inp["global"]["pseudo_dir"], only_scan=False)

if __name__ == "__main__":
    driver_v1("input.json")
