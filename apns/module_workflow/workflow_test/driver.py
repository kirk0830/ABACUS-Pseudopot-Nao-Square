"""version 1 driver is more well-designed, it seperate the whole workflow in a different way from version 0

version 1 considers more about "ITERATE" rather than "FLOW" in version 0.
Therefore version 1 will first do Cartesian direct product on calculation_settings and pseudopotential-
numerical atomic orbital settings. The first iteration layer is caculation_settings, the second iteration-
layer is pseudopotential-numerical atomic orbital settings. The third iteration layer is systems.

Therefore 
"""
def driver_v1(input_file: str):
    """new version of driver"""

    """initialize
    1. create cache directory, if not exist
    2. download structure from Materials Project to cache directory, the cif named as mp-xxx.cif
    3. translate input file to json and modify it
    4. scan valid pseudopotentials and numerical orbitals according to input file
    5. load pseudopotential archive
    6. load numerical orbital archive, not implemented yet
    7. check whether all valid pseudopotentials are compatible with the software

    return:
    1. input: the translated input json, following aspects are modified compared with user-input:  
              1. systems are changed to system with mpids  
              2. pseudopotentials and numerical_orbitals are expanded from list to dict  
              3. default values are set if not explicitly specified
    2. vpspot: valid pseudopotentials for all elements in input file, the first layers keys
                are elements, the second layers keys are "identifiers" of pseudopotentials,
                the third layers keys are "kind", "version", "appendix" and "file", the first
                three keys can concatenate as a "pseudopotential identifier".
    3. vnao: valid numerical orbitals for all elements in input file. The first layers
                keys are elements, the second are "pseudopotential identifier", the third
                are "numerical orbital identifier", the fourth are "type", "rcut", "appendix"
                and "file", the first three keys can concatenate as a "numerical orbital
                identifier".
    4. pspot_arch: pseudopotential archive, a dict whose keys are "identifiers" of pseudopotentials,
                    stored in `pseudo_dir`, values are folders' absolute path.
    5. nao_arch: numerical orbital archive, a dict whose keys are "numerical orbital identifier", stored in
                `nao_dir`, values are folders' absolute path.
    """
    # ---------------------------------------------------------------------------------------------
    # initialize, to get program running environment information, the "environment" includes
    # available pseudopotentials and numerical orbitals (but if return, it seems not good...)
    # does initialize really need to return in principle?
    # An initialization of the whole code, if want to return, should return the runtime information
    # , say all information needed for the code to run.
    # DESIGN: INITIALIZE CAN RETURN RUNTIME INFORMATION
    import apns.module_workflow.workflow_test.initialize as amwinit
    runtime_settings, vupfs, vorbs, upf_arch, orb_arch = amwinit.initialize(input_file)
    # ---------------------------------------------------------------------------------------------
    """iterate
    1. setup_iterables
        - create cross product of system on pseudopotential-numerical atomic orbital settings
        - on calculation settings
        - on extensive settings
    2. iterate"""
    import apns.module_workflow.workflow_test.apns_itertools as amwai
    # calculation_settings is about INPUT
    # extensive_settings is about "outer" loop, "outer" loop usually is about different systems
    # upforb_settings is about "inner" loop
    inner_iter, abacus_inputs, outer_iter = \
        amwai.setup_iterables(system_list=runtime_settings["systems"],
                              pseudopotentials=vupfs,                               # always available
                              numerical_orbitals=vorbs,                             # not always available
                              calculation_settings=runtime_settings["calculation"], # INPUT
                              extensive_settings=runtime_settings["extensive"])     # outer loop
    import apns.module_workflow.workflow_test.iterate as amwi
    folders = amwi.iterate(software=runtime_settings["global"]["software"].lower(),         # with software
                           systems=runtime_settings["systems"],                             # on systems
                           pseudopot_nao_settings=inner_iter,
                           calculation_settings=abacus_inputs,
                           extensive_settings=outer_iter,
                           valid_pseudopotentials=vupfs, valid_numerical_orbitals=vorbs,
                           pspot_archive=upf_arch,       nao_archive=orb_arch,
                           test_mode=False)
    """compress"""
    import time
    import os
    import apns.module_io.compress as amic
    amic.pack(folders, "apns_{}.zip".format(time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())))
    for folder in folders:
        pass
        os.system("rm -rf {}".format(folder))
    """submit?"""
    import apns.module_io.abacustest as amia
    amia.image_information(software=runtime_settings["global"]["software"].lower())
    # sorry I don't connect with abacustest yet
    """citation"""
    import apns.module_io.citation as amicite
    amicite.citation(software=runtime_settings["global"]["software"].lower())

def driver_v0(input_file: str):
    """old version of driver"""
    # Step 0: initialize
    import apns.module_workflow.workflow_test.initialize as init
    init.initialize_cache()
    # Step 1: input -> work_status
    import apns.module_workflow.workflow_test.old.to_work_status as tws
    work_status = tws.to(fname=input_file)
    # Step 2: work_status -> test_status
    import apns.module_workflow.workflow_test.old.to_test_status as tts
    test_status = tts.to(work_status=work_status)
    # Step 3: test_status -> test
    import apns.module_workflow.workflow_test.old.to_test as tt
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
    import apns.module_pseudo.archive as ampua
    ampua.archive(pseudo_dir=inp["global"]["pseudo_dir"], only_scan=False)

if __name__ == "__main__":
    driver_v1("input.json")
