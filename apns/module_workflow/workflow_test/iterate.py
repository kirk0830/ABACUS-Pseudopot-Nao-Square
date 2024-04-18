import apns.module_workflow.identifier as amwi
import apns.module_software.abacus.generation as amsag
import apns.module_software.qespresso.generation as amsqg
import apns.module_structure.basic as amsb
import os

def iterate(**kwargs):
    software = kwargs.get("software", "abacus")
    systems = kwargs.get("systems", {}) # would be like {"formula1": [(fcif, magmom), ...], ...}
    extensive_settings = kwargs.get("extensive_settings", [])
    upforb_bundles = kwargs.get("upforb_bundles", [])
    calculation_settings = kwargs.get("calculation_settings", [])
    upfs = kwargs.get("upfs", {})
    orbs = kwargs.get("orbs", {})
    print(f"""
ITERATIVELY GENERATE TESTS
systems: {systems}
extensive_settings: {extensive_settings}
upforb_bundles: {upforb_bundles}
calculation_settings: {len(calculation_settings)} suites
software: {software}""")
    folders = []
    for ifo, formula in enumerate(systems.keys()): # iterate over structures
        elements = amsb.scan_elements(formula)
        for upforb_fo in upforb_bundles[ifo]: # iterate over pseudopotential and numerical orbital settings
            # spns: system_pseudopot_nao_setting
            for calculation_setting in calculation_settings: # iterate over calculation settings
                for extensive_setting in extensive_settings: # iterate over extensive settings
                    # make folder
                    num_orb = upforb_fo.get("numerical_orbital", [])
                    upforb_fo["numerical_orbital"] = num_orb
                    # create folder for each structure...
                    # in formula, there are several structures [(fcif, magmom), ..., (::AUTOGEN::DIMER, None), ...]
                    for istru in range(len(systems[formula])):
                        first = systems[formula][istru][0]
                        suffix = first.replace("\\", "/").split("/")[-1].split(".")[0]\
                                    if "::AUTOGEN::" not in first else first.split("::")[-1]
                        token = "_".join([formula, suffix]) # token specifies the structure
                        testcase = {"token": token, "fcif": first,
                                    "calculation_setting": calculation_setting,
                                    "extensive_setting": extensive_setting, "upforb_fo": upforb_fo,
                                    "upfs": upfs, "orbs": orbs, "magmom": systems[formula][istru][1],
                                    "elements": elements, "software": software}
                        folder = iterate_structure(**testcase)
                        folders.append(folder)
    return folders

def iterate_structure(**kwargs):
    token, fcif, calculation_setting, extensive_setting, upforb, upfs, orbs, magmom, elements, software = \
        kwargs.get("token", ""), kwargs.get("fcif", ""), kwargs.get("calculation_setting", {}),           \
        kwargs.get("extensive_setting", {}), kwargs.get("upforb_fo", {}), kwargs.get("upfs", {}),         \
        kwargs.get("orbs", {}), kwargs.get("magmom", None), kwargs.get("elements", []),                   \
        kwargs.get("software", "abacus")
        
    print("\nPrepare test for %s"%token.replace("_", " (")+")")
    folder = amwi.folder(system=token,
                         pseudo_nao_identifier=amwi.pseudopot_nao(**upforb).replace(".", "").replace("_", ""),
                         calculation_identifier=amwi.calculation(calculation_setting),
                         extensive_identifier=amwi.extensive(extensive_setting))
    os.makedirs(folder, exist_ok=True)
    # copy pseudopotential
    pseudopotentials = {}
    # upforb_bundle["pseudopotential"] arranges like [pseudopotential_identifier, ...]
    # indiced by index of element
    pspotids = upforb["pseudopotential"] # get pseudopotential identifiers of one combination
    for i in range(len(pspotids)): # iterate over all elements
        pspotid = pspotids[i] # for one element, get its pseudopotential identifier, however, which is this element?
        element = elements[i]
        fpseudo_withpath = upfs[element][pspotid]
        fpseudo = fpseudo_withpath.replace("\\", "/").split("/")[-1]
        os.system(f"cp {fpseudo_withpath} {folder}/{fpseudo}")
        pseudopotentials[element] = fpseudo
    # copy numerical orbital
    numerical_orbitals = {}
    naoids = upforb["numerical_orbital"] # get numerical orbital identifiers of one combination
    for i in range(len(naoids)):
        naoid = naoids[i] + "@" + pspotids[i]
        element = elements[i]
        fnao_withpath = orbs[element][naoid]
        fnao = fnao_withpath.replace("\\", "/").split("/")[-1]
        os.system("cp {} {}/{}".format(fnao_withpath, folder, fnao))
        numerical_orbitals[element] = fnao

    software_setting = {"fcif": fcif, "magmom": magmom,                                                     # define structure
                        "token": token, "target_folder": folder,                                            # define test
                        "calculation_setting": calculation_setting,"extensive_setting": extensive_setting,  # define parameters
                        "pseudopotentials": pseudopotentials, "numerical_orbitals": numerical_orbitals}     # define external files
    if software == "qespresso":
        iterate_qespresso(**software_setting)
    elif software == "abacus":
        iterate_abacus(**software_setting)
    
    return folder

def iterate_abacus(**kwargs) -> None:
    """iterate over all possible combinations of input parameters and generate folders
    especially for ABACUS

    Args:
        `system_with_mpid (str, optional)`: on which structure? system with mp-id. Defaults to "".
        `target_folder (str, optional)`: in which folder? target folder. Defaults to "./".
        `calculation_setting (dict, optional)`: use what calculation setting? Defaults to {}.
        `pseudopotentials (dict, optional)`: which set of pseudopotentials? Defaults to {}.
        `numerical_orbitals (dict, optional)`: which set of numerical orbitals? Defaults to {}.
        `characteristic_length (float, optional)`: how much distortion? Defaults to 0.0.
        `if_crystal (bool, optional)`: crystal or molecule? for crystal, distortion is cell scaling; for molecule, distortion is bond length. Defaults to True.
        `test_mode (bool, optional)`: for test. Defaults to True.
    """
    # define structure
    fcif = kwargs.get("fcif", "")
    magmom = kwargs.get("magmom", None)
    # define test
    token = kwargs.get("token", "")
    if token.count("_") != 1:
        raise ValueError("token should be like 'Er2O3_12345'")
    target_folder = kwargs.get("target_folder", "./")
    # define parameters
    calculation_setting = kwargs.get("calculation_setting", None)
    if calculation_setting is None:
        print("Warning: calculation_setting is empty, use default")
    extensive_setting = kwargs.get("extensive_setting", {})
    # define external files
    pseudopotentials = kwargs.get("pseudopotentials", None)
    if pseudopotentials is None:
        raise ValueError("pseudopotentials is empty")
    numerical_orbitals = kwargs.get("numerical_orbitals", None)
    if numerical_orbitals is None:
        print("Warning: numerical_orbitals is empty, ignore if do pw calculation.")

    # write INPUT
    _input = amsag.INPUT(calculation_setting)
    with open(target_folder + "/INPUT", "w") as f: f.write(_input)
    # write STRU
    numerical_orbitals = None if numerical_orbitals == [] else numerical_orbitals
    structure = token.split("_")[-1].lower()
    stru_param = {"fname": fcif, 
                  "shape": structure,
                  "pseudopotentials": pseudopotentials, "numerical_orbitals": numerical_orbitals, 
                  "cell_scaling": extensive_setting["characteristic_lengths"], 
                  "bond_length": extensive_setting["characteristic_lengths"], 
                  "starting_magnetization": magmom}
    
    if structure in ["dimer", "trimer", "tetramer"]:
        stru, cell = amsag.STRU_Molecule(**stru_param)
    elif structure in ["sc", "bcc", "fcc", "diamond"]:
        stru, cell = amsag.STRU_ACWFRef(**stru_param)
    else:
        stru, cell = amsag.STRU_Pymatgen(**stru_param)
    with open(target_folder + "/STRU", "w") as f: f.write(stru)
    # write KPT
    nks = extensive_setting["nkpoints_in_line"]
    kpt = amsag.KPT(isolated=(nks < 0), cell_parameters=cell)
    kpt = amsag.KLINE(fname=fcif, nkpts_in_line=nks) if nks > 0 else kpt
    with open(target_folder + "/KPT", "w") as f: f.write(kpt)

def iterate_qespresso(**kwargs) -> None:
    """iterate over all possible combinations of input parameters and generate folders
    especially for Quantum Espresso

    Args:
        `system_with_mpid (str, optional)`: on which structure? system with mp-id, or element with shape (e.g. Er_dimer). Defaults to "".
        `target_folder (str, optional)`: in which folder? target folder. Defaults to "./".
        `calculation_setting (dict, optional)`: use what calculation setting? Defaults to {}.
        `pseudopotentials (dict, optional)`: which set of pseudopotentials? Defaults to {}.
        `characteristic_length (float, optional)`: for non isolated, it is cell scaling; for isolated, it is bond length. Defaults to 0.0.
        `isolated (bool, optional)`: crystal or molecule? Defaults to False.
        `test_mode (bool, optional)`: for test. Defaults to True.
    
    Raises:
        `ValueError`: system_with_mpid should be like 'Er2O3_12345'
        `ValueError`: calculation_setting is empty
        `ValueError`: pseudopotentials is empty
    
    Returns:
        `None`
    """
    # define structure
    fcif = kwargs.get("fcif", "")
    magmom = kwargs.get("magmom", None)
    # define test
    system_with_mpid = kwargs.get("token", "")
    if system_with_mpid.count("_") != 1:
        raise ValueError("system_with_mpid should be like 'Er2O3_12345'")
    target_folder = kwargs.get("target_folder", "./")
    # define parameters
    calculation_setting = kwargs.get("calculation_setting", None)
    if calculation_setting is None:
        print("Warning: calculation_setting is empty, use default")
    extensive_setting = kwargs.get("extensive_setting", {})
    # define external files
    pseudopotentials = kwargs.get("pseudopotentials", None)
    if pseudopotentials is None:
        raise ValueError("pseudopotentials is empty")
    
    in_ = ""
    ntype, natom = len(pseudopotentials), 0
    if extensive_setting["nkpoints_in_line"] < 0:
        in_, natom = amsqg._ISOLATED(element=system_with_mpid.split("_")[0], 
                                     shape=system_with_mpid.split("_")[1],
                                     bond_length=extensive_setting["characteristic_lengths"])
    else:
        in_, natom = amsqg._CIF(fname=fcif,
                                cell_scaling=extensive_setting["characteristic_lengths"],
                                constraints=[])
    in_ = amsqg._K_POINTS(fname=fcif, 
                          nkpoints_in_line=extensive_setting["nkpoints_in_line"]) + in_
    in_ = amsqg._ATOMIC_SEPCIES(pseudopotential=pseudopotentials) + in_

    iterate_sections = ["control", "system", "electrons", "ions", "cell"]
    for section in iterate_sections:
        in_ = amsqg._calculation(section=section, ntype=ntype, natom=natom, **calculation_setting) + in_

    with open(target_folder + "/scf.in", "w") as f: f.write(in_)
    return None
