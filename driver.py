import apns.module_workflow.initialize as amwinit
import apns.module_workflow.iterate as amwi
import apns.module_workflow.apns_itertools as amwai

def _driver_(input_file: str):
    """new version of driver"""
    inp, vpspot, vnao, pspot_arch, nao_arch = amwinit.initialize(input_file)
    folders = amwi.iterate(systems=inp["systems"],
                           pseudopot_nao_settings=amwai.system(
                               elements=[amwai.pseudopot_nao(
                                   list(vpspot[element].keys())) for element in inp["systems"]]
                           ),
                           calculation_settings=amwai.calculation(inp["calculation"]),
                           cell_scalings=inp["calculation"]["cell_scaling"],
                           valid_pseudopotentials=vpspot,
                           valid_numerical_orbitals=vnao,
                           pspot_archive=pspot_arch,
                           nao_archive=nao_arch)
    print(folders)

_driver_("input.json")

"""
import apns.module_structure.basic as amsb
import apns.module_pseudo.local_validity_scan as amplvs
import apns.module_io.input_translate as amiwse
import apns.module_software.abacus.generation as amsag
import apns.module_workflow.identifier as amwi
import apns.module_workflow.apns_itertools as amwai
import apns.module_pseudo.upf_archive as ampua
import os

import apns.module_workflow.initialize as amwinit

def driver(input_file: str):

    input = amiwse.inp_translate(fname=input_file, system_with_mpids={"Cr": ["Cr_90"]})

    elements = amsb.scan_elements(input["systems"])
    valid_pseudopotentials = amplvs._svp_(elements, input["pseudopotentials"])
    for element in elements:
        if element not in valid_pseudopotentials.keys():
            raise ValueError("No valid pseudopotential for element {}.".format(element))
        
    pseudopot_nao_settings = []
    if input["calculation"]["basis_type"] == "lcao":
        raise NotImplementedError("lcao calculation is not supported yet.") # NOQA
    else:
        pseudopot_nao_settings = amwai.system(
            elements=[
                amwai.pseudopot_nao(list(valid_pseudopotentials[element].keys())) for element in elements
            ]
        )
    
    pseudopot_arch = ampua.load(input["global"]["pseudo_dir"])

    calculation_settings = amwai.calculation(input["calculation"])

    # to iterate: input["systems"], pseudopot_nao_settings, calculation_settings, input["calculation"]["cell_scaling"]

    for _system in input["systems"]:
        for system_pseudopot_nao_setting in pseudopot_nao_settings:
            for calculation_setting in calculation_settings:
                for cell_scaling in input["calculation"]["cell_scaling"]:
                    # make folder
                    folder = amwi._folder_(_system, amwi.pseudopot_nao(system_pseudopot_nao_setting["pseudopotential"]), amwi.calculation(calculation_setting))
                    folder = amwi.foldr_reduce(folder)
                    os.makedirs(folder, exist_ok=True)
                    # write INPUT
                    _input = amsag.INPUT(calculation_setting)
                    with open(folder + "/INPUT", "w") as f:
                        f.write(_input)
                    # write pseudopotential
                    _elements = amsb.scan_elements(_system)
                    _pnids = system_pseudopot_nao_setting["pseudopotential"]
                    for _pnid in _pnids:
                        for _element in _elements:
                            fpseudo = pseudopot_arch[_pnid] + "/"
                            fpseudo += valid_pseudopotentials[_element][_pnid]["file"]
                            os.system("cp {} {}".format(fpseudo, folder))
                    # write STRU
                    _stru, _cell = amsag._STRU_(
                        fname=amwi.TEMPORARY_FOLDER + "/" + amwi.cif(_system), 
                        pseudopotentials={"Cr": valid_pseudopotentials[_element][_pnid]["file"]},
                        cell_scaling=cell_scaling)
                    with open(folder + "/STRU", "w") as f:
                        f.write(_stru)
                    # write KPT
                    _kpt = amsag.KPT_generation("crystal", cell_parameters=_cell)
                    with open(folder + "/KPT", "w") as f:
                        f.write(_kpt)

"""