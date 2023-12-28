"""Read test_status and, generate template input 
scripts and then sed to get input files"""

from apns.module_structure.basic import scan_elements
import apns.module_software.abacus.generation as abacus
import apns.module_software.qespresso.generation as qe
import os
import apns.module_io.output_pw as op
import apns.module_io.output_lcao as ol
import apns.module_workflow.identifier as id

def _abacus_(test_status: dict, basis_type: str):
    """Iteratively generate template files for each system
    
    Args:
        test_status (dict)
        basis_type (str): "pw" or "lcao"
    Returns:
        None
    """
    for system in test_status["tests"].keys():
        # STRU file
        _cif_name = "mp-" + system.split("_")[-1] + ".cif"
        _elements = scan_elements(system=system)
        _pseudopots = {}
        for _element in _elements:
            _pseudopots[_element] = ""
        if basis_type == "lcao":
            _numer_orbs = {}
            for _element in _elements:
                _numer_orbs[_element] = ""
            _stru, _cell = abacus.STRU_cif(fname=id.TEMPORARY_FOLDER + "/" + _cif_name,
                                           pseudopotentials=_pseudopots,
                                           numerical_orbitals=_numer_orbs,
                                           template=True)
        else:
            _stru, _cell = abacus.STRU_cif(fname=id.TEMPORARY_FOLDER + "/" + _cif_name,
                                           pseudopotentials=_pseudopots,
                                           template=True)
        with open(id.TEMPORARY_FOLDER + "/" + "STRU_"+system, "w") as stru_f:
            stru_f.writelines(_stru)
        # INPUT file
        _input = abacus.INPUT_generation()
        with open(id.TEMPORARY_FOLDER + "/" + "INPUT_"+system, "w") as input_f:
            input_f.writelines(_input)
        # KPT file
        _kpt = abacus.KPT_generation(mode="crystal",
                                     gamma_centered=True,
                                     cell_parameters=_cell)
        with open(id.TEMPORARY_FOLDER + "/" + "KPT_" + system, "w") as kpt_f:
            kpt_f.writelines(_kpt)

def _qespresso_(test_status):
    """Iteratively generate template input scripts of Quantum ESPRESSO
    
    Args:
        test_status (dict)
        
    Returns:
        None
    """
    for system in test_status["tests"].keys():
        _cif_name = "mp-" + system.split("_")[-1] + ".cif"
        fname = qe.cif_to_qespresso(cif_file=id.TEMPORARY_FOLDER + "/" + _cif_name)
        os.rename(fname, id.TEMPORARY_FOLDER + "/" + "qespresso_"+system+".in")
        # this name is required by module_workflow.identifier.qespresso when template = True

def prepare(test_status: dict,
            software: str,
            basis_type: str):
    if software == "ABACUS":
        _abacus_(test_status=test_status,
                 basis_type=basis_type)
    elif software == "qespresso":
        _qespresso_(test_status=test_status)

def to(test_status: dict, 
       software: str = "ABACUS", 
       basis_type: str = "pw",
       functionals: list = ["pbe"],
       cell_scalings: list = [0.0]):

    #ob._mkdir_(id.TEMPORARY_FOLDER)
    # Generate template input scripts
    prepare(test_status=test_status,
            software=software,
            basis_type=basis_type)
    # Generate input files
    if basis_type == "pw":
        op.generate(test_status=test_status,
                    software=software,
                    functionals=functionals,
                    cell_scalings=cell_scalings)
    elif basis_type == "lcao":
        ol.generate(test_status=test_status,
                    functionals=functionals,
                    cell_scalings=cell_scalings)
    else:
        print("Current basis_type: ", basis_type)
        raise ValueError("Not supported basis_type.")
