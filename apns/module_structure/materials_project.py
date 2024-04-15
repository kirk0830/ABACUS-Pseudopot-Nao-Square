from mp_api.client import MPRester
import apns.module_workflow.identifier as amwi
import os
def download(api_key: str, 
             formula: list[str], 
             n_structures: list[int], 
             stru_cache_dir: str = amwi.TEMPORARY_FOLDER, 
             **kwargs):
    """download structures to cif with given formula and number of structures"""
    print("Deprecated param: with_magmom, will always download magnetism information")
    result = {fo: None for fo in formula}
    # result is indiced by:
    # result[formula][impid] = (file, magmom)

    # make sure all filters in kwargs are list
    for key in kwargs.keys():
        if not isinstance(kwargs[key], list):
            kwargs[key] = [kwargs[key]] * len(formula)

    # define degradation behavior of two filters: theoretical and is_stable
    # therefore they two are compulsory
    lst_theoretical = kwargs.get("theoretical", [False] * len(formula))
    lst_is_stable = kwargs.get("is_stable", [True] * len(formula))
    # back-update kwargs
    kwargs.update({"theoretical": lst_theoretical, "is_stable": lst_is_stable})
    # filter now is (theoretical, is_stable), from the most strict to the most loose to be:
    FILTER_DEGRADATION = [(False, True), (False, False), (True, False)]
    # (False, True): only experimentally observed, the most stable
    # (False, False): experimentally observed, not needed to be the most stable
    # (True, False): theoretically predicted is also okay

    print("Establishing connection to Materials Project database...")
    with MPRester(api_key) as mpr:
        for ifo, fo in enumerate(formula):
            # number of structures to download
            ndocs = n_structures[ifo]
            # get filter information for present formula
            filter_ofeach = {_filter: kwargs[_filter][ifo] for _filter in kwargs.keys()}
            theoretical, is_stable = filter_ofeach["theoretical"], filter_ofeach["is_stable"]
            ifilterdeg = FILTER_DEGRADATION.index((theoretical, is_stable))
            while ifilterdeg < len(FILTER_DEGRADATION):
                print("Searching structures for formula {} with filter theoretical = {} and is_stable = {}".format(fo, theoretical, is_stable))
                docs = mpr.materials.summary.search(formula=fo, **filter_ofeach)
                # get mpid information
                mp_ids = [doc.material_id for doc in docs]
                if len(mp_ids) >= ndocs:
                    break
                else:
                    ifilterdeg += 1
                    filter_ofeach.update(dict(zip(["theoretical", "is_stable"], FILTER_DEGRADATION[ifilterdeg])))
            if len(mp_ids) < ndocs:
                raise ValueError("Not enough structure found for formula {}".format(fo))
            # structure information
            fcifs, magmoms = [], []
            for mp_id in mp_ids:
                structure = mpr.get_structure_by_material_id(mp_id)
                fcif = f"{mp_id}.cif"
                structure.to(filename=fcif, fmt="cif")
                os.system(f"mv {fcif} {stru_cache_dir}")
                fcifs.append("{}/{}".format(stru_cache_dir, fcif))
                magmom = mpr.magnetism.search(material_ids=mp_id)[0].magmoms
                magmoms.append(magmom)
            result[fo] = list(zip(fcifs, magmoms))
            
    return result
