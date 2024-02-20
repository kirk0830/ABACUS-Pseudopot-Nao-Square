"""at very first, configure pseudopotentials, sort them and archive, create description.json in each folder"""
import apns.module_pseudo.archive as arch
import json

def _pseudopotential_download_(pseudo_dir: str = "./download/pseudopotentials"):
    """download pseudopotentials from websites provided in module_pseudo/download/resources.json
    NO RECENT PLAN TO IMPLEMENT YET.
    """

def _pseudopotential_archive_(pseudo_dir: str = "./download/pseudopotentials"):
    result = arch.archive(pseudo_dir=pseudo_dir,
                          only_scan=False)
    # statistics
    nelements = len(result)
    npseudos_each_element = [len(result[element]) for element in result.keys()]
    npseudos = sum(npseudos_each_element)
    print("Total number of elements: {}".format(nelements))
    print("Total number of pseudopotentials: {}".format(npseudos))
    print("Number of pseudopotentials for each element: {}".format(npseudos_each_element))
    return result