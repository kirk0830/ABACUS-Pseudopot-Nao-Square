import os

def _mkdir_(folder: str) -> None:
    """create a folder, if already exists, delete first

    Args:
        folder (str): folder name
    """
    if os.path.exists(folder):
        print("Warning: " + folder + " already exists, clean.")
        os.system("rm -rf " + folder)
    os.mkdir(folder)

def _sed_(folder: str,
          file_to_sed: str,
          functional: str,
          ecutwfc: float,
          **kwargs) -> None:
    """adjust input script

    Args:
        folder (str): folder of specific system-test
        functional (str): DFT functional for this test
        ecutwfc (float): ecutwfc for this test
    
    Extension:
        kwargs: to adjust more
    """
    # adjust functional in input file
    os.system("sed -i 's/functional_to_test/%s/g' %s/s"%(functional, folder, file_to_sed))
    # adjust ecutwfc in input file
    os.system("sed -i 's/ecutwfc_to_test/%s/g' %s/s"%(str(ecutwfc), folder, file_to_sed))

    for key in kwargs:
        variable_to_sed = key + "_to_test"
        os.system("sed -i 's/%s/%s/g %s/s' "%(variable_to_sed, kwargs[key], folder, file_to_sed))
