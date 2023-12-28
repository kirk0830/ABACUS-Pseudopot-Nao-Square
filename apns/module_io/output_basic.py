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
          **kwargs) -> None:
    """adjust input script

    Args:
        folder (str): folder of specific system-test
        functional (str): DFT functional for this test
        ecutwfc (float): ecutwfc for this test
    
    Extension:
        kwargs: to adjust more
    """

    for key in kwargs.keys():
        variable_to_sed = key + "_to_test"
        os.system("sed -i 's/%s/%s/g' %s/%s "%(variable_to_sed, str(kwargs[key]), folder, file_to_sed))
