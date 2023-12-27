"""initialize is the module should be called everytime will the program starts"""
import os
import apns.module_workflow.identifier as id

def initialize_cache() -> None:

    print("Current working directory: {}".format(os.getcwd()))
    """change id.TEMPORARY_FOLDER to absolute path"""
    id.TEMPORARY_FOLDER = os.path.join(os.getcwd(), id.TEMPORARY_FOLDER)
    """create cache directory if not exist"""
    if not os.path.exists(id.TEMPORARY_FOLDER):
        os.mkdir(id.TEMPORARY_FOLDER)
    