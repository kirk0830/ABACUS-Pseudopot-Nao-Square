"""this script is for calculating the basis set completeness"""

def collect_jobs(folder: str):
    """luckily, abacustest supports a reuse mode that can automatically
    complete the PW vs LCAO EOS test. The energy minimal can be used to
    measure the basis completeness. However, if the minimal given by any
    of PW-LCAO pair distinct too much, indicates the basis set is far from
    correct
    
    Args:
        folder (str): the folder that contains all the jobs
    

    """
