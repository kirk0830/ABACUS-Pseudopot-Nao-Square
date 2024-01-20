def is_normal_end(fname: str) -> bool:
    """check if the job is normally ended for QE"""
    with open(fname, "r") as f:
        for line in f:
            if line.startswith("JOB DONE."):
                return True
    return False

def fname_setting(job_dir = "", calculation = "scf", suffix = "QE", include_path = False):
    """return the file names (stdout, log and input script) for QE
    The suffix parameter is not used here because QE use it to name the folder"""
    job_dir = job_dir if not job_dir.endswith("/") else job_dir[:-1]
    if include_path:
        return job_dir+"/out.log", job_dir+"/out.log", job_dir+"/%s.in"%calculation
    else:
        return "out.log", "out.log", "%s.in"%calculation