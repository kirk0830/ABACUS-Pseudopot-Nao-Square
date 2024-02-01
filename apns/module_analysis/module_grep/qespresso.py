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

def split_band(line: str) -> list:
    # first split by space:
    result = []
    words = line.split()
    for word in words:
        if(word.count("-")) > 1:
            subwords =  [f"-{x}" for x in word.split("-") if x]
            for subword in subwords:
                result.append(float(subword))
        else:
            result.append(float(word))
    return result

def grep_band(filename, nband = 0):
    """
    parse QE stdout
    INPUT: stdout redirected file name
    OUTPUT: band energies, weight of kpoints
    """
    read_band_section = False
    read_kpoint_info = False
    band_energies = []
    band_index = 0
    line_index = 0
    kpoint_index = 0
    nkpoints = 0
    wk = []
    nelec = 0
    efermi = 0
    with open(filename, 'r') as f:
        line = f.readline()
        line_index += 1
        # if line != EOF
        while(True):
            clean_line = line.strip()
            """
            Read mode switch
            """
            # reach the end, jump out
            if clean_line.startswith("JOB DONE."):
                break
            # the beginning report
            if clean_line.startswith("number of electrons"):
                words = clean_line.split()
                nelec = int(words[-1].split(".")[0])
            # the beginning report
            if clean_line.startswith("number of Kohn-Sham states"):
                words = clean_line.split()
                nband = int(words[-1])
            # band structure result begin
            if clean_line.startswith("End of self-consistent calculation"):
                read_band_section = True
                read_kpoint_info = False
                band_index = 0
                kpoint_index = -1
                band_energies.clear()
            # band structure result end
            if clean_line.startswith("the Fermi energy is"):
                read_band_section = False
                read_kpoint_info = False
                words = clean_line.split()
                efermi = float(words[-2])
                #print("efermi = ", efermi)
            # the beginning report
            if clean_line.startswith("number of k points="):
                words = clean_line.split()
                for word in words:
                    if word.isdigit():
                        nkpoints = int(word)
                        read_band_section = False
                        read_kpoint_info = True
                        break
            """
            Contents parsing of each mode
            """
            # BAND STRUCTURE RESULT PARSE
            if read_band_section and len(clean_line) > 1:
                if (clean_line.find("k =") != -1) and (clean_line.find("PWs") != -1):
                    # if line starts with "k =", it is the beginning of a new kpoint
                    kpoint_index += 1
                    pass
                else:
                    # if line starts with number, it is still of present kpoint
                    if clean_line[0].isdigit() or clean_line[0] == '-':
                        # read band energy and store in band_energies
                        energies_one_line = [float(e) for e in split_band(clean_line)]
                        for e in energies_one_line:
                            band_energies.append((e, wk[kpoint_index]))
                            band_index += 1

            # KPOINT INFO PARSE
            if read_kpoint_info and len(wk) < nkpoints:
                if clean_line.startswith("k("):
                    kpoint_index += 1
                    words = clean_line.split()
                    for index, word in enumerate(words):
                        if word.startswith("wk"):
                            wk.append(float(words[index+2]))
                        
            line = f.readline() # read next line
            line_index += 1
            if (line_index % 100) == 0:
                #print("line_index = ", line_index)
                pass

    return band_energies, nelec, nband, efermi