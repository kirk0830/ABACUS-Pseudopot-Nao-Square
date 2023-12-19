import re

# COMPLETE
def _PP_HEADER_(data: str) -> dict:
    
    result = {
        "data": [],
        "attributes": {}
    }
    for line in data.split("\n"):
        if line.endswith("Version Number"):
            result["attributes"]["version"] = line.split()[0]
        elif line.endswith("Element"):
            result["attributes"]["element"] = line.split()[0]
        elif line.endswith("pseudopotential"):
            result["attributes"]["pseudopotential_type"] = line.split()[0]
        elif line.endswith("Nonlinear Core Correction"):
            result["attributes"]["nlcc"] = True if line.split()[0] == "T" else False
        elif line.endswith("Exchange-Correlation functional"):
            result["attributes"]["xc_functional"] = line.split()[-3]
            result["attributes"]["xc_functional_components"] = line.split()[:-3]
        elif line.endswith("Z valence"):
            result["attributes"]["z_val"] = int(float(line.split()[0]))
        elif line.endswith("Total energy"):
            result["attributes"]["total_energy"] = float(line.split()[0])
        elif line.endswith("Suggested cutoff for wfc and rho"):
            result["attributes"]["ecutwfc"] = float(line.split()[0])
            result["attributes"]["ecutrho"] = float(line.split()[1])
        elif line.endswith("Max angular momentum component"):
            result["attributes"]["l_max"] = int(line.split()[0])
        elif line.endswith("Number of points in mesh"):
            result["attributes"]["mesh"] = int(line.split()[0])
        elif line.endswith("Number of Wavefunctions, Number of Projectors"):
            result["attributes"]["number_of_wfc"] = int(line.split()[0])
            result["attributes"]["number_of_proj"] = int(line.split()[1])

    return result
# COMPLETE
def _PP_INPUTFILE_(data: str) -> dict:
    
    return {
        "data": [],
        "attributes": {}
    }
# COMPLETE
def _PP_INFO_(data: str) -> dict:
    
    possible_attribute_names = ["Author", "Generation date"]
    table_titles = ["nl", "pn", "l", "occ", "rcut", "rcut_us", "epseu"]

    result = {
        "data": [],
        "attributes": {}
    }
    for table_title in table_titles:
        result["attributes"][table_title] = []
    contents = data.split("\n")

    number_of_digit_line = 0

    for line in contents:
        line = line.strip()
        if "version" in line:
            version_number = ""
            words = line.split()
            read_version = False
            for word in words:
                if word == "version":
                    read_version = True
                else:
                    if read_version:
                        version_number += word
        else:
            number_of_attributes = line.count(":")
            print("number_of_attributes: ", number_of_attributes)
            if number_of_attributes:
                pass
            else:
                if line[0].isdigit():
                    number_of_digit_line += 1
                    if number_of_digit_line == 1:
                        if line == "1":
                            result["attributes"]["relativistic"] = "scalar"
                        else:
                            result["attributes"]["relativistic"] = "full"
                    elif number_of_digit_line == 2:
                        result["attributes"]["local_potential_cutoff_radius"] = float(line.split()[0])
                    else:
                        # read data of electronic configuration table...
                        words = line.split()
                        for iw in range(len(words)):
                            data = words[iw]
                            if iw == 1 or iw == 2:
                                data = int(data)
                            elif iw == 0:
                                pass
                            else:
                                data = float(data)
                            result["attributes"][table_titles[iw]].append(data)
    return result

def _PP_BETA_(data: str) -> dict:
    
    return {
        "data": [],
        "attributes": {}
    }

def _PP_DIJ_(data: str) -> dict:
    
    return {
        "data": [],
        "attributes": {}
    }

def _PP_QIJ_(data: str) -> dict:
    
    return {
        "data": [],
        "attributes": {}
    }