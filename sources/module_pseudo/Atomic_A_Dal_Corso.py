from module_pseudo.upf_basic import read_line_according_to_state

# COMPLETE
def _PP_INPUTFILE_(data: str) -> dict:
    """
    convert the <PP_INPUTFILE>["data"][0] to the dict that can be assign to PP_INPUTFILE
    """

    result = {
        "data": [],
        "attributes": {}
    }
    lines = data.split("\n")
    for line in lines:
        if line.count("="):
            attributes = read_line_according_to_state(line, "attribute")
            for key in attributes.keys():
                result["attributes"][key] = attributes[key]
        """
        we omit the reference states definition like:
        3
        4S  1  0  2.00  0.00  2.10  2.10  0.0
        4P  2  1  3.00  0.00  2.10  2.10  0.0
        4D  3  2  0.00  0.10  2.10  2.10  0.0
        """
    return result

def _PP_INFO_(pp_info: dict) -> dict:
    """PP_INFO special treatment for Atomic_A_Dal_Corso kind pseudopotential

    Args:
        pp_info (dict): PP_INFO tag

    Returns:
        dict: PP_INFO tag
    """
    contents = pp_info["data"][0]
    words = [word for word in contents.split() if word != " "]

    search_table = False
    read_table = False
    table = []

    valence_configuration_table = False
    generation_configuration_table = False
    for iw, word in enumerate(words):
        if word == "Author:":
            pp_info["attributes"]["author"] = words[iw+1]
        elif word == "date:" and words[iw-1] == "Generation":
            pp_info["attributes"]["generation_date"] = words[iw+1]
        elif word == "type:" and words[iw-1] == "Pseudopotential":
            pp_info["attributes"]["pseudopotential_type"] = words[iw+1]
        elif word == "Element:":
            pp_info["attributes"]["element"] = words[iw+1]
        elif word == "wavefunctions:":
            pp_info["attributes"]["ecutwfc"] = float(words[iw+1])
        elif word == "density:" and words[iw-1] == "charge":
            pp_info["attributes"]["ecutrho"] = float(words[iw+1])
        elif word.endswith("Relativistic"):
            rel, _ = word.split("-")
            if rel == "Scalar":
                pp_info["attributes"]["relativistic"] = False
            else:
                pp_info["attributes"]["relativistic"] = True
        elif word == "Valence" and words[iw+1] == "configuration:":
            search_table = True
            valence_configuration_table = True
        elif word == "Generation" and words[iw+1] == "configuration:":
            # end of Valence configuration table
            generation_configuration_table = True
            pp_info["attributes"]["valence_configuration"] = table
            table = []
            read_table = False
            search_table = True
        elif word == "Pseudization" and words[iw+1] == "used:":
            # end of Generation configuration table
            pp_info["attributes"]["generation_configuration"] = table
            search_table = False
            read_table = False
        else:
            if search_table:
                if word[0].isdigit():
                    read_table = True
                    search_table = False
                    for i in range(7):
                        table.append([])
                    table[0].append(word)
            elif read_table:
                _len = min([len(column) for column in table])
                for column in table:
                    if _len == len(column):
                        column.append(word)
                        break

    # reorganize two tables
    tables_to_reorganize = ["valence_configuration", "generation_configuration"]
    for i, b in enumerate([valence_configuration_table, generation_configuration_table]):
        if b:
            _table = pp_info["attributes"][tables_to_reorganize[i]]
            pp_info["attributes"][tables_to_reorganize[i]] = {
                "nl": [], "pn": [], "l": [], "occ": [], "Rcut": [], "Rcut_US": [], "E_pseu": []
            }
            for i in range(len(_table[0])):
                for ik, key in enumerate(pp_info["attributes"]["valence_configuration"].keys()):
                    pp_info["attributes"]["valence_configuration"][key].append(_table[ik][i])

    return pp_info
