from apns.module_pseudo.upf_basic import *
import ONCVPSP_D_R_Hamann as ODRH
import GBRV_Vanderbilt as GBRV
import Atomic_A_Dal_Corso as AADC
import ATOMPAW_wentzcovitch as ATOMPAW

# COMPLETE
def pp_info_special_treatment(pp_info: dict, kind: str) -> dict:
    """PP_INFO special treatment for Atomic_A_Dal_Corso kind pseudopotential

    Args:
        pp_info (dict): PP_INFO tag

    Returns:
        dict: PP_INFO tag
    """
    contents = pp_info["data"][0]
    if kind == "ONCVPSP_D_R_Hamann":
        return pp_info
    elif kind == "GBRV_Vanderbilt":
        return GBRV._PP_INFO_(contents)
    elif kind == "Atomic_A_Dal_Corso":
        return AADC._PP_INFO_(contents)
    else:
        raise NotImplementedError
# COMPLETE
def pp_header_special_treatment(pp_header: dict, kind: str) -> dict:

    contents = pp_header["data"][0]
    if kind == "ONCVPSP_D_R_Hamann":
        return pp_header
    elif kind == "GBRV_Vanderbilt":
        return GBRV._PP_HEADER_(contents)
    elif kind == "Atomic_A_Dal_Corso":
        return pp_header
    else:
        raise NotImplementedError
# COMPLETE
def pp_inputfile_special_treatment(pp_inputfile: dict, kind: str) -> dict:

    contents = pp_inputfile["data"][0]
    if kind == "ONCVPSP_D_R_Hamann":
        return ODRH._PP_INPUTFILE_(contents)
    elif kind == "GBRV_Vanderbilt":
        return pp_inputfile
    elif kind == "Atomic_A_Dal_Corso":
        return AADC._PP_INPUTFILE_(contents)
    else:
        raise NotImplementedError
# COMPLETE
def special_treatment(data: dict, tag: str, kind: str) -> dict:
    """special treatment for some tags

    Args:
        data (dict): parsed data
        tag (str): tag name

    Returns:
        dict: parsed data
    """
    if tag == "PP_INFO":
        data = pp_info_special_treatment(data, kind)
    elif tag == "PP_HEADER":
        data = pp_header_special_treatment(data, kind)
    elif tag == "PP_INPUTFILE":
        data = pp_inputfile_special_treatment(data, kind)
    else:
        print("Requested parsed Tag: ", tag)
        raise NotImplementedError
    return data

def basic_parse(fname: str, text_treatment_tags: list = []):
    """parse the pseudopotential file

    Args:
        fname (str): pseudopotential file name
        text_treatment_tags (list): list of tag names to omit, make it possible for compatible with different kinds of pseudopotentials
    Returns:
        dict: parsed contents
        dict: tag connections
    Note:
        ONLY for GGA UPF, compatible pseudopotential kinds: SG15, PD04, DOJO.
    """

    tags, states = filter_to_get_tags(fname)
    connection = recoginize_connection_from_state(tags, states)
    attribute_states = ["enclose_end", "enclose_start", "start", "attributes", "start_end"]
    # attribute_states records states possibly have attributes in it

    _present = ""
    result = {}
    eof = "</UPF>"
    with open(fname, 'r') as f:
        line = f.readline().strip()
        if not line.startswith("<UPF"):
            eof = "</PP_RHOATOM>"
        while line != eof:
            line = f.readline().strip()
            state = recognize_line_state(line)
            
            # tag connection structure determine
            if state.endswith("start"):
                tag = extract_tag(line)
                result[tag] = tag_initialization(tag)
                _present = tag
            elif state.endswith("end") and state != "enclose_end" and state != "start_end":
                _present = connection[_present]["parents"]
                # for enclose_end, there will also be attribute states, read that first, then jump back
                # for start_end, it is not the real end of a tag, so do nothing
            # parse tag contents
            if _present in text_treatment_tags:
                result[_present] = text_treatment(result[_present], line)
            else:
                if state in attribute_states:
                    _key_value_pairs = read_line_according_to_state(line, state)
                    for attribute in _key_value_pairs["attributes"].keys():
                        result[_present]["attributes"][attribute] = _key_value_pairs["attributes"][attribute]

                else:
                    if state == "data":
                        _list = read_line_according_to_state(line, state)["data"]
                        for data in _list:
                            result[_present]["data"].append(data)
                    else:
                        pass
                # after attributes read, jump back
                if state == "enclose_end":
                    _present = connection[_present]["parents"] # pointer jump-back

    return result, connection
# 
def ONCVPSP_D_R_Hamann(fname: str):

    text_treatment_tags = ["PP_INFO", "PP_INPUTFILE"]
    data, connection = basic_parse(fname, text_treatment_tags)

    for tag in text_treatment_tags:
        data[tag] = special_treatment(data[tag], tag, "ONCVPSP_D_R_Hamann")
    return data, connection
#
def GBRV_Vanderbilt(fname: str):

    text_treatment_tags = ["PP_HEADER", "PP_INFO"]
    data, connection = basic_parse(fname, text_treatment_tags)

    for tag in text_treatment_tags:
        data[tag] = special_treatment(data[tag], tag, "GBRV_Vanderbilt")
    return data, connection
#
def Atomic_A_Dal_Corso(fname: str):

    pass

# COMPLETE
def determine_pseudopotential_kind(fname: str) -> str:

    line = "start"
    with open(fname, 'r') as f:
        while line != "</UPF>":
            line = f.readline().strip()
            if line.startswith("Generated using \"atomic\" code by A. Dal Corso"):
                return "Atomic_A_Dal_Corso"
            elif line.startswith("Generated using Vanderbilt code"):
                return "GBRV_Vanderbilt"
            elif line.startswith("ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)"):
                return "ONCVPSP_D_R_Hamann"
    return "unknown"
# COMPLETE
def read(fname: str):

    kind = determine_pseudopotential_kind(fname)
    if kind == "Atomic_A_Dal_Corso":
        return Atomic_A_Dal_Corso(fname)
    elif kind == "GBRV_Vanderbilt":
        return GBRV_Vanderbilt(fname)
    elif kind == "ONCVPSP_D_R_Hamann":
        return ONCVPSP_D_R_Hamann(fname)
    else:
        raise NotImplementedError

if __name__ == "__main__":

    import json

    Directory = 'D:/Documents/GitHub/abacus-develop_(Gitee)/tests/PP_ORB/'
    FileName = 'Al.pz-vbc.UPF'

    result, connection = ONCVPSP_D_R_Hamann(Directory + FileName)
    # save to json
    with open("result.json", 'w') as f:
        json.dump(result, f, indent=4)
    with open("connection.json", 'w') as f:
        json.dump(connection, f, indent=4)