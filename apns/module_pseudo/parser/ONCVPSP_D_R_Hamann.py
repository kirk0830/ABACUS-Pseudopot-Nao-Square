from apns.module_pseudo.upf_basic import read_line_according_to_state, recognize_line_state

# COMPLETE
def _PP_INPUTFILE_(data: str) -> dict:
    
    if data.endswith("\n"):
        data = data[:-1]
    result = {
        "data": [],
        "attributes": {}
    }
    lines = data.split("\n")
    states = [recognize_line_state(line) for line in lines]
    for il, line in enumerate(lines):
        state = states[il]
        if state == "table":
            table = read_line_according_to_state(line, state)["table"]
            result["data"].append(table)
            # initialize the __displacement__ key when a new table is read
            result["attributes"][table] = {"__displacement__": 0}

        if state == "table_titles":
            # if under the same table title, there are more than one table,
            # then the __displacement__ will be added a non-zero value, then
            # the key can be properly matched when importing data
            result["attributes"][table]["__displacement__"] += len(result["attributes"][table].keys()) - 1
            _list = [title.replace(",", "") for title in read_line_according_to_state(line, state)["table_titles"]]

            for title in _list:
                table = result["data"][-1]
                result["attributes"][table][title] = []

        if state == "data":
            _list = read_line_according_to_state(line, state)["data"]
            table = result["data"][-1]
            keys = [key for key in result["attributes"][table].keys() if key != "__displacement__"]
            _displacement = result["attributes"][table]["__displacement__"]
            # the _displacement variable is for the case like this:
            """
            # ATOM AND REFERENCE CONFIGURATION
            # atsym  z   nc   nv     iexc    psfile
            Er 68.00    9    5       4      upf
            #
            #   n    l    f        energy (Ha)
                1    0    2.00
                2    0    2.00
                2    1    6.00
                3    0    2.00
            """
            # there are not only one table under one capitalized title, and the 
            # columns do not necessarily have the same number of data
            for id, data in enumerate(_list):
                key = keys[id + _displacement]
                result["attributes"][table][key].append(data)
            continue
    # delete the key __displacement__
    for table in result["attributes"].keys():
        del result["attributes"][table]["__displacement__"]
    return result
