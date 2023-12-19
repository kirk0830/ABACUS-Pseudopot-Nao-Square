# COMPLETE
def recognize_line_state(line: str) -> str:
    """recognize if one line is a tag

    Args:
        line (str): line

    Returns:
        _type_: table: '# PSEUDOPOTENTIAL AND OPTIMIZATION'
                table_titles: '# atsym  z   nc   nv     iexc    psfile'
                delimiter: '#'
                attributes: 'author="anonymous"'
                data: '0   1.75000  -0.53265    3    8   5.30000'
                enclose_end: 'number_of_proj="6"/>'
                enclose_start: '<PP_HEADER'
                end: '</PP_NLCC>'
                comment: '<!-- END OF HUMAN READABLE SECTION -->'
                comment_start: '<!--'
                comment_end: '-->'
                start: '<PP_R type="real"  size="1358" columns="8">'
                start_end: 'cutoff_radius="    1.8300000000E+00" >'
    """
    # cases line NOT start with '<'
    if not line.startswith("<"):
        if not line.endswith(">"):
            if line.startswith("#"):
                if len(line) > 1:
                    _section_title = False
                    for iw, word in enumerate(line):
                        if word.isupper():
                            if line[iw-1:iw+4] != "(Ha)":
                                _section_title = True
                            break
                    if _section_title:
                        return "table"
                    else:
                        return "table_titles"
                else:
                    return "delimiter"
            else:
                if line.count("=") > 0:
                    return "attributes"
                else:
                    return "data"
        else:
            if line.endswith("/>"):
                return "enclose_end" # number_of_proj="6"/>
            elif line.endswith("-->"):
                return "comment_end" #                               -->
            else:
                return "start_end"
    # cases line start with '<'
    else:
        if not line.endswith(">"):
            return "enclose_start" # <PP_HEADER
        else:
            if line[1] == "/":
                return "end" # </PP_NLCC>
            elif line.startswith("<!--"):
                if line.endswith("-->"):
                    return "comment" # <!--                               -->
                else:
                    return "comment_start" # <!--
            else:
                return "start" # <PP_RHOATOM type="real"  size="1358" columns="4">
# COMPLETE
def convert_data(data: str):
    """convert data from string to float or int, if possible

    Args:
        data (str): data in string

    Raises:
        ValueError: unknown value format

    Returns:
        _type_: float or int
    """
    try:
        _data = float(data)
        ndot = data.count(".")
        if ndot == 1:
            return float(data)
        elif ndot == 0:
            return int(data)
        else:
            print("Read value: ", data)
            raise ValueError("Unexpected value read, please check")
    except ValueError:
        return data
# COMPLETE
def split_attributes(line: str) -> dict:
    """split attributes in one line

    Args:
        line (str): line

    Returns:
        dict: attributes
    """
    attributes = {}
    line = line.replace("/", "").replace(">", "")
    _expressions = [expression for expression in line.split() if not expression.startswith("<")]
    expressions = []
    for _expression in _expressions:
        if "=" in _expression:
            expressions.append(_expression)
        else:
            expressions[-1] += " " + _expression
    for expression in expressions:
        key, value = expression.split("=")
        attributes[key] = convert_data(value.replace('"', '').replace("'", ""))
    return attributes
# COMPLETE
def read_line_according_to_state(line: str, state: str, text: bool = False) -> dict:
    """according to the state of line, parse with different method

    Args:
        line (str): line contents
        state (str): state, defined in function recognize_line_state
        text (bool): for pure text data, will not split, useful for PP_INFO section reading

    Raises:
        ValueError: unknown value format

    Returns:
        dict: parsed contents, key can be 'table', 'table_titles', 'data', 'attributes'
    """
    result = {}
    if state == "table":
        result["table"] = line[2:]
    elif state == "table_titles":
        result["table_titles"] = [title for title in line.split()[1:] if title != "(Ha)"]
    elif state == "delimiter" or state == "end" or state == "comment":
        pass
    elif state == "data":
        if text:
            result["data"] = line
        else:
            result["data"] = [convert_data(data) for data in line.split()]
    elif state == "enclose_end" or state == "enclose_start" or state == "start" or state == "attributes" or state == "start_end":
        result["attributes"] = split_attributes(line)

    return result
# COMPLETE
def extract_tag(line: str) -> str:
    """extract tag from line

    Args:
        line (str): line read

    Returns:
        str: tag
    """
    return line.replace("<", "").replace("/", "").replace(">", "").split()[0]
# COMPLETE
def filter_to_get_tags(fname: str):
    """scan the whole file to get all tags and its state

    Args:
        fname (str): pseudopotential file name

    Returns:
        tuple: tags and corresponding states

    Note:
        there are replication in tags returned, helpful for determining tag connections
    """
    not_tags_tags = ["table", "table_titles", "delimiter", "attributes", "data",
                     "comment", "comment_start", "comment_end"]
    return_tags = []
    return_states = []

    eof = "</UPF>"
    with open(fname, 'r') as f:
        line = f.readline().strip()
        if not line.startswith("<UPF"):
            eof = "</PP_RHOATOM>"
        while line != eof:
            line = f.readline().strip()
            state = recognize_line_state(line)

            if state not in not_tags_tags:
                #print(line, state)
                if state != "start_end":
                    return_states.append(state)
                else:
                    return_states[-1] = "start"
                    continue
                if state != "enclose_end":
                    return_tags.append(extract_tag(line))
                else:
                    return_tags.append(return_tags[-1])
                #print(return_tags[-1], return_states[-1])
    return return_tags, return_states
# COMPLETE
def recoginize_connection_from_state(redundant_tags: list, states: list) -> dict:
    """Generate a dictionary to record the connection between tags like:
    <UPF version="2.0.1">
        <PP_INFO>
            <PP_HEADER>
                <PP_MESH/>
                <PP_RHOATOM/>
                <PP_NLCC/>
            </PP_HEADER>
        </PP_INFO>
    </UPF>
    Args:
        redundant_tags (list): everytime when a tag appears, no matter it is start or end, it will be recorded
        states (list): tag states, defined in function recognize_line_state

    Returns:
        dict: connection between tags
    """

    connection = {}
    if "UPF" not in redundant_tags:
        redundant_tags.insert(0, "UPF")
        redundant_tags.append("UPF")
        states.insert(0, "start")
        states.append("end")
    print(redundant_tags)
    for tag in redundant_tags:
        if tag not in connection.keys():
            connection[tag] = {"childen": [], "parents": ""}

    iteration_depth = 0
    _parent = ""
    for i in range(len(redundant_tags)):
        if states[i].endswith("start") and not states[i].startswith("comment"):
            iteration_depth += 1
            #print(iteration_depth)
            if i == 0:
                _parent = redundant_tags[i]
            else:
                connection[_parent]["childen"].append(redundant_tags[i])
                connection[redundant_tags[i]]["parents"] = _parent
                _parent = redundant_tags[i]
        elif states[i].endswith("end") and not states[i].startswith("comment"):
            iteration_depth -= 1
            _parent = connection[redundant_tags[i]]["parents"]
    assert iteration_depth == 0
    return connection
# COMPLETE
def tag_initialization(tag: str = "") -> dict:
    """initialize a tag

    Args:
        tag (str): tag name

    Returns:
        dict: initialized tag
    """

    return {
        "data": [],
        "attributes": {}
        }
# COMPLETE
def text_treatment(tag: dict, line: str) -> dict:
    """for tags marked as text_treatment, line will directly be appended to the data key

    Args:
        tag (dict): PP_INFO tag
        line (str): line read

    Returns:
        dict: PP_INFO tag
    """
    if len(tag["data"]) == 0:
        tag["data"].append(line + "\n")
    else:
        tag["data"][0] += line + "\n"
    return tag