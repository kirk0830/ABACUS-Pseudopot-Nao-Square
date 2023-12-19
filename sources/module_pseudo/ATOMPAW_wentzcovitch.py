from module_pseudo.upf_basic import read_line_according_to_state

def _PP_INFO_(data: str) -> dict:

    result = {
        "data": [],
        "attributes": {}
    }
    lines = data.split("\n")
    for il, line in enumerate(lines):
        if line.count("="):
            attributes = read_line_according_to_state(line, "attribute")
            for key in attributes.keys():
                result["attributes"][key] = attributes[key]
        else:
            result["data"].append(line)
