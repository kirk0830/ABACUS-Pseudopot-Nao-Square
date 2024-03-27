import re
def is_numeric_data(data):
    """judge if the data line is full of numbers (including scientific notation) separated by spaces, tabs or newlines"""
    if re.match(r"^\s*[+-]?\d+(\.\d*)?([eE][+-]?\d+)?(\s+[+-]?\d+(\.\d*)?([eE][+-]?\d+)?)*\s*$", data):
        return True
    else:
        return False

def decompose_data(data):
    """to decompose all numbers in one line, but need to judge whether int or float"""
    if re.match(r"^\s*[+-]?\d+(\.\d*)([eE][+-]?\d+)?(\s+[+-]?\d+(\.\d*)([eE][+-]?\d+)?)*\s*$", data):
        return [float(x) for x in data.split()] if " " in data else float(data.strip())
    elif re.match(r"^\s*[+-]?\d+(\s+[+-]?\d+)*\s*$", data):
        return [int(x) for x in data.split()] if " " in data else int(data.strip())
    else:
        print(data)
        raise ValueError("data is not numeric")