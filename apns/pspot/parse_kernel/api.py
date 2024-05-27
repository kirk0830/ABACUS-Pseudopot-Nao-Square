"""this file provides interface to several kinds of kernel parsers supporting different pseudopotential formats

1. UPF (XML)
2. psp8
"""

def determine_format(fname: str) -> str:
    """a function for determine the format of pseudopotential how organizing information
    currently this function is not programmed yet"""
    
    return "xml"

import apns.pspot.parse_kernel.xml as ampkx
import apns.pspot.parse_kernel.psp8 as ampp8
def parse(fname: str) -> dict:
    """parse pseudopotential file, presently only support UPF (XML)"""
    if determine_format(fname) == "xml":
        return ampkx.parse(fname)
    if determine_format(fname) == "psp8":
        raise NotImplementedError("psp8 format is not supported yet")
    else:
        raise ValueError("Pseudopotential format not recognized")
    
def convert(fname: str, target_format: str):
    """convert pseudopotential file to target format"""
    format0 = determine_format(fname)
    if format0 == target_format:
        return fname
    elif format0 == "xml":
        if target_format == "psp8":
            ampkx.to_psp8(fname)
        else:
            raise ValueError("Target format not recognized")
    elif format0 == "psp8":
        if target_format == "xml":
            ampp8.to_xml(fname)
        else:
            raise ValueError("Target format not recognized")