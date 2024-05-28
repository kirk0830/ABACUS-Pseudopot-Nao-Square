import re
import apns.pspot.parse_kernel.util as ampku
import uuid
def xml_standardize(fname: str):
    """ Pseudopotential XML file preprocess

    preprocess is relatively hard-coded. There are some cases the UPF file is not standard xml, therefore
    some modifications are needed. Cases considered:
    1. ADC pseudopotential has & symbol at the beginning of line, which is not allowed in xml, replace `&` 
    with `&amp;`
    2. GBRV pseudopotential does not startswith <UPF version="2.0.1">, and not endswith </UPF>, add <UPF version="2.0.1">
    to the beginning of the file and </UPF> to the end of the file
    3. some of ADC pseudopotentials have `CDATA` tag, which is not allowed in xml, remove them
    """
    ftemp = f"{str(uuid.uuid3(uuid.NAMESPACE_DNS, fname))}.xml"
    print(f"Preprocessing {fname}, will write standard XML formatted temporaray file to {ftemp}")
    # because a pseudopotential file would not be large, directly read all lines into memory
    with open(fname, "r") as f:
        lines = f.readlines()
    # it is compulsory for all XML files to start with <UPF version="2.0.1">
    if not lines[0].startswith("<UPF version="): # not expected case, but how?
        if lines[0].strip().startswith("<UPF version="):
            # remove all left whitespaces
            lines[0] = lines[0].strip() + "\n"
        else:
            lines.insert(0, "<UPF version=\"unknown\" comment=\"added to complete xml format\">\n")
    # from the last line, check if the first line startswith `<` is not </UPF>. If not, add </UPF> to the end of the file
    i = -1
    while not lines[i].strip() and i > -len(lines):
        i -= 1
    if not "</UPF>" in lines[i]:
        lines.append("</UPF>\n")
    # write the modified lines back to the file
    # there are more than one tasks remain:
    # for ADC pseudopotential, replace `&` with `&amp;`
    # for PseudoDojo v1.0 pseudopotential, carefully check consistency between opening and closing tags
    # for ADC pseudopotential or ld1 generated, replace eliminate both `<![CDATA[` and `]]>` if they exist
    lines = [line.replace("<![CDATA[", "").replace("]]>", "") for line in lines]
    with open(ftemp, "w") as f:
        # replace the `&` symbol at the beginning of the line with `&amp;`, this is done by xml_syntax_filter
        for _, line in ampku.xml_syntax_filter(lines):
        #for line in lines:
            _match = re.match(r"^([\s]*)(\&[\w]+)([^;]*)(\s*)$", line)
            line = f"{_match.group(1)}{_match.group(2).replace('&', '&amp;')}{_match.group(3)}{_match.group(4)}\n"\
                   if _match else line
            f.write(line)
            if "</UPF>" in line: # there are some pseudopotential files endswith ppgen file, but will crash the xml parser
                break
    return os.path.abspath(ftemp)

def postprocess(parsed: dict):

    for section in parsed:
        """first the data"""
        if parsed[section]["data"] is not None:
            parsed[section]["data"] = parsed[section]["data"].strip()
            if ampku.is_numeric_data(parsed[section]["data"]):
                parsed[section]["data"] = ampku.decompose_data(parsed[section]["data"])
        """then the attributes"""
        if parsed[section]["attrib"] is not None:
            for attrib in parsed[section]["attrib"]:
                parsed[section]["attrib"][attrib] = parsed[section]["attrib"][attrib].strip()
                if ampku.is_numeric_data(parsed[section]["attrib"][attrib]):
                    parsed[section]["attrib"][attrib] = ampku.decompose_data(parsed[section]["attrib"][attrib])
                elif parsed[section]["attrib"][attrib] == "T":
                    parsed[section]["attrib"][attrib] = True
                elif parsed[section]["attrib"][attrib] == "F":
                    parsed[section]["attrib"][attrib] = False
    return parsed

from lxml import etree
def parse(fname: str):
    """Main function of parsing the XML formatted pseudopotential file
    
    Several steps included:
    1. preprocess the file: standardize the pseudopotential file if needed. See annotation of function preprocess for details
    2. parse the XML file using ElementTree. If not properly prepocessed, the file will raise an error
    3. iterate through the tree and return a dictionary
    4. postprocess the dictionary, convert the data to the correct type
    """
    ftemp = xml_standardize(fname)
    tree = etree.parse(ftemp)
    root = tree.getroot()
    # iterate through the tree and return a dictionary
    parsed = {}
    for child in root:
        if isinstance(child.tag, str):
            parsed[child.tag] = {"data": child.text, "attrib": child.attrib}
    #parsed = postprocess(parsed)
    os.remove(ftemp)
    return parsed

def to_psp8(fname: str):
    """convert UPF (XML) format pseudopotential to psp8 format. This code is from cpmd2upf.x
    written by qespresso official"""
    raise NotImplementedError("psp8 format is not supported yet")

import os
import unittest
class TestXML(unittest.TestCase):
    def test_xml_standardize(self):
        contents = """There is neither a beginning tag nor an ending tag\n"""
        with open("test.upf", "w") as f:
            f.write(contents)
        ftemp = xml_standardize("test.upf")
        with open(ftemp, "r") as f:
            lines = f.readlines()
        self.assertEqual(lines[0], "<UPF version=\"unknown\" comment=\"added to complete xml format\">\n")
        self.assertEqual(lines[-1].strip(), "</UPF>")
        os.remove("test.upf")
        os.remove(ftemp)
        
        contents = """<UPF version="2.0.1">
There is a beginning tag but not an ending tag
"""
        with open("test.upf", "w") as f:
            f.write(contents)
        ftemp = xml_standardize("test.upf")
        with open(ftemp, "r") as f:
            lines = f.readlines()
        self.assertEqual(lines[0], "<UPF version=\"2.0.1\">\n")
        self.assertEqual(lines[-1].strip(), "</UPF>")
        os.remove("test.upf")
        os.remove(ftemp)

        contents = """<UPF version="2.0.1">
There is a beginning tag but with wrong ending tag
</!-->
"""
        with open("test.upf", "w") as f:
            f.write(contents)
        ftemp = xml_standardize("test.upf")
        with open(ftemp, "r") as f:
            lines = f.readlines()
        self.assertEqual(lines[0], "<UPF version=\"2.0.1\">\n")
        self.assertEqual(lines[-1].strip(), "</UPF>")
        os.remove("test.upf")
        os.remove(ftemp)

        contents = """      <UPF version="2.0.1">
    The correct case but with accidental indents
    </UPF>
"""
        with open("test.upf", "w") as f:
            f.write(contents)
        ftemp = xml_standardize("test.upf")
        with open(ftemp, "r") as f:
            lines = f.readlines()
        self.assertEqual(lines[0], "<UPF version=\"2.0.1\">\n")
        self.assertEqual(lines[-1].strip(), "</UPF>")
        os.remove("test.upf")
        os.remove(ftemp)

        contents = """<UPF version="2.0.1">
The correct case but with & symbol at the beginning of the line
&This is a test
</UPF>
"""
        with open("test.upf", "w") as f:
            f.write(contents)
        ftemp = xml_standardize("test.upf")
        with open(ftemp, "r") as f:
            lines = f.readlines()
        self.assertEqual(lines[0], "<UPF version=\"2.0.1\">\n")
        self.assertEqual(lines[-1].strip(), "</UPF>")
        self.assertEqual(lines[2], "&amp;This is a test\n")
        os.remove("test.upf")
        os.remove(ftemp)

if __name__ == "__main__":
    unittest.main()
    #xml_standardize("/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/pseudopotentials/ps-library/Ar.pbe-n-rrkjus_psl.1.0.0.UPF")