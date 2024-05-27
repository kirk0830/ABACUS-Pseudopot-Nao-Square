import re
from lxml import etree

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

def xml_syntax_filter(content: str|list[str]):
    """it is found that some pseudopotential may have incorrect xml format, this function is to check the consistency of tags
    and correct if possible. A two-member tuple is returned, the first element is a boolean indicating whether the nearest tag is closed,
    the second element is the corrected line.

    ```python
    with open("file.xml", "r") as f:
        content = f.readlines()
    for is_closed, line in xml_syntax_filter(content):
        print(line)
    ```
    """
    content = content.split("\n") if isinstance(content, str) else content
    # for there are tags like "<PP_HEADER" in the file, we need to keep the buffer
    # then we can combine the buffer with the current line to form a complete tag
    # like "<PP_HEADER .../>", then it will be a single tag
    # or <PP_HEADER to form a complete "opening" tag <PP_HEADER ...>
    buf = ""

    # regular expression     example
    resg    = r"<[^>]+/>"    # <.../>
    reop    = r"<[^/][^>]+>" # <...>
    recls   = r"</[^>]+>"    # </...>
    relicmp = r"<[^>]+"      # <...
    rericmp = r"[^>]+>"      # ...>
    recmt   = r"<!--.*-->"   # <!--...-->
    
    # the stack to store the names of the tags, when a closing tag is found, it should be the same as the nearest opened tag
    # otherwise, it is a mismatch error and correction is needed

    # there is an exception caused by the opening comment: 20240506
    inbuilt_name_stack = []
    for line in content: # loop over all lines
        l = line.strip() # copy the value then do the strip operation
        # XML comments can be add both in one line and across line. However,
        # only in pslibrary 1.0.0 the across-line-comment is encountered...
        # but it is not clever to add more complexity (I mean more regular
        # expressions) to handle this case
        if re.match(recmt, l):
            yield True, line
        # single tag case, with the form <.../>, directly yield it, remember
        # to yield the original line, rather than the stripped one
        elif re.match(resg, l):
            yield True, line
        # opening tag case, with the form <...>. If this form appears, it
        # is possible that there are plenty of data belonging to this tag.
        # To save the tag name into the stack, then yield the original line
        elif re.match(reop, l):
            name = re.match(r"<([^ >]+)", l).group(1)
            inbuilt_name_stack.append(name)
            yield False, line
        # closing tag case, with the form </...>, the end of the possibilty
        # above. If this form appears, then it is the end of the tag, so pop
        # the tag name from the stack, then yield the original line. However,
        # for PseudoDojo v1.0 pseudopotential, there would be mismatch between
        # the opening tag and the closing tag, so need to correct it
        elif re.match(recls, l):
            name = re.match(r"</([^ >]+)", l).group(1)
            yield (True, line) if inbuilt_name_stack[-1] == name \
                else (True, line.replace(name, inbuilt_name_stack[-1]))
            inbuilt_name_stack.pop()
        # the case that a "left tag incomplete", means the tag is not closed but like
        # <PP_HEADER ...
        # comment="PP_HEADER is the most typical case that not-closed tag"
        # additional="but other tag can also be like this, e.g. PSWFC"
        # ... />
        # in this case, will use buffer to store all content until the tag is closed
        # then yield-back the complete tag.
        # however, must be sure the buffer is empty before the next tag
        # EXCEPTION: read <!-- in
        elif re.match(relicmp, l) and l != "<!--":
            assert buf == "", "buffer is not empty"
            buf = l
            continue
        elif re.match(rericmp, l) and l != "-->":
            assert buf != "", "buffer is empty"
            buf = " ".join([buf, l])
            if re.match(resg, buf):
                yield True, buf + "\n"
                buf = ""
            elif re.match(reop, buf):
                name = re.match(r"<([^ >]+)", buf).group(1)
                inbuilt_name_stack.append(name)
                yield False, buf + "\n"
                buf = ""
            else:
                raise ValueError(f"incorrect tag: {buf}")
        # for not-tag case, might be data, so in this case, do not touch
        # it. but it can also be within tag, say incomplete tag.
        else:
            # if buffer is not empty and the current line does not have any tag
            # related information, then it would be the situation within the tag
            # or out of any tags. Fortunately the tag will always start with its
            # name, or say the symbol `<` always appear to be in the same line
            # as the tag name, therefore once it is really within the tag, the
            # buffer will at least has the content of tag name, rather than empty
            if buf != "":
                buf = " ".join([buf, l])
                continue
            # if buffer is empty, then it is not within any tag, just yield the line
            else:
                yield True, line

import os
import unittest
class TestPseudoKernelUtil(unittest.TestCase):
    def test_xml_syntax_filter(self):
        with open("/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/pseudopotentials/NCPP-PD04-PBE/Kr.PD04.PBE.UPF", "r") as f:
            content = f.read()
        with open("unittest.xml", "w") as f:
            for state, line in xml_syntax_filter(content):
                f.write(line + "\n")
        
        # check if unittest.xml has correct XML format
        try:
            tree = etree.parse("unittest.xml")
        except etree.XMLSyntaxError as e:
            self.fail(f"unittest.xml is not a valid XML file: {e}")
        
        #os.remove("unittest.xml")

if __name__ == '__main__':
    unittest.main()