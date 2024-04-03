import re
import xml.etree.ElementTree as ET

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

def xml_tagchecker(content: str|list[str]):
    """it is found that some pseudopotential may have incorrect xml format, this function is to check the consistency of tags
    and correct if possible. A two-member tuple is returned, the first element is a boolean indicating whether the nearest tag is closed,
    the second element is the corrected line.

    ```python
    with open("file.xml", "r") as f:
        content = f.readlines()
    for is_closed, line in xml_tagchecker(content):
        print(line)
    ```
    """
    content = content.split("\n") if isinstance(content, str) else content
    # for there are tags like "<PP_HEADER" in the file, we need to keep the linebuffer
    # then we can combine the linebuffer with the current line to form a complete tag
    # like "<PP_HEADER .../>", then it will be a single tag
    # or <PP_HEADER to form a complete "opening" tag <PP_HEADER ...>
    _buffer = ""
    # in XML, there are only three types of tags: single, opening, closing
    # complete tag: <.../>
    _single, _opening, _closing = r"<[^>]+/>", r"<[^/][^>]+>", r"</[^>]+>"
    # incomplete tag: <... or ...>
    _l_incomplete, _r_incomplete = r"<[^>]+", r"[^>]+>"
    # tags to be ignored
    _comment, _instruction = r"<!--.*-->", r"<\?.*\?>"
    
    # the stack to store the names of the tags, when a closing tag is found, it should be the same as the nearest opened tag
    # otherwise, it is a mismatch error and correction is needed
    names = []
    for line in content:
        # copy the value then do the strip operation
        _line = line.strip()
        # ignore the comment and instruction
        if re.match(_comment, _line) or re.match(_instruction, _line):
            yield True, line
        elif re.match(_single, _line):
            yield True, line
        elif re.match(_opening, _line):
            name = re.match(r"<([^ >]+)", _line).group(1)
            names.append(name)
            yield False, line
        elif re.match(_closing, _line):
            name = re.match(r"</([^ >]+)", _line).group(1)
            yield (True, line) if names[-1] == name else (True, line.replace(name, names[-1]))
            names.pop()
        elif re.match(_l_incomplete, _line):
            assert _buffer == "", "buffer is not empty"
            _buffer = _line
            continue
        elif re.match(_r_incomplete, _line):
            assert _buffer != "", "buffer is empty"
            _buffer = " ".join([_buffer, _line])
            if re.match(_single, _buffer):
                yield True, _buffer
                _buffer = ""
            elif re.match(_opening, _buffer):
                name = re.match(r"<([^ >]+)", _buffer).group(1)
                names.append(name)
                yield False, _buffer
                _buffer = ""
            else:
                raise ValueError(f"incorrect tag: {_buffer}")
        else:
            if _buffer != "":
                _buffer = " ".join([_buffer, _line])
                continue
            else:
                yield True, line

import os
import xml.etree.ElementTree as ET
import unittest
class TestPseudoKernelUtil(unittest.TestCase):
    def test_xml_tagchecker(self):
        content = """<UPF version="2.0.1">
  <PP_INFO>

 This pseudopotential file has been produced using the code
 ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)
 fully-relativistic version 4.0.1 06/20/2107 by D. R. Hamann
 The code is available through a link at URL www.mat-simresearch.com.
 Documentation with the package provides a full discription of the
 input data below.

    <PP_INPUTFILE>
    </PP_INPUTFILE>
  </PP_INFO>
  <!--                               -->
  <!-- END OF HUMAN READABLE SECTION -->
  <!--                               -->
    <PP_HEADER
       generated="Generated using ONCVPSP code by D. R. Hamann"
       author="anonymous"
       date="230825"
       comment=""
       element="Fm"
       pseudo_type="NC"
       relativistic="full"
       is_ultrasoft="F"
       is_paw="F"
       is_coulomb="F"
       has_so="T"
       has_wfc="F"
       has_gipaw="F"
       core_correction="T"
       functional="PBE"
       z_valence="   40.00"
       total_psenergy="  -1.67775776374E+03"
       rho_cutoff="   1.67900000000E+01"
       l_max="3"
       l_local="-1"
       mesh_size="  1680"
       number_of_wfc="11"
       number_of_proj="15"/>
 <PP_MESH>
   <PP_R type="real"  size="1680" columns="8">
    0.0000    0.0100    0.0200    0.0300    0.0400    0.0500    0.0600    0.0700
   </PP_R>
   <PP_RAB type="real"  size="1680" columns="8">
    0.0100    0.0100    0.0100    0.0100    0.0100    0.0100    0.0100    0.0100
   </PP_RAB>
 </PP_MESH>
  <PP_LOCAL type="real"  size="1680" columns="4">
   -1.4586394845E+02   -1.4585673047E+02   -1.4583505967E+02   -1.4579891619E+02
  </PP_LOCAL>
 <PP_NONLOCAL>
   <PP_BETA.1
       type="real"
       size="1680"
       columns="4"
       index="1"
       angular_momentum="0"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -1.3114099939E-07    1.6417536395E-01    3.2681146636E-01    4.8639972411E-01
   </PP_BETA.1>
   <PP_BETA.2
       type="real"
       size="1680"
       columns="4"
       index="2"
       angular_momentum="0"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -3.9236904481E-07    7.4429059472E-02    1.4634569308E-01    2.1331827448E-01
   </PP_BETA.2>
   <PP_BETA.3
       type="real"
       size="1680"
       columns="4"
       index="3"
       angular_momentum="0"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -1.9415543229E-08    3.6726765784E-02    7.3299258560E-02    1.0956696975E-01
   </PP_BETA.3>
   <PP_BETA.4
       type="real"
       size="1680"
       columns="4"
       index="4"
       angular_momentum="1"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -9.1552479486E-10    1.3026340132E-03    5.2085524789E-03    1.1711773561E-02
   </PP_BETA.4>
   <PP_BETA.5
       type="real"
       size="1680"
       columns="4"
       index="5"
       angular_momentum="1"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -3.0334206147E-08    3.7539768897E-03    1.4995576146E-02    3.3663477396E-02
   </PP_BETA.5>
   <PP_BETA.6
       type="real"
       size="1680"
       columns="4"
       index="6"
       angular_momentum="1"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -5.0889758740E-08    4.5180556847E-03    1.8040976428E-02    4.0474528040E-02
   </PP_BETA.6>
   <PP_BETA.7
       type="real"
       size="1680"
       columns="4"
       index="7"
       angular_momentum="1"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -4.2503100914E-08    5.4416393372E-05    2.3669754290E-04    6.0317273181E-04
   </PP_BETA.7>
   <PP_BETA.8
       type="real"
       size="1680"
       columns="4"
       index="8"
       angular_momentum="2"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -7.4937106364E-09    1.5215967842E-05    1.2202885229E-04    4.1352707322E-04
   </PP_BETA.8>
   <PP_BETA.9
       type="real"
       size="1680"
       columns="4"
       index="9"
       angular_momentum="2"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -3.2273620387E-08    1.3526083956E-04    1.0801185018E-03    3.6343162861E-03
   </PP_BETA.9>
   <PP_BETA.10
       type="real"
       size="1680"
       columns="4"
       index="10"
       angular_momentum="2"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
    6.8528415470E-08    1.3043983734E-04    1.0392880039E-03    3.4839232380E-03
   </PP_BETA.10>
   <PP_BETA.11
       type="real"
       size="1680"
       columns="4"
       index="11"
       angular_momentum="2"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -5.2719028205E-09    1.9120204016E-05    1.5310070945E-04    5.1748766384E-04
   </PP_BETA.11>
   <PP_BETA.12
       type="real"
       size="1680"
       columns="4"
       index="12"
       angular_momentum="3"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
    1.5602303824E-08    2.2305985161E-07    3.5955509417E-06    1.8425026772E-05
   </PP_BETA.12>
   <PP_BETA.13
       type="real"
       size="1680"
       columns="4"
       index="13"
       angular_momentum="3"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
    1.5642224010E-08    2.0826409589E-07    3.3592128397E-06    1.7231887613E-05
   </PP_BETA.13>
   <PP_BETA.14
       type="real"
       size="1680"
       columns="4"
       index="14"
       angular_momentum="3"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -2.5487721665E-08    1.3168348068E-07    2.0615713046E-06    1.0057012553E-05
   </PP_BETA.14>
   <PP_BETA.15
       type="real"
       size="1680"
       columns="4"
       index="15"
       angular_momentum="3"
       cutoff_radius_index=" 140"
       cutoff_radius="    1.3900000000E+00" >
   -2.5143515767E-08    1.2640912565E-07    1.9778658079E-06    9.6389866815E-06
   </PP_BETA.15>
   <PP_DIJ type="real"  size=" 225" columns="4">
    4.2884002764E+01    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    5.6439197789E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
   -3.5066915411E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
   -9.3922015413E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    1.4776525958E+01    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    8.7206907180E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
   -2.0826582538E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
   -6.8245389897E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    2.8472373156E+01    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
   -7.7174456147E-01    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
   -4.7621864914E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
   -1.1406988275E+01    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
   -1.0760341953E+01    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
   -3.7754819214E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
   -3.5202548746E+00
   </PP_DIJ>
 </PP_NONLOCAL>
 <PP_PSWFC>
   <PP_CHI.1
       type="real"
       size="1680"
       columns="4"
       index="1"
       occupation=" 2.000"
       pseudo_energy="   -0.3195715086E+02"
       label="5S"
       l="0" >
    3.5525860032E-11    4.3858871451E-02    8.7681181416E-02    1.3143029307E-01
   </PP_CHI.1>
   <PP_CHI.2
       type="real"
       size="1680"
       columns="4"
       index="2"
       occupation=" 4.000"
       pseudo_energy="   -0.1900691083E+02"
       label="5P"
       l="1" >
    6.7716082741E-11    9.7518331765E-04    3.8977802540E-03    8.7589418091E-03
   </PP_CHI.2>
   <PP_CHI.3
       type="real"
       size="1680"
       columns="4"
       index="3"
       occupation=" 2.000"
       pseudo_energy="   -0.2612778243E+02"
       label="5P"
       l="1" >
    3.9618466841E-10    1.4334559282E-03    5.7275838435E-03    1.2863697621E-02
   </PP_CHI.3>
   <PP_CHI.4
       type="real"
       size="1680"
       columns="4"
       index="4"
       occupation=" 6.000"
       pseudo_energy="   -0.9560174160E+01"
       label="5D"
       l="2" >
    1.5655958127E-10    1.9996004375E-05    1.5980039719E-04    5.3838467537E-04
   </PP_CHI.4>
   <PP_CHI.5
       type="real"
       size="1680"
       columns="4"
       index="5"
       occupation=" 4.000"
       pseudo_energy="   -0.1061047229E+02"
       label="5D"
       l="2" >
    1.5956532450E-10    2.1317200240E-05    1.7035958882E-04    5.7396362224E-04
   </PP_CHI.5>
   <PP_CHI.6
       type="real"
       size="1680"
       columns="4"
       index="6"
       occupation=" 6.857"
       pseudo_energy="   -0.2915543856E+00"
       label="5F"
       l="3" >
   -3.7498630145E-10    1.9589796909E-07    3.1313179540E-06    1.5826598365E-05
   </PP_CHI.6>
   <PP_CHI.7
       type="real"
       size="1680"
       columns="4"
       index="7"
       occupation=" 5.143"
       pseudo_energy="   -0.4352125751E+00"
       label="5F"
       l="3" >
   -4.0091289858E-10    2.0700815165E-07    3.3088830052E-06    1.6723854119E-05
   </PP_CHI.7>
   <PP_CHI.8
       type="real"
       size="1680"
       columns="4"
       index="8"
       occupation=" 2.000"
       pseudo_energy="   -0.4387812550E+01"
       label="6S"
       l="0" >
    6.0105032063E-11   -1.9704979474E-02   -3.9394608404E-02   -5.9053316685E-02
   </PP_CHI.8>
   <PP_CHI.9
       type="real"
       size="1680"
       columns="4"
       index="9"
       occupation=" 4.000"
       pseudo_energy="   -0.1594822595E+01"
       label="6P"
       l="1" >
   -1.8824768480E-10   -4.1733287536E-04   -1.6677789679E-03   -3.7466874943E-03
   </PP_CHI.9>
   <PP_CHI.10
       type="real"
       size="1680"
       columns="4"
       index="*"
       occupation=" 2.000"
       pseudo_energy="   -0.2669524562E+01"
       label="6P"
       l="1" >
   -4.4805911759E-10   -7.0718311849E-04   -2.8246880836E-03   -6.3404070224E-03
   </PP_CHI.*>
   <PP_CHI.11
       type="real"
       size="1680"
       columns="4"
       index="*"
       occupation=" 2.000"
       pseudo_energy="   -0.2988599205E+00"
       label="7S"
       l="0" >
    4.8731907398E-11    5.4235207341E-03    1.0843363006E-02    1.6255763261E-02
   </PP_CHI.*>
 </PP_PSWFC>
 <PP_NLCC type="real"  size="1680" columns="4">
    5.0209028196E+01    5.0175470342E+01    5.0074922722E+01    4.9907762540E+01
 </PP_NLCC>
 <PP_RHOATOM type="real"  size="1680" columns="4">
    0.0000000000E+00    4.6922174181E-03    1.8868740833E-02    4.2828369886E-02
 </PP_RHOATOM>
 <PP_SPIN_ORB>
   <PP_RELBETA.1  index="1"  lll="0" jjj="0.5"/>
   <PP_RELBETA.2  index="2"  lll="0" jjj="0.5"/>
   <PP_RELBETA.3  index="3"  lll="0" jjj="0.5"/>
   <PP_RELBETA.4  index="4"  lll="1" jjj="0.5"/>
   <PP_RELBETA.5  index="5"  lll="1" jjj="1.5"/>
   <PP_RELBETA.6  index="6"  lll="1" jjj="0.5"/>
   <PP_RELBETA.7  index="7"  lll="1" jjj="1.5"/>
   <PP_RELBETA.8  index="8"  lll="2" jjj="1.5"/>
   <PP_RELBETA.9  index="9"  lll="2" jjj="2.5"/>
   <PP_RELBETA.10 index="10" lll="2" jjj="1.5"/>
   <PP_RELBETA.11 index="11" lll="2" jjj="2.5"/>
   <PP_RELBETA.12 index="12" lll="3" jjj="2.5"/>
   <PP_RELBETA.13 index="13" lll="3" jjj="3.5"/>
   <PP_RELBETA.14 index="14" lll="3" jjj="2.5"/>
   <PP_RELBETA.15 index="15" lll="3" jjj="3.5"/>
   <PP_RELWFC.1  index="1"  lchi="0" jchi="0.5" nn="1"/>
   <PP_RELWFC.2  index="2"  lchi="1" jchi="1.5" nn="2"/>
   <PP_RELWFC.3  index="3"  lchi="1" jchi="0.5" nn="2"/>
   <PP_RELWFC.4  index="4"  lchi="2" jchi="2.5" nn="3"/>
   <PP_RELWFC.5  index="5"  lchi="2" jchi="1.5" nn="3"/>
   <PP_RELWFC.6  index="6"  lchi="3" jchi="3.5" nn="4"/>
   <PP_RELWFC.7  index="7"  lchi="3" jchi="2.5" nn="4"/>
   <PP_RELWFC.8  index="8"  lchi="0" jchi="0.5" nn="5"/>
   <PP_RELWFC.9  index="9"  lchi="1" jchi="1.5" nn="6"/>
   <PP_RELWFC.10  index="10"  lchi="1" jchi="0.5" nn="6"/>
   <PP_RELWFC.11  index="11"  lchi="0" jchi="0.5" nn="7"/>
 </PP_SPIN_ORB>
</UPF>"""
        
        with open("unittest.xml", "w") as f:
            for state, line in xml_tagchecker(content):
                f.write(line + "\n")
        
        # check if unittest.xml has correct XML format
        try:
            tree = ET.parse("unittest.xml")
        except ET.ParseError:
            self.fail("unittest.xml is not a valid XML file")
        
        os.remove("unittest.xml")

if __name__ == '__main__':
    unittest.main()