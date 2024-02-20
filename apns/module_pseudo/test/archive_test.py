import unittest


if __name__ == '__main__':
    import re
    #unittest.main()
    pattern = r"^([A-Z][a-z]?)(\.pbe-.*-)(.*)(_psl.*)(\.UPF)$"
    
    import os
    for file in os.listdir("./download/pseudopotentials/pslibrary-pbe.0.3.1/PSEUDOPOTENTIALS/"):
        _match = re.match(pattern, file)
        if _match:
            print("element", _match.group(1))
            print("type", _match.group(3))
            print("version", _match.group(4))
    pattern = r"^([A-Z][a-z]?)(\.pbe-.*)(kjpaw|rrkjus)(_psl.*)(\.UPF)$"
    fname = "Xe.pbe-n-kjpaw_psl.0.3.1.UPF"
    _match = re.match(pattern, fname)
    if _match:
        print("element", _match.group(1))
        print("type", _match.group(3))
        print("version", _match.group(4))
    else:
        print("no match")
    fname = "H.pbe-kjpaw_psl.0.1.UPF"
    _match = re.match(pattern, fname)
    if _match:
        print("element", _match.group(1))
        print("type", _match.group(3))
        print("version", _match.group(4))
    else:
        print("no match")
    pattern = r"^([A-Z][a-z]?)(\.pbe\-)([^\-]*)(\-)?(hgh.UPF)$"
    fname = "Al.pbe-hgh.UPF"
    _match = re.match(pattern, fname)
    if _match:
        print("element", _match.group(1))
        print("type", _match.group(3))
        print("version", _match.group(4))
    else:
        print("no match")
    fname = "Be.pbe-s-hgh.UPF"
    _match = re.match(pattern, fname)
    if _match:
        print("element", _match.group(1))
        print("type", _match.group(3))
        print("version", _match.group(4))
    else:
        print("no match")
    pattern = r"^([a-z]{1,2})(\_pbe\_)(.*)(v1.*)(\.uspp\.F\.UPF)$"
    fname = "hf_pbe_v1.uspp.F.UPF"
    _match = re.match(pattern, fname)
    if _match:
        print("element", _match.group(1))
        print("type", _match.group(3))
        print("version", _match.group(4))
    else:
        print("no match")
    fname = "hf_pbe_plus4_v1.uspp.F.UPF"
    _match = re.match(pattern, fname)
    if _match:
        print("element", _match.group(1))
        print("type", _match.group(3))
        print("version", _match.group(4))
    else:
        print("no match")