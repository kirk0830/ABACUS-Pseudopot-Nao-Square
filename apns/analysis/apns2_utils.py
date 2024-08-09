#########
# utils #
#########
def read_apnsjob_desc(fdesc: str):
    import json
    with open(fdesc, 'r') as f:
        desc = json.load(f)

    print("Read APNS job description from file: ", fdesc)
    atom_species_symbols = [as_["symbol"] for as_ in desc["AtomSpecies"]]
    s = ", ".join(atom_species_symbols)
    print(f"Atom species: {s}")
    pps = [as_["pp"] for as_ in desc["AtomSpecies"]]
    s = "\n".join(pps)
    print(f"Pseudopotentials bond with:\n{s}")
    # naos = [as_["nao"] for as_ in desc["AtomSpecies"]]
    # s = ", ".join(naos)
    # print(f"Number of atomic orbitals: {s}")
    cellgen = desc["CellGenerator"]
    print(f"""CellGenerator (where the cell generated from)
identifier: {cellgen["identifier"]}
source: {cellgen["config"]}
""")
    return desc["AtomSpecies"], desc["CellGenerator"]

def handle_hgh(fpp: str):
    import re
    import os
    family, version = "Hartwigsen-Goedecker-Hutter", ""
    appendix = re.match(r"([A-Z][a-z]?\.pbe\-)(.*)(hgh\.UPF)", os.path.basename(fpp)).group(2)
    appendix = "" if appendix is None else appendix
    return family, version, appendix

def handle_pd04(fpp: str):
    import re
    import os
    family, version = "PD04", ""
    appendix = re.match(r"([A-Z][a-z]?)([\d\+\-\_\w]*)(\.PD04\.PBE\.UPF)", os.path.basename(fpp)).group(2)
    if appendix is None or len(appendix) == 0:
        appendix = ""
    elif appendix[0] in ["+", "-"]:
        appendix = appendix[1:]
    return family, version, appendix

def handle_gbrv(fpp: str):
    family, version, appendix = "GBRV", "1.5", ""
    return family, version, appendix

def handle_pd03(fpp: str):
    family, version, appendix = "PD03", "", ""
    return family, version, appendix

def handle_sg15(fpp: str):
    import re
    import os
    family = "SG15"
    match_ = re.match(r"([A-Z][a-z]?(_ONCV_PBE)(_)?(FR)?(\-)(\d\.\d)(\.upf))", os.path.basename(fpp))
    version = match_.group(6)
    appendix = "fr" if match_.group(4) is not None else "sr"
    return family, version, appendix

def handle_psl(fpp: str):
    import re
    import os
    family = "PSlibrary"
    match_ = re.match(r"([A-Z][a-z]?)(\.)(rel-)?(pbe|pz)(-\w+)?(-)(rrkjus|kjpaw|nc)(_psl\.)?([\.\d]+)?(\.UPF)", os.path.basename(fpp))
    version = "0.3.1" if match_.group(9) is None else match_.group(9)
    apps = []
    if match_.group(7): apps.append(match_.group(7).upper())
    if match_.group(3): apps.append("fr")
    if match_.group(5): apps.append(match_.group(5)[1:])
    appendix = ", ".join(apps)
    return family, version, appendix

def handle_pseudo_dojo(fpp: str):
    if "nc-sr-05_pbe_standard_upf" in fpp:
        family, version, appendix = "PseudoDojo", "0.5", "sr"
    elif "nc-fr-04_pbe_standard" in fpp:
        family, version, appendix = "PseudoDojo", "0.4", "fr"
    elif "pbe_s_sr" in fpp:
        family, version, appendix = "PseudoDojo", "0.3", "sr"
    elif "nc-sr-04_pbe_standard_upf" in fpp or "nc-sr-04-3plus_pbe_standard_upf" in fpp:
        family, version, appendix = "PseudoDojo", "0.4", "sr"
    elif "pseudos_ac_she" in fpp:
        family, version = "PseudoDojo", "1.0"
        appendix = "fr" if fpp.endswith("_r.upf") else "sr"
    else:
        raise ValueError(f"Unrecognized pseudopotential file: {fpp}")
    
    return family, version, appendix

def handle_gth(fpp: str):
    import os
    family, version = "Goedecker-Teter-Hutter", ""
    appendix = os.path.basename(fpp).split("_")[-1].split(".")[0]
    return family, version, appendix

def handle_20240723(fpp: str):
    import os
    family, version = "HighPressure", "20240723"
    appendix = "rcut="+os.path.basename(fpp).split("-")[1]
    return family, version, appendix

def convert_fpp_to_ppid(fpp: str):
    """Convert pseudopotential file name to pseudopotential identifier, the one
    more human readable. The conversion is based on the pseudopotential family
    and version. The appendix is also included in the identifier if it is not empty."""
    func_map = {
        "hgh": handle_hgh,
        "NCPP-PD04-PBE": handle_pd04,
        "GBRV_pbe_UPF_v1.5": handle_gbrv,
        "NCPP-PD03-PBE": handle_pd03,
        "sg15_oncv_upf_2020-02-06": handle_sg15,
        "nc-sr-05_pbe_standard_upf": handle_pseudo_dojo,
        "nc-fr-04_pbe_standard": handle_pseudo_dojo,
        "pbe_s_sr": handle_pseudo_dojo,
        "nc-sr-04_pbe_standard_upf": handle_pseudo_dojo,
        "nc-sr-04-3plus_pbe_standard_upf": handle_pseudo_dojo,
        "pseudos_ac_she": handle_pseudo_dojo,
        "gth": handle_gth,
        "psl": handle_psl,
        "high_pressure_oncv_upf": handle_20240723
    }
    print(f"Converting {fpp}")
    for key in func_map:
        if key in fpp:
            family, version, appendix = func_map[key](fpp)
            return f"{family} v{version} ({appendix})".replace("v ", "").replace("()", "")
    print(f"Unrecognized pseudopotential file: {fpp}")
    return fpp

def convert_forb_to_orbid(forb: str):
    import os, re
    forb = os.path.basename(forb)
    rcut_match = re.search(r"\d+(\.\d+)?au", forb)
    ecut_match = re.search(r"\d+(\.\d+)?Ry", forb)
    if not rcut_match or not ecut_match:
        raise ValueError("Invalid format for forbidden string")
    rcut = rcut_match.group(0)
    ecut = ecut_match.group(0)
    conf = forb.split("_")[-1].split(".")[0]
    return f"{rcut}, {ecut} ({conf})"

def cal_dict_diff(desc1: dict, desc2: dict) -> dict:
    """calculate diff between two dict. For the same key, if the value is different,
    record the difference in tuple, the first element is from desc1, the second is from desc2.
    NOTE: if parent key has different value, the child key will not be compared.
    """
    if desc1 == desc2:
        return {}
    if not desc1:
        # if desc1 is None, then desc2 should not have tuple value
        if any([isinstance(v, tuple) for v in desc2.values()]):
            raise ValueError("desc2 should not have tuple value")
        return {k: (None, v) for k, v in desc2.items()} if desc2 else {}
    if not desc2:
        # if desc2 is None, then desc1 should not have tuple value
        if any([isinstance(v, tuple) for v in desc1.values()]):
            raise ValueError("desc1 should not have tuple value")
        return {k: (v, None) for k, v in desc1.items()} if desc1 else {}
    
    # if they are not concurrently None, they should be dict
    if not (isinstance(desc1, dict) and isinstance(desc2, dict)):
        raise ValueError("desc1 and desc2 should be dict")
    # make sure before comparison, the value is not tuple
    if not all([not isinstance(v, tuple) for v in desc1.values()]):
        raise ValueError("desc1 should not have tuple value")
    if not all([not isinstance(v, tuple) for v in desc2.values()]):
        raise ValueError("desc2 should not have tuple value")

    diff = {}
    for k, v in desc1.items():
        v_ = desc2.get(k)
        if isinstance(v, dict):
            _diff = cal_dict_diff(v, v_)
            diff.update({k: _diff}) if _diff else None
        else:
            if v != v_:
                diff[k] = (v, v_)
    for k, v in desc2.items():
        if k not in desc1:
            diff[k] = (None, v)
    return diff

def cal_desc_diff(desc1: dict, desc2: dict) -> dict:
    """more specifically, calculate the difference between two APNS job
    descriptions. Because the AtomSpecies key has value as list of dict,
    therefore make difference for each element in list and return"""
    keys = ["AtomSpecies", "Cell", "DFTParamSet", "CellGenerator"]
    assert set(desc1.keys()) == set(desc2.keys()), f"desc1 and desc2 should have the same keys: {desc1.keys()} != {desc2.keys()}"
    assert set(keys) == set(desc1.keys()), f"desc1 should have keys: {keys} != {desc1.keys()}"
    assert set(keys) == set(desc2.keys()), f"desc2 should have keys: {keys} != {desc2.keys()}"
    # the assertation above is to ensure the correctness of structure
    # it must be satisfied before comparison
    diff = {}
    # first the AtomSpecies
    asdiff = [cal_dict_diff(as1, as2) for as1, as2 in zip(desc1["AtomSpecies"], desc2["AtomSpecies"])]
    diff.update({"AtomSpecies": asdiff}) if not all([not asd for asd in asdiff]) else None
    # other keys
    for key in keys[1:]:
        d = cal_dict_diff(desc1[key], desc2[key])
        diff.update({key: d}) if d else None
    return diff

def stru_rev_map(structure: str, basename: bool = False):
    """export the reverse map from file name to chemical formula"""
    import json, os
    try:
        with open(structure, "r") as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"File not found: {structure}")
        return {}
    except json.JSONDecodeError:
        print(f"Invalid JSON format in file: {structure}")
        return {}
    rev_map_ = {}
    for k, v in data.items():
        for v_ in v:
            f = os.path.basename(v_["file"]) if basename else v_["file"]
            rev_map_[f] = k + f" ({f})"
    return rev_map_

import unittest
class APNS2UtilsTest(unittest.TestCase):
    def test_cal_dict_diff(self):
        desc1 = {
            "a": 1,
            "b": {
                "c": 2,
                "d": 3
            },
            "e": 4
        }
        desc2 = {
            "a": 1,
            "b": {
                "c": 2,
                "d": 3
            },
            "e": 5
        }
        diff = cal_dict_diff(desc1, desc2)
        self.assertEqual(diff, {"e": (4, 5)})
        diff = cal_dict_diff(desc1, None)
        self.assertEqual(diff, {k: (v, None) for k, v in desc1.items()})
        diff = cal_dict_diff(None, desc2)
        self.assertEqual(diff, {k: (None, v) for k, v in desc2.items()})
        diff = cal_dict_diff(None, None)
        self.assertEqual(diff, {})
        diff = cal_dict_diff({}, {})
        self.assertEqual(diff, {})
        diff = cal_dict_diff({}, desc2)
        self.assertEqual(diff, {k: (None, v) for k, v in desc2.items()})
        diff = cal_dict_diff(desc1, {})
        self.assertEqual(diff, {k: (v, None) for k, v in desc1.items()})
        diff = cal_dict_diff(desc1, desc1)
        self.assertEqual(diff, {})
        diff = cal_dict_diff(desc2, desc2)
        self.assertEqual(diff, {})
        diff = cal_dict_diff(desc2, desc1)
        self.assertEqual(diff, {"e": (5, 4)})

        desc2 = {
            "a": 1,
            "b": {
                "c": 3,
                "d": 3
            },
            "e": 4
        }
        diff = cal_dict_diff(desc1, desc2)
        self.assertEqual(diff, {"b": {"c": (2, 3)}})

        desc3 = {
            "a": 1,
            "b": {
                "c": 2,
                "d": 3,
                "e": 4
            }
        } # someone write the e key in b by mistake
        diff = cal_dict_diff(desc1, desc3)
        self.assertEqual(diff, {"b": {"e": (None, 4)}, "e": (4, None)})

        desc4 = {
            "a": 1,
            "b": {
                "c": 2
            },
            "e": 4
        } # someone forget to write the d key in b
        diff = cal_dict_diff(desc1, desc4)
        self.assertEqual(diff, {"b": {"d": (3, None)}})

    def test_cal_desc_diff(self):
        desc1 = {
            "AtomSpecies": [
                {"symbol": "H", "pp": "H.pbe-rrkjus.UPF", "nao": 1},
                {"symbol": "O", "pp": "O.pbe-rrkjus.UPF", "nao": 2}
            ],
            "Cell": {
                "a": 1,
                "b": 2,
                "c": 3
            },
            "DFTParamSet": {
                "basis_type": "pw"
            },
            "CellGenerator": {
                "identifier": "pwscf",
                "config": "pw.x"
            }
        }
        desc2 = {
            "AtomSpecies": [
                {"symbol": "H", "pp": "H.pbe-rrkjus.UPF", "nao": 1},
                {"symbol": "O", "pp": "O.pbe-rrkjus.UPF", "nao": 3}
            ],
            "Cell": {
                "a": 1,
                "b": 2,
                "c": 3
            },
            "DFTParamSet": {
                "basis_type": "pw"
            },
            "CellGenerator": {
                "identifier": "pwscf",
                "config": "pw.x"
            }
        }
        diff = cal_desc_diff(desc1, desc2)
        self.assertEqual(diff, {
            "AtomSpecies": [{}, {"nao": (2, 3)}]
        })

    def test_handle_hgh(self):
        fpp = "H.pbe-hgh.UPF"
        family, version, appendix = handle_hgh(fpp)
        self.assertEqual(family, "Hartwigsen-Goedecker-Hutter")
        self.assertEqual(version, "")
        self.assertEqual(appendix, "")
    
    def test_handle_pd04(self):
        fpp = "H.PD04.PBE.UPF"
        family, version, appendix = handle_pd04(fpp)
        self.assertEqual(family, "PD04")
        self.assertEqual(version, "")
        self.assertEqual(appendix, "")
        fpp = "H-sp.PD04.PBE.UPF"
        family, version, appendix = handle_pd04(fpp)
        self.assertEqual(family, "PD04")
        self.assertEqual(version, "")
        self.assertEqual(appendix, "sp")
        fpp = "Gd3+_f--core-icmod1.PD04.PBE.UPF"
        family, version, appendix = handle_pd04(fpp)
        self.assertEqual(family, "PD04")
        self.assertEqual(version, "")
        self.assertEqual(appendix, "3+_f--core-icmod1")

    def test_handle_gbrv(self):
        fpp = "as_pbe_v1.uspp.F.UPF"
        family, version, appendix = handle_gbrv(fpp)
        self.assertEqual(family, "GBRV")
        self.assertEqual(version, "1.5")
        self.assertEqual(appendix, "")

    def test_handle_pd03(self):
        fpp = "H.PD03.PBE.UPF"
        family, version, appendix = handle_pd03(fpp)
        self.assertEqual(family, "PD03")
        self.assertEqual(version, "")
        self.assertEqual(appendix, "")
    
    def test_handle_sg15(self):
        fpp = "H_ONCV_PBE-1.0.upf"
        family, version, appendix = handle_sg15(fpp)
        self.assertEqual(family, "SG15")
        self.assertEqual(version, "1.0")
        self.assertEqual(appendix, "sr")
        fpp = "H_ONCV_PBE_FR-1.0.upf"
        family, version, appendix = handle_sg15(fpp)
        self.assertEqual(family, "SG15")
        self.assertEqual(version, "1.0")
        self.assertEqual(appendix, "fr")
    
    def test_handle_psl(self):
        fpp = "H.pbe-rrkjus_psl.1.0.0.UPF"
        family, version, appendix = handle_psl(fpp)
        self.assertEqual(family, "PSlibrary")
        self.assertEqual(version, "1.0.0")
        self.assertEqual(appendix, "RRKJUS")
        fpp = "H.rel-pbe-kjpaw_psl.0.1.UPF"
        family, version, appendix = handle_psl(fpp)
        self.assertEqual(family, "PSlibrary")
        self.assertEqual(version, "0.1")
        self.assertEqual(appendix, "KJPAW, fr")
        fpp = "Ac.rel-pbe-spfn-rrkjus_psl.1.0.0.UPF"
        family, version, appendix = handle_psl(fpp)
        self.assertEqual(family, "PSlibrary")
        self.assertEqual(version, "1.0.0")
        self.assertEqual(appendix, "RRKJUS, fr, spfn")


if __name__ == "__main__":
    unittest.main()