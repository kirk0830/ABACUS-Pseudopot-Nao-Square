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

def convert_fpp_to_ppid(fpp: str):
    import re
    import os
    if "hgh" in fpp:
        family, version = "Hartwigsen-Goedecker-Hutter", ""
        appendix = re.match(r"([A-Z][a-z]?\.pbe\-)(.*)(hgh\.UPF)", os.path.basename(fpp)).group(2)
        appendix = "" if appendix is None else appendix
    elif "NCPP-PD04-PBE" in fpp:
        family, version = "PD04", ""
        appendix = re.match(r"([A-Z][a-z]?)([\d\+\-\_\w]*)(\.PD04\.PBE\.UPF)", os.path.basename(fpp)).group(2)
        appendix = "" if appendix is None else appendix[1:] if appendix[0] in ["+", "-"] else appendix
    elif "GBRV_pbe_UPF_v1.5" in fpp:
        family, version, appendix = "GBRV", "1.5", ""
    elif "nc-sr-05_pbe_standard_upf" in fpp:
        family, version, appendix = "PseudoDojo", "0.5", "sr"
    elif "nc-fr-04_pbe_standard" in fpp:
        family, version, appendix = "PseudoDojo", "0.4", "fr"
    elif "pbe_s_sr" in fpp:
        family, version, appendix = "PseudoDojo", "0.3", "sr"
    elif "NCPP-PD03-PBE" in fpp:
        family, version, appendix = "PD03", "", ""
    elif "sg15_oncv_upf_2020-02-06" in fpp:
        family = "SG15"
        match_ = re.match(r"([A-Z][a-z]?(_ONCV_PBE)(_)?(FR)?(\-)(\d\.\d)(\.upf))", os.path.basename(fpp))
        version = match_.group(6)
        appendix = "fr" if match_.group(4) is not None else "sr"
    elif "nc-sr-04_pbe_standard_upf" in fpp:
        family, version, appendix = "PseudoDojo", "0.4", "sr"
    elif "nc-sr-04-3plus_pbe_standard_upf" in fpp:
        family, version, appendix = "PseudoDojo", "0.4", "sr"
    elif "pseudos_ac_she" in fpp:
        family, version = "PseudoDojo", "1.0"
        appendix = "fr" if fpp.endswith("_r.upf") else "sr"
    elif "gth" in fpp:
        family, version = "Goedecker-Teter-Hutter", ""
        appendix = os.path.basename(fpp).split("_")[-1].split(".")[0]
    elif "psl" in fpp:
        match_ = re.match(r"([A-Z][a-z]?)(\.)(rel-)?(pbe|pz)(-\w+)?(-)(rrkjus|kjpaw)(_psl\.)([\.\d]+)(\.UPF)", os.path.basename(fpp))
        family = "PSlibrary"
        version = match_.group(9)
        apps = []
        if match_.group(7): apps.append(match_.group(7).upper())
        if match_.group(3): apps.append("fr")
        if match_.group(5): apps.append(match_.group(5)[1:])
        appendix = ", ".join(apps)
    else:
        raise ValueError(f"Unrecognized pseudopotential file: {fpp}")
    return f"{family} v{version} ({appendix})".replace("v ", "").replace("()", "")

def cal_dict_diff(desc1: dict, desc2: dict) -> dict:
    """calculate diff between two dict. For the same key, if the value is different,
    record the difference in tuple, the first element is from desc1, the second is from desc2.
    NOTE: if parent key has different value, the child key will not be compared.
    """
    if desc1 == desc2:
        return {}
    if not desc1:
        assert all([not isinstance(v, tuple) for v in desc2.values()])
        return {k: (None, v) for k, v in desc2.items()} if desc2 else {}
    if not desc2:
        assert all([not isinstance(v, tuple) for v in desc1.values()])
        return {k: (v, None) for k, v in desc1.items()} if desc1 else {}
    
    assert isinstance(desc1, dict) and isinstance(desc2, dict)
    # make sure before comparison, the value is not tuple
    assert all([not isinstance(v, tuple) for v in desc1.values()]), "desc1 should not have tuple value"
    assert all([not isinstance(v, tuple) for v in desc2.values()]), "desc2 should not have tuple value"

    diff = {}
    for k, v in desc1.items():
        v_ = desc2.get(k, None)
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
    keys = ["AtomSpecies", "Cell", "ParamSet", "CellGenerator"]
    assert set(keys) == set(desc1.keys()) == set(desc2.keys())
    # the assertation above is to ensure the correctness of structure
    # it must be satisfied before comparison
    diff = {}
    # first the AtomSpecies
    asdiff = [cal_dict_diff(as1, as2) for as1, as2 in zip(desc1["AtomSpecies"], desc2["AtomSpecies"])]
    diff.update({"AtomSpeices": asdiff}) if not all([not asd for asd in asdiff]) else None
    # other keys
    for key in keys[1:]:
        d = cal_dict_diff(desc1[key], desc2[key])
        diff.update({key: d}) if d else None
    return diff

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
            "ParamSet": {
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
            "ParamSet": {
                "basis_type": "pw"
            },
            "CellGenerator": {
                "identifier": "pwscf",
                "config": "pw.x"
            }
        }
        diff = cal_desc_diff(desc1, desc2)
        self.assertEqual(diff, {
            "AtomSpeices": [{}, {"nao": (2, 3)}]
        })

if __name__ == "__main__":
    unittest.main()