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