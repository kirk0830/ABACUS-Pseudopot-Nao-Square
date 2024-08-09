"""APNS2 has been released for weeks and APNS1 is not compatible with APNS2.
During the refactor of APNS2 website, the data for estabilishing APNS1 database
should be re-collected. Those jobs are either submitted by abacustest or
Bohrium-Lesbegue. I post the relevant information here for even future use:

job-abacustest-v0.3.42-03ca7c InX band gap pseudopotential tests
job-abacustest-v0.3.43-636969 InX band gap pseudopotential tests
job-abacustest-v0.3.47-d111d6 InX band gap pseudopotential tests
job-abacustest-v0.3.53-3c5d89 InX band gap pseudopotential tests
job-abacustest-v0.3.53-89f65e InX band gap pseudopotential tests

job-abacustest-v0.3.68-8aaad5 Er numerical atomic orbital test, bond length scan test

11549318 job-abacustest-v0.3.76-d072bb ecutwfc test failed: first 3 period
11549321 job-abacustest-v0.3.76-966d98 ecutwfc test failed: 4th period except FeMn
NO FORCES AND STRESS DATA job-abacustest-v0.3.90-6dfe2e ecutwfc test stopped: La series elements except Yb
NO FORCES AND STRESS DATA 11583294 job-abacustest-v0.3.90-c131ee ecutwfc test failed: 6th period: Cs - Bi
11587903 job-abacustest-v0.3.90-b5e66e ecutwfc test success: 5th period: Rb - Xe
11587900 job-abacustest-v0.3.90-bcf834 ecutwfc test failed: 6th period: Cs - Bi
11587867 job-abacustest-v0.3.90-c24ce1 ecutwfc test failed: La series elements except Yb
11588012 job-abacustest-v0.3.90-157582 ecutwfc test success: Lu
11628080 job-abacustest-v0.3.92-c3df9d ecutwfc test success: Fe
11662133 job-abacustest-v0.3.92-712464 ecutwfc test failed: Mn
11668634 job-abacustest-v0.3.92-58bd04 ecutwfc test success: Co
11841284 job-abacustest-v0.3.106-7290da ecutwfc test success: Cr
11845898 job-abacustest-v0.3.106-6fb67b ecutwfc test success: Ac
UNCONVERGED 11854128 job-abacustest-v0.3.107-a6a489 ecutwfc test failed: Pu
UNCONVERGED 11854477 job-abacustest-v0.3.107-4fd227 ecutwfc test stopped: Pu
11854514 job-abacustest-v0.3.107-b63f46 ecutwfc test success: Pu (with O)
11880642 job-abacustest-v0.3.108-7388cb ecutwfc test success: Th
11880707 job-abacustest-v0.3.108-fac27c ecutwfc test success: Pa
11880756 job-abacustest-v0.3.108-55984c ecutwfc test success: U
11880811 job-abacustest-v0.3.108-5d5317 ecutwfc test success: Np
11955072 job-abacustest-v0.3.113-75b0c3 ecutwfc test success: Mn
11955089 job-abacustest-v0.3.113-8a31dc ecutwfc test success: P
11955093 job-abacustest-v0.3.113-ef0906 ecutwfc test success: Se
BOHRIUM FAILED job-abacustest-v0.3.114-484bec ecutwfc test success: Yb
12310698 job-abacustest-v0.3.114-a947c0 ecutwfc test success: Yb

job-abacustest-v0.3.112-affbf8 EOS test success: Mg dojo04sr
job-abacustest-v0.3.112-b42de7 EOS test success: Mg dojo04sr
job-abacustest-v0.3.112-5afa2f EOS test success: Li, Be, B, C, Na, Mg, Al, Si, S
job-abacustest-v0.3.112-2914b8 EOS test success: K - Zn
job-abacustest-v0.3.112-04cd27 EOS test success: Li
job-abacustest-v0.3.112-5d28eb EOS test success: Li
job-abacustest-v0.3.112-54c25d EOS test success: Li, Be, B, C, Na, Mg, Al, Si, S
job-abacustest-v0.3.113-d43a70 EOS test success: Ac, Th, Pu, Pa, U, Np
job-abacustest-v0.3.113-e06c22 EOS test success: Ac, Th, Pu, Pa, U, Np
job-abacustest-v0.3.113-556a09 EOS test success: Ac, Th, Pu, Pa, U, Np
job-abacustest-v0.3.113-07844a EOS test success: Th
job-abacustest-v0.3.113-a715e3 EOS test success: Th
job-abacustest-v0.3.113-658e8c EOS test success: Co
job-abacustest-v0.3.113-ab0125 EOS test success: W
job-abacustest-v0.3.113-75b0c3 ecutwfc test success: Mn
job-abacustest-v0.3.113-8a31dc ecutwfc test success: P
job-abacustest-v0.3.113-ef0906 ecutwfc test success: Se
11961199 job-abacustest-v0.3.113-a46d73 EOS test success: H, He, Li, Be, B
11961535 job-abacustest-v0.3.113-5ac7a5 EOS test success: Hf, Cu, Al
11961537 job-abacustest-v0.3.113-16074d EOS test success: Mg, Mo, Nb
11961545 job-abacustest-v0.3.113-df7db2 EOS test success: Ta, Ti, V, Zr
11961891 job-abacustest-v0.3.113-48852a EOS test success: C, N, O, F, Ne
11962200 job-abacustest-v0.3.113-911e14 EOS test success: Na, P, Si
11962399 job-abacustest-v0.3.113-790df1 EOS test success: S, Cl, Ar, K, Ca
11963248 job-abacustest-v0.3.113-8b6308 EOS test success: Cr, Mn, Sc
11963394 job-abacustest-v0.3.113-d55a54 EOS test success: Fe, Ni, Zn
11963775 job-abacustest-v0.3.113-d9d71c EOS test success: As, Ga, Ge
11963851 job-abacustest-v0.3.113-763187 EOS test success: Se, Kr, Br
11963881 job-abacustest-v0.3.113-138a61 EOS test success: Rb, Sr, Y
11963919 job-abacustest-v0.3.113-be8684 EOS test success: Tc, Ru, Rh
11964051 job-abacustest-v0.3.113-2d6d4d EOS test success: Pd, Ag, Cd
11964236 job-abacustest-v0.3.113-1d7e2e EOS test success: In, Sb, Sn
11964329 job-abacustest-v0.3.113-7d717b EOS test success: Te, I, Xe
job-abacustest-v0.3.113-750676 EOS test success: Fe pslnc031
11972507 EOS test success: Cs, ...
11972582 EOS test success: Ba, ...
11972714 EOS test success: Re, ...
11972721 EOS test success: Pt, ...
11972723
11972782
11973093
11973674
11974609

job-abacustest-v0.3.113-8cb494 stopped: Li3MCl6 systems
job-abacustest-v0.3.113-49aad4 stopped: Li3MCl6 systems
job-abacustest-v0.3.113-208103 stopped: Li3MCl6 systems
job-abacustest-v0.3.113-0ccbf9
job-abacustest-v0.3.113-4110ba
job-abacustest-v0.3.113-c5aba8
"""

def _SG15(contracted: str):
    """the case the contracted name startswith "sg15"."""
    assert len(contracted) in [6, 8], f"not recognized contracted ppname: {contracted}"
    version = "v" + contracted[4] + "." + contracted[5]
    out = f"SG15 {version}"
    if len(contracted) == 8:
        out += " (" + contracted[6:] + ")"
    return out

def _PD03(contracted: str):
    """the case the contracted name startswith "pd03" """
    return "PD03"

def _PD04(contracted: str):
    """the case the contracted name startswith "pd04" """
    out = "PD04"
    if len(contracted) > 4:
        _ap = contracted[4:]
        _ap = "sp-exc" if _ap == "sp-e" else _ap
        out += " (" + _ap + ")"
    return out

def _PseudoDojo(contracted: str):
    """the case the contracted name startswith "dojo" """
    out = "PseudoDojo"
    version = "v" + contracted[4] + "." + contracted[5]
    out += " " + version
    if len(contracted) > 6:
        out += " (" + contracted[6:] + ")"
    return out

def _pslnc(contracted: str):
    """the case the contracted name startswith "pslnc" """
    out = "PSLibrary NC"
    if len(contracted) > 5:
        out += " (" + contracted[5:] + ")"
    return out

def _GTH(contracted: str):
    return "Goedecker-Teter-Hutter (LnPP1)"

def ppfullname(contracted: str):
    """because the APNS1 workflow has been deprecated since the release
    of APNS2, the contracted pseudopotential name now is not under
    maintainment anymore. Note that there are only limited number
    of contracted names, I will list them and convert to full name
    manually"""
    if contracted.startswith("sg15"):
        return _SG15(contracted)
    elif contracted.startswith("pd03"):
        return _PD03(contracted)
    elif contracted.startswith("pd04"):
        return _PD04(contracted)
    elif contracted.startswith("dojo"):
        return _PseudoDojo(contracted)
    elif contracted.startswith("pslnc"):
        return _pslnc(contracted)
    elif contracted.startswith("gth"):
        return _GTH(contracted)
    else:
        raise ValueError(f"not recognized contracted ppname: {contracted}")
    
def apns1zval(elem: str, contracted: str):
    """because the APNS1 workflow has been deprecated since the release
    of APNS2, the contracted pseudopotential name now is not under
    maintainment anymore. Here a full library is pasted from old
    commits. Find zval by element and contracted pseudopotential name"""
    import json
    from apns.pspot.parse import z_valence
    # this is hard-coded path, the pseudo_db.json is the mapping saved in
    # json format recording the encoding strategy of pseudopotential files
    # in APNS1
    with open("/root/abacus-develop/pseudopotentials/pseudo_db.json") as f:
        data = json.load(f).get(elem)
    contracted = "pd04sp-exc" if contracted == "pd04sp-e" else contracted
    fname = "invalid"
    for key, val in data.items():
        the_key = key.replace("_", "").replace(".", "")
        if contracted.startswith("dojo"):
            if not contracted.endswith("sr") and not contracted.endswith("fr") and not contracted.endswith("3plus"):
                contracted += "sr"
        if contracted == "gth":
            contracted = "gthLnPP1"
        if the_key == contracted:
            fname = val
            break
    assert fname != "invalid", f"not found {elem} {contracted}"
    # because APNS1 uses relative path, which is, not robust enough. APNS2 uses absolute path
    fname = fname.replace("./download/", "/root/abacus-develop/")
    return z_valence(fname)