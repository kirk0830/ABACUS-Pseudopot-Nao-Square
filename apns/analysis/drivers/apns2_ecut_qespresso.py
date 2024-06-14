def collect_apnsjob_data(folder: str):
    print("* * * Collect QESPRESSO result * * *".center(100))
    import apns.analysis.postprocess.read_qespresso_out as read_qe
    import apns.analysis.postprocess.read_abacus_out as read_abacus
    import apns.pspot.parse as ppparse
    import os, re, json
    result = {}
    for root, _, files in os.walk(folder):
        if "description.json" in files:
            with open(os.path.join(root, "description.json"), "r") as f:
                desc = json.load(f)
            natom = desc["DFTParamSet"]["system"]["nat"]
            eks = read_qe.read_e(os.path.join(root, "pwscf.out"))
            press = read_qe.read_press(os.path.join(root, "pwscf.out"))
            _, bs = read_qe.read_bs(os.path.join(root, "pwscf.out"))
            if None in [eks, press, bs]:
                print(f"WARNING: Present APNS job is broken: {root}")
                continue
            