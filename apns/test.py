def run(finp: str):
    import apns.new.main as amn
    import json
    import apns.new.citation as amc
    with open(finp, "r") as f:
        inp = json.load(f)
    amn.main(inp)
    amc.citation()