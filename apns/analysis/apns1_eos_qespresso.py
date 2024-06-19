"""
For calibration of ABACUS PW results with QE.
ABACUS input converted to QE input by abacustest workflow
"""
import json
def collect(fjson):
    with open(fjson, 'r') as f:
        data = json.load(f)

    result = {}
    # collect data
    for key, value in data.items():
        element, mpid, pspotid, _ = key.split('_')
        result.setdefault(element, {}).setdefault(mpid, {})\
            .setdefault(pspotid, dict(zip(["volume", "energy"], [[], []])))\
            ["volume"].append(value["volume"])
        result[element][mpid][pspotid]["energy"].append(value["energy"])
    
    # sort data along volume
    for element, v in result.items():
        for mpid, v1 in v.items():
            for pspotid, v2 in v1.items():
                v2["volume"], v2["energy"] = zip(*sorted(zip(v2["volume"], v2["energy"])))
    
    return result

def calculate(data):
    # calculate EOS on data yielded by collect() function
    from apns.analysis.apns2_eos_utils import fit_birch_murnaghan
    result = {}

    for element, v in data.items():
        for mpid, v1 in v.items():
            for pspotid, v2 in v1.items():
                v = data[element][mpid][pspotid]["volume"]
                e = data[element][mpid][pspotid]["energy"]
                bm_fit = fit_birch_murnaghan(v, e, as_dict=True)
                if bm_fit != (None, None, None, None):
                    result.setdefault(element, {}).setdefault(mpid, {})\
                        .setdefault(pspotid, {})["bm_fit"] = bm_fit
                    result[element][mpid][pspotid]["data"] = [v, e]

    return result

import apns.analysis.apns1_eos_abacus as amaddeos
def calculate_delta(data):
    # calculate delta value on data yielded by calculate() function
    result = {}

    for element, v in data.items():
        for mpid, v1 in v.items():
            for pspotid, v2 in v1.items():
                bm_fit = v2["bm_fit"]
                vmin = min(v2["data"][0])
                vmax = max(v2["data"][0])
                delta, _ = amaddeos.cal_delta_wrtacwf(element=element, 
                                                      bmfit=bm_fit, 
                                                      vmin=vmin, vmax=vmax,
                                                      bravis=mpid)
                result.setdefault(element, {}).setdefault(mpid, {})[pspotid] = delta
    
    return result

def main(fjson):
    data = collect(fjson)
    result = calculate(data)
    delta = calculate_delta(result)
    print(delta)
    return delta

if __name__ == "__main__":
    main("./result.json")