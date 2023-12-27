from apns.module_pseudo.upf_import import GBRV_Vanderbilt as bp
import json

directory = 'D:/pseudopotential_test_workflow/abacus_pseudopotential_square/pseudopotentials/resources/'

test_result, test_connection = bp(directory + 'gbrv_15/as_pbe_v1.uspp.F.UPF')

with open("test_result.json", 'w') as f:
    json.dump(test_result, f, indent=4)
with open("test_connection.json", 'w') as f:
    json.dump(test_connection, f, indent=4)

"""
FAILED
"""