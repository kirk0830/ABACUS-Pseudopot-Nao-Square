from apns.module_pseudo.upf_import import ONCVPSP_D_R_Hamann as bp
import json

directory = 'D:/pseudopotential_test_workflow/abacus_pseudopotential_square/pseudopotentials/resources/'

test_result, test_connection = bp(directory + 'dojo_05/As_PBE_dojo_05_NCSR.upf')
with open("test_result.json", 'w') as f:
    json.dump(test_result, f, indent=4)
with open("test_connection.json", 'w') as f:
    json.dump(test_connection, f, indent=4)