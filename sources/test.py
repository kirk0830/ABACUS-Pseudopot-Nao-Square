import module_io.work_status_expand as wse

"""test of work_status_expand.py"""
_structures = {
    "Er2O3": ["Er2O3_2460", "Er2O3_1225560", "Er2O3_679"]
}
work_status = wse.render("input.json", system_with_mpids=_structures)
print(work_status)
# this will generate a dict as work_status

import module_workflow.to_test_status as tts
test_status = tts.to_test_status(work_status=work_status)
print(test_status)
# this will generate a dict as test_status