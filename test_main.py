"""very initial configure"""

"""test of module_pseudo.archive.py"""
import apns.module_pseudo.upf_archive as arch
arch.archive(only_scan=True)

"""routinely run"""

"""test of initialization of work_folder"""
import apns.module_workflow.initialize as init
init.initialize_cache()
"""test of work_status_expand.py"""
_structures = {
    "Er2O3": ["Er2O3_2460", "Er2O3_1225560", "Er2O3_679"]
}
import apns.module_io.work_status_expand as wse
work_status = wse.render("input.json", system_with_mpids=_structures)
print(work_status)
# this will generate a dict as work_status

"""test of to_test_status.py"""
import apns.module_workflow.to_test_status as tts
test_status = tts.to(work_status=work_status)
#print(test_status)
# this will generate a dict as test_status

"""test of to_test.py"""
import apns.module_workflow.to_test as tt
tt.to(test_status=test_status,
      software=work_status["global"]["software"],
      basis_type=work_status["calculation"]["basis_type"],
      functionals=["pbe"])