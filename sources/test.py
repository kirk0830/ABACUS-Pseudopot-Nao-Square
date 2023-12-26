"""very initial configure"""
"""test of module_pseudo.archive.py"""
import module_pseudo.archive as arch
arch.archive(only_scan=True)

"""routinely run"""
"""test of work_status_expand.py"""
_structures = {
    "Er2O3": ["Er2O3_2460", "Er2O3_1225560", "Er2O3_679"]
}
import module_io.work_status_expand as wse
work_status = wse.render("input.json", system_with_mpids=_structures)
print(work_status)
# this will generate a dict as work_status

"""test of to_test_status.py"""
import module_workflow.to_test_status as tts
test_status = tts.to_test_status(work_status=work_status)
#print(test_status)
# this will generate a dict as test_status

"""test of to_test.py"""
import module_workflow.to_test as tt
tt.generate(test_status=test_status,
            software=work_status["global"]["software"],
            basis_type=work_status["calculation"]["basis_type"],
            functionals=["pbe"])