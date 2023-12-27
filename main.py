"""main"""

if __name__ == "__main__":
    # Step 0: initialize
    import apns.module_workflow.initialize as init
    init.initialize_cache()
    # Step 1: input -> work_status
    import apns.module_workflow.to_work_status as tws
    work_status = tws.to(fname="input.json")
    # Step 2: work_status -> test_status
    import apns.module_workflow.to_test_status as tts
    test_status = tts.to(work_status=work_status)
    # Step 3: test_status -> test
    import apns.module_workflow.to_test as tt
    tt.to(test_status=test_status,
          software=work_status["global"]["software"],
          basis_type=work_status["calculation"]["basis_type"],
          functionals=work_status["calculation"]["functionals"])