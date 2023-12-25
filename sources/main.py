"""main"""

if __name__ == "__main__":
    # Step 1: input -> work_status
    import sources.module_workflow.to_work_status as tws
    work_status = tws.to_work_status(fname="input.json", api_key="API_KEY", num_cif=1)
    # Step 2: work_status -> test_status
    import sources.module_workflow.to_test_status as tts
    test_status = tts.to_test_status(work_status=work_status)
    # Step 3: test_status -> test