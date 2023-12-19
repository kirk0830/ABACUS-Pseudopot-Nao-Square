"""main"""

if __name__ == "__main__":
    # Step 1: Generate work_status
    import sources.module_workflow.to_work_status as tws
    work_status = tws.to_work_status(fname="input.json", api_key="API_KEY", num_cif=1)
    # Step 2: Generate test_status
    import sources.module_workflow.to_test_status as tts
    test_status = tts.to_test_status(fname="input.json", work_status=work_status)