import argparse
def initialize():
    parser = argparse.ArgumentParser(description="APNS")
    parser.add_argument("-v", "--version", help="Version of the workflow", default="v1")
    parser.add_argument("-i", "--input", help="input file specifying the workflow", default="input.json")

    input_file = parser.parse_args().input
    version = parser.parse_args().version
    return input_file, version

"""main"""
import apns.module_workflow.driver as amwd
import apns.module_workflow.workflow_test.driver as amwtd
def main():

    input_file, version = initialize()
    if version == "v1":
        driver = amwd.spawn_driver(input_file)
        driver.setup()
        driver.run()
    else:
        amwtd.driver_v0(input_file)

if __name__ == "__main__":
    main()