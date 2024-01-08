"""main"""

def main(input_file: str, version: str = "v1"):
    import apns.module_workflow.driver as amwd
    amwd.driver_v1(input_file) if version == "v1" else amwd.driver_v0(input_file)

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2: # interactive mode
        main("input.json")
    elif len(sys.argv) < 3: # command line mode
        main(sys.argv[1])
    else: # command line mode with version
        main(sys.argv[1], sys.argv[2])