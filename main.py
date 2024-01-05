"""main"""

def main(input_file: str):
    import apns.module_workflow.driver as amwd
    amwd.driver_v1(input_file)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        main("input.json")
    else:
        main(sys.argv[1])