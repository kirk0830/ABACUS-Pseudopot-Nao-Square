if __name__ == "__main__":
    import time
    import argparse
    import json
    parser = argparse.ArgumentParser(description="Orbital generation workflow driver")
    parser.add_argument("-finp", "-i", help="input file for orbital generation workflow")
    args = parser.parse_args()
    print(f"Opening input file: {args.finp}")
    with open(args.finp, "r") as f:
        inp = json.load(f)
    print(inp)
    time.sleep(2)