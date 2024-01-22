def valelec_config(parsed: dict) -> list:

    result = {}

    content = parsed["PP_INFO"]["data"]
    lines = content.split("\n")

    read_valence_config = False
    for line in lines:
        line = line.strip()
        if line.startswith("Valence configuration:"):
            read_valence_config = True
            continue

        if line.startswith("Generation configuration:"):
            break

        if read_valence_config and line[0].isdigit():
            words = line.split()
            if len(words) == 7:
                symbol = words[0][-1]
                if symbol not in result:
                    result[symbol] = []
                result[symbol].append(words[0])
    # then convert to list
    sequence = ["S", "P", "D", "F", "G", "H", "I", "K", "L", "M", "N"]
    result_list = []

    for isym, symbol in enumerate(sequence):
        if symbol in result:
            result_list.append(result[symbol])
        else:
            if isym > len(result):
                break
            else:
                result_list.append([])
    return result_list