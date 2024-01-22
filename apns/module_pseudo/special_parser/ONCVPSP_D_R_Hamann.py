def valelec_config(parsed: dict) -> list:

    result = {}

    zval = parsed["PP_HEADER"]["attrib"]["z_valence"]
    
    content = parsed["PP_INPUTFILE"]["data"]
    lines = [line.strip() for line in content.split("\n")]

    reference_config = []
    read_valence_config = False
    for line in lines:
        if line.startswith("#   n    l    f        energy (Ha)"):
            read_valence_config = True
            continue

        if read_valence_config:
            if line.startswith("#"):
                break
            else:
                if len(line.split()) >= 3:
                    reference_config.append(line)
    # reversely read the reference_config
    sequence = ["S", "P", "D", "F", "G", "H", "I", "K", "L", "M", "N"]
    reference_config.reverse()
    for line in reference_config:
        if zval <= 0:
            break
        else:
            words = line.split()
            index = int(words[1])
            symbol = sequence[index]
            if symbol not in result:
                result[symbol] = []
            result[symbol].append(words[0]+symbol)
            zval -= float(words[2])
    # then convert to list
    result_list = []
    for isym, symbol in enumerate(sequence):
        if symbol in result:
            result_list.append(result[symbol])
        else:
            if isym >= len(result):
                break
            else:
                result_list.append([])
    return result_list