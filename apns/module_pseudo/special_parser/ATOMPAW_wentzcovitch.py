def valelec_config(parsed: dict) -> list:

    result = {}

    contents = parsed["PP_INFO"]["data"]
    lines = [line.strip() for line in contents.split("\n")]
    
    sequence = ["S", "P", "D", "F", "G", "H", "I", "K", "L", "M", "N"]

    nmax = []
    norb = []
    iorb = 0
    read_valence_config = False
    for il, line in enumerate(lines):
        if il == 3: # this is the line specifying how many orbitals considered for s, p, d, f and g
            nmax = [int(x) for x in line.split()]
            norb = [(nmax[i] - i) if nmax[i] > i else 0 for i in range(len(nmax))]
            continue
        if line == "0 0 0":
            read_valence_config = True
            continue
        if read_valence_config:
            if not line[0].isalpha():
                break
            else:
                if line == "c": # core electron
                    pass
                elif line == "v": # valence electron
                    present_l = 0
                    present_n = 1
                    while iorb >= sum(norb[:present_l]):
                        present_l += 1
                    present_n += iorb - sum(norb[:present_l - 1]) + present_l - 1
                    symbol = sequence[present_l - 1]
                    if symbol not in result:
                        result[symbol] = []
                    result[symbol].append(str(present_n)+symbol)
                else:
                    raise ValueError("Unknown line in valelec_config: {}".format(line))
                iorb += 1
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