def table(titles: list, data: list, suptitle: str = None, level: int = 2):

    if suptitle is not None:
        result = "#" + "#" * level + " " + suptitle + "\n"
    assert len(titles) == len(data)
    # presently force all items in data to have the same length
    data_length = 0
    for it, item in enumerate(titles):
        if data_length == 0:
            data_length = len(data[it])
        else:
            assert data_length == len(data[it])
    result += "|" + " | ".join(titles) + "|\n"
    result += "| " + "--- |"*len(titles) + "\n"
    for i in range(data_length):
        result += "| " + " | ".join([str(data[j][i]) for j in range(len(titles))]) + " |\n"
    
    return result

if __name__ == "__main__":
    
    ecutwfc = [10, 20, 30]
    ecutrho = [100, 200, 300]
    data = [ecutwfc, ecutrho]
    titles = ["ecutwfc", "ecutrho"]

    print(table(titles, data, "Convergence test"))
