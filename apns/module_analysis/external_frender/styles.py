def styles_factory(property: str = "color", val: float = -1, ndim: int = None) -> list:
    if property == "color":
        colorpool = ["#2b316f", "#d8006a", "#24b5a5", "#e8cc47", "#005bbd"]
        return [colorpool[i % len(colorpool)] for i in range(ndim)]
    elif property == "marker":
        markerpool = ["o", "s", "D", "v", "^", "<", ">", "p", "P", "*", "h", "H", "+", "x", "X", "|", "_"]
        return [markerpool[i % len(markerpool)] for i in range(ndim)]
    elif property == "markersize":
        return [5 if val < 0 else val] * ndim
    elif property == "linestyle":
        linestylepool = ["-", "--", "-.", ":"]
        return [linestylepool[i % len(linestylepool)] for i in range(ndim)]
    elif property == "linewidth":
        return [1 if val < 0 else val] * ndim
    elif property == "alpha":
        return [1.0 if val < 0 else val] * ndim
    else:
        raise ValueError("Unknown property")
