def scan_elements(system: str|list) -> list:

    if isinstance(system, str):
        pass
    elif isinstance(system, list):
        system = "_".join(system)
    else:
        raise TypeError("system must be str or list")
    
    # 1. remove all digits
    system = ''.join([i for i in system if not i.isdigit()])
    # 2. remove all special characters
    system = ''.join([i for i in system if i.isalpha()])
    # 3. split the system word by uppercase letters
    elements = []
    element = ""
    for index, letter in enumerate(system):
        if letter.isupper():
            if index == 0:
                element += letter
            else:
                elements.append(element)
                element = letter
        else:
            element += letter
    elements.append(element)

    # 4. remove duplicates but keeps order
    elements = list(dict.fromkeys(elements))
    # 5. remove empty strings
    elements = [element for element in elements if element != '']
    return elements

if __name__ == "__main__":
    
    print(scan_elements("Er2O3"))