def scan_elements(system: str) -> list:

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

    # 4. remove duplicates
    elements = list(set(elements))
    # 5. remove empty strings
    elements = [element for element in elements if element != '']
    return elements