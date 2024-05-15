def structure(config: str, characteristic: float):
    if config == 'sc':
        return simple_cubic(characteristic)
    elif config == 'fcc':
        return face_centered_cubic(characteristic)
    elif config == 'bcc':
        return body_centered_cubic(characteristic)
    elif config == 'diamond':
        return diamond(characteristic)
    elif config == 'dimer':
        return dimer(characteristic)
    elif config == 'trimer':
        return trimer(characteristic)
    elif config == 'tetramer':
        return tetramer(characteristic)
    elif config == 'X2Y':
        return X2Y(characteristic)
    elif config == 'X2Y3':
        return X2Y3(characteristic)
    elif config == 'X2Y5':
        return X2Y5(characteristic)
    elif config == 'XY':
        return XY(characteristic)
    elif config == 'XY2':
        return XY2(characteristic)
    elif config == 'XY3':
        return XY3(characteristic)
    else:
        raise NotImplementedError(f'config not supported: {config}')
    
def dimer(bl: float):
    """build dimer with bond length `bl`
    
    Args:
        bl (float): bond length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [20, 20, 20, 90, 90, 90], \
        [0, 0], np.array([[0, 0, 0], [bl, 0, 0]])

def trimer(bl: float):
    """build triangle with bond length `bl`
    
    Args:
        bl (float): bond length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [20, 20, 20, 90, 90, 90], \
        [0, 0, 0], np.array([[0, 0, 0], [bl, 0, 0], [bl/2, bl/2*3**0.5, 0]])

def tetramer(bl: float):
    """build tetrahedron with bond length `bl`
    
    Args:
        bl (float): bond length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [20, 20, 20, 90, 90, 90], \
        [0, 0, 0, 0], np.array([[0, 0, 0], [bl, 0, 0], [bl/2, bl/2*3**0.5, 0], [bl/2, bl/2*3**0.5/3, bl/2*3**0.5]])

def simple_cubic(celldm: float):
    """build simple cubic with characteristic length `celldm`
    
    Args:
        celldm (float): characteristic length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [celldm, celldm, celldm, 90, 90, 90], \
        [0], np.array([[0, 0, 0]])

def body_centered_cubic(celldm: float):
    """build body centered cubic with characteristic length `celldm`
    
    Args:
        celldm (float): characteristic length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [celldm, celldm, celldm, 109.4712206, 109.4712206, 109.4712206], \
        [0], np.array([[0, 0, 0]])

def face_centered_cubic(celldm: float):
    """build face centered cubic with characteristic length `celldm`
    
    Args:
        celldm (float): characteristic length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [celldm, celldm, celldm, 60.0, 60.0, 60.0], \
        [0], np.array([[0, 0, 0]])

def diamond(celldm: float):
    """build diamond with characteristic length `celldm`
    
    Args:
        celldm (float): characteristic length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [celldm, celldm, celldm, 60.0, 60.0, 60.0], \
        [0, 0], np.array([[0, 0, 0], [0.25, 0.25, 0.25]])

def X2Y(celldm: float):
    """build X2Y structure with characteristic length `celldm`
    
    Args:
        celldm (float): characteristic length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [celldm, celldm, celldm, 60, 60, 60], \
        [0]*2 + [1]*1, \
        np.array([[0.75, 0.75, 0.75], [0.25, 0.25, 0.25], [0, 0, 0]])

def X2Y3(celldm: float):
    """build X2Y3 structure with characteristic length `celldm`
    
    Args:
        celldm (float): characteristic length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [celldm, celldm, celldm, 90, 90, 90], \
        [0]*4 + [1]*6, \
        np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], 
                    [0.25, 0.75, 0.75], 
                    [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0], 
                    [0.5, 0, 0.5], [0, 0.5, 0.5]])

def X2Y5(celldm: float):
    """build X2Y5 structure with characteristic length `celldm`
    
    Args:
        celldm (float): characteristic length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [celldm, celldm, celldm, 90, 90, 90], \
        [0]*4 + [1]*10, \
        np.array([[0.75, 0.75, 0.75], [0.25, 0.25, 0.75], [0.25, 0.75, 0.25],
                    [0.75, 0.25, 0.25],
                    [0.25, 0.75, 0.75], [0.75, 0.25, 0.75], [0.75, 0.75, 0.25],
                    [0.25, 0.25, 0.25], [0.5, 0.5, 0], [0.5, 0, 0.5],
                    [0, 0.5, 0.5], [0.5, 0, 0], [0, 0.5, 0],
                    [0, 0, 0.5]])

def XY(celldm: float):
    """build XO structure with characteristic length `celldm`
    
    Args:
        celldm (float): characteristic length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [celldm, celldm, celldm, 60, 60, 60], \
        [0]*1 + [1]*1, \
        np.array([[0, 0, 0], [0.5, 0.5, 0.5]])

def XY2(celldm: float):
    """build XO2 structure with characteristic length `celldm`
    
    Args:
        celldm (float): characteristic length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [celldm, celldm, celldm, 60, 60, 60], \
        [0]*1 + [1]*2, \
        np.array([[0, 0, 0], [0.75, 0.75, 0.75], [0.25, 0.25, 0.25]])

def XY3(celldm: float):
    """build XO3 structure with characteristic length `celldm`
    
    Args:
        celldm (float): characteristic length in Angstrom
    Returns:
        lat (list): lattice parameters, including a, b, c and alpha, beta, gamma
        species_map (list): mapping the index of coords to species, size = n_atoms
        coords (np.ndarray): coordinates of atoms
    """
    import numpy as np
    return [celldm, celldm, celldm, 60, 60, 60], \
        [0]*1 + [1]*3, \
        np.array([[0, 0, 0], 
                    [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]])
