def load(orbital_dir: str):
    """similar with load in module_pseudo/manage.py"""
    pass

def valid_nao(orbital_dir: str, elements: list, nao_settings: dict):
    """similar with valid_pseudo in module_pseudo/manage.py"""
    return {element: {} for element in elements}
    