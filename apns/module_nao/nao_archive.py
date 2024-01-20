"""no matter how these orbitals are stored, should we find them"""

def _scan_(orbital_dir: str, allow_skip: bool = True):
    """scanning numerical orbital is different from pseudopotential. From pseudopotential,
    the kind, version and other information is self-contained, but numerical orbital does
    not have its corresponding pseudopotential information
    
    Therefore archiving numerical orbital requires new name or careful design of how to name
    or store orbitals.
    
    Information of numerical orbitals should be divided into layers like:
    1. pseudopotential, including kind, version and appendix, therefore it would be natural
       to use identifier to name this layer
    2. (to be continue)"""

def archive(orbital_dir: str = "./download/numerical_orbitals/", only_scan: bool = True):

    if not orbital_dir.endswith("/") and not orbital_dir.endswith("\\"):
        orbital_dir += "/"

    if not only_scan:
        _, folders = _scan_(orbital_dir=orbital_dir, allow_skip=True)
