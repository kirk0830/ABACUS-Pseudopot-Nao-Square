"""citation information got from Google Scholar, Web of Science, etc."""

def citation(software: str = "abacus",
             pseudopot_nao: dict = None):
    
    result = "\n"
    result += "="*100
    result += "START OF CITATION INFORMATION\n"
    result += sssp()
    result += qespresso() if software == "qespresso" else abacus()
    result += seekpath()
    result += materials_project()
    result += "END OF CITATION INFORMATION\n"
    result += "="*100
    result += "\n"
    print(result)

def sssp():

    result = "-"*100
    result += "\n"
    result += """
We recommend user to read works done included in Standard Solid State Pseudopotentials (SSSP) library, which should be cited as:
Prandini G, Marrazzo A, Castelli I E, et al. Precision and efficiency in solid-state pseudopotential calculations[J]. npj Computational Materials, 2018, 4(1): 72.
Bosoni E, Beal L, Bercx M, et al. How to verify the precision of density-functional-theory implementations via reproducible and universal workflows[J]. Nature Reviews Physics, 2024, 6(1): 45-58.
"""
    result += "-"*100
    result += "\n"
    return result

def seekpath():

    result = "-"*100
    result += "\n"
    result += """
The crystal symmetry analysis is done by py-seekpath, which should be cited as:
Hinuma Y, Pizzi G, Kumagai Y, et al. Band structure diagram paths based on crystallography[J]. Computational Materials Science, 2017, 128: 140-184.

SPGLIB is the kernel of symmetry analysis, which should be cited as:
Togo A, Tanaka I. Spglib: a software library for crystal symmetry search[J]. arXiv preprint arXiv:1808.01590, 2018.
"""
    result += "-"*100
    result += "\n"
    return result

def materials_project():

    result = "-"*100
    result += "\n"
    result += """
The materials properties are obtained from Materials Project, which should be cited as:
Jain A, Ong S P, Hautier G, et al. Commentary: The Materials Project: A materials genome approach to accelerating materials innovation[J]. APL materials, 2013, 1(1).
"""
    result += "-"*100
    result += "\n"
    return result

def abacus():
    
    result = "-"*100
    result += "\n"
    result += """
The calculation is done by ABACUS, which should be cited as:
Chen M, Guo G C, He L. Systematically improvable optimized atomic basis sets for ab initio calculations[J]. Journal of Physics: Condensed Matter, 2010, 22(44): 445501.
Li P, Liu X, Chen M, et al. Large-scale ab initio simulations based on systematically improvable atomic basis[J]. Computational Materials Science, 2016, 112: 503-517.
"""
    result += "-"*100
    result += "\n"
    return result

def qespresso():

    result = "-"*100
    result += "\n"
    result += """
The calculation is done by Quantum ESPRESSO, which should be cited as:
Giannozzi P, Baroni S, Bonini N, et al. QUANTUM ESPRESSO: a modular and open-source software project for quantum simulations of materials[J]. Journal of physics: Condensed matter, 2009, 21(39): 395502.
Giannozzi P, Andreussi O, Brumme T, et al. Advanced capabilities for materials modelling with Quantum ESPRESSO[J]. Journal of physics: Condensed matter, 2017, 29(46): 465901.
Giannozzi P, Baseggio O, Bonfà P, et al. Quantum ESPRESSO toward the exascale[J]. The Journal of chemical physics, 2020, 152(15).
"""
    result += "-"*100
    result += "\n"
    return result
    
def pseudopotential(kind: list):

    print("%d kinds of pseudopotentials are used in this test, which are cited as:" % len(kind))

def pspot_sg15():

    print("")

def pspot_pd():

    print("")

def pspot_pslibrary():

    print("")

def numerical_orbital():

    print("The numerical orbitals are generated based on methods developed by ABACUS team, which are cited as:")
