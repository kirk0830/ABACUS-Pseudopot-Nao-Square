def citation(software: str = "abacus",
             pseudopot_nao: dict = None):
    
    if pseudopot_nao is None:
        pseudopot_nao = {}
    # software: ABACUS, qespresso
    if software == "abacus":
        abacus()
    elif software == "qespresso":
        qespresso()
    # python packages
    seekpath()
    materials_project()
    

def seekpath():

    print("The crystal symmetry analysis is done by py-seekpath, which is cited as:")
    print("- Y. Hinuma, et al., Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017)")
    print("SPGLIB is the kernel of symmetry analysis, which is cited as:")
    print("- A. Togo, et al., Spglib: a software library for crystal symmetry search, arXiv:1808.01590 (2018)")

def materials_project():

    print("The materials properties are obtained from Materials Project, which is cited as:")
    print("- A. Jain, et al., A high-throughput infrastructure for density functional theory calculations, Comp. Mat. Sci. 50, 2295 (2011)")

def abacus():
    
    print("The calculation is done by ABACUS, which is cited as:")
    print("- C. Mohan, et al., Systematically improvable optimized atomic basis sets for ab initio calculations. Journal of Physics: Condensed Matter 22.44 (2010): 445501.")
    print("- L. Pengfei, et al., Large-scale ab initio simulations based on systematically improvable atomic basis. Computational Materials Science 112 (2016): 503-517.")

def qespresso():

    print("The calculation is done by Quantum ESPRESSO, which is cited as:")
    print("- P. Giannozzi, et al., QUANTUM ESPRESSO: a modular and open-source software project for quantum simulations of materials, J. Phys.: Condens. Matter 21, 395502 (2009)")
    
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
