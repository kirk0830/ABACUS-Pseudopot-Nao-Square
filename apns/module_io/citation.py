def citation(software: str = "ABACUS",
             pseudopot_nao: dict = {}):
    # software: ABACUS, qespresso
    if software == "ABACUS":
        abacus()
    elif software == "qespresso":
        qespresso()
    # python packages
    seekpath()
    materials_project()
    # pseudopotentials
    pseudopotential(kind=list(pseudopot_nao["pseudopotentials"]))
    if "numeric_orbitals" in pseudopot_nao.keys():
        numerical_orbital()

def seekpath():

    print("The crystal symmetry analysis is done by py-seekpath, which is cited as:")
    print("- Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017)")
    print("SPGLIB is the kernel of symmetry analysis, which is cited as:")
    print("- A. Togo, I. Tanaka, Spglib: a software library for crystal symmetry search, arXiv:1808.01590 (2018)")

def materials_project():

    print("The materials properties are obtained from Materials Project, which is cited as:")
    print("- A. Jain, S.P. Ong, G. Hautier, W. Chen, W.D. Richards, S. Dacek, S. Cholia, D. Gunter, D. Skinner, G. Ceder, K.A. Persson, A high-throughput infrastructure for density functional theory calculations, Comp. Mat. Sci. 50, 2295 (2011)")

def abacus():
    
    print("The calculation is done by ABACUS, which is cited as:")

def qespresso():

    print("The calculation is done by Quantum ESPRESSO, which is cited as:")
    print("- P. Giannozzi, S. Baroni, N. Bonini, M. Calandra, R. Car, C. Cavazzoni, D. Ceresoli, G.L. Chiarotti, M. Cococcioni, I. Dabo, A. Dal Corso, S. de Gironcoli, S. Fabris, G. Fratesi, R. Gebauer, U. Gerstmann, C. Gougoussis, A. Kokalj, M. Lazzeri, L. Martin-Samos, N. Marzari, F. Mauri, R. Mazzarello, S. Paolini, A. Pasquarello, L. Paulatto, C. Sbraccia, S. Scandolo, G. Sclauzero, A.P. Seitsonen, A. Smogunov, P. Umari, R.M. Wentzcovitch, QUANTUM ESPRESSO: a modular and open-source software project for quantum simulations of materials, J. Phys.: Condens. Matter 21, 395502 (2009)")
    
def pseudopotential(kind: list = []):

    print("%d kinds of pseudopotentials are used in this test, which are cited as:" % len(kind))

def pspot_sg15():

    print("")

def pspot_pd():

    print("")

def pspot_pslibrary():

    print("")

def numerical_orbital():

    print("The numerical orbitals are generated based on methods developed by ABACUS team, which are cited as:")