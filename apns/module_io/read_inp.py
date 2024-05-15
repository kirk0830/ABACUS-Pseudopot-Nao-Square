"""This is a new version of input read module
Logic is designed as the following:

Test workflow

APNS is mainly the workflow collection for solving the problem "how to perform high-throughput screening over
both pseudopotentials (pp) and numerical atomic orbitals (nao)?"
Based on properties of interest, there are:
- energy (per atom)
- stress and pressure
- band structure (both band structure on kmesh and kline, the former can also be used to calculate band gap)
- density of states (DOS)
- equation of states (EOS)
- cohesive energies
- phonon frequencies
- ...
some of them determine the outermost loop of ABACUS. To impose those varying parameters, some of them requires
change on INPUT (or batch of INPUTs), some STRU, some KPT. Therefore there should be three modules for generating
those three files respectively. With more details, for those properties of interest, requires changes on:
-------------------------------------------------------------------
property         test                   files vary      files fixed
-------------------------------------------------------------------
energy           ecutwfc convergence    INPUT           STRU, KPT
stress/pressure  ecutwfc convergence    INPUT           STRU, KPT
band structure   ecutwfc convergence    INPUT           STRU, KPT   set KPT with care, for KLINE mode, use fixed INPUT
DOS              -                                      INPUT, STRU, KPT
EOS              accurancy              STRU            INPUT, KPT
cohesive energy  ecutwfc convergence    INPUT           STRU, KPT
phonon freq      ecutwfc convergence    INPUT           STRU, KPT
-------------------------------------------------------------------
based on the above table, the pp and nao will always varies. 

Nomenclature: INPUT, STRU and KPT
Indeed, there will be some interference between those three files:
- the KPT file would be overwritten if there is `gamma_only` or `kspacing` keyword defined in INPUT file. 
- latname in INPUT will be used to generate STRU file (but currently I have no interest to implement this)
- 
"""
a = "\t123"
print(a.expandtabs(6))