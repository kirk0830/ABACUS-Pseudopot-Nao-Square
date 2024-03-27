def generate_result_page(test: str = "Pseudopotential", 
                         element: str = "H",
                         pseudopotential_type: str = "norm-conserving",
                         functional: str = "PBE",
                         software: str = "ABACUS",
                         eos: bool = False,
                         cohesive_energy: bool = False,
                         dos: bool = True) -> str:
    """Generate result page for a given test.
    
    Args:
    
    test: str
        Test name. Currently, only "Pseudopotential" and "Pseudopot-Nao" are supported.
    element: str
        Element name.
    pseudopotential_type: str
        Pseudopotential type. Currently, only "norm-conserving" is supported.
    functional: str
        DFT xc functional. Currently, only "PBE" is supported.
    software: str
        DFT software. Currently, only "ABACUS" is supported.
    eos: bool
        Whether to show EOS plot.
    cohesive_energy: bool
        Whether to show cohesive energy plot.
    dos: bool
        Whether to show DOS plot.

    Returns:

    str
        Result page in markdown format.
    """

    if test != "Pseudopotential" and test != "Pseudopot-Nao":
        raise ValueError("test must be either Pseudopotential or Pseudopot-Nao")
    if pseudopotential_type != "norm-conserving" and pseudopotential_type != "ultrasoft":
        raise ValueError("pseudopotential_type must be either norm-conserving or ultrasoft")
    if pseudopotential_type != "norm-conserving":
        raise NotImplementedError("ananlysis for ultrasoft pseudopotential is not implemented yet.")
    if functional != "PBE":
        raise ValueError("functional must be PBE presently")
    if software != "ABACUS":
        raise ValueError("software must be ABACUS presently")
    
    convergence = element + ".svg"

    feos = ""
    if eos:
        feos = element + "_eos.svg"
        eos = """<img src="{}" class="plain-figure">""".format(feos)
    else:
        eos = """No data yet!"""
    feco = ""
    if cohesive_energy:
        feco = element + "_cohesive_energy.svg"
        cohesive_energy = """<img src="{}" class="plain-figure">""".format(feco)
    else:
        cohesive_energy = """No data yet!"""
    fdos = ""
    if dos:
        fdos = element + "_dos.svg"
        dos = """<img src="{}" class="plain-figure">""".format(fdos)
    else:
        dos = """No data yet!"""
    
    template = """---
layout: result
test: {0}
title: {1}
---

<h1>Pseudopotential tests</h1>
<h2>Test information</h2>
<ul>
    <li>element: {1}</li>
    <li>pseudopotential type: {2}</li>
    <li>DFT xc functional: {3}</li>
    <li>software: {4} (version: latest commit)</li>
</ul>
<h2>Test results</h2>
<table>
<tr><td>
<table class="banner-frame">
    <tr>
        <td class="banner-header">Convergence test</td>
    </tr>
    <tr>
        <td class="banner-body">
            <p align="center">
                <img src="{5}" class="plain-figure">
            </p>
        </td>
    </tr>
</table>
</td></tr><tr><td>
<table class="banner-frame">
    <tr>
        <td class="banner-header">Equation of States (EOS)</td>
    </tr>
    <tr>
        <td class="banner-body">
            <p align="center">
                {6}
            </p>
        </td>
    </tr>
</table>
</td></tr><tr><td>
<table class="banner-frame">
    <tr>
        <td class="banner-header">Cohesive energy curve</td>
    </tr>
    <tr>
        <td class="banner-body">
            <p align="center">
                {7}
            </p>
        </td>
    </tr>
</table>
</td></tr><tr><td>
<table class="banner-frame">
    <tr>
        <td class="banner-header">Density of States (DOS)</td>
    </tr>
    <tr>
        <td class="banner-body">
            <p align="center">
                {8}
            </p>
        </td>
    </tr>
</table>
</td></tr>
</table>""".format(test, element, pseudopotential_type, functional, software, convergence, eos, cohesive_energy, dos)

    return template

if __name__ == "__main__":
    print(generate_result_page(test = "Pseudopotential", 
                               element = "H",
                               pseudopotential_type = "norm-conserving",
                               functional = "PBE",
                               software = "ABACUS",
                               eos = False,
                               cohesive_energy = False,
                               dos = False))