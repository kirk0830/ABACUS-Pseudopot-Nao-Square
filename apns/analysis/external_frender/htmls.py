"""This file contains html templates for APNS pages, which need to be dynamically rendered with data from the database"""

"""utils"""
def liquid_header(layout: str, test: str, title: str):
    """Liquid header, a jekyll header for liquid layout"""
    return f"---\nlayout: {layout}\ntest: {test}\ntitle: {title}\n---\n\n"

def banner_frame(title: str, content: str|list):
    """Banner frame"""
    content = [content] if isinstance(content, str) else content
    content = [f"<p align=\"center\">{c}</p>" for c in content]
    result =  f"""<table class="banner-frame">
    <tr>
        <td class="banner-header">{title}</td>
    </tr>
    <tr>
        <td class="banner-body">
"""
    result += "\n".join(content)
    result += """
        </td>
    </tr>
</table>"""
    return result

def header(title: str, level: int = 1):
    """Header"""
    return f"<h{level}>{title}</h{level}>"

def unordered_list(content: list):
    """Unordered list"""
    content = [f"<li>{c}</li>" for c in content]
    result = "<ul>"
    result += "\n".join(content)
    result += "</ul>"
    return result

def ordered_list(content: list):
    """Ordered list"""
    content = [f"<li>{c}</li>" for c in content]
    result = "<ol>"
    result += "\n".join(content)
    result += "</ol>"
    return result

"""Pseudopotentials - [ELEMENT]
This is a html-string-like template for pseudopotential tests,
presently this page is designed for presenting following tests:

1. Convergence test
Convergence test of `ecutwfc` set for ABACUS planewave calculation. For more example please
refer to Standard Solid State Pseudopotential (SSSP) library. The convergence of `ecutwfc` is
crucial for the accuracy of the calculation. If not enough `ecutwfc` is used, the calculation
actually uses a far less complete basis set, and the results are not reliable.

Two factors are considered in the convergence test:
    - energy per atom
    - pressure of whole cell
Thresholds are set as:
      standards:                    accurate    normal      rough
    - energy per atom/meV/atom          1e-3      1e-2       1e-1
    - pressure of whole cell/kbar        0.1       0.5        0.5

About plotting
Line plot, with markers

2. Equation of States (EOS)
Equation of States (EOS) is a test to determine the equation of state of a material. EOS can
be used to determine the ability of pseudopotential to predict the lattice constant and 
structural properties of a material. The test is performed by fitting the energy-volume data
to a Birch-Murnaghan equation of state.
After the "convergence test" is performed, the `ecutwfc` here is fixed to the converged value
at "accurate" level. EOS tests need first relax the whole cell, then for each side, shrink or
stretch the lattice constant so that volume change 2% per step, each side 3 steps.

About plotting
Scatter plot, with line

3. Cohesive energy
Cohesive energy is the energy required to break a material into its constituent atoms. It is
a measure of the strength of the bonds in the material. The cohesive energy is calculated by
the difference between the total energy of the bulk material and the total energy of the
isolated atoms.
This test is performed at each `ecutwfc` value, as many as the "convergence test" performed.

About plotting
Line plot, with markers

4. Density of States (DOS)
Density of States (DOS) is a test to determine the electronic structure of a material. DOS
can be used to determine the band gap of a material, and the electronic properties of a
material. The test is performed by directly integrating the `istate.info` file from ecutwfc
converged calculation, only the converged `ecutwfc` is used.

About plotting
Line plot, with shaded area under fermi level
"""
def pseudopotentials(element: str, 
                     xc_functional: str = "PBE", 
                     software: str = "ABACUS",
                     commit: str = "latest commit",
                     fconv: str = None, 
                     fconvlog: str = None, 
                     fdos: str = None, 
                     feos: str = None, 
                     fecoh: str = None):
    """Pseudopotential test"""
    # Liquid header
    html = liquid_header(layout="result", test="Pseudopotential", title=element)
    # Test information
    html += header("Pseudopotential tests", level=1) + "\n"
    html += header("Test information", level=2) + "\n"
    content = [f"element: {element}", f"pseudopotential type: {element}", f"DFT XC (exchange-correlation) functional: {xc_functional}", f"software: {software} (version: {commit})"]
    html += unordered_list(content)
    # Test results
    html += header("Test results", level=2) + "\n"
    html += "<table>\n"
    # BLOCK: Convergence test
    html += "<tr><td>\n"
    content = [f for f in [fconv, fconvlog] if f is not None]
    content = [f"<img src=\"{f}\" class=\"plain-figure\">" for f in content]
    content = "Planewave ecutwfc convergence test data will be available soon." if len(content) == 0 else content
    html += banner_frame(title="Convergence test", content=content)
    html += "</td></tr>\n"
    # BLOCK: Equation of States (EOS)
    # EOS data can be presented with javascript plotly, for user to interact with
    html += "<tr><td>\n"
    content = f"<img src=\"{feos}\" class=\"plain-figure\">" if feos is not None else "Equation of States (EOS) data will be available soon."
    html += banner_frame(title="Equation of States (EOS)", content=content)
    html += "</td></tr>\n"
    # BLOCK: Cohesive energy
    html += "<tr><td>\n"
    content = f"<img src=\"{fecoh}\" class=\"plain-figure\">" if fecoh is not None else "Cohesive energy data will be available soon."
    html += banner_frame(title="Cohesive energy", content=content)
    html += "</td></tr>\n"
    # BLOCK: DOS
    html += "<tr><td>\n"
    content = f"<img src=\"{fdos}\" class=\"plain-figure\">" if fdos is not None else "Density of States (DOS) data will be available soon."
    html += banner_frame(title="Density of States (DOS)", content=content)
    html += "</td></tr>\n"
    # End of table
    html += "</table>\n"
    return html

"""for another page, Pseudopot-Nao, for pseudopotential-nao bundle test
This is a html-string-like template for pseudopotential-nao bundle tests,
presently this page is designed for presenting following tests:
1. Equation of States (EOS)
Importance of EOS is not reviewed again, EOS test here is for determining whether the 
equation of state calculated by numerical atomic orbital ABACUS lcao calculation agrees
with the ABACUS pw calculation.

About plotting
Line plot, to show many EOS curves

2. Cohesive energy


3. Band structure reproduction
Band structure reproduction test is to determine whether the band structure calculated by
numerical atomic orbital ABACUS lcao calculation agrees with the ABACUS pw calculation.
The specific k-path is generated by SeeKpath Python package, and the band structure is
calculated by both ABACUS pw and lcao calculation.

About plotting
Line plot, to show the band structure
"""
def pseudopot_nao(element: str, xc_functional: str = "PBE", 
                  feos: str = None, fecoh: str = None, fband: str = None):
    """Pseudopotential-Nao test
    Not designed..."""
    return ""