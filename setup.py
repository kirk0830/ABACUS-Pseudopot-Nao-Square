import pathlib
import setuptools


here = pathlib.Path(__file__).parent.resolve()
readme = (here / 'README.md').read_text(encoding='utf-8')

# did not include torch and pyscf here
install_requires=["mp_api",
                   "numpy",
                   "pymatgen", 
                   "seekpath", 
                   "ruamel.yaml<0.18",
                   "urllib3<2.1,>=1.25.4"]

setuptools.setup(
    name="apns",
    author="Yike HUANG",
    author_email="huangyk@aici.ac.cn",
    description="ABACUS-Pseudopot-Nao-Square: for ABACUS calculation with higher efficiency and precision",
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(include=['apns', 'apns.*']),
    classifiers=[
        'Development Status :: 3 - Alpha',
        "Programming Language :: Python :: 3.7",
    ],
    install_requires=install_requires,
    python_requires=">=3.7",
    version="0.0.2"
)