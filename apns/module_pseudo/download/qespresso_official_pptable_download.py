"""WARNING
Please do not make the sleep time too short, otherwise you may be blocked by the server.
Politeness is the key to success.

TO USE:
pseudopotential_kind: str, "hgh", "psl", "original-qe", "fhi-aims"
functional: str, "pbe", "pbesol", "pz"

"""

import bs4
import urllib.request as ur

DOMAIN = {
            "hgh": "hartwigesen-goedecker-hutter-pp", 
            "psl": "ps-library", 
            "original-qe": "original-qe-pp-library", # not recommended
            "fhi-aims": "fhi-pp-from-abinit-web-site" # not recommended
            }
# you must change the header to avoid being blocked by the server, windows 11
HEADER = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) \
                        AppleWebKit/537.36 (KHTML, like Gecko) \
                        Chrome/80.0.3987.149 Safari/537.36'}
# do not send too many requests to the server in a short time
SLEEP = 60 # seconds, this number appears to be safe for not bothering QE server too much

import numpy as np
import time
def sleep(mode: str = "switch-element"):

    mode_factors = {
        "click-download": 1/3, "switch-element": 1, "switch-table": 1
    }
    time_to_sleep = np.random.rand() * SLEEP * mode_factors[mode]
    print("You should not bother too much. Sleep for {} seconds...".format(time_to_sleep))
    time.sleep(time_to_sleep)

def get_available_elements(pseudopotential_kind: str):

    # Get the url
    url = "https://pseudopotentials.quantum-espresso.org/legacy_tables/" + DOMAIN[pseudopotential_kind]
    # Get the html source with specified header but... how?
    html = ur.urlopen(url).read()

    # Parse the html source
    soup = bs4.BeautifulSoup(html, features="html.parser")

    # find all tags whose class is "element_anchor"
    elements = soup.find_all("a", {"class": "element_anchor"})
    elements = [element.text.strip() for element in elements if element.text.strip()[0].isalpha()]

    return elements

import re
def get_available_pseudopotentials(element: str, 
                                   pseudopotential_kind: str, 
                                   functional: str, 
                                   pattern = r"([A-Z][a-z]?)(\.)(rel\-)?([\w]+)(\-)(.*)(.UPF)$"):

    element_page_header = "https://pseudopotentials.quantum-espresso.org/legacy_tables/" + DOMAIN[pseudopotential_kind] # then followed by element name

    url = element_page_header + "/" + element.lower()
    print("Getting available pseudopotentials for %s from %s" % (element, url))
    html = ur.urlopen(url).read()

    soup = bs4.BeautifulSoup(html, features="html.parser")

    pseudopotentials = soup.find_all("a", {"class": "element_anchor"})
    pseudopotentials = [pseudopotential.text.strip() for pseudopotential in pseudopotentials]
    print("Found {} pseudopotentials".format(len(pseudopotentials)))
    for pseudopotential in pseudopotentials:
        print(pseudopotential)
    # filter by functional, the 3rd domain in the file name
    pseudopotentials = [pseudopotential for pseudopotential in pseudopotentials if re.match(pattern, pseudopotential).group(4) == functional]

    return pseudopotentials

import os
def download_pseudopotential(pseudopotentials: list):
    pseudopotential_page_header = "https://pseudopotentials.quantum-espresso.org/upf_files/" # then followed by exact file names
    for pseudopotential in pseudopotentials:
        furl = pseudopotential_page_header + pseudopotential
        print("Downloading {}".format(pseudopotential))
        if os.path.exists(pseudopotential):
            print("{} already exists. Skip".format(pseudopotential))
            continue
        sleep("click-download")
        ur.urlretrieve(furl, pseudopotential)
        print("Done")

def driver(pseudopotential_kind: str,
           functional: str,
           elements: list,
           startsfrom: str = "H",
           destination: str = "./download/pseudopotentials/"):
    """Main function to download pseudopotentials from Quantum ESPRESSO official website pptable
    
    Args:
        pseudopotential_kind (str): "hgh", "psl", "original-qe", "fhi-aims"
        functional (str): "pbe", "pbesol", "pz"
        elements (list, optional): Defaults to []. If empty, all elements available will be downloaded. Otherwise, only elements in the list will be downloaded.
        startsfrom (str, optional): Defaults to "H".
    """
    if elements == []:
        elements = get_available_elements(pseudopotential_kind)

    startsfrom = startsfrom.capitalize()
    elements = elements[elements.index(startsfrom):]

    fnames = []
    for element in elements:
        sleep("switch-element")
        pseudopotentials = get_available_pseudopotentials(element, pseudopotential_kind, functional)
        download_pseudopotential(pseudopotentials)
        fnames += pseudopotentials
    
    # move all downloaded pseudopotentials to destination
    path_backup = os.getcwd()
    if not os.path.exists(destination):
        print("Warning: destination {} does not exist. Create it.".format(destination))
        os.makedirs(destination)
    
    destination = destination + "/" if not destination.endswith("/") else destination
    destination += DOMAIN[pseudopotential_kind] + "/"
    os.makedirs(destination, exist_ok=True)
    
    for fname in fnames:
        os.rename(fname, destination + fname)
    os.chdir(path_backup)
    
if __name__ == "__main__":

    """example1: download all pseudopotentials for all elements of PBE functional of HGH kind, but download is interrupted, so restart from element At"""
    #driver(pseudopotential_kind="hgh", functional="pbe", startsfrom="At")
    """example2: download all pseudopotentials for H of PBE functional of pslibrary kind"""
    #driver(pseudopotential_kind="psl", functional="pbe", elements=["H"])