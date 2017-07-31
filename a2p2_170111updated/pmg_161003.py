from pymatgen import *

import fnmatch
import re
import os

def pmg_converter():
    rule = re.compile(fnmatch.translate('*.cif'), re.IGNORECASE)               # case-insensitive glob
    cif_files =  [name for name in os.listdir('.') if rule.match(name)]        # run in '.' : current directory

    for i in cif_files:
        a = Structure.from_file(i)
        a.to(filename="POSCAR_"+i[:-4]+".vasp")

pmg_converter()
