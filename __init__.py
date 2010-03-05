import sys
import os

#Add Pyprop location to path
MODULE_PATH = os.path.dirname(os.path.realpath(__file__))
PypropLocation = os.path.realpath("%s/pyprop" % MODULE_PATH)
sys.path.append(PypropLocation)
import pyprop

#Contains function and class name from package that might appear in config
#files. To be passed to Pyprop for resolving during loading of said files.
from utils import ProjectNamespace, RegisterProjectNamespace

import libpotential
from libpotential import *
for key in libpotential.__dict__.iterkeys():
	if not key.startswith("__"):
		RegisterProjectNamespace(eval(key))

__all__  = ["example", "absorber", "oneelectrontools", "libpotential", "utils", "laserfunctions", "analysis"]
