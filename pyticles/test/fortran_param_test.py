""" This is a simpler one. """

import sys
import fkernel
import splib
import feos
import numpy as np
import sphforce3d
from time import time
from forces import PairwiseForce

adash = 4.0
bdash = 0.5
kbdash = 1000.0

# Are the parameters set by simple assignment?
feos.eos.adash = adash
feos.eos.bdash = bdash
feos.eos.kbdash = kbdash
feos.eos.eos_print_vdw_params()

feos.eos.eos_set_vdw_params(adash,bdash,kbdash)
feos.eos.eos_print_vdw_params()

# Yes they are set. This is because they are defined in the module.
# When they are imported, as for sphforce3d, which gets these parameters
# from eos, they are not mutable.

#sphforce3d.sphforce3d.core_sigma = sigma
#sphforce3d.sphforce3d.core_rcoef = rcoef
#sphforce3d.sphforce3d.cgrad = cgrad
#sphforce3d.sphforce3d.thermalk = thermalk

sphforce3d.sphforce3d.adash = adash
sphforce3d.sphforce3d.bdash = bdash
sphforce3d.sphforce3d.kbdash = kbdash
sphforce3d.sphforce3d.print_vdw_params()

sphforce3d.sphforce3d.set_vdw_params(adash,bdash,kbdash)
sphforce3d.sphforce3d.print_vdw_params()

# The moral of the story, for me, is to avoid imports in fortran
# modules I want to wrap from python, even if this means passing a
# a few parameters as arguments.


