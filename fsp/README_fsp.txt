===================
Fortran Soure Files
===================

run.py

makefile -- run make. 
testall.py -- run all tests


Makefile Notes
--------------
The style of the makefile is a little old, but it does the job. Prefer
manually specifying each filename for clarity and so that dependencies are
clear rather than implicit in the order of source file listing.

make clean will remove all object files and other junk so you can start
from scratch.


Run control programs
--------------------
fsph.f90 -- Main module for running 2D simulations
global.f90 -- Configuration variables used by numerous modules.
fsph3d.f90 -- A test module for 3D. Could be extended to run a full sim.


Sub-Program Modules
-------------------
Modules that use neighbour list derived type. Modules that use the particle
derived type. These modules implement moderately complex algorithms.

density.f90
adaptive_timestep.f90
boundary.f90



Stand-alone modules with only primitive argument subroutines
------------------------------------------------------------
These modules do not use structures such as the particle or neighbour list.

art_viscosity.f90 -- Artificial viscosity.
collision.f90 -- Hard collisions.
core_potential.f90 -- Repulsive potential for reducing particle clumping.
eos.f90 -- Equations of state.
splib.f90 -- Smooth particle library.
art_viscosity.f90
boundary.f90


Test Modules
------------
These exist as short programs to test other modules, not in the formal unit
test sense but in the more informal 'this should run' kind of test.

test_sph.f90 -- this is pretty much the run control script fsph, but simpler
test_force3d.f90 -- force3d is a large and complicated module. This simply
tests that it can get out the starting blocks.



Python wrappers for stand-alone modules
---------------------------------------
Run these scripts to build and test the python bindings. Python bindings only
exist for fortran code with primitive argument subroutines.

wrapfort.py -- build the fortran
wrapspect.py -- working tests, plots etc
wraptest.py -- compliance tests
wrapcalc.py -- there may be a test case needing to be migrated to wraptest


Manual
------
In the manual subdirectory is a latex document with detailed specifications
of the code and its implementation. While in hindsight I'd prefer this to
have been in restructured text though a conversion seems unlikely. The
manual is purely based on the fortran code (with some discussion of the run
scripts), and doesn't cover the python wrapping, which is described in this
readme file.


Python wrapping notes
=====================
These are notes on how to make f2py wrappers.

If you want variables written to and the values returned, always
use a 

  !f2py intent(in,out,overwrite) :: u

This is important! You MUST ensure that all subroutines intended to return
values have an intent statement in the variable declaration.

Also, how it is passed to fortran is very, very important.
And a bit fiddly. And seems to require trial and error.

    #u = np.reshape(p.u[0:n].copy(),(n,1))#,order='F')
    u = p.u[0:n].reshape((n))#,order='F')

Design notes
------------
Where possible have constants declared at the module level.
Because this is a global namespace when imported in fortran,
ensure that the name is not likely to clash (e.g., K,g are bad)

GOTCHAS
=======
When arguments are variable shape, you cannot pass

np.array(1.0) 

rather you need to pass

np.array([1.0])

You cannot change the value of a module level variable in an imported
module directly. That is, if the module called from python imports another
module, that module's variable cannot be changed by simple python statements.

Miscellaneous
=============
scraps.f90 -- code fragments that may come in handy one day?


