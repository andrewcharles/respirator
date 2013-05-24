#! /bin/sh

# Wraps certain modules of my Fortran smooth particle code for
# use by python.

# ==========
# REFERENCES
# ==========
# http://cens.ioc.ee/projects/f2py2e/usersguide/index.html#the-smart-way
# http://www.astro.uni-bonn.de/~bertoldi/boa/html/7_4Adding_Fortran90_code.html 

# The REPORT_ON_ARRAY_COPY flag reports when ... well, when a array is copied.
#--g3-numpy is another flag to try in the future with f2py upgrade
#-DF2PY_REPORT_ATEXIT ?

HOSTNAME=`hostname`
if [ $HOSTNAME = "tango.vpac.org" ] ; then
    F2PY_PYF_FLAGS='--overwrite-signature'
	F2PYFLAGS='--fcompiler=pg'
else
    F2PY_PYF_FLAGS='--overwrite-signature'
	F2PYFLAGS='--fcompiler=gnu95 --noarch' # -DF2PY_REPORT_ON_ARRAY_COPY=1000'
fi

f2py -m fkernel -h kernel.pyf $F2PY_PYF_FLAGS kernel.f90
f2py -c kernel.pyf $F2PYFLAGS kernel.f90

f2py -m fpairsep -h fpairsep.pyf $F2PY_PYF_FLAGS pairsep.f90
f2py -c fpairsep.pyf $F2PYFLAGS pairsep.f90

f2py -m art_viscosity -h art_viscosity.pyf $F2PY_PYF_FLAGS art_viscosity.f90 
f2py -c art_viscosity.pyf $F2PYFLAGS art_viscosity.f90

f2py -m collision -h collision.pyf $F2PY_PYF_FLAGS collision.f90
f2py -c collision.pyf $F2PYFLAGS collision.f90

f2py -m feos -h eos.pyf $F2PY_PYF_FLAGS eos.f90
f2py -c eos.pyf $F2PYFLAGS eos.f90

f2py -m fsplib -h splib.pyf $F2PY_PYF_FLAGS splib.f90
f2py -c splib.pyf $F2PYFLAGS splib.f90 kernel.o

f2py -m art_viscosity -h art_viscosity.pyf $F2PY_PYF_FLAGS art_viscosity.f90
f2py -c art_viscosity.pyf $F2PYFLAGS art_viscosity.f90

f2py -m core_potential -h core_potential.pyf $F2PY_PYF_FLAGS core_potential.f90
f2py -c core_potential.pyf $F2PYFLAGS core_potential.f90

f2py -m sphforce3d -h sphforce3d.pyf $F2PY_PYF_FLAGS sphforce3d.f90
f2py -c sphforce3d.pyf $F2PYFLAGS sphforce3d.f90 art_viscosity.o core_potential.o splib.o eos.o kernel.o


# No longer used
#f2py -m fsph_global -h global.pyf --overwrite-signature global.f90
#f2py -c global.pyf $F2PYFLAGS global.f90

