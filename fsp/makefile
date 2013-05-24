# Makefile for Fortran Smooth particle hydrodynamics program
# Andrew Charles August 2009, June 2010
#
# References:
# http://mrbook.org/tutorials/make/
# http://www.gnu.org/software/make/manual/make.html
# Based on mmd2 makefile by P.Daivis
#
# INSTRUCTIONS
# 1. You must export the environment variable HOSTNAME

ifeq ($(HOSTNAME),LOULI)
	F90C = gfortran
	F90FLAGS = -Wall -fbounds-check -fimplicit-none -fPIC -O3
else ifeq ($(HOSTNAME),spinode)
	F90C = gfortran
	F90FLAGS = -Wall -fbounds-check -fimplicit-none \
	-Wsurprising -Waliasing -fPIC
else ifeq ($(HOSTNAME),tango.vpac.org)
	F90C = pgf95
	F90FLAGS = -O3 -fPIC
	#DEBUGFLAGS = -heap-arrays -check all -warn all -g
else ifeq ($(HOSTNAME),edda)
	F90C = xlf90 -q64 -qsuffix=f=f90
else ifeq ($(HOSTNAME),ho-aifs-ws1440)
	F90C = /home/acharles/local/bin/gfortran
	F90FLAGS = -Wall -fbounds-check -fPIC -LLIBDIR=/home/acharles/local/lib
else ifeq ($(HOSTNAME),acharles-netbook)
	F90C = gfortran
	F90FLAGS = -Wall -fbounds-check -fimplicit-none \
	 -Wsurprising -Waliasing -fPIC
else ifeq ($(HOSTNAME),doctor)
	F90C = gfortran
	F90FLAGS = -Wall -fbounds-check -fimplicit-none \
	 -Wsurprising -Waliasing -fPIC
else
	F90C = gfortran
	F90FLAGS = -Wall -fbounds-check -fimplicit-none \
	 -Wsurprising -Waliasing -fPIC  -arch x86_64 -arch i686
endif

SRCDIR = ./src

SOURCES =  \
	global.f90 \
	art_viscosity.f90 \
	collision.f90 \
	kernel.f90 \
	splib.f90 \
	eos.f90 \
	reader.f90 \
	simulation_box.f90 \
	neighbour_list.f90 \
	particle.f90 \
	adaptive_timestep.f90 \
	core_potential.f90 \
	system_properties.f90 \
	thermostat.f90 \
	writer.f90 \
	density.f90 \
	boundary.f90 \
	sphforce.f90 \
	sphforce_ab.f90 \
	sphstep.f90
	
#test_sph.f90

OBJECTS = $(SOURCES:.f90=.o)
PROG = fsph

# This part of the makefile compactly builds all fortran source
# files into object files using a 'suffix rule'.
# Suffix rules are now obsolete because they are confusing

# Define suffixes.
.SUFFIXES: .f90 .o

# Rule to create object files from FORTRAN source.
.f90.o:
	${F90C} ${F90FLAGS} -c $< 

# Rule for making fsph
${PROG}: $(OBJECTS)
	${F90C} ${F90FLAGS} $@.f90 -o $@ $(OBJECTS)

# Target to make the main program
main: ${SOURCES} ${PROG}

#all: main test 3d pymod doc
all: main test 3d pymod

debug: $(OBJECTS)
	${F90C} ${F90FLAGS} $(DEBUGFLAGS) -o ${PROG} $(OBJECTS)

test: $(OBJECTS) 
	${F90C} $(F90FLAGS) test_sph.f90 -o testsph $(OBJECTS)

pairsep:
	${F90C} $(F90FLAGS) -c pairsep.f90

3d:
	${F90C} $(F90FLAGS) -c global.f90
	${F90C} $(F90FLAGS) -c eos.f90
	${F90C} $(F90FLAGS) -c kernel.f90
	${F90C} $(F90FLAGS) -c art_viscosity.f90
	${F90C} $(F90FLAGS) -c core_potential.f90
	${F90C} $(F90FLAGS) -c splib.f90
	${F90C} $(F90FLAGS) -c sphforce3d.f90
	${F90C} $(F90FLAGS) -c test_force3d.f90
	${F90C} $(F90FLAGS) -o sp3d sphforce3d.o splib.o eos.o kernel.o art_viscosity.o core_potential.o global.o test_force3d.o
	${F90C} $(F90FLAGS) -c time_force3d.f90
	${F90C} $(F90FLAGS) -o t3d time_force3d.o sphforce3d.o splib.o eos.o kernel.o art_viscosity.o core_potential.o global.o

doc:
	./gen_doc.sh

pymod: main 3d
	./wrapfort.sh

#Delete the object files
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.so
	rm -f neighbour.00*
	rm -f snap.00*
	rm -f sphstate.00*
	rm -f gofr.0*
	rm -f properties.output
	rm -f list
	rm -fr test_cases/test_output/*
	rm -fr test_cases/test_report/*
	rm -fr manual/*.log
	rm -fr manual/*.toc
	rm -fr manual/*.brf
	rm -fr manual/*.pdf
	rm -fr manual/*.lof
	rm -fr manual/*.log
	rm -fr manual/*.lot
	rm -fr manual/*.out
	rm -fr manual/*.idx
	rm -fr manual/*.aux
	rm -f *.pyf
	rm -f *.pyc
	rm -f testsph
	rm -f sp3d
	rm -f t3d

put:
	rsync -avtub --exclude "*.o" --exclude "*.mod" --exclude "fsph" --exclude "*.svn*" --exclude "makefile" ~/masters/active/fsph/ acharles@edda.vpac.org:~/fsph/
