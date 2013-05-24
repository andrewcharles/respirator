# Makefile for Pyticles
# Andrew Charles August 2009, June 2010
export PATH:=/home/acharles/rsp/pyticles:$(PATH)

cython:
	$(PYTICLES_BASE)/build.sh
	./build.sh

all: cython

#Delete the object files
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.so
	rm -f *.pyf
	rm -f *.pyc
	rm -fr build/*

