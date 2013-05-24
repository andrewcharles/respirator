#! /bin/bash

# Sets up environment variables for RSP
# The point of this is to have directories for data writing, 
# image writing and source code control and not have to specify paths manually
# in multiple different codes and scripts.

export HOSTNAME=`hostname`
echo Hello $USER, I am intialising an environment for you.


export PHDIR=/Users/acharles/phd/rsp
if [ $HOSTNAME = "acharles-netbook" ] ; then
	echo 'Running on '$HOSTNAME
	export PHDIR=/home/acharles/rsp
elif [ $HOSTNAME == 'spinode' ]; then
	echo 'Running on '$HOSTNAME
	export PHDIR=/Users/acharles/phd/rsp
	export TEXPATH=/usr/local/texlive/2007/bin/i386-darwin
	NCPATH=/opt/netcdf4/bin 
	path=($VSP_BASE $TEXPATH $path)
	alias epython='/opt/epd-6.2-2/bin/ipython'
	export PATH=$PATH:/opt/epd-6.2-2/bin
    PYTHON=/Library/Frameworks/Python.framework/Versions/Current/bin/python
    export PYTHONPATH=$PYTHONPATH:$PYTHON/lib/python2.7/site-packages
elif [ $HOSTNAME == 'cyclonic.bom.gov.au' ]; then
	export PHDIR=/Users/acharles/phd/rsp
	NCPATH=/opt/netcdf4/bin 
    PYTHON=/Library/Frameworks/Python.framework/Versions/2.7
elif [[ $HOSTNAME == tango* ]]; then
	echo 'Running on '$HOSTNAME
	export PHDIR=/home/acharles/rsp
elif [[ $HOSTNAME == LOULI ]]; then
	echo 'Running on '$HOSTNAME
	export PHDIR=/home/ac/rsp
elif [[ $HOSTNAME == doctor ]]; then
	echo 'Running on '$HOSTNAME
	export PHDIR=/home/acharles/rsp
    export PYTHONPATH=/home/acharles/rsp/extern/lib64/python2.7/site-packages:$PYTHONPATH
elif [[ $HOSTNAME == ANDROLISK ]]; then
	echo 'Running on '$HOSTNAME
	export PHDIR=/home/$USER/rsp
elif [[ $HOSTNAME == ho-aifs-ws1440 ]]; then
    export CC=/home/acharles/local/bin/gcc
    export CXX=/home/acharles/local/bin/g++
	export PHDIR=/data/acharles/rsp
    export PATH=/data/acharles/rsp/extern/cmake-2.8.5/bin:${PATH}
    source /home/acharles/pythonpath7.sh
    #export EXTERN_BASE=$PHDIR/extern/lib/python2.6/site-packages
    export PYTHONPATH=/home/acharles/local/python/2.7.2/lib/python2.7/site-packages:$PYTHONPATH
fi

export RSP_BASE=`pwd`
# Set the paths for fortran sph run file
export NPY_NO_DEPRECATED_API=True
#export EXTERN_BASE=$PHDIR/extern/lib/python2.6/site-packages
export FSP_BASE=$RSP_BASE/fsp
export FSP_RUN=$RSP_BASE/data/fsp_runs
export VSP_BASE=$RSP_BASE/vsp
export VSPL_BASE=$RSP_BASE/splib
export VSPL_BASE=$VSPL_BASE:$RSP_BASE/vsp/scripts
export VSPL_BASE=$VSPL_BASE:$RSP_BASE/vsp/fsp_post
export VSPL_BASE=$VSPL_BASE:$RSP_BASE/vsp/vdwlib
export SPDATA=$RSP_BASE/data
export PYTICLES_BASE=$RSP_BASE/pyticles
export PATH=$PYTICLES_BASE:$FSP_BASE:$VSP_BASE:$PATH
export PYTHONPATH=$EXTERN_BASE:${RSP_BASE}:${VSP_BASE}:${VSPL_BASE}:${FSP_BASE}:${PYTICLES_BASE}:${PYTHONPATH}
export PYTHONPATH=/Users/acharles/local/lib/python2.7/site-packages:${PYTHONPATH}



