#! /bin/bash

cd fsp
make all
cd ../pyticles
./build.sh
cd ../splib
./build.sh


