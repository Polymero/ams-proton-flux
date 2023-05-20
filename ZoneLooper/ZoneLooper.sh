#!/bin/bash

export ROOTSYS=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.20.08/x86_64-centos7-gcc48-opt
source $ROOTSYS/bin/thisroot.sh

root -b -q '/afs/cern.ch/user/s/svenenda/public/ams-proton-flux/ZoneLooper/ZoneLooper.C('$1')'

