#!/bin/bash

if [ ${#1} -eq 0 ] || [ ${#2} -eq 0 || ${#3} -eq 0 ] ; then
    echo "./SelectAllX.sh <WhichX> <Process> <Trigger>"
    exit
fi

WhichX=$1
#echo $WhichX
Trigger=$3
#echo $Trigger

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc530
source $VO_CMS_SW_DIR/cmsset_default.sh
cd $VO_CMS_SW_DIR/$SCRAM_ARCH/cms/cmssw/CMSSW_8_0_32/src
eval `scramv1 runtime -sh`

cd /cms/ldap_home/mambroza/DrellYan2016/SelectedX/

Process=$2
root -l -q -b MakeSelectedX.C'("'$WhichX'", "'$Process'", "'$Trigger'")'&

exit
