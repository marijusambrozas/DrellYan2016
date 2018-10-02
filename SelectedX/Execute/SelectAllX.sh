#!/bin/bash

if [ ${#1} -eq 0 ] || [ ${#2} -eq 0 ] ; then
    echo "./SelectAllX_B.sh <WhichX> <Process> <Trigger>"
    exit
fi

WhichX=$1
#echo $WhichX
Trigger=$3
#echo $Trigger

source /cvmfs/cms.cern.ch/cmsset_default.sh
scram p CMSSW CMSSW_8_0_6&
wait %1
cd CMSSW_8_0_6/src
eval `scram runtime -sh`
cd /cms/ldap_home/mambroza/DrellYan2016/SelectedX/


Process=$2
root -l -q -b MakeSelectedX.C'("'$WhichX'", "'$Process'", "'$Trigger'")'
if [[ $Process == *'QCDEM'* ]] ; then
    root -l -q -b MakeSelectedX.C'("QCDfail", "", "'$Trigger'")'&
    wait %1
    root -l -q -b MakeSelectedX.C'("QCDmerge", "", "")'
fi