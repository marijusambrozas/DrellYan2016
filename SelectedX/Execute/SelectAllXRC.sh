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
eval `scram runtime -sh`&
wait %1
cd /cms/ldap_home/mambroza/DrellYan2016/SelectedX/

#voms-proxy-init --voms cms

Process=$2
root -l -q -b MakeSelectedX.C'("'$WhichX'", "'$Process'", "'$Trigger'", kTRUE)'&
wait%1
if [[ $Process == *'QCDEM_120to170'* ]] ; then
    root -l -q -b MakeSelectedX.C'("QCDfail", "", "'$Trigger'")'&
    wait %1
    root -l -q -b MakeSelectedX.C'("QCDmerge", "", "")'
fi

exit
