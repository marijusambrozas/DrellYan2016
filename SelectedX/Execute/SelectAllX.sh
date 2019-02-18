#!/bin/bash

if [ ${#1} -eq 0 ] || [ ${#2} -eq 0 ] ; then
    echo "./SelectAllX_B.sh <WhichX> <Process> <Trigger>"
    exit
fi

WhichX=$1
#echo $WhichX
Trigger=$3
#echo $Trigger

#source /cvmfs/cms.cern.ch/cmsset_default.sh
#source /cms/ldap_home/mambroza/DrellYan2016/SelectedX/Execute/cmsset_default.sh&
#wait %1
#scram p CMSSW CMSSW_8_0_6&
#wait %1
#cd CMSSW_8_4_0/src
#eval `scram runtime -sh`&
#wait %1
/cms/ldap_home/mambroza/DrellYan2016/SelectedX/Execute/setup.sh

cd /cms/ldap_home/mambroza/DrellYan2016/SelectedX/
#ls
#voms-proxy-init --voms cms

Process=$2
root -l -q -b MakeSelectedX.C'("'$WhichX'", "'$Process'", "'$Trigger'")'&
wait %1
if [[ $Process == *'QCDEM_120to170'* ]] ; then
    root -l -q -b MakeSelectedX.C'("QCDfail", "", "'$Trigger'")'&
    wait %1
    root -l -q -b MakeSelectedX.C'("QCDmerge", "", "")'
fi

exit
