#!/bin/bash

cd /cms/ldap_home/mambroza/DrellYan2016/SelectedX/etc/RoccoR
rm *.d
rm *.so
rm *_cxx
rm *.pcm
cd -

cd /cms/ldap_home/mambroza/DrellYan2016/SelectedX/Execute/JetSelection_ele

condor_submit Submit_WJets_ext2v5.jds
sleep 30
condor_submit Submit_DY_50to100.jds
sleep 30
condor_submit Submit_ttbar.jds
sleep 30
condor_submit Submit_DoubleEG_B.jds
sleep 30
condor_submit Submit_DoubleEG_G.jds
sleep 30
condor_submit Submit_DoubleEG_H.jds
sleep 30
condor_submit Submit_DoubleEG_C.jds
sleep 30
condor_submit Submit_DoubleEG_D.jds
sleep 30
condor_submit Submit_DoubleEG_E.jds
sleep 30
condor_submit Submit_DoubleEG_F.jds
sleep 30
condor_submit Submit_DY_10to50.jds
sleep 30
condor_submit Submit_WJets.jds
sleep 30
condor_submit Submit_ttbar_700to1000.jds
sleep 30
condor_submit Submit_ttbar_1000toInf.jds
sleep 30
condor_submit Submit_DY_100toInf.jds
sleep 30
condor_submit Submit_tW.jds
sleep 30
condor_submit Submit_tbarW.jds
sleep 30
condor_submit Submit_WW.jds
sleep 30
condor_submit Submit_WZ.jds
sleep 30
condor_submit Submit_ZZ.jds
sleep 30
condor_submit Submit_QCDEM_20to30.jds
sleep 30
condor_submit Submit_QCDEM_30to50.jds
sleep 30
condor_submit Submit_QCDEM_50to80.jds
sleep 30
condor_submit Submit_QCDEM_80to120.jds
sleep 30
condor_submit Submit_QCDEM_120to170.jds
sleep 30
condor_submit Submit_QCDEM_170to300.jds
sleep 30
condor_submit Submit_QCDEM_300toInf.jds

cd -
