#!/bin/bash
cd /cms/ldap_home/mambroza/DrellYan2016/SelectedX/Execute/FRSelection_mu

#	condor_submit Submit_WJets_ext2v5.jds
#	sleep 30
#	condor_submit Submit_DY_50to100.jds
#	sleep 30
#	condor_submit Submit_ttbar.jds
#	sleep 30
condor_submit Submit_SingleMuon_B.jds
sleep 30
#condor_submit Submit_SingleMuon_G.jds
#sleep 30
condor_submit Submit_SingleMuon_H.jds
sleep 30
#condor_submit Submit_SingleMuon_C.jds
#sleep 30
#condor_submit Submit_SingleMuon_D.jds
#sleep 30
condor_submit Submit_SingleMuon_E.jds
sleep 30
#condor_submit Submit_SingleMuon_F.jds
#sleep 30
#condor_submit Submit_DY_10to50.jds
#sleep 30
#condor_submit Submit_WJets.jds
#sleep 30
#condor_submit Submit_ttbar_700to1000.jds
#sleep 30
#condor_submit Submit_ttbar_1000toInf.jds
#sleep 30
#condor_submit Submit_DY_100toInf.jds
#sleep 30
#condor_submit Submit_tW.jds
#sleep 30
#condor_submit Submit_tbarW.jds
#sleep 30
#condor_submit Submit_WW.jds
#sleep 30
#condor_submit Submit_WZ.jds
#sleep 30
#condor_submit Submit_ZZ.jds
#sleep 30
#condor_submit Submit_QCDMu_15to20.jds
#sleep 30
#condor_submit Submit_QCDMu_20to30.jds
#sleep 30
#condor_submit Submit_QCDMu_30to50.jds
#sleep 30
#condor_submit Submit_QCDMu_50to80.jds
#sleep 30
#condor_submit Submit_QCDMu_80to120.jds
#sleep 30
#condor_submit Submit_QCDMu_120to170.jds
#sleep 30
#condor_submit Submit_QCDMu_170to300.jds
#sleep 30
#condor_submit Submit_QCDMu_300to470.jds
#sleep 30
#condor_submit Submit_QCDMu_470to600.jds
#sleep 30
#condor_submit Submit_QCDMu_600to800.jds
#sleep 30
#condor_submit Submit_QCDMu_800to1000.jds
#sleep 30
#condor_submit Submit_QCDMu_1000toInf.jds
#
cd -
