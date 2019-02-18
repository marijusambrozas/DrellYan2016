export SCRAM_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

## ROOT thisroot.sh PATH
ROOT_PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_20
cd $ROOT_PATH/src

## cmsenv
eval `scramv1 runtime -sh`
