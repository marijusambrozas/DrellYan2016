#!/bin/bash

if [ $KP_ANALYZER_PATH ]; then
    echo "KP_ANALYZER_PATH is already defined: use a clean shell!"
    return 1
fi

# -- analyzer path (do not insert "/" in the end of the path)-- #
export KP_ANALYZER_PATH=$(pwd)
export KP_INCLUDE_PATH=$KP_ANALYZER_PATH/Include
export KP_ROOTFILE_PATH=$KP_ANALYZER_PATH/Results_ROOTFiles

# -- root setup: only works from root6 -- #
export ROOT_INCLUDE_PATH=${KP_ANALYZER_PATH}:${ROOT_INCLUDE_PATH}

# -- Directory for official style plots -- #
export KP_PLOT_PATH=$KP_ANALYZER_PATH/Local/Plots

# -- ntuple path -- #
export KP_DATA_PATH=""
if [ $HOSTNAME == "tamsa2.snu.ac.kr" ]; # -- 147.47.242.67 -- # 
then 
	# KP_DATA_PATH="/data5/DrellYan2016"
	# KP_DATA_PATH="/data5/DYntuple/v2.0"
	KP_DATA_PATH="/data9/DATA/DYntuple/v2.0"

	echo "cmsenv under CMSSW_8_0_26 (for ROOT6 setup) ..."

	# -- cmssw setup (for ROOT6 & compatible with RooUnfold in tamsa2) -- #
	export SCRAM_ARCH=slc6_amd64_gcc530
	export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
	source $VO_CMS_SW_DIR/cmsset_default.sh

	cd /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_26
	eval `scramv1 runtime -sh` # -- equivalent to cmsenv (cmsenv doesn't work. why?) -- #
	cd $KP_ANALYZER_PATH
		
elif [ $HOSTNAME == "cms.snu.ac.kr" ]; then
	KP_DATA_PATH=""

elif [ $HOSTNAME == "muon" ]; then
	KP_DATA_PATH="/scratch/kplee/DrellYan2016"
	echo "Muon server: please set [scl enable devtoolset-2 bash] to use ROOT6"
fi

echo "================ environment ================"
echo "KP_ANALYZER_PATH:" ${KP_ANALYZER_PATH}
echo "KP_INCLUDE_PATH:" ${KP_INCLUDE_PATH}
echo "KP_ROOTFILE_PATH:" ${KP_ROOTFILE_PATH}

echo "KP_DATA_PATH:" ${KP_DATA_PATH}
if [ -z $KP_DATA_PATH ]; then
    echo "     [WARNING]: ntuples are not available in this machine"
fi

echo "KP_PLOT_PATH:" ${KP_PLOT_PATH}
echo "============================================="
echo "setup is finished. Welcome :)"
