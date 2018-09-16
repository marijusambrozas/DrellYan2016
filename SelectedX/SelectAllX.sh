#!/bin/bash

if [ $# -eq 0 ] ; then
	echo "./SelectAllX.sh <DileptonChannel>"
	echo "(Dilepton channel - 'EE', 'MuMu' or 'EMu')"
	exit
elif [ $1 != 'EE' ] && [ $1 != 'ee' ] && [ $1 != 'MuMu' ] && [ $1 != 'mumu' ] && [ $1 != 'MUMU' ] && [ $1 != 'EMU' ] && [ $1 != 'EMu' ] && [ $1 != 'emu' ] ; then
	echo "./SelectAllX.sh <DileptonChannel>"
	echo "(Dilepton channel - 'EE', 'MuMu' or 'EMu')"
	exit
fi

echo "******************************************** SelectAllX ********************************************"

sessionName="SelectAll$1"
macroLocation="~/DrellYan2016/SelectedX/"
#Uncoment the line below and coment the line above if you want to specify different address for ROOT macro OR just change the value.
#macroLocation=$2

tmux new -s $sessionName -d
echo "tmux session name is $sessionName."
sleep 1
tmux send-keys -t $sessionName "cd ~/" enter
# CHECK IF THIS IS REALLY HOW IT SHOULD BE DONE
tmux send-keys -t $sessionName ". ~/.bash_profile" enter
tmux send-keys -t $sessionName "cd $macroLocation" enter
tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"DY$1 10to200\\\")" enter

echo "Waiting 30 seconds in case macro needs to be compiled..."
sleep 30
echo "Please wait ~10 more seconds. Opening more tmux windows..."

tmux new-window -n "window1"
sleep 1
tmux send-keys -t $sessionName "cd ~/" enter
tmux send-keys -t $sessionName ". ~/.bash_profile" enter
tmux send-keys -t $sessionName "cd $macroLocation" enter
tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"DY$1 200toInf\\\")" enter

tmux new-window -n "window2"
sleep 1
tmux send-keys -t $sessionName "cd ~/" enter
tmux send-keys -t $sessionName ". ~/.bash_profile" enter
tmux send-keys -t $sessionName "cd $macroLocation" enter
tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"DYTauTau_full\\\")" enter

tmux new-window -n "window3"
sleep 1
tmux send-keys -t $sessionName "cd ~/" enter
tmux send-keys -t $sessionName ". ~/.bash_profile" enter
tmux send-keys -t $sessionName "cd $macroLocation" enter
tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"ttbar_full\\\")" enter

tmux new-window -n "window4"
sleep 1
tmux send-keys -t $sessionName "cd ~/" enter
tmux send-keys -t $sessionName ". ~/.bash_profile" enter
tmux send-keys -t $sessionName "cd $macroLocation" enter
tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"VVnST\\\")" enter

tmux new-window -n "window5"
sleep 1
tmux send-keys -t $sessionName "cd ~/" enter
tmux send-keys -t $sessionName ". ~/.bash_profile" enter
tmux send-keys -t $sessionName "cd $macroLocation" enter
tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"WJets\\\")" enter

tmux new-window -n "window6"
sleep 1
tmux send-keys -t $sessionName "cd ~/" enter
tmux send-keys -t $sessionName ". ~/.bash_profile" enter
tmux send-keys -t $sessionName "cd $macroLocation" enter
tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"ttbar_full\\\")" enter

if [ $1 == 'EE' ] || [ $1 == 'ee' ] ; then
	tmux new-window -n "window7"
	sleep 1
	tmux send-keys -t $sessionName "cd ~/" enter
	tmux send-keys -t $sessionName ". ~/.bash_profile" enter
	tmux send-keys -t $sessionName "cd $macroLocation" enter
	tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"QCDEM_full\\\")" enter

	tmux new-window -n "window8"
	sleep 1
	tmux send-keys -t $sessionName "cd ~/" enter
	tmux send-keys -t $sessionName ". ~/.bash_profile" enter
	tmux send-keys -t $sessionName "cd $macroLocation" enter
	tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"DoubleEG_full\\\")" enter

	tmux new-window -n "window9"
	sleep 1
	tmux send-keys -t $sessionName "cd ~/" enter
	tmux send-keys -t $sessionName ". ~/.bash_profile" enter
	tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"SingleElectron_full\\\")" enter

	tmux new-window -n "window10"
	sleep 1
	tmux send-keys -t $sessionName "cd ~/" enter
	tmux send-keys -t $sessionName ". ~/.bash_profile" enter
	tmux send-keys -t $sessionName "cd $macroLocation" enter
	tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"QCDfail\\\")&" enter
	tmux send-keys -t $sessionName "wait %1" enter
	tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"QCDmerge\\\")" enter
else
	tmux new-window -n "window7"
	sleep 1
	tmux send-keys -t $sessionName "cd ~/" enter
	tmux send-keys -t $sessionName ". ~/.bash_profile" enter
	tmux send-keys -t $sessionName "cd $macroLocation" enter
	tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"QCDMu_full\\\")" enter

	tmux new-window -n "window8"
	sleep 1
	tmux send-keys -t $sessionName "cd ~/" enter
	tmux send-keys -t $sessionName ". ~/.bash_profile" enter
	tmux send-keys -t $sessionName "cd $macroLocation" enter
	tmux send-keys -t $sessionName "root -l -q -b MakeSelectedX.C+(\\\"$1\\\", \\\"SingleMuon_full\\\")" enter
fi

echo "All the ROOT macros are set up and running."
echo "Type in 'tmux a -t $sessionName' to attach the tmux session."
echo "You can check running tmux sessions with 'tmux ls'."
echo "When scripts finish running you can kill the tmux session with 'tmux kill-session -t $sessionName'."
echo "****************************************************************************************************"
