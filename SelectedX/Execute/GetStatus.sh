#!/bin/bash

echo "**************************** GetStatus ****************************"
tail -n 1 ../Status/*
echo "*******************************************************************"

#if [ $# -eq 0 ] ; then
#	echo "./GetStatus.sh <tmux-session-name>"
#	echo "(You can find name with 'tmux ls')"
#	exit
#fi
#
#echo "*********************************** GetStatus ***********************************"
#
#stop=0
#line=1
#
#tmux list-windows -t $1 -F '#I'  |   
#  while read w; do tmux list-panes -F '#P' -t $w | 
#     while read p; do
#	echo -n  "${w}.${p}"
#	z=$(tmux capture-pane -p -t "${w}.${p}" | tail -n $line)
#	while [ -z "$z" ] || [[ $z == *"]$" ]] ; do
#		((line++))
#		z=$(tmux capture-pane -p -t "${w}.${p}" | tail -n $line | head -1)
#	done
#	echo $z
#     done 
#  done
#
#echo "*********************************************************************************"
