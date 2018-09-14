#!/bin/bash

if [ $# -eq 0 ] ; then
	echo "./GetStatus.sh <tmux-session-name>"
	echo "(You can find name with 'tmux ls')"
	exit
fi

echo "*********************************** GetStatus ***********************************"

tmux list-windows -t $1 -F '#I'  |   
  while read w; do tmux list-panes -F '#P' -t $w | 
     while read p; do echo -n  "${w}.${p}" ; tmux capture-pane -p -t "${w}.${p}" | 
        tail -n 2 | head -1
     done 
  done

echo "*********************************************************************************"
