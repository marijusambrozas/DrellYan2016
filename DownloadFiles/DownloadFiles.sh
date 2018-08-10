#!/bin/bash

eos=/usr/bin/xrdcp

server="cms-xrdr.sdfarm.kr:1094"
path="/xrd/store/user/dpai/"

if [ ${#2} -eq 0 ] ; then
    echo "./DownloadFiles.sh destPath filenames_file(s)"
    exit
fi

destPath=$1
shift 1

for f in `cat $@` ; do
    echo "File to be downoladed:"
    echo $f
#    d=`echo $f | awk -F'/' '{print $6$7"/"$8"/"$11}'`
    d=`echo $f | awk -F'/' '{print $9}'`
#    d=${d/_v2p3/v2p3}
#    dpath1=`echo $d | awk -F'/' '{print $1}'`
#    dpath2=`echo $d | awk -F'/' '{print $1"/"$2}'`
    echo "Will be downloaded to:"
    echo "${destPath}$d"
#    echo "     $dpath1 , $dpath2"

#    mkdir ${destPath}/$dpath1 2>/dev/null
#    mkdir ${destPath}/$dpath2 2>/dev/null
    echo "Start download?"
    read isOk
    if [[ "$isOk" == "y" || "$isOk" == "yes" ]] ; then
         ${eos} root://$server//$f ${destPath}$d
    elif [[ "$isOk" == "n" || "$isOk" == "no" ]] ; then
	echo "Uh.. Ok."
    else
        echo "Type \"yes\" or \"no\""
    fi
done
