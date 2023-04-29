#!/bin/bash
set -euo pipefail

PREFIX=""
if [[ $# -gt 0 ]]; then
	TURBORVB=$1
	OUT=$2
        if [[ $# -gt 2 ]]; then
	        PREFIX=$3
	fi
fi


#cd vmc
echo " VMC speed test " 
echo "" 
echo " dir = test/vmc" 
echo "" 

$PREFIX $TURBORVB < datasvmc.ds > $OUT

grep "Total time with no initialization =" $OUT

echo "End VMC speed test " 
echo "" 


