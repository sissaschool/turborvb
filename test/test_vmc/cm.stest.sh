#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	TURBORVB=$1
	OUT=$2
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


