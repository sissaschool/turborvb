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
echo " lrdmc speed test " 
echo "" 
echo " dir = test/lrdmc_new" 
echo "" 

$PREFIX $TURBORVB < datasfn.ds > $OUT

grep "Total time with no initialization" $OUT
grep "Total time with no measures" $OUT

echo "" 
echo "End VMC speed test " 
echo "" 


