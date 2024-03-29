#!/bin/bash
set -euo pipefail

PREFIX=""
if [[ $# -gt 0 ]]; then
	PREP=$1
	OUT=$2
        if [[ $# -gt 2 ]]; then
	        PREFIX=$3
	fi
fi

echo ""
echo "DFT PBC Twist complex speed test"
echo ""

$PREFIX $PREP < prep.ds > $OUT

grep "Total initialization time (sec.)" $OUT
grep "Total loading time matrices (sec.)" $OUT
grep "Total self-consistent cycle time (sec.)" $OUT 

echo "DFT PBC Twist complex speed test end"
echo ""
