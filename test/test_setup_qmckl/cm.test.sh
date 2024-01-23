#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	TURBO=$1
	READF=$2
	IN=$3
	OUT=$4
	REF_FORT21=$5
	ROUND_OFF=$6
else
	source ../settings.sh
fi

if [ ! -f "$TURBO" ]; then
    echo "Executable $TURBO does not exists"
    exit 1
fi

$TURBO < $IN

