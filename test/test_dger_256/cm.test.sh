#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	TZGEMM=$1
	OUT=$2
else
	source ../settings.sh
fi

if [ ! -f "$TZGEMM" ]; then
    echo "Executable $TZGEMM does not exists"
    exit 1
fi

echo " Test zgeru_ wrapper " 
cat s | $TZGEMM > $OUT
[ $? -eq 0 ] && echo "Run without non-zero exit code" || exit 1

grep OK $OUT
