#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	TDGEMM=$1
	OUT=$2
else
	source ../settings.sh
fi

if [ ! -f "$TDGEMM" ]; then
    echo "Executable $TDGEMM does not exists"
    exit 1
fi

echo " Test dgemm_ wrapper " 
cat s | $TDGEMM > $OUT
[ $? -eq 0 ] && echo "Run without non-zero exit code" || exit 1

grep OK $OUT
