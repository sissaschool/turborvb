#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	ZSKTRS=$1
	UPLO=$2
	OUT=$3
else
	source ../settings.sh
fi

if [ ! -f "$ZSKTRS" ]; then
    echo "Executable $ZSKTRS does not exists"
    exit 1
fi

echo " Test zsktrs wrapper " 
cat s$UPLO | $ZSKTRS > $OUT
[ $? -eq 0 ] && echo "Run without non-zero exit code" || exit 1

grep OK $OUT
