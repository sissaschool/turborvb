#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	ZSKTRI=$1
	UPLO=$2
	OUT=$3
else
	source ../settings.sh
fi

if [ ! -f "$ZSKTRI" ]; then
    echo "Executable $ZSKTRI does not exists"
    exit 1
fi

echo " Test zsktri wrapper " 
cat s$UPLO | $ZSKTRI > $OUT
[ $? -eq 0 ] && echo "Run without non-zero exit code" || exit 1

grep OK $OUT
