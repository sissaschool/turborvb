#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	UPWINVP=$1
	OUT=$2
else
	source ../settings.sh
fi

if [ ! -f "$UPWINVP" ]; then
    echo "Executable $UPWINVP does not exists"
    exit 1
fi

echo " Test upwinvp wrapper " 
echo 1 | $UPWINVP > $OUT
[ $? -eq 0 ] && echo "Run without non-zero exit code" || exit 1

grep OK $OUT
