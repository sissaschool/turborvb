#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	UPWINVP_COMPLEX=$1
	OUT=$2
else
	source ../settings.sh
fi

if [ ! -f "$UPWINVP_COMPLEX" ]; then
    echo "Executable $UPWINVP_COMPLEX does not exists"
    exit 1
fi

echo " Test upwinvp_complex wrapper " 
echo 1 | $UPWINVP_COMPLEX > $OUT
[ $? -eq 0 ] && echo "Run without non-zero exit code" || exit 1

grep OK $OUT
