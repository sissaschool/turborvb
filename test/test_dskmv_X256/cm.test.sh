#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	TDGEMM=$1
	UPLO=$2
	OUT=$3
else
	source ../settings.sh
fi

if [ ! -f "$TDGEMM" ]; then
    echo "Executable $TDGEMM does not exists"
    exit 1
fi

echo " Test zskmv_ wrapper " 
cat s$UPLO | $TDGEMM > $OUT
[ $? -eq 0 ] && echo "Run without non-zero exit code" || exit 1

grep OK $OUT
