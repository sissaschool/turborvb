#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	CHECKDER=$1
else
	source ../settings.sh
fi

echo " Checking AAD derivatives vs finite difference ones (yes_sparse no)"
sed "s/yes_sparse=.true./yes_sparse=.false./g" datasder.d | $CHECKDER > out_der
if [ $? -ne 0 ]
then
    exit 1
fi

if [ $(grep -c ERR out_der) -gt 0 ]; then
  echo "    Errors in output:"
  grep ERR $OUT 
  exit 1
fi
