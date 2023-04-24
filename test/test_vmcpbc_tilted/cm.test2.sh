#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	CHECKDER=$1
else
	source ../settings.sh
fi

echo " Checking AAD derivatives vs finite difference ones "
$CHECKDER < datasder.d > out_der
if [ $? -ne 0 ]
then
    exit 1
fi

if [ $(grep -c ERR out_der) -gt 0 ]; then
  echo "    Errors in output:"
  grep ERR out_der
  exit 1
fi
