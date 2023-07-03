#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	CHECKDER=$1
else
	source ../settings.sh
fi

echo " Checking AAD derivatives vs finite difference ones "
$CHECKDER < datasder.d > out_der
echo `grep ERR out_der`

if [ $(grep -c ERR out_der) -gt 0 ]; then
  echo "    Errors in output:"
  echo `grep ERR $OUT`
  exit 1
fi
