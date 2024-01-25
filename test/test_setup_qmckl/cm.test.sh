#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	TURBO=$1
	READF=$2
	IN=$3
	OUT=$4
	REF_FORT21=$5
	ROUND_OFF=$6
else
	source ../settings.sh
fi

if [ ! -f "$TURBO" ]; then
    echo "Executable $TURBO does not exists"
    exit 1
fi

cat $IN | $TURBO | tee > $OUT

[ $? -eq 0 ] && echo " Run without non-zero exit code" || exit 1

if [ $(grep -c ERR $OUT) -gt 0 ]; then
  echo "    Errors in output:"
  echo `grep ERR $OUT`
  exit 1
fi

#check energies
echo "  Calculate local energies:'0 1 1 1' | readf.x" 
echo "0 1 1 1" | $READF >& /dev/null

echo "    Rounds off values in fort.21 < 10**-${ROUND_OFF}".
cat fort.21 | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f  %.*f  %.*f\n", ROUND_OFF, $1, ROUND_OFF, $2, ROUND_OFF, $3)}' > fort.21_roundoff
cat ${REF_FORT21} | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f  %.*f  %.*f\n", ROUND_OFF, $1, ROUND_OFF, $2, ROUND_OFF, $3)}' > REFERENCE_fortXXI_roundoff
#
#diff fort.21
echo "    Compares fort.21_roundoff and REFERENCE_fortXXI_roundoff."
echo "    If you do not see any "diff" here, they are consistent."

if [ $(diff fort.21_roundoff REFERENCE_fortXXI_roundoff | wc -l) -gt 0 ]; then
   diff fort.21_roundoff REFERENCE_fortXXI_roundoff
   exit 1
fi
