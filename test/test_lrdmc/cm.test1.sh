#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	TURBORVB=$1
	READF=$2
	FORCEFN=$3
	IN=$4
	OUT=$5
	TRUEOUT=$6
	REF_FORT21=$7
	ROUND_OFF=$8
        if [[ $# -gt 8 ]]; then
		PREFIX=$9
	else
		PREFIX=""
	fi
else
	source ../settings.sh
fi

if [ ! -f "$TURBORVB" ]; then
    echo "Executable $TURBORVB does not exists"
    exit 1
fi

echo " LRDMC without 3-4body  TEST " 
echo " dir=test_lrdmc" 
$PREFIX $TURBORVB < $IN > $OUT
[ $? -eq 0 ] && echo " Run without non-zero exit code" || exit 1

echo "  Calculate local energies:'0 1 1 1' | readf.x" 
echo "0 1 1 1" | $READF >& /dev/null

if [ $(grep -c ERR $OUT) -gt 0 ]; then
  echo "    Errors in output:"
  grep ERR $OUT 
  exit 1
fi

#check energies
echo "    Rounds off values in fort.21 < 10**-${ROUND_OFF}".
cat fort.21 | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f  %.*f  %.*f\n", ROUND_OFF, $1, ROUND_OFF, $2, ROUND_OFF, $3)}' > fort.21_roundoff
cat ${REF_FORT21} | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f  %.*f  %.*f\n", ROUND_OFF, $1, ROUND_OFF, $2, ROUND_OFF, $3)}' > REFERENCE_fortXXI_roundoff

#diff fort.21
echo "    Compares fort.21_roundoff and REFERENCE_fortXXI_roundoff."
echo "    If you do not see any "diff" here, they are consistent."

if [ $(diff fort.21_roundoff REFERENCE_fortXXI_roundoff | wc -l) -gt 0 ]; then
   diff fort.21_roundoff REFERENCE_fortXXI_roundoff
   exit 1
elif [ $FORCEFN == NA ]; then
   exit 0
else
   exit 0
fi

