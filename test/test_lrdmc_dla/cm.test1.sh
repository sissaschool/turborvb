#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	TURBORVB=$1
	READF=$2
	FORCEFN=$3
	OUT=$4
	TRUEOUT=$5
	REF_FORT21=$6
	ROUND_OFF=$7
        if [[ $# -gt 7 ]]; then
		PREFIX=$8
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

echo " LRDMC with dla" 
echo " dir=test_lrdmc_dla" 
$PREFIX $TURBORVB < datasfn.d > $OUT
[ $? -eq 0 ] && echo " Run with non-zero exit code" || exit 1

echo "  Calculate local energies:'0 1 1 1' | readf.x" 
echo "0 1 1 1" | $READF >& /dev/null
[ $? -eq 0 ] && echo " Execution of readf with non-zero exit code" || exit 1

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

