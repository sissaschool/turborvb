#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	PREP=$1
	OUT=$2
	TRUEOUT=$3
	FORT10=$4
	ROUND_OFF=$5
        if [[ $# -gt 5 ]]; then
		PREFIX=$6
	else
		PREFIX=""
	fi
else
	source ../settings.sh
fi

if [ ! -f "$PREP" ]; then
    echo "Executable $PREP does not exists"
    exit 1
fi

echo " DFT-PBC-TWIST-COMPLEX TEST " 
echo " dir=test_dft_pbc_twist_complex" 
$PREFIX $PREP < prep.d > $OUT
[ $? -eq 0 ] && echo " Run without non-zero exit code" || exit 1
DFT_ene_ref=`grep "Final variational DFT" ${TRUEOUT} | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f\n", ROUND_OFF, $7)}'`
DFT_ene=`grep "Final variational DFT" $OUT | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f\n", ROUND_OFF, $7)}'`
dft_ene_diff=`echo "scale=${ROUND_OFF}; ${DFT_ene_ref} - ${DFT_ene}" | bc -l`
echo "   -DFT energy = ${DFT_ene}"
echo "   -DFT energy (ref) = ${DFT_ene_ref}"
echo "   -The diff = ${dft_ene_diff}"
echo "   -If the diff is finite, there is something wrong."
echo `grep ERR $OUT`

[ -z "${dft_ene_diff}" ] && exit 1

#re="^[+-][0-9]*.[0-9]*$"
#if ! [[ ${dft_ene_diff} =~ $re ]] ; then
#  exit 1
#fi

exit_code=`echo "$dft_ene_diff>0" | bc -l`
exit ${exit_code}
