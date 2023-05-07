#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	TURBORVB=$1
	READF=$2
	FORCEVMC=$3
	IN=$4
	OUT=$5
	OUT_FORCEVMC=$6
	TRUEOUT=$7
	TRUEOUT_FORCEMVC=$8
	REF_FORT21=$9
	REF_FORCEVMC=${10}
	ROUND_OFF=${11}
        if [[ $# -gt 11 ]]; then
		PREFIX=${12}
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

#cd vmc
echo " VMC and VMC forces TEST " 
echo " dir=test_vmc" 
$PREFIX $TURBORVB < $IN > $OUT
[ $? -eq 0 ] && echo " Run without non-zero exit code" || exit 1

#check energies
echo "  Calculate local energies:'0 1 1 1' | readf.x" 
echo "0 1 1 1" | $READF >& /dev/null

if [ $(grep -c ERR $OUT) -gt 0 ]; then
  echo "    Errors in output:"
  echo `grep ERR $OUT`
  exit 1
fi

echo "    Rounds off values in fort.21 < 10**-${ROUND_OFF}".
cat fort.21 | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f  %.*f  %.*f\n", ROUND_OFF, $1, ROUND_OFF, $2, ROUND_OFF, $3)}' > fort.21_roundoff
cat ${REF_FORT21} | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f  %.*f  %.*f\n", ROUND_OFF, $1, ROUND_OFF, $2, ROUND_OFF, $3)}' > REFERENCE_fortXXI_roundoff

#diff fort.21
echo "    Compares fort.21_roundoff and REFERENCE_fortXXI_roundoff."
echo "    If you do not see any "diff" here, they are consistent."

if [ $(diff fort.21_roundoff REFERENCE_fortXXI_roundoff | wc -l) -gt 0 ]; then
   diff fort.21_roundoff REFERENCE_fortXXI_roundoff
   exit 1
elif [ $FORCEVMC == NA ]; then
   exit 0
fi

#calc forces
echo "  Calculate forces:forcevmc.sh 1 1 1 " 
$FORCEVMC 1 1 1 > $OUT_FORCEVMC
[ $? -eq 0 ] && echo " Run without non-zero exit code" || exit 1

#check force
echo "    Rounds off values in forces_vmc.dat < 10**-${ROUND_OFF}".
cat forces_vmc.dat | tr ' ' '\n' |  sed  '/^$/d' > tmp1
cat ${REF_FORCEVMC} | tr ' ' '\n' |  sed  '/^$/d' > tmp2
#sdiff -s tmp1 tmp2 | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f\n", ROUND_OFF, $1)}' > forces_vmc_roundoff
#sdiff -s tmp1 tmp2 | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f\n", ROUND_OFF, $3)}' > forces_REFERENCE_roundoff
diff -y tmp1 tmp2 | echo `grep \|` | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f\n", ROUND_OFF, $1)}' > forces_vmc_roundoff
diff -y tmp1 tmp2 | echo `grep \|` | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f\n", ROUND_OFF, $3)}' > forces_REFERENCE_roundoff

#diff forces
echo "    Compares forces_vmc_roundoff and forces_REFERENCE_roundoff."
echo "    If you do not see any "diff" here, they are consistent."
if [ $(diff forces_vmc_roundoff forces_REFERENCE_roundoff | wc -l) -gt 0 ]; then
  diff forces_vmc_roundoff forces_REFERENCE_roundoff | wc -l
  exit 1
else
  exit 0
fi

