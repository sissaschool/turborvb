#!/bin/bash
set -euo pipefail

# Set default values
TURBORVB=turborvb-serial.x
READF=readf.x
FORCEVMC=NA
INPUT_TURBORVB=datas.input
OUTPUT_TURBORVB=out
OUTPUT_FORCEVMC=out_forcevmc
OUTPUT_REF=out_true
OUTPUT_FORCEVMC_REF=out_forcevmc_true
FORT21_REF=fort.21_true
FORCEVMC_REF=force_vmc_true
ROUND_OFF=6
EXECUTE_PREFIX=""

if [[ $# -gt 0 ]]; then
  for ARG in "$@"; do
    if [[ $ARG =~ ^([^=]+)=([^\']*)$ ]]; then
  
      KEY="${BASH_REMATCH[1]}"
      VALUE="${BASH_REMATCH[2]}"
  
      eval "${KEY}='${VALUE}'"
    else
      echo "Uknown argument ${ARG}"
      echo "Parameters has to be in key value pairs: $0 [key1=value1] [key2=value2] ..."
      exit 1
    fi
  done
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
$EXECUTE_PREFIX $TURBORVB < $INPUT_TURBORVB > $OUTPUT_TURBORVB
[ $? -eq 0 ] && echo " Run without non-zero exit code" || exit 1

#check energies
echo "  Calculate local energies:'0 1 1 1' | readf.x" 
echo "0 1 1 1" | $READF >& /dev/null

if [ $(grep -c ERR $OUTPUT_TURBORVB) -gt 0 ]; then
  echo "    Errors in output:"
  echo `grep ERR $OUTPUT_TURBORVB`
  exit 1
fi

echo "    Rounds off values in fort.21 < 10**-${ROUND_OFF}".
cat fort.21 | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f  %.*f  %.*f\n", ROUND_OFF, $1, ROUND_OFF, $2, ROUND_OFF, $3)}' > fort.21_roundoff
cat ${FORT21_REF} | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f  %.*f  %.*f\n", ROUND_OFF, $1, ROUND_OFF, $2, ROUND_OFF, $3)}' > REFERENCE_fortXXI_roundoff

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

