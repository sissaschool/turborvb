#!/bin/bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
	TURBORVB=$1
	READF=$2
	FORCEVMC=$3
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

realpath() {
  case "$1" in /*) ;; *) printf '%s/' "$PWD";; esac; echo "$1"
}

TURBORVB_abs=`realpath $TURBORVB`
READF_abs=`realpath $READF`

if [ ! -f "$TURBORVB_abs" ]; then
    echo "Executable $TURBORVB_abs does not exists"
    exit 1
fi

if [ ! -f "$READF_abs" ]; then
    echo "Executable $READF_abs does not exists"
    exit 1
fi

#cd vmc
echo " VMC TEST for 8 crystals at k=(0.50 0.50 0.50)" 
echo " dir=test_vmc_pbc_trim_8_crystals_insulator"

exit_code_arr=()

debug_root=`pwd`

vmc_dir_list=`ls $debug_root/results`

for vmc_dir_prefix in $vmc_dir_list
do
	vmc_dir=$debug_root/results/$vmc_dir_prefix
	ref_dir=$vmc_dir/vmc-HF-workflow
	debug_dir=$vmc_dir/vmc-HF-debug
        cp $ref_dir/fort.10 $debug_dir/fort.10
        cp $ref_dir/pseudo.dat $debug_dir/pseudo.dat

	cd $debug_dir
	label=`basename $vmc_dir`
	echo " =${label}="
	$PREFIX $TURBORVB_abs < datasvmc.d > $OUT
	[ $? -eq 0 ] && echo " Run without non-zero exit code" || exit 1
	
	#check energies
	echo "  Calculate local energies:'0 1 1 1' | readf.x" 
	echo "0 1 1 1" | $READF_abs >& /dev/null
	
	if [ $(grep -c ERR $OUT) -gt 0 ]; then
	echo "    Errors in output:"
	grep ERR $OUT 
	exit 1
	fi
	
	echo "    Rounds off values in fort.21 < 10**-${ROUND_OFF}".
	cat fort.21 | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f  %.*f  %.*f\n", ROUND_OFF, $1, ROUND_OFF, $2, ROUND_OFF, $3)}' > fort.21_roundoff
	cat ${REF_FORT21} | awk -v ROUND_OFF=${ROUND_OFF} '{printf("%.*f  %.*f  %.*f\n", ROUND_OFF, $1, ROUND_OFF, $2, ROUND_OFF, $3)}' > REFERENCE_fortXXI_roundoff
	
	#diff fort.21
	echo "    Compares fort.21_roundoff and REFERENCE_fortXXI_roundoff."
	echo "    If you do not see any "diff" here, they are consistent."
	echo ""
	if [ $(diff fort.21_roundoff REFERENCE_fortXXI_roundoff | wc -l) -gt 0 ]; then
		diff fort.21_roundoff REFERENCE_fortXXI_roundoff
		exit_code_arr+=(1)
	elif [ $FORCEVMC == NA ]; then
		exit_code_arr+=(0)
	else
		exit_code_arr+=(0)
	fi

done

for exit_code in ${exit_code_arr[@]}; do
	if [[ ${exit_code} = 1 ]]; then
		exit 1
	fi
done
exit 0
