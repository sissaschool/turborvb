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

realpath() {
  case "$1" in /*) ;; *) printf '%s/' "$PWD";; esac; echo "$1"
}

PREP_abs=`realpath $PREP`

if [ ! -f "$PREP_abs" ]; then
    echo "Executable $PREP_abs does not exists"
    exit 1
fi

echo " DFT-OPEN TEST for 38 molecules " 
echo " dir=test_dft_open_38_molecules" 

exit_code_arr=()

debug_root=`pwd`

prep_dir_list=`ls $debug_root/results`

for prep_dir_prefix in $prep_dir_list
do
	prep_dir=$debug_root/results/$prep_dir_prefix
        ref_dir=$prep_dir/dft-LDA-workflow
	debug_dir=$prep_dir/dft-LDA-debug
        cp $ref_dir/fort.10 $debug_dir/fort.10
        cp $ref_dir/pseudo.dat $debug_dir/pseudo.dat

	cd $debug_dir
	label=`basename $prep_dir`
	echo " =${label}="
	$PREFIX $PREP_abs < prep.d > $OUT
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
	exit_code=`echo "$dft_ene_diff != 0" | bc -l`
	exit_code_arr+=($exit_code)
	
	cd $debug_root
done

for exit_code in ${exit_code_arr[@]}; do
	if [[ ${exit_code} = 1 ]]; then
		exit 1
	fi
done
exit 0

