#!/bin/bash 

TURBORVB="../../bin/turborvb-serial.x"
READF="../../bin/readf.x"
FORT21=REFERENCE_fortXXI
OUT=out_true.o

realpath() {
  case "$1" in /*) ;; *) printf '%s/' "$PWD";; esac; echo "$1"
}

TURBORVB_abs=`realpath $TURBORVB`
READF_abs=`realpath $READF`

debug_root=`pwd`

vmc_dir_list=`ls $debug_root/results`

for vmc_dir_prefix in $vmc_dir_list
do
	vmc_dir=$debug_root/results/$vmc_dir_prefix
	cd $vmc_dir
	# mkdir debug dir
	ref_dir=$vmc_dir/vmc-HF-workflow
	debug_dir=$vmc_dir/vmc-HF-debug
	mkdir -p $debug_dir
	
	# copy fort.10 and vmc.input
	cp $ref_dir/fort.10 $debug_dir/fort.10
	cp $ref_dir/pseudo.dat $debug_dir/pseudo.dat
	cp $ref_dir/datasvmc_0.input $debug_dir/datasvmc.d
	
	# replace grid and Lbox
	new_ngen=10
	grep -E ^\\s\+ngen $debug_dir/datasvmc.d | awk -F "=" '{print $2}' > tmp; read ngen < tmp; rm tmp
	grep -E -n ^\\s\+ngen $debug_dir/datasvmc.d  | sed -e 's/:.*//g' > tmp; read lngen < tmp; rm tmp
	sed -i -e "${lngen} s/${ngen}/${new_ngen}/g" $debug_dir/datasvmc.d
	
	# run VMC
	cd $debug_dir
	$TURBORVB_abs < datasvmc.d > $OUT
	echo "0 1 1 1" | $READF_abs >& /dev/null 
	cp fort.21 ${FORT21}

	cd $debug_root
done

