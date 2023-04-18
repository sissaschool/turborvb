#!/bin/bash 

PREP=../../bin/prep-serial.x
FORT10=REFERENCE_fort.10_new
OUT=out_true.o

realpath() {
  case "$1" in /*) ;; *) printf '%s/' "$PWD";; esac; echo "$1"
}

PREP_abs=`realpath $PREP`

debug_root=`pwd`

prep_dir_list=`ls $debug_root/results`

for prep_dir_prefix in $prep_dir_list
do
	prep_dir=$debug_root/results/$prep_dir_prefix
	cd $prep_dir
	# mkdir debug dir
	ref_dir=$prep_dir/dft-LDA-workflow
	debug_dir=$prep_dir/dft-LDA-debug
	mkdir -p $debug_dir
	
	# copy fort.10 and prep.input
	cp $ref_dir/fort.10 $debug_dir/fort.10
	cp $ref_dir/pseudo.dat $debug_dir/pseudo.dat
	cp $ref_dir/prep.input $debug_dir/prep.d
	
	# replace grid and Lbox
	m_a=6; m_n=15
	grep -E ^\\s\+ax $debug_dir/prep.d | awk -F "=" '{print $2}' > tmp; read ax < tmp; rm tmp
	grep -E ^\\s\+ay $debug_dir/prep.d | awk -F "=" '{print $2}' > tmp; read ay < tmp; rm tmp
	grep -E ^\\s\+az $debug_dir/prep.d | awk -F "=" '{print $2}' > tmp; read az < tmp; rm tmp 
	grep -E ^\\s\+nx $debug_dir/prep.d | awk -F "=" '{print $2}' > tmp; read nx < tmp; rm tmp
	grep -E ^\\s\+ny $debug_dir/prep.d | awk -F "=" '{print $2}' > tmp; read ny < tmp; rm tmp
	grep -E ^\\s\+nz $debug_dir/prep.d | awk -F "=" '{print $2}' > tmp; read nz < tmp; rm tmp
	grep -E ^\\s\+epsdft $debug_dir/prep.d | awk -F "=" '{print $2}' > tmp; read epsdft < tmp; rm tmp

	new_ax=`echo "scale=0; $ax * ${m_a}" | bc -l`
	new_ay=`echo "scale=0; $ay * ${m_a}" | bc -l`
	new_az=`echo "scale=0; $az * ${m_a}" | bc -l`
	
	new_nx=`echo "scale=0; $nx / ${m_n}" | bc -l`
	new_ny=`echo "scale=0; $ny / ${m_n}" | bc -l`
	new_nz=`echo "scale=0; $nz / ${m_n}" | bc -l`
	
	grep -E -n ^\\s\+ax $debug_dir/prep.d  | sed -e 's/:.*//g' > tmp; read lax < tmp; rm tmp
	grep -E -n ^\\s\+ay $debug_dir/prep.d  | sed -e 's/:.*//g' > tmp; read lay < tmp; rm tmp
	grep -E -n ^\\s\+az $debug_dir/prep.d  | sed -e 's/:.*//g' > tmp; read laz < tmp; rm tmp
	grep -E -n ^\\s\+nx $debug_dir/prep.d  | sed -e 's/:.*//g' > tmp; read lnx < tmp; rm tmp
	grep -E -n ^\\s\+ny $debug_dir/prep.d  | sed -e 's/:.*//g' > tmp; read lny < tmp; rm tmp
	grep -E -n ^\\s\+nz $debug_dir/prep.d  | sed -e 's/:.*//g' > tmp; read lnz < tmp; rm tmp
	grep -E -n ^\\s\+epsdft $debug_dir/prep.d  | sed -e 's/:.*//g' > tmp; read lepsdft < tmp; rm tmp

	sed -i -e "${lax} s/${ax}/${new_ax}/g" $debug_dir/prep.d
	sed -i -e "${lay} s/${ay}/${new_ay}/g" $debug_dir/prep.d
	sed -i -e "${laz} s/${az}/${new_az}/g" $debug_dir/prep.d
	sed -i -e "${lnx} s/${nx}/${new_nx}/g" $debug_dir/prep.d
	sed -i -e "${lny} s/${ny}/${new_ny}/g" $debug_dir/prep.d
	sed -i -e "${lnz} s/${nz}/${new_nz}/g" $debug_dir/prep.d
	sed -i -e "${lepsdft} s/${epsdft}/0.001/g" $debug_dir/prep.d

	# run DFT
	cd $debug_dir
	$PREP_abs < prep.d > $OUT
	cp fort.10_new ${FORT10}

	cd $debug_root
done

