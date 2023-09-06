#!/bin/bash

# Copyright (C) 2022 TurboRVB group
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# Kosuke Nakano created on 24 Jan. 2019
# Kosuke Nakano mmodified on 20 Feb. 2019 to avoid reading
# fort.21 repeatedly. When you put a negative number for
# Initial bin for averages $2, the new function works.
# The flag_read_once in forcevmc.sh is set "True" in this case.
#
# updated by Kosuke Nakano on 6 Nov. 2020 to remove pulay_ratio
# option. Indeed, one can omit pulay_ratio input, in this case,
# pulay_ratio is set 1 automatically.
#
# how to use
# forcevmc_kpoints_parallel.sh  bins  initials  pulay_ratio  number_of_processors
# or
# forcevmc_kpoints_parallel.sh  bins  initials  number_of_processors


start_time=`date +%s`

BASEDIR=$(dirname "$0")

# setting for parallel calculations
if [ ! -z $4 ] ; then
    scale_pulay=$3
    num_cpu_max=`expr $4`
else
    scale_pulay=1
    num_cpu_max=`expr $3`
fi

root_dir=`pwd`

echo "max cpu num = ${num_cpu_max}"

######################## MAIN SCRIPT ########################

awk ' NR == 1 {print $7}' parminimized.d  > c1val
read nvar < c1val
awk ' NR == 1 {print $8}' parminimized.d  > c1val
read nfor < c1val
rm -f c1val
awk ' NR == 1 {print $11}' parminimized.d  > c1val
read ipcr < c1val
awk ' NR == 1 {print $1}' parminimized.d  > c1val
read iesinv < c1val
awk ' NR == 1 {print $2}' parminimized.d  > c1val
read iesm < c1val
awk ' NR == 1 {print $3}' parminimized.d  > c1val
read iesd < c1val
awk ' NR == 1 {print $4}' parminimized.d  > c1val
read iesfree < c1val
awk ' NR == 1 {print $5}' parminimized.d  > c1val
read iessw < c1val

# check if k-points are present
if [ -f "kp_info.dat" ]  ; then
	awk ' NR == 1 {print $1}' kp_info.dat  > c1val
	read nkps < c1val
else
	nkps=1
    echo " 1 " > kp_info.dat
    echo "#" >> kp_info.dat
    echo " 1 0.0 0.0 0.0 1 " >> kp_info.dat
fi

ndimj=$( echo "$iesfree+$iesd+$iesm+$iesinv" | bc -l)
ndims=$( echo "$iesfree+$iesd+$iesm+$iesinv+$iessw" | bc -l)
ntot=$( echo "$iesfree+$iesd+$iesm+$iesinv+$iessw+$nvar" | bc -l)

#define ipc
if [ "$ipcr" -gt 0 ] ;then
	ipc=$( echo "$ipcr" | bc -l)
else
	ipc=$( echo "(-1)*$ipcr" | bc -l)
fi

echo "number of k-points: $nkps"
echo "number of ionic forces $nvar"
echo "number of VMC parameters  $nfor"
echo "number of Jastrow VMC parameters  $ndimj"
echo "Bin length   $1"
echo "Initial bin for averages   $2"

if [ $ntot -gt 0 ] ; then
	rm -f forces_vmc.dat
fi

touch forces_vmc.dat

# clean the previous trials
ls tmp_dir_K* >/dev/null 2>&1
if [ $? -eq 0 ]; then
	rm -r tmp_dir_K*
fi

ls fort12_tmp >/dev/null 2>&1
if [ $? -eq 0 ]; then
	rm -r fort12_tmp
fi

# loop conter for kpoints
counter_for_parallel_calc=1

# make temp directory
mkdir fort12_tmp
python ${BASEDIR}/split_fort.12.py ./fort.12 ./fort.12 ./kp_info.dat
mv fort.12_K* fort12_tmp/

for i in `seq 1 $nkps`

do

	mkdir tmp_dir_K${i}
	cd tmp_dir_K${i}

	# create copy
    cp ../c1val ./c1val
	cp ../parminimized.d ./parminimized.d
    echo "     1" > ./kp_info.dat
    sed -n '2p' ../kp_info.dat >> ./kp_info.dat
    sed -n ''$((2+$i))'p' ../kp_info.dat >> ./kp_info.dat
    sed -n ''$((2+$nkps+1))'p' ../kp_info.dat >> ./kp_info.dat
    sed -n ''$((2+$nkps+1+$i))'p' ../kp_info.dat >> ./kp_info.dat
    mv ../fort12_tmp/fort.12_K${i} ./fort.12

	touch forces_vmc.dat

	# calc forces
    ${BASEDIR}/forcevmc.sh $1 $2 ${scale_pulay} 1 1>out_K${i}.o 2>err_K${i}.o & # run background
	echo "  calculating for K${i}..., pid = $!, counter = ${counter_for_parallel_calc}"
        counter_for_parallel_calc=`expr ${counter_for_parallel_calc} + 1`

	# waiting
	if [ ${counter_for_parallel_calc} -gt ${num_cpu_max} ] ; then
		counter_for_parallel_calc=1
		echo -n "  waiting for the calculations..."
		wait
		echo "done."
	fi

	cd ../

done

# wait for the all calculations
echo -n "  waiting for all the calculations..."
wait
echo "done."

# post procedures

echo -n "  collecting the results..."

for i in `seq 1 $nkps`

do
	cd tmp_dir_K${i}

        mv pip0.d ../pip0_K$i

        if [ $ntot -gt 0 ] ; then
			mv forces_vmc.dat ../forces_K$i
        fi

        cd ../

done

echo "done."

python ${BASEDIR}/energy.py  > pip0.d

if [ $ntot -gt 0 ]; then
	python ${BASEDIR}/force.py  > forces_vmc.dat
fi

cd $root_dir
rm -r fort12_tmp
rm -r tmp_dir_K*

end_time=`date +%s`
time=$((end_time - start_time))
echo "--------------------------------------------------------------------------"
echo "The elapse time of forcevmc_kpoints_parallel.sh is $time [sec] by $num_cpu_max cores"
echo "--------------------------------------------------------------------------"
exit
