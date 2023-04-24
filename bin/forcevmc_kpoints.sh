#!/bin/sh

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

# Sandro Sorella created on 29th Jun. 2016.
# Kosuke Nakano modified on 20 Feb. 2019 to avoid reading
# fort.21 repeatedly. When you put a negative number for
# Initial bin for averages $2, the new function works.
# The flag_read_once in forcevmc.sh is set "True" in this case.


start_time=`date +%s`

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

if [ -f "kp_info.dat" ]; then
    awk ' NR == 1 {print $1}' kp_info.dat  > c1val
    read nkps < c1val
else
    nkps=1
fi

ndimj=$( echo "$iesfree+$iesd+$iesm+$iesinv" | bc -l)
ndims=$( echo "$iesfree+$iesd+$iesm+$iesinv+$iessw" | bc -l)
ntot=$( echo "$iesfree+$iesd+$iesm+$iesinv+$iessw+$nvar" | bc -l)

#define ipc
if [ "$ipcr" -gt 0 ]; then
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

if [ ! -z $3 ]
then

scale_pulay=$3

else

scale_pulay=1

fi



if [ $ntot -gt 0 ]; then
    rm -f forces_vmc.dat
fi

touch forces_vmc.dat

for i in `seq 1 $nkps`
do
    forcevmc.sh $1 $2 ${scale_pulay} $i

    wkp=$(awk -v ind="${i}" ' NR == ind+2 {print $5} ' kp_info.dat)
    #echo " $i $wkp " >> prova

    mv pip0.d pip0_K$i

    if [ $ntot -gt 0 ]; then
        mv forces_vmc.dat forces_K$i
    fi
done

python $(which energy.py)  > pip0.d

if [ $ntot -gt 0 ]; then
    python $(which force.py)  > forces_vmc.dat
fi

end_time=`date +%s`
time=$((end_time - start_time))
echo "--------------------------------------------------------"
echo "The elapse time of forcevmc_kpoints.sh is $time [sec]"
echo "--------------------------------------------------------"
exit
