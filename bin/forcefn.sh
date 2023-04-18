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

# Sandro Sorella created on 14th Dec. 2007.
# Kosuke Nakano modified on 26 Mar. 2019 to speed up.
# It is possible to avoid from repeatedly reading fort.21.
# when you put a negative number for Initial bin for averages $3
# The flag_read_once is set "True" in this case.


start_time=`date +%s`

BASEDIR=$(dirname "$0")

# read number of parameters from "parminimized.d"
awk ' NR == 1 {print $7}' parminimized.d  > c1val
read nvar < c1val
awk ' NR == 1 {print $8}' parminimized.d  > c1val
read nfor < c1val
awk ' NR == 1 {print $10}' parminimized.d  > c1val
read niese < c1val
awk ' NR == 1 {print $9}' parminimized.d  > c1val
read nisfix < c1val
awk ' NR == 1 {print $11}' parminimized.d  > c1val
read ipcr < c1val

# check if k-points are present
if [ -f "kp_info.dat" ]; then
    awk ' NR == 1 {print $1}' kp_info.dat  > c1val
    read nkps < c1val
else
    nkps=1
fi

rm -f c1val

#define ipc
if [ "$ipcr" -gt 0 ]; then
    ipc=$( echo "$ipcr" | bc -l)
else
    ipc=$( echo "(-1)*$ipcr" | bc -l)
fi

ntot=$( echo "$nfor+$nvar" | bc -l)


echo "number of k-points: $nkps"
echo "number of ionic forces $nvar"
echo "number of VMC parameters  $nfor"

echo "Bin length   $1"
echo "Correcting factors   $2"
echo "Initial bin for averages   $3"

if [ ! -z $4 ]
then

scale_pulay=$4

else

scale_pulay=1

fi


echo " ScalePulay    ${scale_pulay} "

if [ $# -gt 4 ]; then
    echo " Averages refers to K value # $5 "
fi

if [ $3 -ge 0 ]; then
    flag_read_once=False # old version
else
    flag_read_once=True
    num_fort12_columns=`python ${BASEDIR}/check_fort12_columns_length.py`
    num_all_columns=`expr ${num_fort12_columns} - 4`  #  skip wbuf(j), pip, wtot(j), eskip(1), see readffn.x
    echo "num_all_columns=$num_all_columns"

    # read all values here
    echo "  reading all values in fort.12...."
    echo " $2 $1 $3 ${num_all_columns} " | readffn.x $5
    mv fort.21 fort.21_master

    # make each fort.21_col
    mkdir fort.21_tmp
    python ${BASEDIR}/read_columns_fort21.py fort.21_master ./fort.21_tmp/fort.21_col
fi

if [ $ntot -gt 0 ]; then
    rm -f forces_fn.dat
    touch forces_fn.dat
fi

if [ ${flag_read_once} = True ]; then
    # just copy the corresponding column
    cp ./fort.21_tmp/fort.21_col_1 ./fort.21
    mv ./fort.21 fort.21.1
else
    echo " $2 $1 $3 1 " | readffn.x $5
    mv fort.21 fort.21.1
fi

if [ "$ipc" -gt 1 ]; then
    if [ ${flag_read_once} = True ]; then
        # just copy the corresponding column
        cp ./fort.21_tmp/fort.21_col_2 ./fort.21
        mv ./fort.21 fort.21.1i
    else
        echo " $2 $1 $3 2 " | readffn.x $5
        mv fort.21 fort.21.1i
    fi
fi

if [ "$nvar" -gt 0 ]; then
for i in `seq 1 $nvar`
do

    n=$( echo "2 + 3*($i-1)+$nfor+$niese" | bc -l)

    if [ ${flag_read_once} = True ]; then
        # just copy the corresponding one
        cp ./fort.21_tmp/fort.21_col_$n ./fort.21
        mv fort.21 fort.21.2
    else
        echo " $2 $1 $3 $n " | readffn.x $5
        mv fort.21 fort.21.2
    fi

    n=$( echo "1 + 3*($i-1)+$nfor+$niese" | bc -l)

    if [ ${flag_read_once} = True ]; then
        # just copy the corresponding one
        cp ./fort.21_tmp/fort.21_col_$n ./fort.21
        mv fort.21 fort.21.3
    else
        echo " $2 $1 $3 $n " | readffn.x $5
        mv fort.21 fort.21.3
    fi

    n=$( echo "3 + 3*($i-1)+$nfor+$niese" | bc -l)

    if [ ${flag_read_once} = True ]; then
        # just copy the corresponding one
        cp ./fort.21_tmp/fort.21_col_$n ./fort.21
        mv fort.21 fort.21.4
    else
        echo " $2 $1 $3 $n " | readffn.x $5
        mv fort.21 fort.21.4
    fi

    echo " Force component $i " >> forces_fn.dat
    echo "${scale_pulay}" | corrforza.x >> forces_fn.dat

done
fi

if [ "$nfor" -gt 0 ]; then
echo " Scalepulay ${scale_pulay} "

for i in `seq 1 $nfor`
do

    n=$( echo "$niese + $i" | bc -l)

    if [ ${flag_read_once} = True ]; then
        # just copy the corresponding one
        cp ./fort.21_tmp/fort.21_col_$n ./fort.21
        mv fort.21 fort.21.2
    else
        echo " $2 $1 $3 $n " | readffn.x $5
        mv fort.21 fort.21.2
    fi

    n=$( echo "1 + $i+3*$nvar+$nfor+$niese+$nisfix-1" | bc -l)

    if [ ${flag_read_once} = True ]; then
        # just copy the corresponding one
        cp ./fort.21_tmp/fort.21_col_$n ./fort.21
        mv fort.21 fort.21.3
    else
        echo " $2 $1 $3 $n " | readffn.x $5
        mv fort.21 fort.21.3
    fi

    echo "Parameter number $i" >> forces_fn.dat
    corrforzap.x >> forces_fn.dat

done
fi

n=$( echo "1+$niese+$nfor+3*$nvar" | bc -l)

if [ ${flag_read_once} = True ]; then
    # just copy the corresponding one
    cp ./fort.21_tmp/fort.21_col_$n ./fort.21
    cp fort.21 fort.22
    cp fort.21.1 fort.21
else
    echo " $2 $1 $3 $n " | readffn.x $5
    cp fort.21 fort.22
    cp fort.21.1 fort.21
fi

echo "$1" | corrvar.x > pip0_fn.d
    echo " $2 $1 $3 1 " | readf.x $5

    tail -1 fort.20 >  c1val
    awk ' NR == 1 {print "  Energy (ave) =",$1,$2}' c1val >> pip0_fn.d
    rm -f c1val
#echo "$1" | corrvarvar.x >> pip0_fn.d

# necessary?
#if [ ${flag_read_once} = True ]; then
    :
#else
#fi

ls fort.21_tmp >/dev/null 2>&1
if [ $? -eq 0 ]; then
	\rm -r fort.21_tmp
fi

end_time=`date +%s`
time=$((end_time - start_time))
echo "-------------------------------------------------"
echo "The elapse time of forcefn.sh is $time [sec]"
echo "-------------------------------------------------"

exit
