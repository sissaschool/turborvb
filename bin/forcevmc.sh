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

# Sandro Sorella created on 14th Dec. 2007.
# Kosuke Nakano modified on 20 Feb. 2019 to speed up.
# It is possible to avoid from repeatedly reading fort.21.
# when you put a negative number for Initial bin for averages $2
# The flag_read_once is set "True" in this case.


start_time=`date +%s`

BASEDIR=$(dirname "$0")

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

rm -f c1val

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



if [ $# -gt 3 ]; then
    echo " Averages refers to K value # $4 "
fi

if [ $2 -ge 0 ]; then
    flag_read_once=False
else
    flag_read_once=True
    num_fort12_columns=`python3 ${BASEDIR}/check_fort12_columns_length.py`
    num_all_columns=`expr ${num_fort12_columns} - 3`  # wbuf(j), wtot(j), and eskip(1). see readf.x

    # read all values here
    echo "  reading all values in fort.12...."
    echo " 0 $1 $2 ${num_all_columns} " | $BASEDIR/readf.x $4
    mv fort.21 fort.21_master

    # make each fort.21_col
    mkdir fort.21_tmp
    python3 ${BASEDIR}/read_columns_fort21.py fort.21_master ./fort.21_tmp/fort.21_col
fi

if [ $ntot -gt 0 ]; then
    rm -f forces_vmc.dat
    touch forces_vmc.dat
fi

# estimate energy, variance and error bars with boostrap resampling technique
n=$( echo "1+$ipc+$nfor+3*$nvar" | bc -l)

if [ ${flag_read_once} = True ]; then
    # just copy the corresponding column
    cp ./fort.21_tmp/fort.21_col_$n ./fort.22
    cp ./fort.21_tmp/fort.21_col_1 ./fort.21
    cp ./fort.21 fort.21.1

else
    echo " 0 $1 $2 $n " | $BASEDIR/readf.x $4
    #average_variance $nkps $1 $2 $n
    cp fort.21 fort.22
    echo " 0 $1 $2 1 " | $BASEDIR/readf.x $4
    cp fort.21 fort.21.1
fi

echo "$1" | $BASEDIR/corrvar.x > pip0.d

# energy imaginary part
# average $nkps $1 $2 1 "fort.21.1"

if [ "$ipc" -gt 1 ]; then
    # average $nkps $1 $2 2 "fort.21.1i"
    if [ ${flag_read_once} = True ]; then
        # just copy the corresponding column
        cp ./fort.21_tmp/fort.21_col_2 ./fort.21
        cp ./fort.21 fort.21.1i
    else
        echo " 0 $1 $2 2 " | $BASEDIR/readf.x $4
        cp fort.21 fort.21.1i
    fi
fi

if [ "$nvar" -gt 0 ]; then
    echo " Scalepulay ${scale_pulay} "

    for i in `seq 1 $nvar`
    do
        n=$( echo "$ipc+2 + 3*($i-1)+$nfor" | bc -l)

        if [ ${flag_read_once} = True ]; then
            # just copy the corresponding one
            cp ./fort.21_tmp/fort.21_col_$n ./fort.21.2

        else
            echo "0 $1 $2 $n" | $BASEDIR/readf.x $4
            cp fort.21 fort.21.2
        fi


        n=$( echo "$ipc+1 + 3*($i-1)+$nfor" | bc -l)
        #average $nkps $1 $2 $n "fort.21.3"
        if [ ${flag_read_once} = True ]; then
            # just copy the corresponding one
            cp ./fort.21_tmp/fort.21_col_$n ./fort.21.3

        else
            echo "0 $1 $2 $n" | $BASEDIR/readf.x $4
            cp fort.21 fort.21.3
        fi

        n=$( echo "$ipc+3+ 3*($i-1)+$nfor" | bc -l)
        #average $nkps $1 $2 $n "fort.21.4"
        if [ ${flag_read_once} = True ]; then
            # just copy the corresponding one
            cp ./fort.21_tmp/fort.21_col_$n ./fort.21
            cp fort.21 fort.21.4
        else
            echo  "0 $1 $2 $n" | $BASEDIR/readf.x $4
            cp fort.21 fort.21.4
        fi

        echo " Force component $i " >> forces_vmc.dat

        echo "${scale_pulay} " | $BASEDIR/corrforza.x  >> forces_vmc.dat

    done

else
    nvar=0

fi  # [ "$nvar" -gt 0 ]; then


if [ "$nfor" -gt 0 ]; then

    if [ "$ipcr" -gt 1 ]; then
        ndone=$ndimj
    else
        ndone=$nfor
    fi  #  [ "$ipcr" -gt 1 ]; then


    if [ "$ndone" -gt 0 ]; then

        for i in `seq 1 $ndone`
        do
            n=$( echo "$ipc + $i" | bc -l)
            #average $nkps $1 $2 $n "fort.21.2"
            if [ ${flag_read_once} = True ]; then
                # just copy the corresponding one
                cp ./fort.21_tmp/fort.21_col_$n ./fort.21.2

            else
                echo "0 $1 $2 $n" | $BASEDIR/readf.x $4
                cp fort.21 fort.21.2
            fi

            n=$( echo "$ipc+1+ $i+3*$nvar+$nfor" | bc -l)
            #average $nkps $1 $2 $n "fort.21.3"
            if [ ${flag_read_once} = True ]; then
                # just copy the corresponding one
                cp ./fort.21_tmp/fort.21_col_$n ./fort.21
                cp fort.21 fort.21.3
            else
                echo "0 $1 $2 $n" | $BASEDIR/readf.x $4
                cp fort.21 fort.21.3
            fi

            echo "Parameter number $i" >> forces_vmc.dat
            $BASEDIR/corrforzap.x  >> forces_vmc.dat
        done

    fi #  [ "$ndone" -gt 0 ]; then

    if [ "$ipcr" -gt 1 ]; then

        ndone=$( echo "($nfor-$ndimj)/2" | bc )
        ninit=$( echo "$ndimj+$ipc" | bc )


        if [ "$ndone" -gt 0 ]; then

            for i in `seq 1 $ndone`
            do

                n=$( echo "$ninit + 2*$i-1" | bc -l)
                #average $nkps $1 $2 $n "fort.21.2"
                if [ ${flag_read_once} = True ]; then
                    # just copy the corresponding one
                    cp ./fort.21_tmp/fort.21_col_$n ./fort.21.2
                else
                    echo "0 $1 $2 $n" | $BASEDIR/readf.x $4
                    cp fort.21 fort.21.2
                fi

                n=$( echo "$ninit+2*$i+3*$nvar+$nfor" | bc -l)
                #average $nkps $1 $2 $n "fort.21.3"
                if [ ${flag_read_once} = True ]; then
                    # just copy the corresponding one
                    cp ./fort.21_tmp/fort.21_col_$n ./fort.21.3
                else
                    echo "0 $1 $2 $n" | $BASEDIR/readf.x $4
                    cp fort.21 fort.21.3
                fi

                n=$( echo "$ninit + 2*$i" | bc -l)
                #average $nkps $1 $2 $n "fort.21.4"
                if [ ${flag_read_once} = True ]; then
                    # just copy the corresponding one
                    cp ./fort.21_tmp/fort.21_col_$n ./fort.21.4
                else
                    echo "0 $1 $2 $n" | $BASEDIR/readf.x $4
                    cp fort.21 fort.21.4
                fi

                n=$( echo "$ninit+1+2*$i+3*$nvar+$nfor" | bc -l)
                #average $nkps $1 $2 $n "fort.21.5"
                if [ ${flag_read_once} = True ]; then
                    # just copy the corresponding one
                    cp ./fort.21_tmp/fort.21_col_$n ./fort.21
                    cp fort.21 fort.21.5
                else
                    echo "0 $1 $2 $n" | $BASEDIR/readf.x $4
                    cp fort.21 fort.21.5
                fi

                n=$( echo "2*$i-1+$ninit-$ipc" | bc -l)
                np=$( echo "2*$i+$ninit-$ipc" | bc -l)

                echo " Complex Parameter number Real/Imag:  $n $np" >> forces_vmc.dat
                $BASEDIR/corrforzap_complex.x  >> forces_vmc.dat

            done

        fi   # [ "$ndone" -gt 0 ]; then

    fi   #  [ "$ipcr" -gt 1 ]; then

else
    nfor=0

fi    # [ "$nfor" -gt 0 ]; then

ls fort.21_tmp >/dev/null 2>&1
if [ $? -eq 0 ]; then
	\rm -r fort.21_tmp
fi

end_time=`date +%s`
time=$((end_time - start_time))
echo "-------------------------------------------------"
echo "The elapse time of forcevmc.sh is $time [sec]"
echo "-------------------------------------------------"
exit
