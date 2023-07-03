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

# Andrea Zen created on 14th Jan. 2016.

# Script that reads one output file of a TurboRVB optimization and
# plot the ratiodyn per MD step. Use in the following way:
#    plot_ratiodyn.sh output
# where output is the output file.
# This script produces some auxiliary files:
#    _plot.gpl , _XXplotDataXX*
# that can be removed after execution.

### To use in CINECA, uncomment the folowing line
#module load grep

if [ $# -gt "0" ]; then
	INPUT=$( echo $@ )
	nout=$( echo $INPUT | wc -w )
echo "n_args: " $nout
echo "args:" $INPUT
echo
echo "output files: "
	for (( i=1, j=0 ; i<=$nout ; i++ )); do
		out=$( echo $INPUT | awk -v campo=$i '{ print $(campo) }' );
		if [ -e $out ]; then
			(( j++ ))
			outfile[$j]=$out
			echo "outfile[ $j ] = ${outfile[$j]} --- nsteps MD: ` cat ${outfile[$j]} | grep "Ratio dyn =" | wc -l ` "
		fi
	done
	nfiles=$j
        BIN="100"
        if [ ! -e $out ]; then
                if [ $out -gt 0 -a $out -le 1000 ]; then
                        BIN=$out
                fi
                if [ $out -lt 0 ]; then
                        PLOT="0"
                fi
        fi


	if [ ` cat ${outfile[@]} | grep "Ratio dyn =" | wc -l ` -gt "0" ]; then
		cat ${outfile[@]} | grep "Ratio dyn =" | awk '{ if (NF==5) {
			print NR, $4, $5; }
		else {
			print NR, $4, $4;} }'
	fi > _XXplotDataXX

	echo "set title \" plot of $outfile \" " > _plot.gpl

YMAX=$( cat _XXplotDataXX | howmuch.sh 2 | awk '{ print 3*$10 }' )
echo $YMAX

cat << EOF >> _plot.gpl
set xlabel "MD step"
set ylabel " Ratio dyn "
set xrange [*:*]
set yrange [*: ${YMAX} ]
plot "_XXplotDataXX" u 1:2 noti  w l
replot "< running_averages.sh  _XXplotDataXX  $BIN " u 1:2:3 ti "running av. bin=$BIN" w ye
pause -1

EOF

###	echo
###	cat _XXplotDataXX | awk '{ cont++; ratiodymav+=$2 }END{ print "Average Ratio dyn on ",cont," MD steps is ",ratiodymav/cont;}'

	echo
		cat ${outfile[@]} | grep "Ratio dyn =" | howmuch.sh 4 | awk '{ print "Average Ratio dyn = ", $10, "+/-", $12; }'
	echo

if [ $PLOT ]; then
        echo "no plot"
else
        gnuplot _plot.gpl
fi

# rimuovo file temporanei creati
rm -f _XXplotDataXX* _plot.gpl

else
	echo
	echo "Script that reads one output file of a TurboRVB optimization and" 
	echo "plot the ratiodyn per MD step. Use in the following way:"
	echo "   plot_ratiodyn.sh output"
	echo "where output is the output file."
	echo "This script produces some auxiliary files:"
	echo "   _plot.gpl , _XXplotDataXX*  "
	echo "that can be removed after execution. "
	echo
	echo "To use in CINECA, edit the file and uncomment the line:"
	echo "#module load gnuplot"
	echo
#	echo "ERROR: file < $outfile > DO NOT EXIST!"
	echo
fi
