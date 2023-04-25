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

# Andrea Zen created on 1st Jul. 2015.

# Script that reads the output files of a TurboRVB optimization and
# plot the energy per step. Use in the following way:
#    multiplot_Energy.sh out_1 out_2 ... out_n bin
# where out_? are the output files, and bin is the number of
# bins for the running average (it is optional; if it is not
# defined it is 10 by default).
# This script produces some auxiliary files:
#    _plot.gpl , _XXplotDataXX*
# that can be removed after execution.

if [ $# -gt "0" ]; then
	INPUT=$( echo $@ )
	nout=$( echo $INPUT | wc -w )
#echo $nout
#echo $INPUT
	for (( i=1, j=0 ; i<=$nout ; i++ )); do
		out=$( echo $INPUT | awk -v campo=$i '{ print $(campo) }' );
		if [ -e $out ]; then
			(( j++ ))
			outfile[$j]=$out
			echo "outfile[ $j ] = ${outfile[$j]} --- nsteps: ` cat ${outfile[$j]} | grep "New Energy =" | wc -l ` "
		fi
	done
	nfiles=$j
	BIN="10"
	if [ ! -e $out ]; then
		if [ $out -gt 0 -a $out -le 500 ]; then
			BIN=$out
		fi
		if [ $out -lt 0 ]; then
			PLOT="0"
		fi
	fi

#echo ${outfile[@]}
#echo $BIN


	cat ${outfile[@]} | grep "New Energy =" | awk '{ print NR, $4, $5; }' > _XXplotDataXX

	out1=$( cat _XXplotDataXX | wc -l )
	INIZ=$( head -n 1 _XXplotDataXX | awk '{ printf("%.4f(%.4f)\n", $2, $3); }' )
	FINE=$( tail -n 1 _XXplotDataXX | awk '{ printf("%.4f(%.4f)\n", $2, $3); }' )
	DIFFOPT=$( echo " $INIZ $FINE " | sed 's/(/ /g' | sed 's/)/ /g' | awk '{ printf("%.4f(%.4f)\n",$3-$1,sqrt($2^2+$4^2)); }' )
	out2=$( cat _XXplotDataXX | head -n ${BIN} | howmuch.sh 2 | awk '{ printf("%.4f(%.4f)\n", $10, $12); }' )
	out3=$( cat _XXplotDataXX | tail -n ${BIN} | howmuch.sh 2 | awk '{ printf("%.4f(%.4f)\n", $10, $12); }' )
	out4=$( cat _XXplotDataXX | tail -n ${BIN} | howmuch.sh 3 | awk '{ print $10 }' )
	echo "* ${outfile[@]}: $out1 | $INIZ I(${BIN}):$out2 | $FINE F(${BIN}):$out3 | F-I: $DIFFOPT "
#		echo

echo "set title \" plot of ${outfile[@]} \" " > _plot.gpl

cat << EOF >> _plot.gpl
set xlabel "optimization step"
set ylabel "New Energy"
set xrange [*:*]
set yrange [*:*]
plot "< running_averages.sh  _XXplotDataXX  $BIN " u 1:2:3 ti "running av. bin=$BIN" w ye
replot "_XXplotDataXX" u 1:2:3 ti "${outfile[@]}"  w yerr
pause -1
EOF

if [ $PLOT ]; then
	echo "no plot"
else
	gnuplot _plot.gpl
fi

# rimuovo file temporanei creati
rm -f _XXplotDataXX* _plot.gpl

else
	echo
	echo "Script that reads the output files of a TurboRVB optimization and"
	echo "plot the energy per step. Use in the following way:"
	echo "   plot_Energy.sh out_1 out_2 ... out_n bin   "
	echo "where out_? are the output files, and bin is the number of "
	echo "bins for the running average (it is optional; if it is not "
	echo "defined it is 10 by default). "
	echo
	echo "This script produces some auxiliary files: "
	echo "   _plot.gpl , _XXplotDataXX* "
	echo "that can be removed after execution. "
	echo
	echo "To use in CINECA, edit the file and uncomment the line:"
	echo "#module load gnuplot"
	echo
	echo
fi

