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

# Andrea Zen created on 18th Oct. 2011.

# Script that reads one output file of a TurboRVB optimization and
# plot the devmax values per step. Use in the following way:
#    plot_devmax.sh output
# where output is the output file.
# This script produces some auxiliary files:
#    _plot.gpl , _XXplotDataXX*
# that can be removed after execution.

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
			echo "outfile[ $j ] = ${outfile[$j]} --- nsteps: ` cat ${outfile[$j]} | grep "devmax par Normal" | wc -l ` "
		fi
	done
	nfiles=$j

	cat ${outfile[@]} | grep "devmax par Normal"   | awk '{ print NR, $5; }' > _XX1
	cat ${outfile[@]} | grep "devmax SR step"      | awk '{ print     $5; }' > _XX2
	cat ${outfile[@]} | grep "devmax force"        | awk '{ print     $4; }' > _XX3
	cat ${outfile[@]} | grep "devmax par ions"     | awk '{ print NR, $5; }' > _XXplotDataIonXX

	paste _XX1 _XX2 _XX3 > _XXplotDataXX
	rm -f _XX1 _XX2 _XX3

	NparI=$(wc -l _XXplotDataIonXX | awk '{print $1}')
	if [ "$NparI" -gt "0" ]; then
		echo "Structural optimization!"
	else
		rm _XXplotDataIonXX
	fi

	echo "set title \" plot of $outfile \" " > _plot.gpl

cat << EOF >> _plot.gpl
set xlabel "optimization step"
set ylabel "devmax"
set xrange [*:*]
set yrange [*:*]
plot "_XXplotDataXX" u 1:2 ti "devmax  par Normal" w l , "" u 1:3 ti "devmax SR step" w l , "" u 1:4 ti "devmax force" w l
repl 4 ti "threshold <4" w l
EOF

	if [ -e _XXplotDataIonXX ]; then
cat << EOF >> _plot.gpl
repl "_XXplotDataIonXX" u 1:2 ti "devmax par ions" w l
EOF
	fi

cat << EOF >> _plot.gpl
pause -1
EOF


	gnuplot _plot.gpl


else
	echo
	echo "Script that reads one output file of a TurboRVB optimization and"
	echo "plot the devmax values per step. Use in the following way:"
	echo "   plot_devmax.sh output"
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
