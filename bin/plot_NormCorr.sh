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
# plot the Norm Correction per step. Use in the following way:
#    plot_NormCorr.sh output
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
			echo "outfile[ $j ] = ${outfile[$j]} --- nsteps: ` cat ${outfile[$j]} | grep "Norm correction" | wc -l ` "
		fi
	done
	nfiles=$j

	if [ ` cat ${outfile[@]} | grep "Warning decelerated by" | wc -l ` -gt "0" ]; then
		cat ${outfile[@]} | grep "Norm correction" -A 1 | awk '{
			if (/Norm/) {
				step++; normc=$4; dec=1.; scrivi=1; }
			else {
				if (/Warning/) dec=$4;
				if (scrivi==1) scrivi=2; }
			if (scrivi==2) {
				epsi=dec*normc;
				print step, normc, epsi; scrivi=0; }
		}'
	else
		cat ${outfile[@]} | grep "Norm correction" | awk '{ if (NF==5) {
			print NR, $4, $5; }
		else {
			print NR, $4, $4;} }'
	fi > _XXplotDataXX

	echo "set title \" plot of $outfile \" " > _plot.gpl

cat << EOF >> _plot.gpl
set xlabel "optimization step"
set ylabel " | D wf | / | wf | "
set xrange [*:*]
set yrange [*:*]
plot "_XXplotDataXX" u 1:2 ti "Norm Correction" w l , "" u 1:3 ti "epsi (effective)" w l
pause -1

EOF

	echo
	cat _XXplotDataXX | awk '{ if ($2==0.) {cont++} }END{ print "Number Rejected move increasing energy: ",cont,"/",NR," = "100.*cont/NR,"%";}'

	echo
	if [ ` cat ${outfile[@]} | grep "Warning decelerated by" | wc -l ` -gt "0" ]; then
		cat ${outfile[@]} | grep "Warning decelerated by" | wc -l | awk '{ if ($1>0) print "Number decelerated steps: ",$1; }'
		cat ${outfile[@]} | grep "Warning decelerated by" | howmuch.sh 4 | awk '{ print "Average deceleration: ", $10, "+/-", $12; }'
	fi
	echo

	gnuplot _plot.gpl

else
	echo 
	echo "Script that reads one output file of a TurboRVB optimization and"
	echo "plot the Norm Correction per step. Use in the following way:"
	echo "   plot_NormCorr.sh output"
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
