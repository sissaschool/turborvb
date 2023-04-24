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

# Andrea Tirelli created on 28th Aug. 2020.

n_proc=$(cat $1 | grep "Number of mpi proc"| awk '{print $6}')
nweight=$(cat $1 | grep "nweight before" | awk '{print $4}')

output_fname=${1}_lr_mpi_proc_${n_proc}_nweight_${nweight}
cat $1 | grep "before adjust" | awk '{print $4}' >${output_fname}
if [[ ($# -ne 1) ]]; then
	echo 'set terminal "png"
	set output "'${output_fname}'.png"
	plot "'${output_fname}'" with lines' >gnu_command
	gnuplot gnu_command
	if [[ ($# -eq 2) ]]; then
		mv ${output_fname}.png $2
	fi
else
	echo 'plot "'${output_fname}'" with lines' >gnu_command
	gnuplot --persist gnu_command

fi

rm ${output_fname}
rm gnu_command
