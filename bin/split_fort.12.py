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

# Kosuke Nakano created on 22nd Feb. 2019.
# program:split_fort.12.py
# purpose:split_fort.12.py

# python modules
import sys
import numpy as np
from scipy.io import FortranFile


def main():

    fort12_input_name = str(sys.argv[1])
    fort12_output_name = str(sys.argv[2])
    kp_info = str(sys.argv[3])

    # fort12_input_name = "fort.12"
    # fort12_output_name = "fort.12"
    # kp_info = 'kp_info.dat'

    f = FortranFile(fort12_input_name, "r")
    a = f.read_reals(dtype="float64")
    column_length = len(a)
    f.close()

    # start reading fort.12
    head = ("head", "<i")
    tail = ("tail", "<i")

    dt = np.dtype([head, ("a", "<{}d".format(column_length)), tail])
    fd = open(fort12_input_name, "r")

    fort12_original = np.fromfile(fd, dtype=dt, count=-1)

    fd.close()

    with open(kp_info, "r") as kpt_file:
        num_of_kpt = int(kpt_file.readline().split()[0])

    for k_num in range(num_of_kpt):
        fort12_original[k_num::num_of_kpt].tofile(
            fort12_output_name + "_K{}".format(k_num + 1)
        )


if __name__ == "__main__":
    main()
