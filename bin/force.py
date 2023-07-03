#!/usr/bin/env python3

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

# Sandro Sorella created on 30th Jun. 2016.

import numpy as np
import re

w = []
# read kp_info.dat
with open("kp_info.dat") as infile:
    line = infile.readline()
    while re.search("up spin electrons", line, re.IGNORECASE) is None:
        line = infile.readline()
    line = infile.readline()
    while re.search("down spin electrons", line, re.IGNORECASE) is None:
        w.append(float(line.split()[4]))
        line = infile.readline()

f = []
ferr = []
var = []
varerr = []

for i in range(len(w)):
    pipfile = "forces_K" + str(i + 1)
    with open(pipfile) as infile:
        fx = []
        fxerr = []
        for line in infile:
            if re.search("Force   = ", line, re.IGNORECASE) is not None:
                fx.append(float(line.split()[2]))
                fxerr.append(float(line.split()[3]))
    f.append(fx)
    ferr.append(fxerr)

w = np.array(w)
f = np.array(f)
ferr = np.array(ferr)
ftot = (w * f.T).sum(axis=1)

ferrtot = (w * ferr.T) ** 2
ferrtot = ferrtot.sum(axis=1)
ferrtot = ferrtot**0.5

for i in range(len(ftot)):
    print("Total Force[", i + 1, "] = ", ftot[i], " +/- ", ferrtot[i])
