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

# Sandro Sorella created on 18th Nov. 2019.

import re
from math import sqrt

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

en = []
enerr = []
var = []
varerr = []
time = []
timeerr = []


for i in range(len(w)):
    pipfile = "pip0_K" + str(i + 1)
    with open(pipfile) as infile:
        lines = infile.readlines()
        # print(lines)

        for line in lines:
            # search energy (old)
            # if (re.search('[\s]+Energy[\s]*=[\s]*', line, re.IGNORECASE)):
            # 	en.append(float(line.split()[2]))
            # 	enerr.append(float(line.split()[3]))

            # search energy (new one)
            if re.search(
                "[\s]+Energy[\s]+\(ave\)[\s]*=[\s]*", line, re.IGNORECASE
            ):
                en.append(float(line.split()[3]))
                enerr.append(float(line.split()[4]))

            # search variance
            if re.search(
                "[\s]+Variance[\s]+square[\s]*=[\s]*", line, re.IGNORECASE
            ):
                # line = infile.readline()
                # print(line)
                var.append(float(line.split()[3]))
                varerr.append(float(line.split()[4]))

            # search corr. time
            if re.search(
                "[\s]+Est.*[\s]+corr.*[\s]+time[\s]*=[\s]*",
                line,
                re.IGNORECASE,
            ):
                # print(line)
                time.append(float(line.split()[4]))
                timeerr.append(float(line.split()[5]))


entot = 0.0
vartot = 0.0
timetot = 0.0
for i in range(len(w)):
    entot += w[i] * en[i]
    vartot += w[i] * var[i]
    timetot += w[i] * time[i]

enerrtot = 0.0
varrtot = 0.0
timertot = 0.0
for i in range(len(w)):
    den = enerr[i]
    enerrtot += (w[i] * den) ** 2
    den = varerr[i]
    varrtot += (w[i] * den) ** 2
    den = timeerr[i]
    timertot += (w[i] * den) ** 2
enerrtot = sqrt(enerrtot)
varrtot = sqrt(varrtot)
timertot = sqrt(timertot) / len(w)
timetot = timetot / len(w)

print("Total Energy = ", entot, " +/- ", enerrtot)
print("Variance = ", vartot, " +/- ", varrtot)
print("Correlation  Time (unit nbra) = ", timetot, " +/- ", timertot)
