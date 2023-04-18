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
        line = infile.readline()
        while re.search("Energy =", line, re.IGNORECASE) is None:
            line = infile.readline()
        en.append(float(line.split()[2]))
        enerr.append(float(line.split()[3]))
        while re.search("Variance square = ", line, re.IGNORECASE) is None:
            line = infile.readline()
        var.append(float(line.split()[3]))
        varerr.append(float(line.split()[4]))
        while re.search("Est. corr. time  = ", line, re.IGNORECASE) is None:
            line = infile.readline()
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
