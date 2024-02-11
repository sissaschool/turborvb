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
#
# Otto Kohulák created on 5th Nov. 2021.

# copy in the present directory makefuns
cp ../../src/c_adjoint_forward/makefun.f90 .
cp ../../src/c_adjoint_backward/makefun_b.f90 .
cp ../../src/c_adjoint_forward/makefun_pbc.f90 .
cp ../../src/c_adjoint_backward/makefun_pbc_b.f90 .
cp ../../src/c_adjoint_forward/makefun_bump.f90 .
cp ../../src/c_adjoint_backward/makefun_bump_b.f90 .

# process the mkefuns
python3 preprocess.py -c makefun0:[i0=0,indtmin=0,indtm=0]        -o makefun0.f90        makefun.f90
python3 preprocess.py -c makefun0_b:[i0=0,indtmin=0,indtm=0]      -o makefun0_b.f90      makefun_b.f90
python3 preprocess.py -c makefun0_pbc:[i0=0,indtmin=0,indtm=0]    -o makefun0_pbc.f90    makefun_pbc.f90
python3 preprocess.py -c makefun0_pbc_b:[i0=0,indtmin=0,indtm=0]  -o makefun0_pbc_b.f90  makefun_pbc_b.f90
python3 preprocess.py -c makefun0_bump:[i0=0,indtmin=0,indtm=0]   -o makefun0_bump.f90   makefun_bump.f90
python3 preprocess.py -c makefun0_bump_b:[i0=0,indtmin=0,indtm=0] -o makefun0_bump_b.f90 makefun_bump_b.f90

# copy makefuns0
cp makefun0.f90        ../../src/c_adjoint_forward/makefun0.f90
cp makefun0_b.f90      ../../src/c_adjoint_backward/
cp makefun0_pbc.f90    ../../src/c_adjoint_forward/makefun0_pbc.f90
cp makefun0_pbc_b.f90  ../../src/c_adjoint_backward/
cp makefun0_bump.f90   ../../src/c_adjoint_forward/makefun0_bump.f90
cp makefun0_bump_b.f90 ../../src/c_adjoint_backward/
