#!/bin/bash 

TURBORVB="../../bin/turborvb-serial.x"
READF="../../bin/readf.x"
FORT21=REFERENCE_fortXXI
POSITION_OUT=position_true.dat
VELOCITY_OUT=velocity_true.dat
OUT=out_true.o

debug_root=`pwd`

cd $debug_root
cp fort.10_ori fort.10
$TURBORVB < datasmin.d > $OUT
echo "0 1 1 1" | $READF >& /dev/null 
cp fort.21 ${FORT21}
cp position.dat ${POSITION_OUT}
cp velocity.dat ${VELOCITY_OUT}
cd $debug_root
