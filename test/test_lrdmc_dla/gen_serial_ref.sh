#!/bin/bash 

TURBORVB="../../bin/turborvb-serial.x"
READF="../../bin/readf.x"
FORT21=REFERENCE_fortXXI
OUT=out_true.o

debug_root=`pwd`

cd $debug_root
$TURBORVB < datasfn.d > $OUT
echo "0 1 1 1" | $READF >& /dev/null 
cp fort.21 ${FORT21}
cd $debug_root
