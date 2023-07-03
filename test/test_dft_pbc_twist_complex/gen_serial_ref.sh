#!/bin/bash 

PREP=../../bin/prep-serial.x
FORT10=REFERENCE_fort.10_new
OUT=out_true.o

debug_root=`pwd`

cd $debug_root
$PREP < prep.d > $OUT
cp fort.10_new ${FORT10}
grep ERR $OUT
cd $debug_root
