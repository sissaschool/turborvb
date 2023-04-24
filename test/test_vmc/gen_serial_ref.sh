#!/bin/bash 

TURBORVB="../../bin/turborvb-serial.x"
READF="../../bin/readf.x"
FORT21=REFERENCE_fortXXI
REF_FORCEVMC=forces_REFERENCE
OUT=out_true.o
OUT_FORCE=out_true_forcevmc.o

debug_root=`pwd`

cd $debug_root
$TURBORVB < datasvmc.d > $OUT
echo "0 1 1 1" | $READF >& /dev/null
cp fort.21 ${FORT21}
$FORCEVMC 1 1 1 > $OUT_FORCE
cp forces_vmc.dat $REF_FORCEVMC
cd $debug_root

