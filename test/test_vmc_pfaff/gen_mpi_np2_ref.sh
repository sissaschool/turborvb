#!/bin/bash 

TURBORVB="mpirun -np 2 ../../bin/turborvb-mpi.x"
READF="../../bin/readf.x"
FORT21=REFERENCE_fortXXI_mpi_np2
OUT=out_true_mpi_np2.o

debug_root=`pwd`

cd $debug_root
$TURBORVB < datasvmc.d > $OUT
echo "0 1 1 1" | $READF >& /dev/null 
cp fort.21 ${FORT21}
cd $debug_root
