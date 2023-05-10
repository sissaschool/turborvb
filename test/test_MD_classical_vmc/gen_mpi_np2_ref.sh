#!/bin/bash 

TURBORVB="mpirun -np 2 ../../bin/turborvb-mpi.x"
READF="../../bin/readf.x"
FORT21=REFERENCE_fortXXI_mpi_np2
OUT=out_true_mpi_np2.o
POSITION_OUT=position_true_mpi_np2.dat
VELOCITY_OUT=velocity_true_mpi_np2.dat
debug_root=`pwd`

cd $debug_root
cp fort.10_ori fort.10
$TURBORVB < datasmin.d > $OUT
echo "0 1 1 1" | $READF >& /dev/null 
cp fort.21 ${FORT21}
cp position.dat ${POSITION_OUT}
cp velocity.dat ${VELOCITY_OUT}
cd $debug_root
