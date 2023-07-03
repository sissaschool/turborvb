#!/bin/bash 

echo "#########################################################################"
echo "#                                                                        "
echo "#  TurboRVB Debugging tool                                               "
echo "#                                                                        "
echo "#  How to use:                                                           "
echo "#    After compiling TurboRVB, you can test the compiled codes by        "
echo "#                                                                        "
echo "#    ./jobcheck-mpi-np2.sh                                               "
echo "#                                                                        "
echo "#    Created by Andrea Zen.                                              "
echo "#    Modified by Kosuke Nakano on 10 May 2020.                           "
echo "#    Modified by Kosuke Nakano on 5 July 2020.                           "
echo "#    Modified by Kosuke Nakano on 9 Apr. 2023.                           "
echo "#                                                                        "
echo "#    If you do not see any "diff" and "ERR" in QMC outputs               "
echo "#    and there is no difference in DFT energies,                         "
echo "#    the codes are correctly compiled and modified.                      "
echo "#                                                                        "
echo "#########################################################################"

#################################################################################
PREP="../../bin/prep-mpi.x"
TURBORVB="../../bin/turborvb-mpi.x"
CHECKDER="../../bin/testad.x"
READF="../../bin/readf.x"
FORCEVMC="../../bin/forcevmc.sh"
ROUND_OFF=6
OUT=out_mpi_np2.o
OUT_FORCEVMC=out_forcevmc_mpi_np2.o
OUTDER=out_der
PREFIX="mpirun -np 2"
##################################################################################

debug_root=`pwd`
turborvb_root=$(cd $(dirname $0) && cd ../ && pwd)
cd $debug_root
echo " TurboRVB root = $turborvb_root"
echo " PREP binary = $PREP"
echo " TURBORVB binary = $TURBORVB"
echo " Test dir = $debug_root" 
echo " Warning: This debugging tool rounds off < 10**-${ROUND_OFF}".
echo " Warning: You can set ROUND_OFF in this script.".
echo " Debugging starts:"

exit_scores_arr=()

echo "################################################################################" 

cd test_dft_open
./cm.test.sh $PREP $OUT out_true.o REFERENCE_fort.10_new $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

#cd test_dft_open_38_molecules
#./cm.test.sh $PREP $OUT out_true.o REFERENCE_fort.10_new $ROUND_OFF "$PREFIX"
#exit_code=$?
#exit_scores_arr+=($exit_code)
#if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
#echo "################################################################################"
#cd $debug_root

cd test_dft_pbc_gamma
./cm.test.sh $PREP $OUT out_true.o REFERENCE_fort.10_new $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_dft_pbc_twist_real
./cm.test.sh $PREP $OUT out_true.o REFERENCE_fort.10_new $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_dft_pbc_twist_complex
./cm.test.sh $PREP $OUT out_true.o REFERENCE_fort.10_new $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_vmc
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT $OUT_FORCEVMC out_true_mpi_np2.o out_true_forcevmc_mpi_np2.o REFERENCE_fortXXI_mpi_np2 forces_REFERENCE_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
./cm.test2.sh $CHECKDER
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_vmc_open_100_molecules
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################"
cd $debug_root

cd test_lrdmc
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
./cm.test2.sh $CHECKDER
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_lrdmc_new
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
./cm.test2.sh $CHECKDER
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_vmcpbc_tilted
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
./cm.test2.sh $CHECKDER
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
./cm.test3.sh $CHECKDER
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_lrdmc_tilted
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_vmc_pfaff
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
./cm.test2.sh $CHECKDER
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
./cm.test3.sh $CHECKDER
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_lrdmc_pfaff
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
./cm.test2.sh $CHECKDER
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_vmc_pbc_gamma_8_crystals_insulator
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################"
cd $debug_root

cd test_vmc_pbc_complex_8_crystals_insulator
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################"
cd $debug_root

cd test_vmc_pbc_trim_8_crystals_insulator
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################"
cd $debug_root

cd test_vmcopt_jastrow_1-3b_linear_method
./cm.test1.sh $TURBORVB $READF NA $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################"
cd $debug_root

cd test_vmcopt_jastrow_1-3b_stochastic_reconfiguration
./cm.test1.sh $TURBORVB $READF NA $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################"
cd $debug_root

cd test_vmcopt_agp_matrix_linear_method
./cm.test1.sh $TURBORVB $READF NA $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################"
cd $debug_root

cd test_vmcopt_agp_matrix_stochastic_reconfiguration
./cm.test1.sh $TURBORVB $READF NA $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################"
cd $debug_root

cd test_vmcopt_sd_MOs_linear_method
./cm.test1.sh $TURBORVB $READF NA $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################"
cd $debug_root

cd test_lrdmc_dla
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_lrdmc_dltm
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_lrdmc_la
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_lrdmc_tmove
./cm.test1.sh $TURBORVB $READF $FORCEVMC $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################" 
cd $debug_root

cd test_MD_classical_vmc
./cm.test1.sh $TURBORVB $READF NA $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################"
cd $debug_root

cd test_MD_quantum_vmc
./cm.test1.sh $TURBORVB $READF NA $OUT out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 $ROUND_OFF "$PREFIX"
exit_code=$?
exit_scores_arr+=($exit_code)
if [ $exit_code -ne 0 ]; then echo " Warn: Job failed."; fi
echo "################################################################################"
cd $debug_root

success=0
failure=0
for exit_score in ${exit_scores_arr[@]}
do
if [ $exit_score -eq 0 ]; then
    success=`expr $success + 1`
else
    failure=`expr $failure + 1`
fi
done

echo ""
echo "################################################################################"
echo "Results summary"
echo "  All the tests: ${#exit_scores_arr[@]}"
echo "  Success: ${success}"
echo "  Failures: ${failure}"
echo "################################################################################"

echo ""

if [ $failure -gt 0 ]; then
    echo "Some tests failed. Please compile TurboRVB again with proper settings :-("
else
    echo "All the tests completed. Enjoy Turbo :-)"
fi
