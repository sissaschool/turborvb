# Tests of compilation
We recommend you should run `jobcheck_serial.sh` in this directory every time you compile `TurboRVB` on your machine. If you compile `TurboRVB` with the modern `CMake`, you can use the `ctest` command.

<!--
How to extract the following list.
ls | grep test | awk '{print "- " $0}'
-->

<!--
How to generate reference files.
PREP="../../bin/prep-serial.x"
TURBORVB="../../bin/turborvb-serial.x"
READF="../../bin/readf.x"
FORCEVMC="../../bin/forcevmc.sh"

=DFT check=
$PREP < prep.d > out_true.o
cp fort.10_new REFERENCE_fort.10_new

=VMC check=
$TURBORVB < datasvmc.d > $OUT
echo "0 1 1 1" | $READF >& /dev/null 
cp fort.21 REFERENCE_fortXXI

=LRDMC check=
$TURBORVB < datasfn.d > $OUT
echo "0 1 1 1" | $READF >& /dev/null 
cp fort.21 REFERENCE_fortXXI

-->

# The list of the tests

- test_dft_open
- test_dft_pbc_gamma
- test_dft_pbc_twist_complex
- test_dft_pbc_twist_real
- test_dgemm_NN128
- test_dgemm_NT128
- test_dgemm_TT128
- test_dgemm_b_NN128
- test_dgemm_b_NT128
- test_dgemm_b_TT128
- test_dgemv_N256
- test_dgemv_T256
- test_dger2_256
- test_dger_256
- test_dgetrf_128
- test_dgetrfi_128
- test_dskmv_X256
- test_dtrsm_LRT128
- test_dtrsm_ULN128
- test_lrdmc
- test_lrdmc_new
- test_lrdmc_pfaff
- test_lrdmc_tilted
- test_offload_code_assessment
- test_offload_code_forcycle
- test_offload_data
- test_offload_if
- test_offload_pointer_trancription
- test_openmp_reduction
- test_vmc
- test_vmc_pfaff
- test_vmcpbc_tilted
- test_zgemm_CT128
- test_zgemm_NN128
- test_zgemm_NT128
- test_zgemm_b_NN128
- test_zgemm_b_NT128
- test_zgemm_b_TT128
- test_zgemv_C256
- test_zgemv_N256
- test_zgemv_T256
- test_zger2_256
- test_zgeru_256
- test_zgetrf_128
- test_zgetrfi_128
- test_zskmv_X256
- test_ztrsm_LLC128
- test_ztrsm_URN128
- testderiv
- testmp
- testmpi