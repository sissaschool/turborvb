
from pyscf_wrapper import Pyscf_wrapper

# input variables
pyscf_chkfile="pyscf.chk"
structure_file="F.xyz"

# input variables
omp_num_threads=128
charge=0
spin=1
basis="ccecp-ccpv6z"
ecp="ccecp"
scf_method="HF"
dft_xc="NA"
solver_newton=False
MP2_flag=False
CCSD_flag=False
pyscf_output="out.pyscf"
twist_average=False
exp_to_discard=0.0
kpt=[0.0, 0.0, 0.0]
kpt_grid=[1, 1, 1]
smearing_method="fermi"
smearing_sigma=0.0

pyscf_calc=Pyscf_wrapper(
                        structure_file=structure_file,
                        chkfile=pyscf_chkfile,
                        )

pyscf_calc.run_pyscf(
                  omp_num_threads=omp_num_threads,
                  charge=charge,
                  spin=spin,
                  basis=basis,
                  ecp=ecp,
                  scf_method=scf_method,
                  dft_xc=dft_xc,
                  solver_newton=solver_newton,
                  MP2_flag=MP2_flag,
                  CCSD_flag=CCSD_flag,
                  pyscf_output=pyscf_output,
                  twist_average=twist_average,
                  exp_to_discard=exp_to_discard,
                  kpt=kpt,
                  kpt_grid=kpt_grid,
                  smearing_method=smearing_method,
                  smearing_sigma=smearing_sigma
                  )
                