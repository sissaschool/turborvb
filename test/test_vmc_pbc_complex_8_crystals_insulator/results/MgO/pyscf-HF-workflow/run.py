
from pyscf_wrapper import Pyscf_wrapper

# input variables
pyscf_chkfile="pyscf.chk"
structure_file="1000053.cif"

# input variables
omp_num_threads=72
init_guess="minao"
cell_precision=1e-08
multigrid_fftdf=False
level_shift_factor=0.0
charge=0
spin=0
basis="ccecp-ccpvdz"
ecp="ccecp"
scf_method="HF"
dft_xc="NA"
solver_newton=False
MP2_flag=False
CCSD_flag=False
pyscf_output="out.pyscf"
twist_average=False
exp_to_discard=0.1
kpt=[0.25, 0.25, 0.25]
kpt_grid=[1, 1, 1]
smearing_method="fermi"
smearing_sigma=0.0

pyscf_calc=Pyscf_wrapper(
                        structure_file=structure_file,
                        chkfile=pyscf_chkfile,
                        )

pyscf_calc.run_pyscf(
                  omp_num_threads=omp_num_threads,
                  init_guess=init_guess,
                  cell_precision=cell_precision,
                  multigrid_fftdf=multigrid_fftdf,
                  level_shift_factor=level_shift_factor,
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
                