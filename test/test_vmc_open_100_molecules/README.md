# Tests for 38 molecules

The molecules included in this directory were used for checking the consistency between PySCF v.2.0.1 and TurboRVB (VMC w/o Jastrow) calculations. Thus, the fort.10s and thier output files should be correct and be the reference files for future testing.

Each directory contains

- pyscf-HF-workflow  ....   PySCF calcualtion
- vmc-HF-workflow    ....   VMC w/o Jastrow calculation
- vmc-HF-debug       ....   VMC w/o Jastrow calcualtion

Only the third directory is used for the debug, where a very small ngen (MCMC steps) used for reducing the computational cost. Thus, the obtained energy in vmc-HF-debug/ is not consistent with those in pyscf-HF-workflow/ and vmc-HF-workflow/.