# Tests for 38 molecules

The molecules included in this directory were used for checking the consistency between PySCF v.2.0.1 and TurboRVB-DFT(prep) calculations. Thus, the fort.10s and thier output files should be correct and be the reference files for future testing.

Each directory contains

- pyscf-LDA-workflow  ....   PySCF calcualtion
- dft-LDA-workflow    ....   prep calculation
- dft-LDA-debug       ....   prep calcualtion

Only the third directory is used for the debug, where a very coarse grid and a small box are used for reducing the computational cost. Thus, the obtained energy in dft-LDA-debug/ is not consistent with those in dft-LDA-workflow/ and pyscf-LDA-workflow/.
