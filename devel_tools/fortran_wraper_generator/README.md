# Fortran GPU libs C Wrapper generator

This set of scripts generate `fortran.c`, where Fortran-C wrapper is stored for GPU library APIs. To generate `fortran.c` one has to execute `generate_wrapper.sh` script. The file will be stored in correct directory. To add interface other APIs please add entries in `datapack.py` or in `basics.py`.

## Requirements

The script were tested with:

 - python=3.9.2
 - csnake=0.3.5
