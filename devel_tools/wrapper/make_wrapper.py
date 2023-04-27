"""
This script generates fortran wrapper for cudalib C API

Ot(t)o Kohu\'ak
"""

import csnake as cs
from includer import add_this

cw = cs.CodeWriter()

from basics import ( add_complex,
                     default_includes,
                     add_notice,
                   )

add_notice(cw)

default_includes(cw)
add_complex(cw)

b_cublas = True
b_cusolver = True

if b_cublas:
    from datapack import cublas

    cw.start_if_def(f"_CUBLAS")
    cw.include("<cublas.h>")
    cw.end_if_def()

    for entry in cublas:
        add_this(cw, entry)

if b_cublas:
    from datapack import cusolver

    cw.start_if_def(f"_CUSOLVER")
    cw.include('"cusolverDn.h"')
    cw.end_if_def()

    cw.start_if_def(f"_CUSOLVER")
    from basics import add_cudasync
    from basics import ( add_cusolver_handle_init,
                         add_cusolver_handle_destroy,
                       )

    add_cusolver_handle_init(cw)
    add_cusolver_handle_destroy(cw)
    add_cudasync(cw)
    cw.end_if_def()

    for entry in cusolver:
        add_this(cw, entry)

print(cw)
