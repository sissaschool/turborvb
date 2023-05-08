# Copyright (C) 2022 TurboRVB group
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import csnake as cs

def add_license(cw):
    cw.start_comment()
    cw.add("""
  Copyright (C) 2022 TurboRVB group

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
           """)
    cw.end_comment()

def add_notice(cw):
    cw.start_comment()
    cw.add("""

  This wrapper is generated by make_wrapper.py script

  Do not edit, regenerate it instead.

           """)
    cw.end_comment()

def add_complex(cw):
    compl = cs.Struct("complex", typedef = False)
    compl.add_variable(cs.Variable("re", "double"))
    compl.add_variable(cs.Variable("im", "double"))
    cw.add_struct(compl)

def add_cudasync(cw):
    fun = cs.Function("cudasync_", "void")
    fun.add_code("cudaDeviceSynchronize();")
    cw.add_function_definition(fun)

def add_cusolver_handle_init(cw):
    fun = cs.Function("cusolver_handle_init_",
                      "void",
                      arguments = [ cs.Variable("handle",
                                                "long int *")])
    fun.add_code(("cusolverDnHandle_t *h;",
                  "h = (cusolverDnHandle_t*) malloc(sizeof(cusolverDnHandle_t));",
                  "cusolverDnCreate(h);",
                  "*handle = (long int) h;",
                 ))
    cw.add_function_definition(fun)

def add_cusolver_handle_destroy(cw):
    fun = cs.Function("cusolver_handle_destroy_",
                      "void",
                      arguments = [ cs.Variable("handle",
                                                "long int *")])
    fun.add_code(("cusolverDnHandle_t *h = (cusolverDnHandle_t *) *handle;",
                  "cusolverDnDestroy(*h);",
                  "free(h);",
                 ))
    cw.add_function_definition(fun)

def default_includes(cw):
    cw.include("<stdlib.h>")
    cw.include("<stddef.h>")
    cw.include("<ctype.h>")