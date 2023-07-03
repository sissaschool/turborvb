r"""
BEWARE!

This script should work with vanilla python 3.7,
no special libraries are allowed

Author Ot(t)o Kohul\'{a}k
"""

import os
import sys
from turbo_linter import turbo_linter

exclude = ["devel_tools"]

try:
    exclude.extend(os.environ["TURBO_LINTER_EXCLUDE"].split(":"))
except KeyError:
    pass

exclude = [os.path.abspath(x) for x in exclude]

acc = 0
for d, dn, fn in os.walk(os.getcwd() + "/" + sys.argv[1]):
    breakme = False
    for excl in exclude:
        if os.path.samefile(excl, os.path.commonpath([d, excl])):
            breakme = True
            break
    if breakme:
        continue
    for f in fn:
        if f.endswith(".f90"):
            acc += turbo_linter(d + "/" + f)
        if f.endswith(".f"):
            acc += turbo_linter(d + "/" + f)
exit(acc)
