#!/usr/bin/env python3
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
#
# Otto Kohulák created on 5th Nov. 2021.

import argparse
import sys
import re

try:
    from icecream import ic
except:
    pass

class routine_keeper:

    title_endings = [ " *(?:[A-z0-9_])*",
                      " *::*",
                      " *\*",
                      " *,",
                      " *\("]

    def __init__(self):
        self.clear()

    def __add__(self, line):
        line = line.rstrip()
        if not self.head_close:
            if routine_keeper.filter(line):
                self.header.append(line)
                return
        if len(line) == 0: return
        self.head_close = True
        self.lines.append(line)
        return

    def clear(self):
        self.lines = []
        self.header = []
        self.head_close = False

    @classmethod
    def filter_head(cls, line):
        if re.match("^ *subroutine ", line.lower()):
            return True
        return False

    @classmethod
    def filter_title(cls, line):
        if re.match("^ *!", line.lower()):
            return True
        if re.match("^ *\&", line.lower()):
            return True
        return False

    @classmethod
    def filter_use(cls, line):
        if re.match("^ *use (?:[A-z0-9_])*", line.lower()):
            return True
        return False

    @classmethod
    def filter_implicitnone(cls, line):
        if re.match("^ *implicit *none", line.lower()):
            return True
        return False

    @classmethod
    def filter_type(cls, line):
        l = line.lower()
        if re.match("^ *type *\(.*\)", l):
            return True
        return False

    @classmethod
    def filter_logical(cls, line):
        l = line.lower()
        for te in cls.title_endings:
            if re.match("^ *logical" + te, l):
                return True
        return False

    @classmethod
    def filter_character(cls, line):
        l = line.lower()
        for te in cls.title_endings:
            if re.match("^ *character" + te, l):
                return True
        return False

    @classmethod
    def filter_double(cls, line):
        l = line.lower()
        for te in cls.title_endings:
            if re.match("^ *double precision" + te, l):
                return True
        return False

    @classmethod
    def filter_complex(cls, line):
        l = line.lower()
        for te in cls.title_endings:
            if re.match("^ *complex" + te, l):
                return True
        return False

    @classmethod
    def filter_real(cls, line):
        l = line.lower()
        for te in cls.title_endings:
            if re.match("^ *real" + te, l):
                return True
        return False

    @classmethod
    def filter_integer(cls, line):
        l = line.lower()
        for te in cls.title_endings:
            if re.match("^ *integer" + te, l):
                return True
        return False

    @classmethod
    def filter(cls, line):
        methods = [ name for name in dir(cls) if callable(getattr(cls, name)) ]
        methods = [ name for name in methods if name.startswith("filter_") ]
        for method in methods:
            if getattr(cls, method)(line):
                return True
        return False

    def get_lines(self, case):
        for l in self.lines:
            for key, value in case.items():
                key = key.strip()
                r = f"([^A-z0-9_])({key})([^A-z0-9_]|$)"
                rs = f"\g<1>{value}\g<3>"
                l = re.sub(r, rs, l)
                l = re.sub(r, rs, l)
            yield l

    def process_header(self, cases, name):
        header = "\n".join(self.header)
        rex = "(?i)(?:.*)SUBROUTINE *([A-z0-9_]+) *\(([\s\S]*?)\)([\s\S]*)"
        m = re.match(rex, header)
        name = m.group(1)
        variables = m.group(2)
        bodyofhead = m.group(3)
        rex = "& *\n *&"
        bodyofhead = re.sub(rex, "", bodyofhead)
        variables = re.sub(rex, "", variables)

        #process vars:

        variables = ",".join([ x.strip() for x in variables.split(",") if x.strip() not in cases ])

        for key, value in cases.items():
            bodyofhead = re.sub("(\(.*?)(" + key.strip() + ")(.*?\))", f"\g<1>{value}\g<3>", bodyofhead)
            bodyofhead = re.sub("(.*)(, *" + key.strip() + ") *([\n, ].*)", "\g<1>\g<3>", bodyofhead)
            bodyofhead = re.sub("(.*)(" + key.strip() + ") *,([\n ].*)", "\g<1>\g<3>", bodyofhead)

        header = f"SUBROUTINE {name} ({variables})\n{bodyofhead}"
        self.header = header.split("\n")

def p(infile, outfile, cases, name):
    message = """

!###########################################################
!#                                                         #
!#     This is code was generated by preprocess script     #
!#                                                         #
!###########################################################

    """

    entry_re = "(^ *! *PPC) *\{(.*)\}"

    start = False
    lk = routine_keeper()

    with open(infile, "r") as fi:
        with open(outfile, "w") as fo:
            fo.write(message)
            lines = [ l for l in fi ]
            for ii, line in enumerate(lines):
                rex = "(?i)^( *subroutine *)([A-z0-9_]+)( *\()"
                match = re.match(rex, line.lower())
                if match:
                    start = True
                    line = str(re.sub(rex, f"\g<1>{name}\g<3>", line))
                if start:
                    if re.match("^ *end *subroutine", line.lower()):
                        line = f"END SUBROUTINE {name}"
                        start = False
                    elif re.match("^ *end *$", line.lower()):
                        start = False
                    else:
                        lk + line
                        continue
                if len(lk.header) > 0:
                    lk.process_header(cases, name)
                    for l in lk.header:
                        fo.write(l + "\n")
                    for ii, cc in enumerate([cases]):
                        for l in lk.get_lines(cc):
                            fo.write(l.replace("\n", "") + "\n")
                    lk.clear()
                fo.write(line)
def main():

    parser = argparse.ArgumentParser(description='Preproccessng')
    parser.add_argument('input', metavar='INPUT', type=str)
    parser.add_argument('-o', metavar='OUTPUT', type=str, default = None)
    parser.add_argument('-c', metavar='CASES', type=str, default = None)

    args = parser.parse_args()

    infile = args.input

    if args.c is not None:
        rex = "([0-9A-z_]+):\[(.*)\]"
        match = re.match(rex, args.c)
        cases = []
        if match:
            name = match.group(1)
            cc = match.group(2)
            cases = { x.split("=")[0] : x.split("=")[1] for x in cc.split(",") }
            output = f"{name}.f90"
            if args.o is not None:
                output = args.o
            p(infile, output, cases, name)
            sys.exit(0)
        else:
            print("Bad case string")
            sys.exit(1)
    else:
        print("Please provide case string")


if __name__ == "__main__":
    main()
