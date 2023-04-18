r"""
BEWARE!

This script should work with vanilla python 3.7,
no special libraries are allowed

Author Ot(t)o Kohul\'{a}k
"""

import sys
import re

TYPE_FREEFORM = 1
TYPE_FIXEDFORM = 2


def check_reading(line, read):
    if re.match("^[ ]*[!c*][ ]*TL[ ]+on", line):
        return True
    if re.match("^[ ]*[!c*][ ]*TL[ ]+off", line):
        return False
    return read


def get_name(string):
    if not isinstance(string, str):
        return "Unamed"
    ret = re.search("[ ]*name:[ ]+([A-z ]+)", string)
    return ret.groups(1)[0]


class turbo_linter:
    def __new__(cls, filename):

        result = 0

        def rule_opener(cls, rules):
            for r in rules:
                yield getattr(cls, r)

        rules = [
            f
            for f in dir(cls)
            if callable(getattr(cls, f)) and f.startswith("rule_")
        ]

        try:
            form = {"f": TYPE_FIXEDFORM, "f90": TYPE_FREEFORM}[
                filename.split(".")[-1].lower()
            ]
        except KeyError:
            print(
                f"Uknown extension of the source file {filename}",
                file=sys.stderr,
            )
            exit(1)

        with open(filename, "r") as fhandle:
            read = True
            for linenumber, line in enumerate(fhandle):
                read = check_reading(line, read)
                if not read:
                    continue
                for rule in rule_opener(cls, rules):
                    if not rule(line, form):
                        name = get_name(rule.__doc__)
                        index = rule.__name__[5:]

                        # Check if it is a number
                        int(index)

                        print(
                            f" -- Rule #{index} {name.join(['`']*2):>13} for file {filename} at line {linenumber+1}!",
                            file=sys.stderr,
                        )
                        result = 1

        return result

    @classmethod
    def get_ruledoc(cls, val):
        val = f"{int(val):04d}"
        try:
            return getattr(cls, f"rule_{val}").__doc__
        except AttributeError:
            return None
    
    """
    @classmethod
    def rule_0001(cls, line, t):
        '''
        This rule checks the lengths of the lines

        <= 132 words for the freeform
        <=  72 words for the fixed form

        name: length check
        '''
        if t == TYPE_FREEFORM:
            if not len(line) <= 132:
                return False
        if t == TYPE_FIXEDFORM:
            if not len(line) <= 72:
                return False
        return True
    """
    
    @classmethod
    def rule_0009(cls, line, t):
        """
        This rule checks trailing white spaces

        name: Trailing white space
        """
        if re.search(" $", line, re.IGNORECASE):
            return False
        return True

    @classmethod
    def rule_0010(cls, line, t):
        """
        This rule checks implicit nones

        In freeform they has to non-capital, whereas in fixed form
        has to be capital. Only one space is allowed in between
        the words

        name: Implicit none
        """
        if re.match("^\s*implicit\s+none\s*$", line, re.IGNORECASE):
            if t == TYPE_FREEFORM:
                if not re.search("^\s*implicit\s+none\s*$", line):
                    return False
            if t == TYPE_FIXEDFORM:
                if not re.search("^\s*IMPLICIT\s+NONE\s*$", line):
                    return False
        return True

    @classmethod
    def rule_0050(cls, line, t):
        """
        This rule checks end dos

        In freeform they has to non-capital, whereas in fixed form
        has to be capital. Exactly one space is allowed in between
        the words

        name: End do
        """
        if re.match("^\s*end\s*do\s*$", line, re.IGNORECASE):
            if t == TYPE_FREEFORM:
                if not re.search("^\s*end\sdo\s*$", line):
                    return False
            if t == TYPE_FIXEDFORM:
                if not re.search("^\s*END\sDO\s*$", line):
                    return False
        return True

    @classmethod
    def rule_0051(cls, line, t):
        """
        This rule checks end dos

        In freeform they has to non-capital, whereas in fixed form
        has to be capital. Exactly one space is allowed in between
        the words

        name: End if
        """
        if re.match("^\s*end\s*if\s*$", line, re.IGNORECASE):
            if t == TYPE_FREEFORM:
                if not re.search("^\s*end\sif\s*$", line):
                    return False
            if t == TYPE_FIXEDFORM:
                if not re.search("^\s*END\sIF\s*$", line):
                    return False
        return True

    @classmethod
    def rule_0052(cls, line, t):
        """
        This rule checks dos

        In freeform they has to non-capital, whereas in fixed form
        has to be capital. Exactly one space is allowed in between
        the words

        name: Do
        """
        if re.match("^[ ]*do[ ]", line, re.IGNORECASE):
            if t == TYPE_FREEFORM:
                if not re.search("^[ ]*do[ ]", line):
                    return False
            if t == TYPE_FIXEDFORM:
                if not re.search("^[ ]*DO[ ]", line):
                    return False
        return True

    @classmethod
    def rule_0053(cls, line, t):
        """
        This rule checks omparision operators

        The list of checked operators: lt, le, gt, ge, ne, eq

        In freeform they has to non-capital, whereas in fixed form
        has to be capital.

        name: Comparion operator
        """
        operators = ["lt", "le", "gt", "ge", "ne", "eq"]

        for op in operators:
            if re.match(f".*\.{op}\..*", line, re.IGNORECASE):
                if t == TYPE_FREEFORM:
                    if not re.search(f"\.{op}\.", line):
                        return False
                if t == TYPE_FIXEDFORM:
                    if not re.search(f"\.{op.upper()}\.", line):
                        return False
        return True


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Not enough arguments", file=sys.stderr)
        exit(1)

    if len(sys.argv) == 3:
        if sys.argv[1] == "-e":
            message = turbo_linter.get_ruledoc(sys.argv[2])
            if message:
                print(message)
                exit(0)
            print("Uknown rule", file=sys.stderr)
            exit(-1)

    filename = sys.argv[1]

    exit(turbo_linter(filename))
