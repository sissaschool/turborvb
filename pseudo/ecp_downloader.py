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
# Kosuke Nakano, created on 25th Apr. 2023.

# python modules
import os
import shutil
import string
import linecache
import subprocess
import re
import glob
import tempfile
import pathlib
import random
import git

# Logger
from logging import getLogger, StreamHandler, Formatter

# pymatgen
from pymatgen.core.periodic_table import Element
from pymatgen.core.periodic_table import ElementBase

logger = getLogger("pseudo-downloader")

turbo_pp_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "pseudo"))
turbo_pp_bfd_gammes_format_dir = os.path.join(turbo_pp_dir, "bfd_gammes_format")
turbo_pp_ccecp_gammes_format_dir = os.path.join(turbo_pp_dir, "ccecp_gammes_format")
turborvb_bin_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "bin")
)
tmp_dir = os.path.abspath(os.path.join(turbo_pp_dir, "pp_temp"))

chemical_symbols = [
    # 0
    "X",
    # 1
    "H",
    "He",
    # 2
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    # 3
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    # 4
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    # 5
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    # 6
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    # 7
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]


# ccECP PP downloader
class ccECP:
    # PP Downloader, downloading ccECP pseudo potential information from "https://github.com/QMCPACK/pseudopotentiallibrary.git"
    C_URL = "https://github.com/QMCPACK/pseudopotentiallibrary.git"

    def __init__(
        self,
        pseudo_potential_output_dir=turbo_pp_ccecp_gammes_format_dir,
    ):
        self.pseudo_potential_output_dir = pseudo_potential_output_dir

    def to_file(self, element_list=None):
        os.makedirs(self.pseudo_potential_output_dir, exist_ok=True)

        tempdir = pathlib.Path(tempfile.mkdtemp())
        git.Repo.clone_from(self.C_URL, tempdir)
        # logger.debug(list((tempdir/"recipes").glob("**/*")))
        for p in (tempdir / "recipes").glob("**/*"):
            if re.match("[A-z]{1,2}\.ccECP\.gamess", p.name):
                element = str(p.parent.parent.name)
                typ = str(p.parent.name)
                typ_core = typ.replace("ccECP", "")
                logger.debug(typ)
                logger.debug(typ_core)

                if element_list is not None:
                    if element not in element_list:
                        continue
                logger.debug(p.name)
                logger.debug(p.name.replace(".ccECP.gamess", ""))
                with open(p, "r") as fhandle:
                    with open(
                        os.path.join(
                            self.pseudo_potential_output_dir,
                            f"{p.name.replace('.ccECP.gamess', '_ccECP'+typ_core+'.ecp')}",
                        ),
                        "w",
                    ) as fhandle_out:
                        fhandle_out.write(fhandle.read())

    def all_to_file(self, sleep_time=1):
        logger.info(f"Downloading all ccECP files to {os.path.abspath(turbo_pp_dir)}")
        os.makedirs(tmp_dir, exist_ok=True)
        self.to_file(element_list=chemical_symbols)
        ccecp_pp_list = glob.glob(
            os.path.join(turbo_pp_ccecp_gammes_format_dir, "*.ecp")
        )
        pp = Pseudopotentials.parse_pseudopotential_from_gamess_format_files(
            files=ccecp_pp_list
        )
        pp.set_cutoffs()
        pp.write_pseudopotential_turborvb_file(prefix="ccECP")
        shutil.rmtree(tmp_dir)


# ccECP PP downloader
class BFD:
    # PP Downloader, downloading BFD pseudo potential information from "https://github.com/TREX-CoE/BFD-ECP.git"
    C_URL = "https://github.com/TREX-CoE/BFD-ECP.git"

    def __init__(
        self,
        pseudo_potential_output_dir=turbo_pp_bfd_gammes_format_dir,
    ):
        self.pseudo_potential_output_dir = pseudo_potential_output_dir

    def to_file(self, element_list=None):
        os.makedirs(self.pseudo_potential_output_dir, exist_ok=True)

        tempdir = pathlib.Path(tempfile.mkdtemp())
        git.Repo.clone_from(self.C_URL, tempdir)
        # logger.debug(list((tempdir/"recipes").glob("**/*")))
        logger.info(os.path.join(tempdir, "ECP", "GAMESS"))
        for p in (tempdir / "ECP" / "GAMESS").glob("*"):
            if p.name in chemical_symbols:
                element = str(p.name)
                if element_list is not None:
                    if element not in element_list:
                        continue
                with open(p, "r") as fhandle:
                    with open(
                        os.path.join(
                            self.pseudo_potential_output_dir,
                            f"{element+'_BFD.ecp'}",
                        ),
                        "w",
                    ) as fhandle_out:
                        fhandle_out.write(fhandle.read())

    def all_to_file(self, sleep_time=1):
        logger.info(f"Downloading all BFD files to {os.path.abspath(turbo_pp_dir)}")
        os.makedirs(tmp_dir, exist_ok=True)
        self.to_file(element_list=chemical_symbols)
        bfd_pp_list = glob.glob(os.path.join(turbo_pp_bfd_gammes_format_dir, "*.ecp"))
        pp = Pseudopotentials.parse_pseudopotential_from_gamess_format_files(
            files=bfd_pp_list
        )
        pp.set_cutoffs()
        pp.write_pseudopotential_turborvb_file(prefix="BFD")
        shutil.rmtree(tmp_dir)


# PP class
class Pseudopotentials:
    def __init__(
        self,
        max_ang_mom_plus_1=[],
        z_core=[],
        cutoff=[],
        nucleus_index=[],
        element_list=[],
        ang_mom=[],
        exponent=[],
        coefficient=[],
        power=[],
    ):

        self.element_list = element_list
        self.max_ang_mom_plus_1 = max_ang_mom_plus_1
        self.z_core = z_core
        self.cutoff = cutoff
        self.nucleus_index = nucleus_index
        self.ang_mom = ang_mom
        self.exponent = exponent
        self.coefficient = coefficient
        self.power = power

        logger.debug(self.element_list)
        logger.debug(self.max_ang_mom_plus_1)
        logger.debug(self.z_core)
        logger.debug(self.cutoff)
        logger.debug(self.nucleus_index)
        logger.debug(self.ang_mom)
        logger.debug(self.exponent)
        logger.debug(self.coefficient)
        logger.debug(self.power)

        # assertion!!
        assert len(set(self.nucleus_index)) == len(self.max_ang_mom_plus_1)
        assert len(set(self.nucleus_index)) == len(self.z_core)
        assert len(set(self.nucleus_index)) == len(self.cutoff)
        assert len(self.ang_mom) == len(self.exponent)
        assert len(self.ang_mom) == len(self.coefficient)
        assert len(self.ang_mom) == len(self.power)

    @property
    def nuclei_num(self):
        return len(set(self.nucleus_index))

    @property
    def ecp_num(self):
        return len(self.ang_mom)

    @staticmethod
    def run(binary, input_name=None, output_name="out.o"):
        sys_env = os.environ.copy()
        if input_name is None:
            cmd = f"{binary} > {output_name}"
        else:
            cmd = f"{binary} < {input_name} > {output_name}"
        logger.debug(f"Execute command(s): {cmd}")

        subprocess.check_call(cmd, shell=True, env=sys_env)

    @staticmethod
    def return_element_symbol(atomic_number):
        atomic_number = int(float(atomic_number))
        return str(ElementBase.from_Z(atomic_number))

    @staticmethod
    def return_atomic_number(element):
        element = str(element)
        return Element(element).number

    @staticmethod
    def pygrep_lineno(file, keyword):
        cmd = f"grep '{keyword}' -n {file} | cut -d ':' -f 1"
        sys_env = os.environ.copy()
        try:
            lineno = int(subprocess.check_output(cmd, shell=True, env=sys_env)) - 1
        except ValueError:
            cmd = f"grep -c '' {file}"
            lineno = int(subprocess.check_output(cmd, shell=True, env=sys_env))

        return lineno

    @staticmethod
    def pygetline(
        filename, lineno, clearcache=True
    ):  # clearchache should be true!! as a default. # reasons for bugs.
        # logger.debug("get line!")
        line = linecache.getline(filename=filename, lineno=lineno + 1)
        if clearcache:
            linecache.clearcache()
        return line

    def write_pseudopotential_turborvb_file(self, prefix="ecp"):

        for i, nuclei_i in enumerate(set(self.nucleus_index)):
            output = []
            output.append("ECP\n")
            max_ang_mom_plus_1 = self.max_ang_mom_plus_1[i]
            element = self.element_list[i]  # not used.
            z_core = self.z_core[i]  # not used.
            cutoff = self.cutoff[i]
            output.append(
                "{:14d}  {:14f}  {:14d}\n".format(1, cutoff, max_ang_mom_plus_1 + 1)
            )

            nindex_list = [i for i, x in enumerate(self.nucleus_index) if x == nuclei_i]
            ang_mom_n = [self.ang_mom[n] for n in nindex_list]
            exponent_n = [self.exponent[n] for n in nindex_list]
            coefficient_n = [self.coefficient[n] for n in nindex_list]
            power_n = [self.power[n] for n in nindex_list]

            logger.debug(nindex_list)
            logger.debug(ang_mom_n)
            logger.debug(exponent_n)
            logger.debug(coefficient_n)
            logger.debug(power_n)

            for l in range(max_ang_mom_plus_1 + 1):
                output.append("{:14d} ".format(ang_mom_n.count(l)))
            output.append("\n")

            logger.debug(max_ang_mom_plus_1 + 1)
            for l in range(max_ang_mom_plus_1 + 1):
                lindex_list = [i for i, x in enumerate(ang_mom_n) if x == l]
                exponent_l = [exponent_n[n] for n in lindex_list]
                coefficient_l = [coefficient_n[n] for n in lindex_list]
                power_l = [power_n[n] for n in lindex_list]

                logger.debug(lindex_list)
                logger.debug(exponent_l)
                logger.debug(coefficient_l)
                logger.debug(power_l)

                for exp, coeff, pow in zip(exponent_l, coefficient_l, power_l):
                    output.append(
                        "{:14f}  {:14f}  {:14f}\n".format(coeff, pow + 2, exp)
                    )

            file = os.path.join(
                turbo_pp_dir,
                f"Z{z_core}_atomnumber{self.return_atomic_number(element=element.split('-')[0])}.{prefix}",
            )
            with open(file, "w") as f:
                f.writelines(output)

    def set_cutoffs(self, tollerance=0.00001):
        logger.debug(self.cutoff)
        new_cutoff = []
        current_dir = os.getcwd()

        try:
            num_string = 15
            rand_string = "".join(
                random.choices(string.ascii_letters + string.digits, k=num_string)
            )
            tmp_rand_dir = os.path.join(tmp_dir, rand_string)
            os.makedirs(tmp_rand_dir, exist_ok=True)
            os.chdir(tmp_rand_dir)
            for i, nuclei_i in enumerate(set(self.nucleus_index)):

                file = "pseudo.dat"
                with open(file, "x") as f:
                    output = []
                    output.append("ECP\n")

                    max_ang_mom_plus_1 = self.max_ang_mom_plus_1[i]
                    element = self.element_list[i]  # not used.
                    z_core = self.z_core[i]  # not used.
                    cutoff = self.cutoff[i]
                    output.append(
                        "{:14d}  {:14f}  {:14d}\n".format(
                            1, cutoff, max_ang_mom_plus_1 + 1
                        )
                    )

                    nindex_list = [
                        i for i, x in enumerate(self.nucleus_index) if x == nuclei_i
                    ]
                    ang_mom_n = [self.ang_mom[n] for n in nindex_list]
                    exponent_n = [self.exponent[n] for n in nindex_list]
                    coefficient_n = [self.coefficient[n] for n in nindex_list]
                    power_n = [self.power[n] for n in nindex_list]
                    logger.debug(nindex_list)

                    for l in range(max_ang_mom_plus_1 + 1):
                        output.append("{:14d} ".format(ang_mom_n.count(l)))
                    output.append("\n")

                    for l in range(max_ang_mom_plus_1 + 1):
                        lindex_list = [i for i, x in enumerate(ang_mom_n) if x == l]
                        exponent_l = [exponent_n[n] for n in lindex_list]
                        coefficient_l = [coefficient_n[n] for n in lindex_list]
                        power_l = [power_n[n] for n in lindex_list]

                        for exp, coeff, pow in zip(exponent_l, coefficient_l, power_l):
                            output.append(
                                "{:14f}  {:14f}  {:14f}\n".format(coeff, pow + 2, exp)
                            )

                    f.writelines(output)

                out_pp = "out_pp"
                cmd = (
                    f"echo {tollerance} | {os.path.join(turborvb_bin_root, 'pseudo.x')}"
                )
                logger.info(f"cmd={cmd}")
                self.run(binary=cmd, output_name=out_pp)
                lineno = self.pygrep_lineno(
                    file=out_pp, keyword="Suggested cut-off pseudo"
                )
                suggested_cutoff = float(
                    self.pygetline(
                        filename=out_pp, lineno=lineno, clearcache=True
                    ).split()[4]
                )
                new_cutoff.append(suggested_cutoff)
                logger.info(
                    f"suggested_cutoff for ECP index = {i} is rc = {suggested_cutoff} Bohr"
                )

                # clean files
                os.remove(file)
                os.remove(out_pp)
                os.remove("pseudo.plot")

            self.cutoff = new_cutoff

        finally:
            shutil.rmtree(tmp_rand_dir)
            os.chdir(current_dir)

    @classmethod
    def parse_pseudopotential_from_gamess_format_texts(cls, texts):
        max_ang_mom_plus_1 = []
        z_core = []
        cutoff = []
        nucleus_index = []
        element_list = []
        ang_mom = []
        exponent = []
        coefficient = []
        power = []

        for nuc_i, text in enumerate(texts):
            if text is None:
                continue

            pseudo_potential = (
                cls.Pseudopotential.parse_pseudopotential_from_gamess_format_text(text)
            )

            # storing
            max_ang_mom_plus_1.append(pseudo_potential.max_ang_mom_plus_1)
            z_core.append(pseudo_potential.z_core)
            cutoff.append(pseudo_potential.cutoff)
            nucleus_index += [nuc_i] * pseudo_potential.ecp_num
            element_list.append(pseudo_potential.element)

            logger.debug(pseudo_potential.max_ang_mom_plus_1)

            ang_mom += pseudo_potential.ang_mom
            exponent += pseudo_potential.exponent
            coefficient += pseudo_potential.coefficient
            power += pseudo_potential.power

        return cls(
            max_ang_mom_plus_1=max_ang_mom_plus_1,
            z_core=z_core,
            cutoff=cutoff,
            nucleus_index=nucleus_index,
            element_list=element_list,
            ang_mom=ang_mom,
            exponent=exponent,
            coefficient=coefficient,
            power=power,
        )

    @classmethod
    def parse_pseudopotential_from_gamess_format_files(cls, files):

        texts = []
        for file in files:
            if file is None:
                texts.append(None)
            else:
                with open(file, "r") as f:
                    text = f.readlines()
                texts.append("".join(text))

        return cls.parse_pseudopotential_from_gamess_format_texts(texts=texts)

    class Pseudopotential:
        def __init__(
            self,
            max_ang_mom_plus_1=0,
            z_core=0,
            element=None,
            cutoff=0.0,
            ang_mom=[],
            exponent=[],
            coefficient=[],
            power=[],
        ):

            self.element = element
            self.max_ang_mom_plus_1 = max_ang_mom_plus_1
            self.z_core = z_core
            self.cutoff = cutoff
            self.ang_mom = ang_mom
            self.exponent = exponent
            self.coefficient = coefficient
            self.power = power

            logger.debug(self.max_ang_mom_plus_1)
            logger.debug(self.z_core)
            logger.debug(self.cutoff)
            logger.debug(self.ang_mom)
            logger.debug(self.exponent)
            logger.debug(self.coefficient)
            logger.debug(self.power)

            # assertion!!
            assert len(self.ang_mom) == len(self.exponent)
            assert len(self.ang_mom) == len(self.coefficient)
            assert len(self.ang_mom) == len(self.power)

        @property
        def ecp_num(self):
            return len(self.ang_mom)

        @classmethod
        def parse_pseudopotential_from_gamess_format_file(cls, file):
            with open(file, "r") as f:
                text = f.readlines()
            text = "".join(text)
            return cls.parse_pseudopotential_from_gamess_format_text(text=text)

        @classmethod
        def parse_pseudopotential_from_gamess_format_text(cls, text):
            """ "
            http://myweb.liu.edu/~nmatsuna/gamess/input/ECP.html
            -card 1-    PNAME, PTYPE, IZCORE, LMAX+1
            IZCORE is the number of core electrons to be removed.
            Obviously IZCORE must be an even number, or in other
            words, all core orbitals being removed must be
            completely occupied.
            LMAX+1 is the one higher than the maximum angular momentum
            occupied in the core orbitals being removed:
            to remove s,p,d,f core orbitals  (LMAX=0,1,2,3)
            we use p,d,f,g core potentials (LMAX+1=1,2,3,4).
            LMAX+1 is not permitted to exceed 4.
            -card 2-    NGPOT
            NGPOT is the number of Gaussians in this part of the
            fit to the local effective potential.
            -card 3-    CLP,NLP,ZLP   (repeat this card NGPOT times)
            CLP is the coefficient of this Gaussian in the potential.
            NLP is the power of r for this Gaussian, 0 <= NLP <= 2.
            ZLP is the exponent of this Gaussian.
            The potential U(LMAX+1) is given first, followed by
            difference potentials U(L)-U(LMAX+1), for L=0,LMAX.

            """
            lines = text.split("\n")
            c = [
                line for line in lines if re.match(r"^\s*\S+.*", line)
            ]  # remove blank lines

            ang_mom = []
            exponent_list = []
            coefficient_list = []
            power_list = []
            element, _, z_core, lmax_plus_1 = c[0].split()
            z_core = int(z_core)
            max_ang_mom_plus_1 = int(lmax_plus_1)
            lmax_plus_1_num_components = int(c[1].split()[0])
            cutoff = 0.0
            # the first component is the potential U(LMAX+1) is given first,
            # followed by difference potentials U(L)-U(LMAX+1), for L=0,LMAX.

            stride = 2
            for _ in range(lmax_plus_1_num_components):
                logger.debug(c[stride].split())
                coefficient, power, exponent = c[stride].split()
                coefficient = float(coefficient)
                power = int(power)
                exponent = float(exponent)
                ang_mom.append(max_ang_mom_plus_1)
                exponent_list.append(exponent)
                coefficient_list.append(coefficient)
                power_list.append(power - 2)
                stride += 1

            for l in range(max_ang_mom_plus_1):
                num_component = int(c[stride].split()[0])
                stride += 1
                for _ in range(num_component):
                    coefficient, power, exponent = c[stride].split()
                    exponent = float(exponent)
                    power = int(power)
                    coefficient = float(coefficient)
                    ang_mom.append(l)
                    exponent_list.append(exponent)
                    coefficient_list.append(coefficient)
                    power_list.append(power - 2)
                    stride += 1

            logger.debug(ang_mom)
            logger.debug(exponent_list)
            logger.debug(coefficient_list)
            logger.debug(power_list)

            return cls(
                max_ang_mom_plus_1=max_ang_mom_plus_1,
                z_core=z_core,
                element=element,
                cutoff=cutoff,
                ang_mom=ang_mom,
                exponent=exponent_list,
                coefficient=coefficient_list,
                power=power_list,
            )


if __name__ == "__main__":
    logger = getLogger("pseudo-downloader")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("INFO")
    handler_format = Formatter("%(name)s - %(levelname)s - %(lineno)d - %(message)s")
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # set argparse (helper)
    import argparse

    parser = argparse.ArgumentParser(
        description="Download ccECP and BFD ECPs and convert them to the TurboRVB format. How to use: python ecp_downloader.py"
    )
    args = parser.parse_args()

    # ccECP
    ccecp = ccECP()
    ccecp.all_to_file()

    # BFD
    bfd = BFD()
    bfd.all_to_file()
