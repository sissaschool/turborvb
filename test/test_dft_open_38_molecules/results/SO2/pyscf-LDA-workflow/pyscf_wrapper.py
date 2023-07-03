#!/usr/bin/env python
# coding: utf-8

# # pySCF -> TREX-IO (Water molecule)

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)

# load python packages
import numpy as np
import os, sys
import pandas as pd
import scipy
import numpy
import itertools

# load ASE modules
from ase.io import read
from ase.units import Bohr
# load pyscf packages
from pyscf import gto, scf, mp, tools
from pyscf.pbc import gto as gto_pbc
from pyscf.pbc import dft as pbcdft
from pyscf.pbc import scf as pbcscf
from pyscf.scf.chkfile import dump_scf

class Pyscf_wrapper():

    def __init__(self,
            structure_file,
            chkfile="pyscf.chk",
    ):

        # read structural information from the file using ASE.
        logger.info(f"structure file = {structure_file}")

        self.structure_file=structure_file
        self.atom=read(self.structure_file)
        self.pbc_flag=self.atom.get_cell().any()
        if self.pbc_flag:
            logger.info("Periodic System!!")
            #raise NotImplementedError("PBC implementation is in progress.")
        else:
            logger.info("Molecule")
        self.chkfile = chkfile

    def get_mo_index(self, mo_occ_thr = 1.0e-3):
        mf=scf.chkfile.load(self.chkfile, "scf")
        mo_occ=mf['mo_occ']
        mo_index=[]
        logger.debug(f"mo_occ_thr={mo_occ_thr}")
        for i in range(len(mo_occ)):
            mo_index.append(i)
            if mo_occ[i] >= mo_occ_thr and mo_occ[i + 1] < mo_occ_thr:
                logger.info(f"mo_occ < {mo_occ_thr} is 1-{i + 1}")
                return mo_index
        return mo_index

    def run_pyscf(self,
                      omp_num_threads=1,
                      charge=0,
                      spin=0,
                      basis="ccecp-ccpvtz",
                      ecp='ccecp',
                      scf_method="HF",  # HF or DFT
                      dft_xc="LDA_X,LDA_C_PZ",
                      solver_newton=False,
                      MP2_flag=False,
                      CCSD_flag=False,
                      pyscf_output="out_pyscf",
                      max_cycle = 200,
                      symmetry=False,
                      twist_average=False,
                      exp_to_discard=0.10,
                      kpt=[0.0,0.0,0.0], # scaled_kpts!! i.e., crystal coord.
                      kpt_grid=[1,1,1],
                      smearing_method="fermi",
                      smearing_sigma=0.00 # Ha
                      ):

        if self.pbc_flag:
            cell = gto_pbc.M()
            cell.from_ase(self.atom)

            cell.verbose = 5
            cell.output = pyscf_output
            cell.charge = charge
            cell.spin = spin
            cell.symmetry = False
            a = cell.a
            cell.a = np.array([a[0], a[1], a[2]])  # otherwise, we cannot dump a
            # basis set
            cell.basis = basis
            if exp_to_discard != 0.0:
                cell.exp_to_discard = exp_to_discard
                logger.info(f"exp_to_discard={exp_to_discard}")

            # define ecp or pseudo
            if ecp is not None: cell.ecp = ecp

            cell.build(cart=False)

            # electron num
            electron_up_num, electron_dn_num = cell.nelec

            # calc type setting
            logger.info(f"scf_method = {scf_method}")  # HF/DFT

            # check the system
            if electron_up_num == 0 or electron_dn_num == 0:
                one_up_electron_flag = True
            else:
                one_up_electron_flag = False

            if scf_method == "HF":
                # RHF calculation
                if cell.spin == 0:
                    logger.info("HF kernel=RHF")
                    if twist_average:
                        logger.info("twist_average=True")
                        kpt_grid = cell.make_kpts(kpt_grid)
                        mf = pbcscf.khf.KRHF(cell, kpt_grid)

                    else:
                        logger.info("twist_average=False")
                        logger.info(f"kpt={kpt}")
                        logger.info(f'abs kpt = {cell.get_abs_kpts(scaled_kpts=[kpt])[0]}')
                        mf = pbcscf.hf.RHF(cell, kpt=cell.get_abs_kpts(scaled_kpts=[kpt])[0])

                else:
                    logger.info("HF kernel=ROHF")
                    if twist_average:
                        logger.info("twist_average=True")
                        kpt_grid = cell.make_kpts(kpt_grid)
                        mf = pbcscf.krohf.KROHF(cell, kpt_grid)

                    else:
                        logger.info("twist_average=False")
                        logger.info(f"kpt={kpt}")
                        logger.info(f'abs kpt = {cell.get_abs_kpts(scaled_kpts=[kpt])[0]}')
                        mf = pbcscf.rohf.ROHF(cell, kpt=cell.get_abs_kpts(scaled_kpts=[kpt])[0])

                # smearing
                if smearing_sigma != 0.0:
                    logger.info("Smearing is added!")
                    logger.info(f"smearing_sigma={smearing_sigma}, smearing_method={smearing_method}")
                    mf = scf.addons.smearing_(mf, sigma=smearing_sigma, method=smearing_method)

                #newton solver
                if solver_newton:
                    logger.info("solver = newton")
                    mf = mf.newton()

                mf.chkfile = self.chkfile

            elif scf_method == "DFT":
                # DFT calculation
                if cell.spin == 0:
                    logger.info("DFT kernel=RKS")
                    if twist_average:
                        logger.info("twist_average=True")
                        kpt_grid = cell.make_kpts(kpt_grid)
                        mf = pbcdft.krks.KRKS(cell, kpt_grid)
                    else:
                        logger.info("twist_average=False")
                        logger.info(f"kpt={kpt}")
                        logger.info(f'abs kpt = {cell.get_abs_kpts(scaled_kpts=[kpt])[0]}')
                        mf = pbcdft.rks.RKS(cell, kpt=cell.get_abs_kpts(scaled_kpts=[kpt])[0])

                else:
                    logger.info("DFT kernel=ROKS")
                    if twist_average:
                        logger.info("twist_average=True")
                        kpt_grid = cell.make_kpts(kpt_grid)
                        mf = pbcdft.kroks.KROKS(cell, kpt_grid)
                    else:
                        logger.info("twist_average=False")
                        logger.info(f"kpt={kpt}")
                        logger.info(f'abs kpt = {cell.get_abs_kpts(scaled_kpts=[kpt])[0]}')
                        mf = pbcdft.roks.ROKS(cell, kpt=cell.get_abs_kpts(scaled_kpts=[kpt])[0])

                # smearing
                if smearing_sigma != 0.0:
                    logger.info("Smearing is added!")
                    logger.info(f"smearing_sigma={smearing_sigma}, smearing_method={smearing_method}")
                    mf = scf.addons.smearing_(mf, sigma=smearing_sigma, method=smearing_method)

                #newton solver
                if solver_newton:
                    logger.info("solver = newton")
                    mf = mf.newton()

                mf.chkfile = self.chkfile
                mf.xc = dft_xc
            else:
                raise NotImplementedError

            # HF/DFT energy
            total_energy = mf.kernel()
            logger.info(f"Total HF/DFT energy = {total_energy}")
            logger.info("HF/DFT calculation is done.")

            # MP2 part
            if MP2_flag:
                logger.error("MP2 is not implemented for PBC cases")
                raise NotImplementedError

            if CCSD_flag:
                logger.error("CCSD is not implemented for PBC cases")
                raise NotImplementedError

        else: # open system
            chemical_symbols=self.atom.get_chemical_symbols()
            positions=self.atom.get_positions()
            mol_string=""
            for chemical_symbol, position in zip(chemical_symbols, positions):
                mol_string+="{:s} {:.10f} {:.10f} {:.10f} \n".format(chemical_symbol, position[0], position[1], position[2])
            # build a molecule
            mol = gto.Mole()
            mol.atom = mol_string
            mol.verbose = 5
            mol.output = pyscf_output
            mol.unit = 'A' # angstrom
            mol.charge = charge
            mol.spin = spin
            mol.symmetry = symmetry
            #print(mol_string)

            # basis set
            mol.basis = basis

            # define ecp
            if ecp is not None: mol.ecp = ecp

            # molecular build
            mol.build(cart=False)  # cart = False => use spherical basis!!

            #electron num
            electron_up_num, electron_dn_num = mol.nelec

            assert mol.spin == spin
            assert mol.charge == charge
            assert mol.cart == False # this converter supports only the spherical basis.

            # calc type setting
            logger.info(f"scf_method = {scf_method}")  # HF/DFT

            # check the system
            if electron_up_num ==0 or electron_dn_num==0:
                one_up_electron_flag = True
            else:
                one_up_electron_flag = False

            if scf_method == "HF":
                # HF calculation
                if mol.spin == 0:
                    logger.info("HF kernel=RHF")
                    mf = scf.RHF(mol)
                    mf.max_cycle = max_cycle
                else:
                    logger.info("HF kernel=ROHF")
                    #ROHF
                    mf = scf.ROHF(mol)
                    mf.max_cycle = max_cycle

                # smearing
                if smearing_sigma != 0.0:
                    logger.info("Smearing is added!")
                    logger.info(f"smearing_sigma={smearing_sigma}, smearing_method={smearing_method}")
                    mf = scf.addons.smearing_(mf, sigma=smearing_sigma, method=smearing_method)

                # newton solver
                if solver_newton:
                    logger.info("solver = newton")
                    mf = mf.newton()

                mf.chkfile = self.chkfile

            elif scf_method == "DFT":
                # DFT calculation
                if mol.spin == 0:
                    logger.info("DFT kernel=RKS")
                    mf = scf.KS(mol).density_fit()
                    mf.max_cycle = max_cycle

                else:
                    logger.info("DFT kernel=ROKS")
                    #ROKS
                    mf = scf.ROKS(mol)
                    mf.max_cycle = max_cycle

                # smearing
                if smearing_sigma != 0.0:
                    logger.info("Smearing is added!")
                    logger.info(f"smearing_sigma={smearing_sigma}, smearing_method={smearing_method}")
                    mf = scf.addons.smearing_(mf, sigma=smearing_sigma, method=smearing_method)

                # newton solver
                if solver_newton:
                    logger.info("solver = newton")
                    mf = mf.newton()

                mf.chkfile = self.chkfile
                mf.xc = dft_xc

            else:
                raise NotImplementedError

            # Molecular Orbitals and occupations
            logger.info("MOs-HF/DFT")
            #print(mf.mo_coeff)  # HF/DFT coeff
            #print(mf.mo_occ)  # HF/DFT occ
            #print(mf.mo_energy)  # HF/DFT energy
            # Notice!! The mo_i-th molecular orbital is NOT mo_coeff[mo_i], but mo_coeff[:,mo_i] !!

            # HF/DFT energy
            total_energy = mf.kernel()
            logger.info(f"Total HF/DFT energy = {total_energy}")
            logger.info("HF/DFT calculation is done.")

            # MP2 part

            if MP2_flag:
                # MP2 calculation
                logger.info("MP2_flag is True.")
                pt = mp.MP2(mf)
                mp2_E, t2 = pt.kernel(mf.mo_energy, mf.mo_coeff)

                # construct the one body density matrix
                rdm1 = pt.make_rdm1()

                # diagonalize to yield the NOs and NO occupation #s
                no_occ, no = scipy.linalg.eigh(rdm1)
                no_occ = no_occ[::-1]
                no = no[:, ::-1]

                # atomic orbital representation of the NO
                no_coeff = mf.mo_coeff.dot(no)

                # Molecular orbital and occupations
                # overwrite the HF/DFT ones!!!
                logger.warning("The HF/DFT MOs and occ. are overwritten!!")
                mf.mo_coeff = no_coeff # natural coeff
                mf.mo_occ = no_occ # natural orbital
                mf.mo_energy = [0.0] * len(mf.mo_energy) # orbital energy is not defined. So, they are set to 0.0

                logger.debug("MOs-MP2")
                logger.debug(mf.mo_coeff)  # HF/DFT coeff
                logger.debug(mf.mo_occ)  # HF/DFT occ
                logger.debug(mf.mo_energy)  # HF/DFT energy

                total_energy+=mp2_E
                logger.info(f"MP2 correlated energy={mp2_E}")
                logger.info(f"Total MP2 energy = {total_energy}")
                logger.info("MP2 calculation is done.")

            elif CCSD_flag:
                # CCSD calculation
                logger.info("CCSD_flag is True.")
                mycc = cc.CCSD(mf)
                mycc.kernel()

                # construct the one body density matrix
                rdm1 = mycc.make_rdm1()

                # diagonalize to yield the NOs and NO occupation #s
                no_occ, no = scipy.linalg.eigh(rdm1)
                no_occ = no_occ[::-1]
                no = no[:, ::-1]

                # atomic orbital representation of the NO
                no_coeff = mf.mo_coeff.dot(no)

                # Molecular orbital and occupations
                # overwrite the HF/DFT ones!!!
                logger.warning("The HF/DFT MOs and occ. are overwritten!!")
                mf.mo_coeff = no_coeff  # natural coeff
                mf.mo_occ = no_occ  # natural orbital
                mf.mo_energy = [0.0] * len(mf.mo_energy)  # orbital energy is not defined. So, they are set to 0.0

                logger.debug("MOs-CCSD")
                logger.debug(mf.mo_coeff)  # HF/DFT coeff
                logger.debug(mf.mo_occ)  # HF/DFT occ
                logger.debug(mf.mo_energy)  # HF/DFT energy

                logger.info(f"Total CCSD energy = {mycc.e_tot}")
                logger.info("CCSD calculation is done.")

            # this is a special treatment for one-valence electron case!!
            # because PYSCF seems do UKS calc. for one-valence electron case
            # here the alpha ones are taken.
            # only DFT case!!
            if one_up_electron_flag and scf_method == "DFT":
                mf.mo_coeff=mf.mo_coeff[0]
                mf.mo_occ = mf.mo_occ[0]
                mf.mo_energy = mf.mo_energy[0]
                logger.debug(mf.mo_coeff)
                logger.debug(mf.mo_occ)
                logger.debug(mf.mo_energy)

            # dump checkfile explicitly!!
            dump_scf(mol=mol,
                     chkfile=self.chkfile,
                     e_tot=total_energy,
                     mo_energy=mf.mo_energy,
                     mo_coeff=mf.mo_coeff,
                     mo_occ=mf.mo_occ,
                     overwrite_mol=True
                     )

        logger.info("PySCF calculation is done.")

        logger.debug(mf.mo_coeff)
        logger.debug(mf.mo_occ)
        logger.debug(mf.mo_energy)

    @property
    def pyscf_total_energy(self):
        mf = scf.chkfile.load(self.chkfile, "scf")
        return mf['e_tot']

def cli():
    import argparse
    from logging import config, getLogger, StreamHandler, Formatter

    log_level = "INFO"
    logger = getLogger("pyscf-trexio")
    logger.setLevel(log_level)
    stream_handler = StreamHandler()
    stream_handler.setLevel(log_level)
    # handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    handler_format = Formatter('%(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # define the parser
    parser = argparse.ArgumentParser(epilog='pyscf calculation',
                                     usage='python pyscf_wrapper.py structure, options...',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--structure_file', help=f'structure_file', type=str, required=True)
    parser.add_argument('-c', '--chkfile', help=f'chkfile', type=str, default="pyscf.chk")
    parser.add_argument('-o', '--pyscf_output', help=f'pyscf_output', type=str, default="out_pyscf")

    parser.add_argument('--charge', help=f'charge', type=int, default=0)
    parser.add_argument('--spin', help=f'spin', type=int, default=0)
    parser.add_argument('--basis', help=f'basis', type=str, default='ccecp-pvtz')
    parser.add_argument('--ecp', help=f'ecp', type=str, default='ccecp')
    parser.add_argument('--scf_method', help=f'scf_method (HF or DFT)', type=str, default="DFT", choices=["DFT", "HF"])
    parser.add_argument('--dft_xc', help=f'dft_xc', type=str, default="LDA_X,LDA_C_PZ")
    parser.add_argument('--max_cycle', help=f'max_cycle', type=int, default=200)

    parser.add_argument('--solver_newton', help=f'solver_newton', action="store_true", default=False)
    parser.add_argument('--MP2_flag', help=f'MP2_flag', action="store_true", default=False)
    parser.add_argument('--CCSD_flag', help=f'CCSD_flag', action="store_true", default=False)

    parser.add_argument('--symmetry', help=f'symmetry', action="store_true", default=False)
    parser.add_argument('--twist_average', help=f'twist_average', action="store_true", default=False)

    parser.add_argument('--exp_to_discard', help=f'exp_to_discard', type=float, default=0.10)
    parser.add_argument('--kpt', help=f'kpt', type=float, nargs = 3, default = [0.0, 0.0, 0.0])
    parser.add_argument('--kpt_grid', help=f'kpt_grid', type=int, nargs=3, default=[1, 1, 1])


    # parse the input values
    args = parser.parse_args()
    parsed_parameter_dict = vars(args)

    pyscf_wrapper=Pyscf_wrapper(
        structure_file=args.structure_file,
        chkfile=args.chkfile
    )

    pyscf_wrapper.run_pyscf(
    charge = args.charge,
    spin = args.spin,
    basis = args.basis,
    ecp = args.ecp,
    scf_method = args.scf_method,  # HF or DFT
    dft_xc = args.dft_xc,
    solver_newton = args.solver_newton,
    MP2_flag = args.MP2_flag,
    CCSD_flag = args.CCSD_flag,
    pyscf_output = args.pyscf_output,
    max_cycle = args.max_cycle,
    symmetry = args.symmetry,
    twist_average = args.twist_average,
    exp_to_discard = args.exp_to_discard,
    kpt = args.kpt,
    kpt_grid = args.kpt_grid
    )


if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    cli()