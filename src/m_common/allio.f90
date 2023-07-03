! Copyright (C) 2022 TurboRVB group
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.

module allio
    use constants
    use cell
    use Ewald
    use types
    use kpoints_mod
    use io_m, only: lchlen
    ! by E. Coccia (22/11/10)
    use extpot, only: ext_pot, link_atom, mm_restr, write_rwalk
    ! by E. Coccia (23/12/10)
    use van_der_waals, only: vdw
    ! by E. Coccia (9/5/11)
    use link_atoms
    ! for sub communicator by Y. Luo (17/10/14)
    use sub_comm
    ! for mpiio by Y. Luo (24/2/15)
    use mpiio, only: file_obj
    use dielectric
    implicit none

    integer i_main, nw, max_sparse_choice, iseed, ngenc, ngg, ngn, ngs, d, nmax_ion &
   &, ngen, nscra, Lz, iout, indvic, nbra, ng, iopt, nmatr, nshell_det           &
   &, iscramax, iscrapip, iscraipsip, iscramax_c, testderiv                &
   &, nwr, iread, ireadr, nwm, indt, indtupt, nbrar, npm, np, npr, npbra, npmn    &
   &, npbrar, np3, nmp, nmpr, irank, nfatr, indcor, indt4, indt4j      &
   &, npp, nfat, itry, ncore, nbram, dimjas, iesupmax_c                      &
   &, nmat, npsov, npsovr, np3r, ncg, ncgdim, ncgread                        &
   &, inext, ibinit, iend, iendr, irec, writescratch, countscra              &
   &, nelm, nelup, neldo, nelup_mat, nel, nel_mat, indspin, ifz, jmax          &
   &, nwnel, nwdw, nwrep, info, nwfree, nrep, nws                            &
   &, iesfree, iesinv, ieskin, iesking, nwkin, iessw, nwsw, indc, indc0, iesup  &
   &, nwup, iesdr, iesdrr, iesfreer, iesswr, ieskinr, ieskinrp, iesupr, nmats  &
   &, iseedr, npf, iesupind, iesmind, typedyncell, ieskint, nwking, iesswr_eagp &
   &, nwnp, np3m, icount, ireadmin, ireadminr, iessw0, iesup_read            &
   &, nelnw, nelnw0, jold, niesd, ngentry, numcost                          &
   &, Ltab, Lztab, Ltabb, indold, indnn, js, jtype2                   &
   &, nelorb, nshell, nshelldo, nrnel, nelkel, dimee, dimei, iesconv          &
   &, nelnion, nion, indksij, k_ion, i_ion, nelorbj                         &
   &, indtj_nw, indtabj_nw, indtabbj_nw, indkj_nw, indksj_nw, indksij_nw    &
   &, indwwj_nw, indwupj_nw, indwdoj_nw, indaupj_nw, j_nw, j_nws, irej       &
   &, nshellj, indwwjj, indwwjj_nw, nel2wtj, indbar, indjbar, indbar_nw      &
   &, indjbar_nw, nelorbp, nelorbpp, indjbarsz, indjbarsz_nw        &
   &, npar3body, iesswind, nelcol, nelcol_c, inddsw, indn, iflag, nelcolh     &
   &, irankdet, pippo, cpippo, indpippo, dimpippo, ccpippo                 &
   &, nnozeromax, skipforce, pippoc, ndone, lbin, nbinr, nbin        &
   &, nbinread, iesfreesz, ixr, iyr, nbinmax, stepcg, nbindim                &
   &, ieskindim, npdim, ndimpdim, lwork, indfix, repf, kinf, idyn      &
   &, iskipdyn, iskipdynr, idynu, icost, indmax, indj, kp_ion, indfree, npar, initpar     &
   &, nparsw, initparsw, indsw, indinv, endinv, nparinv, initparinv, iflagk   &
   &, maxall, molecular, molecularj, nshell_max, nshellj_max, maxiesup      &
   &, ieskingdim, maxnpar, developer, molopt, nparshellmax, optbra         &
   &, unreliable, iflagerrall, ieskinold, nel2wtfn, indwwjfn     &
   &, nel2barfn, indbarfn, indbarfn_nw, indwwjfn_nw, powermin, npower&
   &, powerminsz, npowersz, ieskinr_pos                           &
   &, indtjr, Lztabr, indtjr_nw, commopt_mpi, commcolopt_mpi, rankopt, nprocopt&
   &, row_comm, col_comm, row_id, col_id, commcov_mpi, nproccov&
   &, commsr_mpi, nprocsr, ranksr, mcol, commrep_mpi, nprocrep, rankrep&
   &, nrep_bead, mcol_rep, rankcolrep, commcolrep_mpi, nproccolrep, commcolsr_mpi&
   &, nproc_diag, nprocu, ref_atom, epscuttyper, ndim_detmat, delay_changeparr&
   &, iesdelay, min_block, max_ortho, maxiter_changeparr, prep&
   &, comm_col, comm_raw, rankraw, rankcol, nbra_cyrus, dim_cyrus, nw_max&
   &, npar_eagp, dim_upwf, dim_ratiovar, dim_uptabtot&
   &, nnozero_eagp, max_target, max_targetsr, nbra_cyrus_read

    type(mpi_sub_comm) :: sub_comm_diag

    integer occ, nnozero, nnozeroj, occtot, kl, freqcheck, nmaxder&
   &, occj, occtotj, indp, iflagnorm, nnozerojr, novar, parcutg, equil_steps&
   &, firstmol, nmolfn, yesfast, iesup_atom, dimtranspip, lastmol, nelorb_at
    integer nel2up, nel2upt, nel2dot, nel2wt, nel2bar, indtot              &
   &, nel2jbar, nel2jbarsz, iesuptouched, vjutouched                      &
   &, indtj, indtabj, indtabbj, indkj, indwupj, indwdoj, indaupj, iskip       &
   &, indwwj, iesupskip, indksj, nmatb, nmatbr                    &
   &, nshellr, nshelljr, nelorbh, nelorbjh, nelorbjh2, irstart       &
   &, nx, ny, nz, nbufd, nmol, nmolmin, nmolmax, nel3&
   &, ndiff, ndiffdim, indndiff, ifreqdump, iimax, iijmax&
   &, nmolmaxw, typereg, niont, niong, epscuttype, ncg_adr, kp0, nbead, neldomax
    real*8 pot_aasunel, gamma, etry, etryr, rata, time, wbra, lambda, eps_umrigar&
   &, ener, ratio, wbra1, wbra2                                           &
   &, ratior(2), ratiodet, sumdiff, timep, timepp, timescra, ratiorn(2)            &
   &, psign, wbran, wcort, timemc, timeopt, timeinit                        &
   &, veff, veffright, beta, temp, kappar                                  &
   &, cost, cost1, cost2, parr, econtnew, jacold                    &
   &, parcutr, epst, epsi, psiav, psisav, pdiag, ttry, psirav, countav, countt  &
   &, avreweight, tleft, tbra, tbrar, counttot, countreg, countcut           &
   &, parcut, costpassed, ttry1, ttry2, pdiag1, pdiag2                      &
   &, wcorwt, epsdgel, eps_dyn5, epstion, epsdgm, memtot                    &
   &, wtotf(2), tpar, tparf, epstl, epstlu, rsignr                          &
   &, parbest, wtotf1, wtotf2                                            &
   &, rweight, hopfraction, ration(2), nontr                                 &
   &, tstep, bcost, ccost, minz, maxz, minzj, maxzj, ukwald                   &
   &, f, fb, vpotint, diffint, costexp, theta_reg                         &
   &, epscut, epstlrat, enermin, varmin, spsi, ratioreg, epscutu     &
   &, fmax, parcutmin, parcutpar, parcute, nacc, naccpseudo         &
   &, npow, normcorr                           &
   &, costwn, signflip, acclarge, weightall, tmes, rsr         &
   &, rmax, rmaxj, rmaxinv, epscutur, epstlur, time_main, time_ratiovar, time_uptabtot&
   &, timewf, timepip, ax, ay, az, zmax, tstepfn, tion, tcell, weight_loc, power &
   &, Klrdmc, epscutdmc, epstldmc, time_meas, time_branch, cutreg, cutweight&
   &, powerwarp, epsrem_contr, tolcg, smoothcut  &
   &, pressclass, dcellclass(3), scalepulay, shift, alat2v, maxtime, smearing&
   &, scale_mass, scale_one_body, scaleeloc, minjonetwobody&
   &, timings(11), timingsb(11), parr_min, parr_max, epsvar, tave_cyrus&
   &, tcount_cyrus, zmin, l0_kousuke, core_pseudo, count_zerowf, count_allwf
    ! added by Andrea Tirelli
    logical tpar_increased, stop_increasing_tpar, yes_cutweight, yes_scemama&
   &, yes_scemama_open, yes_sparse, yes_sparse_choose, yes_dgelscut
    integer counter_very_unstable_tpar, tpar_unstble_stop, &
        inc_tpar_frequency, counter_unstable_tpar, counter_unstable_energy, counter_unstable_err, len_tpar_stable_list
    real(8) min_running_ave, min_running_std, min_running_ave_energy, min_running_std_energy, tpar_max, &
        divide_tpar, multiply_tpar, n_sigmas_tpar
    !end added Andrea Tirelli

    integer*8 :: handle
    integer*4 ldworkspace, lzworkspace, dev_Info(1)
    real*8, allocatable, dimension(:) :: dev_dgetrf_workspace
    complex*16, allocatable, dimension(:) :: dev_zgetrf_workspace
    real*8, allocatable, dimension(:, :) :: dev_dgetri_workspace
    complex*16, allocatable, dimension(:, :) :: dev_zgetri_workspace

    real*8 deriv1, deriv2, deriv3, pulay1, pulay2, pulay3, scalecell(3)&
   &, epsder, scale_grad, norm_corr

    integer icdiff, itest, itestr, jn                                  &
 &, nrest, iesm, isfix, nwfix, itestrr, itestr3, itestr4                   &
 &, iesd, iese, ieser, npsamax, nwm2, nprest, nindt, nwdim, nintpsa          &
 &, ix, iy, iboot, xj, itestrfn, npsar, indberry, true_wagner, nelsquare

!      integer, private :: nshelljmax,indref,indrefg

    real*8 wback, costw&
   &, enerc, jacnew, gradold, gradbarold, friction                         &
   &, delta0, delta0q, delta0k, dt, scalecov, maxdev_dyn, pressfixed&
   &, versoralat(3, 12), epsbas, lepsbas, weight_moroni ! cutoff for the periodic basis set with k-points
    integer ntry, nwinv, nweight, nweightr, nmore_force
    real*8 alat, alat2, sigma_new, sigma_true(2), costa                 &
&, ener_true(2), plat(3), dstep(3)
    integer ndim, ndimp, ndimj, ndimjp, ndims, ndimsp, ndimiesup, movedion

    real(8), dimension(:, :), allocatable :: rion, rion_fast, rionsav     &
   &, efenergy, gradpsi, gradpsibar, efpress                     &
   &, gradpsiold, gradpsibarold, angle, fk, reduce, reducel, rcne            &
   &, warpmat, projmat_c
    real(8), dimension(:, :, :, :), allocatable :: jastrowall_ee, queue_cyrus

    real(8), dimension(:, :, :), allocatable :: rmunew, sov, ef, efp, ivic, jastrowall_ei, rmunewb, agp
    integer, dimension(:, :, :), allocatable :: npip
    real(8), dimension(:, :), allocatable :: kel, keln, dists_kel          &
   &, velion, diffkin, rknew, rknewb

    real(8), dimension(:, :), allocatable :: iond_cart, vpotreg, kdyn, vpotsav_ee

    ! winv vector for complex DFT
    ! variables required for complex total energy
    complex(8) enerc_c, ener_true_c, ener_c, enermin_c

    real(8), dimension(:), allocatable :: duprold                     &
   &, vjurold                                                          &
   &, eig, vj, dup, dsw, zetar, zetar_fast, zetaq, zetamin, distmin            &
   &, dek, dekg, vju, ddw, ddwsz, vjur, dupr, winv, ainv, vpot, enertrue, diffuse &
   &, ainvup, ainvdo, winvup, winvdo, scale, scalpar, scalej                 &
   &, scalejsz, winvj, cnorm, psinew, tmu, winvbar, winvjbar, detmat, jasmat   &
   &, jasmatsz, cnorm_nw                               &
   &, psibar, ainvupb, err, force &
   &, dist, iond, rmu, r, winvs, winvsj, ainvs, dists, wint            &
   &, alphavar, berry_exp, wcorw, wintw                    &
   &, alphab, diagfn, psip, psip_reweight, econf, econfh, wconfn, econfion &
   &, factorsr, vcut, etot, wsto, zeta, tabpip, table, tabler           &
   &, diag, tcore, tcost, gradtot, gradtotbar, costz, costz3                 &
   &, rcarto, rcart, dx_old, dx_new, jasnew_ei                             &
   &, fkav, winvjbarsz, okav, skdiag, cov, cov_old                          &
   &, jasnew_ee, atom_number, atom_number_fast           &
   &, winvbarn, enerint, enerintw, wconfsav, wtot         &
   &, winvbarfn, winvfn, t_cyrus, first_moment, second_moment

    real maxoutput
    real(4), dimension(:), allocatable :: econfw, wconfw
    real(8), dimension(:), allocatable :: psilnw, bufscra, v_adr ! too much sensitive

    integer, dimension(:), allocatable :: ipsip, mult, kion             &
   &, nparam, multj, nparamj, kionj, ioccj, ioptorbj, nozero, nozerodet       &
   &, nozeroj, ioptorb, ioccup, ioccdo, naccm                      &
   &, jbra, jbraw, jbraj, jbrajsz, jbradet, jbradetn                        &
   &, jbrajn, jbraiesup, jbraiesm, whereiesup, whereiesm, vjutouch          &
   &, itouch, icore&
   &, jbraiesup_sav, jbraiesm_sav, nozerojder, ipip_adr, first_cyrus
    logical, dimension(:), allocatable :: sjbradet, slaterorb_read

    type(ion_comp), dimension(:), allocatable :: ion_table

    integer, allocatable :: indpar_tab(:), indorb_tab(:)&
         &, indshell_tab(:), indparj_tab(:), indorbj_tab(:), indshellj_tab(:)&
         &, adr_nion(:), ind_nion(:), adrj_nion(:), indj_nion(:), addr_occ(:)&
         &, type_atom(:), pointvj(:, :)

    integer npsa, lmax, istart, indteff, lzeff
    character(3) pseudoname
    character(60) pseudofile
    integer, dimension(:, :), allocatable :: nparpshell, jpseudo        &
         &, indtm
    integer, dimension(:), allocatable :: kindion, pshell
    real(8), dimension(:, :), allocatable :: parshell, legendre         &
         &, versor, prefactor, enert, agpn
    real(8), dimension(:), allocatable :: rcutoff, wpseudo, pseudolocal &
         &, wintpseudo
    integer npseudopar, npseudoparn, ion, nintpseudo
    real(8) coeff
    real(8) beta_learning
    ! variables for complex/quantum algorithms
    real(8), dimension(:), allocatable :: psiln, psidetln, kdyn_eig
    real(8), dimension(:), allocatable :: psisn
    complex(8), dimension(:), allocatable :: psip_c
    complex(8) ratio_c, ratior_c, ratiodet_c ! complex wave function ratio

    real(8) psioverpsi_new, psioverpsi_old

    real(8) cutoff_p ! cutoff for bump orbitals (abs(Lbox) = 2). Default = 9.d0

!  Variables used only in scalapack (__SCALAPACK) but allocated always for compatibility
!  (irrelevant memory overhead)
    integer :: me_blacs = 0 ! BLACS processor index starting from 0
    integer :: np_blacs = 1 ! BLACS number of processor
    integer :: world_cntx = 0 ! BLACS context of all processor
    integer :: ortho_cntx = 0 ! BLACS context for ortho_comm
    integer :: me_ortho(2) = 0 ! coordinates of the processors
    integer :: me_ortho1 = 0 ! task id for the ortho group
    integer :: np_ortho(2) = 1 ! size of the processor grid used in ortho
    integer :: np_ortho1 = 1 ! size of the ortho group
    integer :: ortho_comm = 0 ! communicator used for fast and memory saving ortho
    integer :: ortho_comm_id = 0 ! id of the ortho_comm
    integer :: leg_ortho = 1 ! the distance in the father communicator
    ! of two neighbour processors in ortho_comm

    integer rank, nproc, nprocr, ist, ien, id1, istm                          &
         &, in1, ierr, skipreshuff
    character(14) chara
    character(14) charaq
    character(60 + lchlen) wherescratch
    character(lchlen) errmsg

    integer nshell_c, occ_c, iesup_c, nnozero_c, nelorb_c, contraction, ikshift     &
         &, nshellj_c, occj_c, npar3body_c, nnozeroj_c, nelorbj_c, contractionj   &
         &, maxparam, maxioccmult, nelorbmax, indocc, maxshell      &
         &, maxparamj, maxshellj, nelorbmaxj, ll, mm, iesupr_2, iesupr_c           &
         &, npar3bodyr_c, occ_tmp, npar3body_fill, npar3body_2, nmolmat&
         &, nmolmatw
    !     parameter(maxioccmult=10)
    integer, dimension(:), allocatable :: mult_c, nparam_c, ioptorb_c   &
         &, kion_c, nozero_c, ioccup_c, occshell                                &
         &, iesuptrans, multranspip, iesuptransb                               &
         &, multj_c, nparamj_c, ioptorbj_c                                     &
         &, kionj_c, nozeroj_c, ioccj_c, occshellj                              &
         &, iesuptransj, multranspipj, iesuptransbj                            &
         &, kiontot, kiontotj, ioptorbja                                       &
         &, typeorb
    type(array_int), allocatable :: transpip(:), transpip_sav(:)
    type(array_int), allocatable :: transpipj(:), transpipj_sav(:)
    type(nkgrid), allocatable :: kgrid(:)
    type(nkgrid), allocatable :: kgrid_atom(:)

    integer, dimension(:, :), allocatable :: mu_touch
    integer, dimension(:, :), allocatable :: muj_touch, adrlambda
    real(8), dimension(:), allocatable :: dup_c, scale_c, detmat_c, projm&
   &, detmat_proj
    real(8), dimension(:), allocatable :: vju_c, scalej_c, scalejsz_c   &
   &, jasmat_c, jasmatsz_c, rpar
    real(8), dimension(:, :), allocatable :: jas_invariant, eagp_pfaff, eagp_pfaffb
    real(8), dimension(:, :), allocatable :: mu_tmp, mu_c
    real(8), dimension(:, :), allocatable :: muj_tmp, muj_c, mat_adr
    logical, dimension(:), allocatable :: orbcost, orbcostn, orbcostl&
   &, yescut, orbps, singdet, allowed_par
    logical, dimension(:, :), allocatable :: allowcontr, allowcontrj

    logical pseudologic, iessz, iesgros, iescost, fncont, iesbra, yeszj &
   &, molyes, moljyes, iesrandoma, iesrandoml, pseudorandom         &
   &, iespbc, fnloc, yesivic, yesnleft&
   &, rejweight, orthoyes, yeszagp, yesdetmat, symmagp, yesdft &
   &, yesfmubar, membig, membigcpu, membigr, onebodysz, yespress, iescostd, twobodyoff&
   &, printoverlap, yesdetmatc, oldscra, symiesup, iesdtwobodyoff, iesdonebodyoff &
   &, warp, yesbump, yeslbox, yespulay, stopdyn, allfit&
   &, yescutjas, yescutdet, defparcutg, detc_proj, yesread10&
   &, stepcg_recount, write_cov, gramyes, fixpar, symmetrize_agp, yesprimitive&
   &, yesQuantum, yesavsr, yesavcov, yesavopt, yeswrite10, yesperiodize, yesturboq& ! main options for the quantum algorithm
   &, yes_complex, yesupel, yessecond, yes_crystal, yes_crystalj, test_aad, eqcellab&
   &, eqcellac, eqcellbc, forcecomplex, oldscaling, ldynsecond, add_onebody2det&
   &, yes_hermite, allowed_averagek, yes_correct, yes_real, srcomplex, killcut&
   &, change_epscut, change_tstep, better_dmc, yesalfe, safelrdmc, yesrootc&
   &, addrognoso, changelambda, cleanrognoso, fixa, fixb, fixc, forcesymm&
   &, signalnoise& ! main options for complex algorithm
   &, real_contracted, gauge_fixing, yesmin_read&
   &, noopt_onebody, real_agp, softcusp, scalermax, yeswritebead, yes_hessc&
   &, no_sjbra, manyfort10, shift_origin, shiftx, shifty, shiftz&
   &, double_mesh, change_parr&
   &, default_epsdgel, read_molecul, hybyes, pfaffup, k6gen, noblocking, add_diff&
   &, lrdmc_der, lrdmc_nonodes, nosingledet, enforce_detailb, nowrite12&
   &, yes_fastbranch, flush_write, yes_adams, only_molecular, add_offmol, novec_loop1
    ! yesupel=.true.    --> means dealing with spin up electrons, now one can choose a different phase for up/down spin

    ! test_aad=.true. --> used only when you are using the program testadc.x for testing AAD derivatives, otherwise set to .false.

    ! yes_crystal=.true. --> use Crystal basis set defined as: \phi_k(r)=\sum_R \phi(r-R_a-R)*exp(ikR). Use the input variable
!   yes_crystalj=.true. --> use Crystal basis  also for the Jastrow
    !                        epsbas to set the cutoff on the sum over the direct lattice vectors.

    ! eqcellab,eqcellac,eqcellbc --> force the cell parameters (a,b),(a,c),(b,c),(a,b,c) respectively to be equal when using
    !                                          typedyncell>0 for cell relaxation.

    ! control io behavior, Y. Luo (15/2/15)
    character(len=80) :: disk_io
    integer :: io_level

    ! MPI-IO file objects
    type(file_obj) :: kelcont, quantcont, details_SP, details_DP

    ! added by K. Nakano for automatic gutta cavat lapidem (to adjust tpar)
    ! parameters read from an input file.
    logical :: change_tpar, use_stable_tpar
    ! internal variables (not read from an input file)
    integer :: inc_counter_tpar, dec_counter_tpar, tpar_buffer_len, len_shorter_buffer, times_tpar_decreased
    real(8) cut_sigma
    logical, dimension(:), allocatable :: tpar_buffer_filled
    real(8), dimension(:), allocatable :: energy_list, error_energy_list
    real(8), dimension(:), allocatable :: tpar_stable_list

! ************ EWALD SUMS *******************!
    integer kmax2, xi, n_body_on, yesmin, yesminr
    real(8), dimension(:, :), allocatable ::                           &
   &rmusin, rmucos
    real(8) selfsum, LBox, LBoxj, Linv, rs, ris(5)                              &
   &, derEVp, errEVp, derEVnopulay, errEVnopulay
    real(8), dimension(:), allocatable :: p_pulay, derEV
    integer add_pulay
! ************ EWALD SUMS *******************!
!       derEV  = derivative of the energy respect to Volume
!       errEV =  error of derivative of the energy respect to Volume

!************* CAFFAREL FORCES *************
! orbderiv = orbital derivatives respect to nuclei
! derpot   = derivatives of potential energy
! dercaf = laplacian of caffarel Q

!************* CAFFAREL FORCES *************
!       Definition default values
!        rank=0
!        nproc=1
!        nw=1

    namelist /simulation/ itestr4, iopt, ngen, nscra, nbra, iseedr, nw, kSq &
    &, kappar, freqcheck, membig, membigcpu, developer, yesfast, maxtime, nproc_diag, disk_io, ip_reshuff &
    &, compute_bands, double_mesh, min_block, max_target, max_targetsr, dielectric_ratio, dielectric_length &
    &, case_diel, neigh, novec_loop1, yes_sparse, yes_sparse_choose, max_sparse_choice

    namelist /pseudo/ nintpsa, npsamax, pseudorandom

    namelist /readio/ ncore, np3, np, iread, writescratch, wherescratch, unreliable, ifreqdump, nowrite12, flush_write

    namelist /vmc/ tstep, hopfraction, epscut, epstlrat, epscuttype, alat2v, shift, change_epscut, change_tstep &
    &, epsvar, theta_reg, true_wagner, cutweight, nbra_cyrus, typereg, npow

    namelist /dmclrdmc/ etry, npow, tbra, gamma, plat, alat2, alat        &
    &, tstepfn, Klrdmc, optbra, parcutg, novar, epscutdmc, typereg&
    &, epstldmc, rejweight, cutreg, cutweight, better_dmc, yesalfe&
    &, safelrdmc, changelambda, noblocking, add_diff, nbra_cyrus, lrdmc_der&
    &, lrdmc_nonodes, enforce_detailb, iesrandoma, zmin, yes_fastbranch&
    &, l0_kousuke, nw_max, true_wagner, weight_moroni

    namelist /optimization/ tpar, nfat, iboot, nweight, nmore_force, epsi, eps_dyn5 &
    &, epsdgel, kl, idyn, nbinr, npbra, ncg, minz, maxz, minzj, maxzj, parr, parcute&
    &, parcut, parcutmin, parcutpar, tion, tcell, molopt, epstion&
    &, onebodysz, twobodyoff, iesdtwobodyoff, iesdonebodyoff, tolcg, minjonetwobody&
    &, symiesup, yescutjas, yescutdet, fixpar, symmetrize_agp&
    &, yesquantum, nbead, yeswrite10, oldscaling, srcomplex&
    &, power, signalnoise, gauge_fixing, beta_learning&
    &, noopt_onebody, scalermax, yeswritebead, yesread10, change_parr&
    &, parr_max, parr_min, delay_changeparr, maxiter_changeparr, k6gen, max_ortho, prep &
    &, change_tpar, inc_tpar_frequency, use_stable_tpar &
    &, eps_umrigar, yes_adams, divide_tpar, multiply_tpar, tpar_buffer_len, tpar_max &
    &, cut_sigma, n_sigmas_tpar, len_tpar_stable_list, yes_dgelscut

    namelist /parameters/ ieser, iesinv, iesm, iesd, isfix, iesfree &
    &, iessw, iesup, ieskin, yespress, warp, powerwarp&
    &, add_pulay, yespulay, typedyncell, scalepulay, ext_pot, vdw, link_atom&
    &, mm_restr, write_rwalk, yes_correct&
    &, yesavopt, yesavsr, yesavcov, nrep_bead, yesperiodize, yes_kpoints, epsbas&
    &, yeszj, yeszagp, decoupled_run, scaleeloc, cutoff_p, fixa, fixb, fixc&
    &, real_contracted, real_agp, no_sjbra, pressfixed, read_molecul, epsder, yes_scemama, yes_scemama_open

    namelist /fitpar/ nparinv, initparinv, rmaxinv, npar, initpar, rmaxj &
    &, nparsw, initparsw, rmax, npower, powermin, npowersz, powerminsz, allfit

    namelist /dynamic/ temp, friction, delta0, delta0q, delta0k, scalecov&
    &, iskipdyn, maxdev_dyn, stepcg_recount, write_cov, normcorr &
    &, yesturboq, yessecond, smoothcut, killcut, scale_mass, eqcellab &
    &, eqcellac, eqcellbc, yesrootc, addrognoso, cleanrognoso

    namelist /unused/ rsignr, beta, testderiv

    namelist /molecul/ epsdgm, nx, ny, nz, nbufd, ax, ay, az, nmolmin, smearing&
    &, nmolmax, weight_loc, orthoyes, epsrem_contr, nmolmaxw, gramyes&
    &, add_onebody2det, shift_origin, shiftx, shifty, shiftz
    ! by E. Coccia (22/11(10): logical flag for the external potential
    ! by E. Coccia (3/2/11): logical flag for the vdw term of the external potential
    !namelist /pot_ext/  ext_pot, vdw, link_atom
    ! by E. Coccia (23/12/10): namelist for writing the density during the vmc run

    namelist /link/ calpha

contains

    ! reduce allio to speed up compilation
    ! read and write fort.10 subroutines are in fort10_io.f90
    ! read and write fort.11 subroutines are in fort11_io.f90
    ! read datas* control files by read_datas.f90
    ! read pseudo potential in read_pseudo.f90
    ! default allocation and deallocation in memOP.f90
    ! write output by writeoutput.f90

    subroutine scontract_genj(nelorbh, nelorb_c         &
         &, detmat, detmat_c, mu_c, psip)
        implicit none
        integer nelorbh, nelorb, nelorb_c, i
        real*8 detmat(2*nelorbh, 2*nelorbh), detmat_c(2*nelorb_c, 2*nelorb_c)&
             &, mu_c(nelorbh, *), psip(nelorbh, *)
!#ifdef __CASO
!    nprocu=nprocopt
!#else
!    nprocu=1
!#endif
        detmat = 0.d0
        call dgemm_my('N', 'N', nelorbh, nelorb_c, nelorb_c, 1.d0, mu_c, nelorbh  &
             &, detmat_c, 2*nelorb_c, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
        call dgemm_my('N', 'T', nelorbh, nelorbh, nelorb_c, 1.d0, psip, nelorbh   &
             &, mu_c, nelorbh, 0.d0, detmat, 2*nelorbh, nprocu, rankopt, commopt_mpi)
!   down-down
        call dgemm_my('N', 'N', nelorbh, nelorb_c, nelorb_c, 1.d0, mu_c, nelorbh  &
       &, detmat_c(nelorb_c + 1, nelorb_c + 1), 2*nelorb_c, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
        call dgemm_my('N', 'T', nelorbh, nelorbh, nelorb_c, 1.d0, psip, nelorbh   &
       &, mu_c, nelorbh, 0.d0, detmat(nelorbh + 1, nelorbh + 1), 2*nelorbh, nprocu, rankopt, commopt_mpi)

!   down-up
        call dgemm_my('N', 'N', nelorbh, nelorb_c, nelorb_c, 1.d0, mu_c, nelorbh  &
       &, detmat_c(nelorb_c + 1, 1), 2*nelorb_c, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
        call dgemm_my('N', 'T', nelorbh, nelorbh, nelorb_c, 1.d0, psip, nelorbh   &
       &, mu_c, nelorbh, 0.d0, detmat(nelorbh + 1, 1), 2*nelorbh, nprocu, rankopt, commopt_mpi)
!   up-down
        call dgemm_my('N', 'N', nelorbh, nelorb_c, nelorb_c, 1.d0, mu_c, nelorbh  &
       &, detmat_c(1, nelorb_c + 1), 2*nelorb_c, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
        call dgemm_my('N', 'T', nelorbh, nelorbh, nelorb_c, 1.d0, psip, nelorbh   &
       &, mu_c, nelorbh, 0.d0, detmat(1, nelorbh + 1), 2*nelorbh, nprocu, rankopt, commopt_mpi)
    end subroutine scontract_genj

    subroutine scontract_mat_jas(nelorbh, nelorb, nelcol, nelorb_c         &
         &, nelcol_c, detmat, detmat_c, mu_c, psip)
        implicit none
        integer nelorbh, nelorb, nelorb_c, nelcol, nelcol_c, i
        real*8 detmat(nelorb, nelcol), detmat_c(nelorb_c, *)              &
             &, mu_c(nelorbh, *), psip(nelorbh, *)
! WARNING  IT REFERS ONLY TO jASTROW SO IPF AND UNPAIRED ARE NOT USEd
!#ifdef __CASO
!    nprocu=nprocopt
!#else
!    nprocu=1
!#endif
        detmat = 0.d0
        call dgemm_my('N', 'N', nelorbh, nelcol_c, nelorb_c, 1.d0, mu_c, nelorbh  &
             &, detmat_c, nelorb_c, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
        call dgemm_my('N', 'T', nelorbh, nelorbh, nelorb_c, 1.d0, psip, nelorbh   &
             &, mu_c, nelorbh, 0.d0, detmat, nelorb, nprocu, rankopt, commopt_mpi)
    end subroutine scontract_mat_jas

    subroutine scontract_mat_det(nelorbh, nelorb, nelcol, nelorb_c         &
         &, nelcol_c, detmat, detmat_c, mu_c, psip)
        implicit none
        integer nelorbh, nelorb, nelorb_c, nelcol, nelcol_c, i
        real*8 detmat(ipc*ipf*nelorb, nelcol), detmat_c(ipc*nelorb_c, *)              &
             &, mu_c(ipc*ipf*nelorbh, *), psip(ipf*ipc*nelorbh, *)
!#ifdef __CASO
!    nprocu=nprocopt
!#else
!    nprocu=1
!#endif
        if (ipc .eq. 2) then
            detmat = 0.d0
            call zgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, (1.d0, 0.d0), mu_c, ipf*nelorbh  &
                 &, detmat_c, nelorb_c, (0.d0, 0.d0), psip, ipf*nelorbh, nprocu, rankopt, commopt_mpi)
            call zgemm_my('N', 'T', ipf*nelorbh, ipf*nelorbh, nelorb_c, (1.d0, 0.d0), psip, ipf*nelorbh   &
                 &, mu_c, ipf*nelorbh, (0.d0, 0.d0), detmat, ipf*nelorb, nprocu, rankopt, commopt_mpi)
            if (nelcol_c .gt. nelorb_c) then
                do i = nelorb_c + 1, nelcol_c
                    call zcopy(ipf*nelorbh, psip(1, i), 1, detmat(1, ipf*nelorb + i - nelorb_c), 1)
                end do
            end if
        else
            detmat = 0.d0
            call dgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, 1.d0, mu_c, ipf*nelorbh&
                 &, detmat_c, nelorb_c, 0.d0, psip, ipf*nelorbh, nprocu, rankopt, commopt_mpi)
            call dgemm_my('N', 'T', ipf*nelorbh, ipf*nelorbh, nelorb_c, 1.d0, psip, ipf*nelorbh&
                    &, mu_c, ipf*nelorbh, 0.d0, detmat, ipf*nelorb, nprocu, rankopt, commopt_mpi)
            if (nelcol_c .gt. nelorb_c) then
                do i = nelorb_c + 1, nelcol_c
                    call dcopy(ipf*nelorbh, psip(1, i), 1, detmat(1, nelorb*ipf + i - nelorb_c), 1)
                end do
            end if
        end if
    end subroutine scontract_mat_det

    subroutine update_kgrid
        implicit none
        integer i, j, ii, jj, kk, ll, count1, count2, count1j, count2j, indpar, indparp, kboundi&
       &, ind, kboundi_max, iii, jjj, kkk, k, maxdim
        integer, dimension(:, :), allocatable :: kshell_map
        real*8 kbound, map_tmp(3, 27), max_rejected, cost_tilted
        logical, external :: slaterorb
        logical not_found
        integer, dimension(:, :), allocatable :: kpip_sav
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization of direct lattice vectors for building the
! periodic basis set. The definition is the same employed in the
! Crystal DFT code. This basis can be used for both real and complex
! wave functions. In the case of complex wave function it can be used
! for an open system too.
! Use the keyword "PBC_C" in the first line of the wave function to enforce this option.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (allocated(kgrid)) deallocate (kgrid)
        if (allocated(kgrid_atom)) deallocate (kgrid_atom)
        ikshift = nshell
        if (abs(LBox) .eq. 3.d0) then
            ! determine the maximum number of direct lattice vectors
            if (yes_crystalj) then
                allocate (kshell_map(3, nshell + nshellj))
            else
                allocate (kshell_map(3, nshell))
            end if
            kshell_map = 0
            indpar = 0
            do i = 1, nshell
                indparp = indpar + 1
                if (lepsbas .gt. 0.d0) then
                    do j = 1, 3
                        if (cellscale(j) .ne. 0.d0) then
                            if (slaterorb(ioptorb(i))) then ! STO
                                kbound = 0.5d0 + (lepsbas/dupr(indparp)/metric_min)/cellscale(j)
                            else ! GTO
                                kbound = 0.5d0 + (dsqrt(lepsbas/dupr(indparp))/metric_min)/cellscale(j)
                            end if
                            kshell_map(j, i) = kbound + 1
                        end if
                    end do
                else
                    do j = 1, 3
                        kshell_map(j, i) = 0
                    end do
                end if
                indpar = indpar + nparam(i)
            end do
            if (yes_crystalj) then
                indpar = 0
                do i = 1, nshellj
                    indparp = indpar + 1
                    if (lepsbas .gt. 0.d0 .and. ioptorbj(i) .ne. 200) then
                        do j = 1, 3
                            if (cellscale(j) .ne. 0.d0) then
                                if (slaterorb(ioptorbj(i))) then ! STO
                                    kbound = 0.5d0 + (lepsbas/vjur(indparp)/metric_min)/cellscale(j)
                                else ! GTO
                                    kbound = 0.5d0 + (dsqrt(lepsbas/vjur(indparp))/metric_min)/cellscale(j)
                                end if
                                kshell_map(j, i + ikshift) = kbound + 1
                            end if
                        end do
                    else
                        do j = 1, 3
                            kshell_map(j, i + ikshift) = 0
                        end do
                    end if
                    indpar = indpar + nparamj(i)
                end do
            end if

            ! build the lattice vector map needed to compute the basis set
            if (yes_crystalj) then
                allocate (kgrid(nshell + nshellj))
                allocate (kgrid_atom(2*nion))
            else
                allocate (kgrid(nshell))
                allocate (kgrid_atom(nion))
            end if

            if (yes2d .or. yes1d) kshell_map(3, :) = 0
            if (yes1d) kshell_map(2, :) = 0

            if (iespbc) then
            do ii = 1, nion
                kboundi_max = 0
                do jj = adr_nion(ii), adr_nion(ii + 1) - 1
                    i = ind_nion(jj)
!    do i=1,nshell
                    kboundi = (2*kshell_map(1, i) + 1)*(2*kshell_map(2, i) + 1)*(2*kshell_map(3, i) + 1)
                    kboundi_max = max(kboundi_max, kboundi)
                    allocate (kgrid(i)%kpip(3, kboundi))
                    kgrid(i)%kpip(:, :) = 0
                end do
                allocate (kgrid_atom(ii)%kpip(3, kboundi_max))
                kgrid(ii)%kpip(:, :) = 0
            end do
            if (yes_crystalj) then
            do ii = 1, nion
                kboundi_max = 0
                do jj = adrj_nion(ii), adrj_nion(ii + 1) - 1
                    i = indj_nion(jj) + nshell
                    kboundi = (2*kshell_map(1, i) + 1)*(2*kshell_map(2, i) + 1)*(2*kshell_map(3, i) + 1)
                    kboundi_max = max(kboundi_max, kboundi)
                    allocate (kgrid(i)%kpip(3, kboundi))
                    kgrid(i)%kpip(:, :) = 0
                end do
                allocate (kgrid_atom(ii + nion)%kpip(3, kboundi_max))
                kgrid(ii + nion)%kpip(:, :) = 0
            end do
            end if
            else
            do i = 1, nshell
                allocate (kgrid(i)%kpip(3, 1))
                kgrid(i)%kpip(:, :) = 0
            end do
            end if

            max_rejected = lepsbas
            count1 = 0
            count2 = 0

            if (iespbc) then
                indpar = 0
                do iii = 1, nion
                    kgrid_atom(iii)%dimshell = 0
                    do jjj = adr_nion(iii), adr_nion(iii + 1) - 1
                        i = ind_nion(jjj)
!       do i=1,nshell
                        indparp = indpar + 1
                        kgrid(i)%dimshell = 0
                        do kk = -kshell_map(3, i), kshell_map(3, i)
                            do jj = -kshell_map(2, i), kshell_map(2, i)
                                do ii = -kshell_map(1, i), kshell_map(1, i)
                                    ! optimize the number of vectors in the summation
                                    ! to be inside a sphere.
                                    count2 = count2 + 1
!In order to estimate a lower bound of exp(- Z_basis  |z|)
!for each ii,jj,kk it is computed  the minimum distance of a point z from origin
!    dist  = [ (x+ii Lx,y+jj Ly, z+kk Lz, metric (x+ii Lx, y+jj Ly, z+kk Lz)]
!  with the  condition |x|<Lx/2 |y|<Ly/2 |z|<Lz/2
! The minimum is typically at the boundaries x=+/-Lx/2,y=+/-Ly/2,z=+/-Lz/2
! apart cases when ii=0 or jj=0 or kk=0, when,  under some conditions
!  depending on the metric, x=0 or y=0 or z=0.

                                    map_tmp(1, 1) = cellscale(1)*ii
                                    map_tmp(2, 1) = cellscale(2)*jj
                                    map_tmp(3, 1) = cellscale(3)*kk

                                    call prep_map(map_tmp, cellscale)

! NB to avoid roundoff problems ii,jj,kk--> min(abs(ii),1) =0 ,1 only.

                                    cost_tilted = norm_metric(map_tmp(1, 1), metric)**2
                                    do ll = 2, 27
                                        cost_tilted = min(cost_tilted, norm_metric(map_tmp(1, ll), metric)**2)
                                    end do

                                    if (slaterorb(ioptorb(i))) then
                                        kbound = dupr(indparp)*dsqrt(cost_tilted)
                                    else
                                        kbound = dupr(indparp)*cost_tilted
                                    end if
                                    if (kbound .lt. lepsbas .or. lepsbas .le. 0.d0) then
                                        count1 = count1 + 1
                                        kgrid(i)%dimshell = kgrid(i)%dimshell + 1
                                        kgrid(i)%kpip(1, kgrid(i)%dimshell) = ii
                                        kgrid(i)%kpip(2, kgrid(i)%dimshell) = jj
                                        kgrid(i)%kpip(3, kgrid(i)%dimshell) = kk

                                        not_found = .true.
                                        do kkk = 1, kgrid_atom(iii)%dimshell
                                            if (kgrid_atom(iii)%kpip(1, kkk) .eq. ii .and.&
                                           & kgrid_atom(iii)%kpip(2, kkk) .eq. jj .and.&
                                           & kgrid_atom(iii)%kpip(3, kkk) .eq. kk) not_found = .false.
                                        end do
                                        if (not_found) then
                                            kgrid_atom(iii)%dimshell = kgrid_atom(iii)%dimshell + 1
                                            kgrid_atom(iii)%kpip(1, kgrid_atom(iii)%dimshell) = ii
                                            kgrid_atom(iii)%kpip(2, kgrid_atom(iii)%dimshell) = jj
                                            kgrid_atom(iii)%kpip(3, kgrid_atom(iii)%dimshell) = kk
                                        end if
                                    else
                                        if (kbound .lt. max_rejected .or. max_rejected .eq. lepsbas)&
                                        &max_rejected = kbound
                                    end if
                                end do
                            end do
                        end do
                        indpar = indpar + nparam(i)
                    end do ! end ii
                end do ! end ion
!       Shrink the memory allocated
                maxdim = 0
                do i = 1, nion
                    maxdim = max(maxdim, kgrid_atom(i)%dimshell)
                end do
                allocate (kpip_sav(3, maxdim))
                do i = 1, nion
                    do j = 1, kgrid_atom(i)%dimshell
                        kpip_sav(:, j) = kgrid_atom(i)%kpip(:, j)
                    end do
                    deallocate (kgrid_atom(i)%kpip)
                    allocate (kgrid_atom(i)%kpip(3, kgrid_atom(i)%dimshell))
                    do j = 1, kgrid_atom(i)%dimshell
                        kgrid_atom(i)%kpip(:, j) = kpip_sav(:, j)
                    end do
                end do
                do i = 1, nshell
                    do j = 1, kgrid(i)%dimshell
                        kpip_sav(:, j) = kgrid(i)%kpip(:, j)
                    end do
                    deallocate (kgrid(i)%kpip)
                    allocate (kgrid(i)%kpip(3, kgrid(i)%dimshell))
                    do j = 1, kgrid(i)%dimshell
                        kgrid(i)%kpip(:, j) = kpip_sav(:, j)
                    end do
                end do
                deallocate (kpip_sav)
                count1j = 0
                count2j = 0
                if (yes_crystalj) then
                    indpar = 0
                    do iii = 1, nion
                        kgrid_atom(iii + nion)%dimshell = 0
                        do jjj = adrj_nion(iii), adrj_nion(iii + 1) - 1
!       do i=1,nshellj
                            i = indj_nion(jjj)
                            indparp = indpar + 1
                            kgrid(i + ikshift)%dimshell = 0
                            do kk = -kshell_map(3, i + ikshift), kshell_map(3, i + ikshift)
                                do jj = -kshell_map(2, i + ikshift), kshell_map(2, i + ikshift)
                                    do ii = -kshell_map(1, i + ikshift), kshell_map(1, i + ikshift)
                                        ! optimize the number of vectors in the summation
                                        ! to be inside a sphere.
                                        count2j = count2j + 1
!  for each ii,jj,kk compute  the minimum distance of a point z from origin
!    dist  = [ (x+ii Lx,y+jj Ly, z+kk Lz, metric (x+ii Lx, y+jj Ly, z+kk Lz)]
!  with the  condition |x|<Lx/2 |y|<Ly/2 |z|<Lz/2
! The minimum is typically at the boundaries x=+/-Lx/2,y=+/-Ly/2,z=+/-Lz/2
! apart cases where ii=0 or jj=0 or kk=0, when,  under some conditions
!  depending on the metric, x=0 or y=0 or z=0.

                                        map_tmp(1, 1) = cellscale(1)*ii
                                        map_tmp(2, 1) = cellscale(2)*jj
                                        map_tmp(3, 1) = cellscale(3)*kk

                                        call prep_map(map_tmp, cellscale)

                                        cost_tilted = norm_metric(map_tmp(1, 1), metric)**2
                                        do ll = 2, 27
                                            cost_tilted = min(cost_tilted, norm_metric(map_tmp(1, ll), metric)**2)
                                        end do

                                        if (ioptorbj(i) .ne. 200) then
                                        if (slaterorb(ioptorbj(i))) then
                                            kbound = vjur(indparp)*dsqrt(cost_tilted)
                                        else
                                            kbound = vjur(indparp)*cost_tilted
                                        end if
                                        else
                                        kbound = 0.d0
                                        end if
                                        if (kbound .lt. lepsbas .or. lepsbas .le. 0.d0) then
                                            count1j = count1j + 1
                                            kgrid(i + ikshift)%dimshell = kgrid(i + ikshift)%dimshell + 1
                                            kgrid(i + ikshift)%kpip(1, kgrid(i + ikshift)%dimshell) = ii
                                            kgrid(i + ikshift)%kpip(2, kgrid(i + ikshift)%dimshell) = jj
                                            kgrid(i + ikshift)%kpip(3, kgrid(i + ikshift)%dimshell) = kk
                                            not_found = .true.
                                            do kkk = 1, kgrid_atom(iii + nion)%dimshell
                                                if (kgrid_atom(iii + nion)%kpip(1, kkk) .eq. ii .and.&
                                               & kgrid_atom(iii + nion)%kpip(2, kkk) .eq. jj .and.&
                                         & kgrid_atom(iii + nion)%kpip(3, kkk) .eq. kk) not_found = .false.
                                            end do
                                            if (not_found) then
                                                kgrid_atom(iii + nion)%dimshell = kgrid_atom(iii + nion)%dimshell + 1
                                                kgrid_atom(iii + nion)%kpip(1, kgrid_atom(iii + nion)%dimshell) = ii
                                                kgrid_atom(iii + nion)%kpip(2, kgrid_atom(iii + nion)%dimshell) = jj
                                                kgrid_atom(iii + nion)%kpip(3, kgrid_atom(iii + nion)%dimshell) = kk
                                            end if
                                        else
                                            if (kbound .lt. max_rejected .or.&
                                                & max_rejected .eq. lepsbas .and.&
                                                & ioptorbj(i) .ne. 200) &
                                                &max_rejected = kbound
                                        end if
                                    end do
                                end do
                            end do
                            indpar = indpar + nparamj(i)
                        end do ! ii
                    end do ! nion
!       Shrink the memory allocated
                    maxdim = 0
                    do i = 1, nion
                        maxdim = max(maxdim, kgrid_atom(i + nion)%dimshell)
                    end do
                    allocate (kpip_sav(3, maxdim))
                    do i = 1, nion
                        do j = 1, kgrid_atom(i + nion)%dimshell
                            kpip_sav(:, j) = kgrid_atom(i + nion)%kpip(:, j)
                        end do
                        deallocate (kgrid_atom(i + nion)%kpip)
                        allocate (kgrid_atom(i + nion)%kpip(3, kgrid_atom(i + nion)%dimshell))
                        do j = 1, kgrid_atom(i + nion)%dimshell
                            kgrid_atom(i + nion)%kpip(:, j) = kpip_sav(:, j)
                        end do
                    end do
                    do i = 1, nshellj
                        do j = 1, kgrid(i + ikshift)%dimshell
                            kpip_sav(:, j) = kgrid(i + ikshift)%kpip(:, j)
                        end do
                        deallocate (kgrid(i + ikshift)%kpip)
                        allocate (kgrid(i + ikshift)%kpip(3, kgrid(i + ikshift)%dimshell))
                        do j = 1, kgrid(i + ikshift)%dimshell
                            kgrid(i + ikshift)%kpip(:, j) = kpip_sav(:, j)
                        end do
                    end do
                    deallocate (kpip_sav)
                end if
            else
                count1 = 0
                do i = 1, nshell
                    kgrid(i)%dimshell = 1
                    kgrid(i)%kpip(1, kgrid(i)%dimshell) = 0
                    kgrid(i)%kpip(2, kgrid(i)%dimshell) = 0
                    kgrid(i)%kpip(3, kgrid(i)%dimshell) = 0
                end do
                cellscale = 0.d0 ! initialize to be sure the irrelevant vector cellscale
            end if
            deallocate (kshell_map)
            if (iespbc) then
            do ii = 1, nion
                do jj = adr_nion(ii), adr_nion(ii + 1) - 1
                    i = ind_nion(jj)
                    allocate (kgrid(i)%tobedone(kgrid_atom(ii)%dimshell))
                    kgrid(i)%tobedone(:) = .false.
                end do
            end do
            if (yes_crystalj) then
            do ii = 1, nion
                do jj = adrj_nion(ii), adrj_nion(ii + 1) - 1
                    i = indj_nion(jj) + nshell
                    allocate (kgrid(i)%tobedone(kgrid_atom(ii + nion)%dimshell))
                    kgrid(i)%tobedone(:) = .false.
                end do
            end do
            end if
            do iii = 1, nion
                do jjj = adr_nion(iii), adr_nion(iii + 1) - 1
                    i = ind_nion(jjj)
                    do j = 1, kgrid_atom(iii)%dimshell
                        do k = 1, kgrid(i)%dimshell
                        if (kgrid_atom(iii)%kpip(1, j) .eq. kgrid(i)%kpip(1, k)&
                        &.and. kgrid_atom(iii)%kpip(2, j) .eq. kgrid(i)%kpip(2, k)&
                        &.and. kgrid_atom(iii)%kpip(3, j) .eq. kgrid(i)%kpip(3, k)) then
                            kgrid(i)%tobedone(j) = .true.
                        end if
                        end do
                    end do
                end do
            end do
            if (yes_crystalj) then
                do iii = 1, nion
                    do jjj = adrj_nion(iii), adrj_nion(iii + 1) - 1
                        i = indj_nion(jjj)
                        do j = 1, kgrid_atom(iii + nion)%dimshell
                            do k = 1, kgrid(i + ikshift)%dimshell
                                if (kgrid_atom(iii + nion)%kpip(1, j) .eq. kgrid(i + ikshift)%kpip(1, k)&
                               &.and. kgrid_atom(iii + nion)%kpip(2, j) .eq. kgrid(i + ikshift)%kpip(2, k)&
                               &.and. kgrid_atom(iii + nion)%kpip(3, j) .eq. kgrid(i + ikshift)%kpip(3, k)) then
                                    kgrid(i + ikshift)%tobedone(j) = .true.
                                end if
                            end do
                        end do
                    end do
                end do
            end if
            end if
            if (rank .eq. 0) then
                write (6, *) 'Warning: updated kgrid considering Det/Jastrow:', count1, count1j
                write (6, *) 'Warning: lowest wf discarded=', exp(-max_rejected)
#ifdef __DEBUG
                write (6, *) ' Grid considered /ion'
                do i = 1, nion
                    write (6, *) ' ion # ', i, 'dimension =', kgrid_atom(i)%dimshell
                    do j = 1, kgrid_atom(i)%dimshell
                        write (6, *) j, kgrid_atom(i)%kpip(1, j), kgrid_atom(i)%kpip(2, j), kgrid_atom(i)%kpip(3, j)
                    end do
                end do
                if (yes_crystalj) then
                    write (6, *) ' Grid considered Jastrow/ion '
                    do i = nion + 1, 2*nion
                        write (6, *) ' ion # ', i - nion, 'dimension =', kgrid_atom(i)%dimshell
                        do j = 1, kgrid_atom(i)%dimshell
                            write (6, *) j, kgrid_atom(i)%kpip(1, j), kgrid_atom(i)%kpip(2, j), kgrid_atom(i)%kpip(3, j)
                        end do
                    end do
                end if
                write (6, *) ' Grid considered '
                do i = 1, nshell
                    write (6, *) ' Shell # ', i, 'dimension =', kgrid(i)%dimshell
                    do j = 1, kgrid(i)%dimshell
                        write (6, *) j, kgrid(i)%kpip(1, j), kgrid(i)%kpip(2, j), kgrid(i)%kpip(3, j)
                    end do
                end do
                if (yes_crystalj) then
                    write (6, *) ' Grid considered Jastrow '
                    do i = 1, nshellj
                        write (6, *) ' Shell # ', i, 'dimension =', kgrid(i + ikshift)%dimshell
                        do j = 1, kgrid(i + ikshift)%dimshell
                            write (6, *) j, kgrid(i + ikshift)%kpip(1, j), &
                                kgrid(i + ikshift)%kpip(2, j), kgrid(i + ikshift)%kpip(3, j)
                        end do
                    end do
                end if
#endif
            end if
!    kgrid%kpip is no longer needed
            do i = 1, nshell
                if (allocated(kgrid(i)%kpip)) deallocate (kgrid(i)%kpip)
            end do
            if (yes_crystalj) then
            do i = nshell + 1, nshell + nshellj
                if (allocated(kgrid(i)%kpip)) deallocate (kgrid(i)%kpip)
            end do
            end if
        end if
    end subroutine update_kgrid

    function norm_metric(r, metric)
        implicit none
        real*8 norm_metric, r(3), metric(3, 3)
        norm_metric = metric(1, 1)*r(1)*r(1) + metric(2, 2)*r(2)*r(2) + metric(3, 3)*r(3)*r(3)&
                & + 2.d0*(metric(1, 2)*r(1)*r(2) + metric(1, 3)*r(1)*r(3) + metric(2, 3)*r(2)*r(3))
        norm_metric = dsqrt(max(norm_metric, 0.d0)) ! protection from roundoff
        return
    end

end module allio

subroutine prep_map(map_tmp, cellscale)
    implicit none
    real*8 map_tmp(3, 3, 3, 3), cellscale(3), cellhalf(3, 3)
    integer i, j, k
! Here we should  find  the possible argmin of the metric.
! For each coordinate there are two possibilities either the minimum is at  the
! boundary +/- L/2 or the minimum is at the current position (when is zero).
! Thus we end up with  27 possibilities 3 for each coordinate including
! the input  map(:,1,1,1).
    cellhalf(:, 1) = 0.d0
    cellhalf(:, 2) = -cellscale(:)/2.d0
    cellhalf(:, 3) = cellscale(:)/2.d0

    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                map_tmp(1, i, j, k) = map_tmp(1, 1, 1, 1) + cellhalf(1, i)
                map_tmp(2, i, j, k) = map_tmp(2, 1, 1, 1) + cellhalf(2, j)
                map_tmp(3, i, j, k) = map_tmp(3, 1, 1, 1) + cellhalf(3, k)
            end do
        end do
    end do
    return
end
