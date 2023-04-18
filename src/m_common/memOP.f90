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

subroutine default_allocate
    use allio
    implicit none
#ifdef  PARALLEL
    include 'mpif.h'
#endif
    yes_sparse_choose = .false.
    yes_sparse = .false.
    novec_loop1 = .true.
    iflagerrall = 0
    writescratch = 0
    yesdft = .true. ! no QMC arrays are allocated by default
    yes_scemama = .false.
    lrdmc_der = .false.
    epscut = 0.d0
    epscuttype = 0
    neigh = 1
    membig = .true.
    yeszj = .false.
    yeszagp = .false.
    epsilon0 = 1.d0
    Cgauss = 1.d0
    case_diel = 0
    test_aad = .false.
    epsbas = 1d-6
    kp_type = -1
    indt = 0
    yesfast = 2
    only_molecular = .false.
    manyfort10 = .false.
    allowed_averagek = .true.
    skip_equivalence = .true.
    symmetrize_agp = .false.
    molopt = 0
    yesmin = 0
    yesmin_read = .false.
#ifdef  PARALLEL
    rankopt = rank
    commopt_mpi = MPI_COMM_WORLD
    nprocopt = nproc
    ! default values for k-points/beads communicators for the tools
    ! needed by the new version of read_fort10
    rankrep = rank
    commrep_mpi = MPI_COMM_WORLD
    nprocrep = nproc
    rankcolrep = 0
    commcolrep_mpi = MPI_COMM_WORLD
    nproccolrep = 1
#else
    rank = 0
    nproc = 1
    rankopt = 0
    commopt_mpi = 0
    nprocopt = 1
    rankrep = 0
    commrep_mpi = 0
    nprocrep = 1
    rankcolrep = 0
    commcolrep_mpi = 0
    nproccolrep = 1
#endif
    ip4 = 4 ! gradient and laplacian in upnewwf
    nw = 1
    nws = 1
    nwm = 1
    in1 = 1
    detc_proj = .false.
    cellderiv = .false.
    add_onebody2det = .false.
    ieskindim = 1
    ieskingdim = 1
    nbindim = 1
    ncgdim = 1
    npdim = 1

    nbinmax = 1

    !        no optimization considered
    np = 0
    npm = 1
    npmn = 1

    nweight = 0
    ibinit = 0
    nfat = 0

    kl = 0

    ! make the simplest vmc allocation
    itestr4 = 2
    itestr3 = 2
    itestrr = 2
    itestr = 2
    itest = 2

    !          No parameters to optimize nor to compute
    ieser = 0
    iese = 0
    isfix = 0
    iesinv = 0
    iesm = 0
    iesd = 0
    iesfree = 0
    iessw = 0
    iesup = 0
    ieskin = 0
    iesking = 0
    ieskint = 0
    typedyncell = -2
    yespress = .false.
    oldscaling = .false.

    !        no use of pseudo is assumed
    nintpsa = 0
    npsa = 0
    npsar = 0
    npsamax = 0
    !        minimal output cnsidered
    iread = 1

end subroutine default_allocate

subroutine deallocate_all
    use allio
    implicit none
    if (allocated(indtm)) deallocate (indtm)
    if (allocated(naccm)) deallocate (naccm, jbra, jbraw, tcost)
    if (allocated(ivic)) deallocate (pseudolocal, ivic)
    if (allocated(angle)) deallocate (prefactor, angle)
    if (allocated(dup_c)) deallocate (dup_c, scale_c&
            &, multranspip, iesuptransb)
    if (allocated(dup)) deallocate (dup)
    if (allocated(vju_c)) deallocate (vju_c, vju, scalej_c, multranspipj   &
            &, iesuptransbj)
    if (allocated(scalejsz_c)) deallocate (scalejsz_c)
    if (allocated(vj)) deallocate (vj, dsw, ddw, ddwsz, dek, dekg, scalpar)

    !###################### INTEGER ARRAY #############################
    if (allocated(nozero_c)) deallocate (nozero_c, jbradet, sjbradet    &
            &, jbradetn, mult_c, kion_c, ioptorb_c, nparam_c)
    if (allocated(nozeroj_c)) deallocate (nozeroj_c, jbraj, jbrajn        &
            &, multj_c, kionj_c, ioptorbj_c, nparamj_c)
    if (allocated(jbrajsz)) deallocate (jbrajsz)
    if (allocated(itouch)) deallocate (itouch)

    ! ################## READ ALL NEC. for input ######################
    if (allocated(wconfw)) deallocate (wconfw, econfw, psilnw            &
            &, rcarto, rcart, dx_old, dx_new                                 &
            &, kel, keln, zetar, rion      &
            &, rionsav, gradpsibar, gradpsi, gradtot                      &
            &, gradtotbar, gradpsiold, gradpsibarold                              &
            &, jastrowall_ee, jasnew_ee, jastrowall_ei, jasnew_ei                  &
            &, rcne)
    if (allocated(ion_table)) deallocate (ion_table)
    !     do i=1,ieskinr
    !     deallocate(ion_table(i)%ion,ion_table(i)%comp)
    !     enddo
    !     deallocate(ion_table)
    !    endif
    if (allocated(atom_number)) deallocate (atom_number, type_atom, pointvj)
    if (allocated(zetar_fast)) deallocate (zetar_fast, atom_number_fast, rion_fast)
    if (allocated(rmucos)) deallocate (rmucos, rmusin)
    if (allocated(ipsip)) deallocate (ipsip)
    if (allocated(dupr)) deallocate (dupr)
    if (allocated(mult)) deallocate (mult, kion, ioptorb, nparam)
    if (allocated(slaterorb_read)) deallocate (slaterorb_read)
    if (allocated(scale)) deallocate (scale)
    if (allocated(nozero)) deallocate (nozero)
    if (allocated(nozerodet)) deallocate (nozerodet)
    if (allocated(vjur)) deallocate (vjur)
    if (allocated(multj)) deallocate (multj, kionj, ioptorbj, nparamj)
    if (allocated(scalej)) deallocate (scalej, nozeroj, nozerojder)

    if (allocated(scalejsz)) deallocate (scalejsz)
    if (allocated(ioccup_c)) deallocate (ioccup_c)
    if (allocated(ioccj_c)) deallocate (ioccj_c)

    if (allocated(ioccup)) deallocate (ioccup, ioccdo)
    if (allocated(ioccj)) deallocate (ioccj)
    if (allocated(iesuptrans)) deallocate (iesuptrans)
    if (allocated(iesuptransj)) deallocate (iesuptransj)

    if (allocated(mu_c)) deallocate (mu_c)
    if (allocated(transpip)) deallocate (transpip)

    if (allocated(muj_c)) deallocate (muj_c)
    if (allocated(transpipj)) deallocate (transpipj)

    if (allocated(iond)) deallocate (winv, iond, rmu, r, econf, econfh&
            &, wconfn, econfion, factorsr, vcut, wsto, zeta, tabpip, table, tabler&
            &, diag, ainv, psiln, psidetln, psisn, vpot, vpotreg, enertrue, diffuse, winvdo &
            &, winvup, winvj, winvbar, winvjbar, wint, dist, tmu, diagfn, enert, berry_exp&
            &, psip, winvsj, psinew, ainvup, ainvupb, ainvs, ainvdo, enerint  &
            &, dists, detmat, dists_kel, jasmat, winvs, cnorm, yescut, diffkin&
            &, wconfsav, iond_cart)

    if (allocated(cnorm_nw)) deallocate (cnorm_nw)

    if (allocated(orbcost)) deallocate (winvjbarsz, orbcost, jasmatsz)
    if (allocated(projm)) deallocate (projm)
    if (allocated(jasmat_c)) deallocate (jasmat_c)
    if (allocated(jasmatsz_c)) deallocate (jasmatsz_c)
    if (allocated(jbraiesup)) deallocate (jbraiesup, jbraiesup_sav)
    if (allocated(jbraiesm)) deallocate (jbraiesm, jbraiesm_sav)
    if (allocated(tcore)) deallocate (icore, tcore)
    if (allocated(jpseudo)) deallocate (nparpshell, jpseudo, kindion     &
            &, pshell, parshell, rcutoff)
    if (allocated(wpseudo))                                           &
            &deallocate (wpseudo, legendre, versor, wintpseudo)
    if (allocated(costz)) deallocate (costz, costz3, zetaq)
    if (allocated(whereiesup)) deallocate (whereiesup)
    if (allocated(whereiesm)) deallocate (whereiesm)
    if (allocated(orbcostn)) deallocate (orbcostn)
    if (allocated(rpar)) deallocate (rpar, adrlambda&
            &, typeorb, jas_invariant, orbps)
    if (allocated(kiontot)) deallocate (kiontot)
    if (allocated(kiontotj)) deallocate (kiontotj)
    if (allocated(ioptorbja)) deallocate (ioptorbja)
    if (allocated(orbcostl)) deallocate (orbcostl)
    if (allocated(cov)) deallocate (cov)
    if (allocated(cov_old)) deallocate (cov_old)
    if (allocated(sov)) deallocate (sov)
    if (allocated(reduce)) then
        deallocate (reduce, velion, ef, efenergy, err, force, efp&
                &, fk, fkav, okav, skdiag, efpress, allowed_par)
    end if
    if (allocated(first_moment)) deallocate (first_moment, second_moment)
    if (allocated(reducel)) deallocate (reducel, v_adr, mat_adr, ipip_adr)
    if (allocated(bufscra)) deallocate (bufscra)
    if (iespbc) then

        if (allocated(q)) then
            deallocate (sum_q_cos_gr, sum_q_sin_gr &
                    &, q, s_cord, qphase, q_cos_gr, q_sin_gr, g, gg&
                    &, factor, ind_g, sum_q_cos, sum_q_sin)
            if (allocated(phsfac1)) deallocate (phsfac1, phsfac2, phsfac3)
        elseif (allocated(sum_q_cos_gr)) then
            deallocate (sum_q_cos_gr, sum_q_sin_gr)
        end if
        if (allocated(cosphase)) deallocate (cosphase, sinphase)
        if (allocated(sinphaseb)) deallocate (sinphaseb, cosphaseb)

    end if
    if (allocated(allowcontr)) deallocate (allowcontr)
    if (allocated(allowcontrj)) deallocate (allowcontrj)
    if (allocated(indpar_tab)) deallocate (indpar_tab, indorb_tab&
            &, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab, adr_nion, ind_nion&
            &, adrj_nion, indj_nion)
    if (allocated(cellscalep)) deallocate (cellscalep)
    if (allocated(rphasep)) deallocate (rphasep)
    if (allocated(s2rp)) deallocate (s2rp)
    if (allocated(x_neigh)) deallocate (x_neigh)
    if (allocated(dist_shift)) deallocate (dist_shift)
    if (allocated(disto_shift)) deallocate (disto_shift)
    if (allocated(distreg_shift)) deallocate (distreg_shift)
    if (allocated(distrego_shift)) deallocate (distrego_shift)
    if (allocated(zetar_fast))&
            &deallocate (zetar_fast, rion_fast, atom_number_fast)
    if (allocated(winvfn)) deallocate (winvfn, winvbarfn)
    if (allocated(vpotsav_ee)) deallocate (vpotsav_ee)
    if (allocated(warpmat)) deallocate (warpmat)
    if (allocated(detmat_proj)) deallocate (detmat_proj)
    if (allocated(projmat_c)) deallocate (projmat_c)
    if (allocated(kgrid)) deallocate (kgrid)
    if (allocated(kgrid_atom)) deallocate (kgrid_atom)
    if (allocated(npip)) deallocate (npip)
    if (allocated(rmunew)) deallocate (rmunew)
    if (allocated(rknew)) deallocate (rknew)
    if (allocated(rmunewb)) deallocate (rmunewb)
    if (allocated(rknewb)) deallocate (rknewb)
    if (allocated(derEV)) deallocate (derEV)
    if (allocated(p_pulay)) deallocate (p_pulay)
    if (allocated(singdet)) deallocate (singdet)
    if (allocated(detmat_c)) deallocate (detmat_c)
    if (allocated(duprold)) deallocate (duprold)
    if (allocated(zetamin)) deallocate (zetamin, distmin)
    if (allocated(eagp_pfaff)) deallocate (eagp_pfaff)
    if (allocated(agp)) deallocate (agp)
    if (allocated(agpn)) deallocate (agpn)
    if (allocated(etot)) deallocate (etot)
end subroutine deallocate_all
