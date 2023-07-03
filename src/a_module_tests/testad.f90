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

program main
    use Constants
    use Cell
    use Ewald
    use allio
    use convertmod
    use IO_m
    implicit none
    integer dim_jasmat, nelorbjmax, indtmax, nshelljmax, iscramaxb, nk_cell, min_cell
    real*8 eloc(2), logpsi(2), press(3), press_pulay(3), scale_spsi
    real*8 drand1, enercont, jacobian, dnrm2, mapping&
            &, costwnn, costexpn, coeff_nw, enerold, costfact, ttrysav&
            &, ratiowconfn, wconfnsav, ttryt, enercuto, psioldinp, voffpseudo&
            &, reweight_fn, srforce, wforce, srforcew, wforcew, costmpi, costsr &
            &, identity, timemcp, timeoptp, vr, vl, celldmsav, costprsav, elocb(2), logpsib(2) &
            &, cost0, costp, costm, epsad, kel_sav, elocp(2), elocm(2), logpsip(2), logpsim(2), pseudolocalb&
            &, scalepulayw, dt4, celldm_write(3), rs_write, aadinp(2)
    integer nleft, i, j, k, ind, ii, jj, kk, imin, imax, imaxn, stodim&
            &, iese_eff, np3p3, maxrank, maxdimeig, nmollimit, firstmolt, new_threads&
            &, iflagpip, indshell, dimfk, perbin, indopen, indopen3, yesfastj, iesdvj
    logical noproj, flagmu, flagcont, forceyes&
            &, flag, iesmeas, dir_exist, nn, yescomm, someparameter, someparameterdet
    integer :: iargc, isdistp, iese_eff_c, iscramaxold, ieskinion, Ltot
    complex(8) voffpseudo_c
    integer ILAENV
    external ILAENV
    real*8 cclock, cellscaleb(3), sr2b(9), cellscaleo(12), cellscalelb(3), cellscalen(12), sr2lb(9)
    real*8, dimension(:, :), allocatable :: kelsav
    real*8, dimension(:), allocatable :: detmat_sav, work, zagp_imag
    character(lchlen) :: path, scratchpath
    character(lchlen + 20) :: ranseedfilename
    real*8, dimension(:, :), allocatable :: kelind
    integer, dimension(:), allocatable :: jbrasymiesup
    logical yeseloc, dofast, yesprojm, yeswrite, acc_dyn, yesdodet, yesforce, loc_old
    ! %%%% memory required for reverse mode

    real*8, allocatable :: kelb(:, :), kelindb(:, :), rionb(:, :)&
            &, tabpipb(:), winvupb(:), winvdob(:), ainvb(:), ainvupbb(:), ainvdobb(:), psipb(:)&
            &, distb(:), rb(:), rmub(:), iond_cartb(:), winvb(:), winvjb(:), winvbarb(:)&
            &, winvjbarszb(:), winvjbarb(:), prefactorb(:), wpseudob(:), legendreb(:)&
            &, rmucosb(:), rmusinb(:), tmub(:), lapnum(:)&
            &, ivicb(:), tabpipsav(:), zout(:, :), zoutb(:, :), zoutsav(:, :), kelindlb(:, :)&
            &, rionlb(:, :), rionlab(:, :), forcedw(:, :), pulaydw(:, :), vjb(:), vjurb(:), duprb(:)&
            &, savetab(:), fbead(:, :), mass_ion(:, :), rionall(:, :, :), rion_write(:, :)&
            &, vjlb(:), vjurlb(:), duprlb(:), detmatb(:), detmat_cb(:), mu_cb(:), jasmatb(:), jasmatszb(:)&
            &, jasmatlb(:), jasmatszlb(:), mu_clb(:), projmlb(:), detmatlb(:), detmat_clb(:)&
            &, jasmat_cb(:), jasmat_clb(:), jasmatsz_cb(:), jasmatsz_clb(:), muj_cb(:), muj_clb(:)&
            &, eagp_pfafflb(:, :)

    real*8 grad1(3), grad1p(3), grad1m(3), grad2, grad2p, grad2m, jastrow&
            &, jastrowp, jastrowm, jee, rc(3), grad2num
    real*8, external :: jastrow_ee
    integer ispin

    !     definition functions
    real ran
    character(20) str
    character(60) name_tool

    !   AAA    Lines to be added just after all definitions of variables.

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

        !          Input the name of the file exactly as it is in /doc
        name_tool = 'testio'
        call help_online(name_tool)

        stop
    end if
    !    AAA   end lines to be added

    rank = 0
    nproc = 1

    call Initializeall
    !         indt4=0
    !         indt4j=0

    write (6, *) ' indt4 indt4j =', indt4, indt4j, firstmol, nelorb_at, nelorb_c, nmolfn
    !     test_aad = .true. ! remove some useless operations in compute_eloc_logpsi
    !     itestrfn=-6
    test_aad = .false.
    indwwjfn = 1
    indbarfn = 1
    iesrandoml = .false.
    pseudologic = .false.
    yesprojm = .true.
    write (6, *) ' yesprojm= ', yesprojm

    if (allocated(winvfn)) deallocate (winvfn, winvbarfn)

    if (yes_complex) then
        allocate (winvfn(2*nmolfn*(indt + ip4)*nel), winvbarfn(nmolfn*2*nel_mat))
    else
        allocate (winvfn(nmolfn*(indt + ip4)*nel), winvbarfn(nmolfn*nel_mat))
    end if
    nelsquare = nel*nel
    allocate (vpotsav_ee(nelsquare, 1))
    winvbarfn = 0.d0
    winvfn = 0.d0
    yeszagp = .true.
    yeszj = .true.
    yesforce = .true.
    yesdodet = .true.
    if (contractionj .eq. 0) then
        yesfastj = 0
    else
        yesfastj = 1
    end if

    !         firstmol=1
    !         nmolfn=nelorb_c

    allocate (kelindlb(3, nel), rionlb(3, nion))

    nelorbjmax = max(nelorbj, 1)
    neldomax = max(neldo, 1)
    indtmax = max(indt, 1)
    nshelljmax = max(nshellj, 1)

    allocate (kelind(3, nel))

    kelind(:, 1:nel) = kel(:, 1:nel)
    !  ALLOCATE all required for reverse mode.
    if (yesfast .ne. 0) then
        iscramax = max(iscramax, ipc*nmolfn*(nelup + ipf*nelorbh))
    else
        iscramax = max(iscramax, ipc*ipf*ipf*nelorbh*nelorbh)
    end if
    iscramax = max(iscramax, nelorbjh*nel)

    deallocate (psip)
    allocate (psip(iscramax))

    isdistp = 20*indt + 20
    nelorbjmax = max(nelorbj, 1)
    indtmax = max(indt, 1)
    nshelljmax = max(nshellj, 1)

    iscramaxb = isdistp
    iscramaxb = max(iscramaxb, (nion*(nion + 1))/2)

    !  ALLOCATE all required for reverse mode.
    if (iessz) then
        iscramaxb = max(iscramaxb, 2*max(nel, indt + 5) + 2*max(nelorbjh, 1) + 20*(indt + 1)&
                & + max(nelorbjh, 1)*(2*indt + 11))
    else
        iscramaxb = max(iscramaxb, 2*max(nel, indt + 5) + max(ipj*nelorbjh, 1) + 20*(indt + 1)&
                & + max(nelorbjh, 1)*(2*indt + 11))
    end if
    !  task9
    if (iessz) then
        isdistp = 27*(indt + 1)*nshelljmax + nelorbjmax*(indt + 6) &
                  + 1 + (7 + indt)*nelorbjmax + 2*max(nel, indt + 5) + nel + nion
    else
        isdistp = 27*(indt + 1)*nshelljmax + nelorbjmax*(indt + 6) &
                  + 1 + (ipj + 5 + indt)*nelorbjmax + 2*max(nel, indt + 5) + nel + nion
    end if
    iscramaxb = max(iscramaxb, isdistp)

    !   NB in TASK8_B reverse mode psip is not used

    iscramaxb = max(iscramaxb, nparshellmax) ! task2 nparshellmax defined  after read_pseudo -->

    iscramaxb = max(iscramaxb, 2*ipc*nelup_mat*nelorb_c) ! task5 upper bound

    !   NB in TASK8_B reverse mode psip is not used

    iscramaxb = max(iscramaxb, nparshellmax) ! task2 nparshellmax defined  after read_pseudo -->
    ! read_fort10 after read_pseudo

    iscramaxb = max(iscramaxb, 2*ipc*nelup_mat*nelorb_c) ! task5 upper bound

    iscramaxb = max(iscramaxb, 2*nelorbj_c*nel) ! task5 Jastrow
    iscramaxb = max(iscramaxb, (nion*(nion - 1))/2, ipc*nelup*nelup)
    iscramaxb = max(iscramaxb, ipc*(indt + 5)*nelorbh + nelorbh*(indt + 6) + 27*(indt + 1)*nshell + nelorbh) ! task10

    if (iscramaxb .lt. iscramax) iscramaxb = iscramax
    write (6, *) ' iscramaxb/iscramax  used ', iscramaxb, iscramax
!   iscramaxb=10000000

    if (iscramaxb .gt. iscramax) then
        iscramax = iscramaxb
        deallocate (psip)
        allocate (psip(iscramax))
    end if

    movedions = 1

    write (6, *) ' Optimal block size for dgetri ='&
            &, nelup*ILAENV(1, 'DGETRI', ' ', nelup, -1, -1, -1), nbdgetri

    allocate (kelb(3, nel*(indt + 1)), kelindb(3, nel), rionb(3, nion)&
            &, tabpipb(nel*(indt + 4)), winvupb(ipc*nelup*(indt + 4)), winvdob(ipc*neldomax*(indt + 4))&
            &, ainvb(ipc*nelup_mat*nelup_mat), ainvupbb(ipc*nelup*nelorbh), ainvdobb(ipc*neldomax*nelorbh)&
            &, psipb(iscramaxb), distb(nel*nion), rb(nion*(indt + 1)), rmub(3*nion*(indt + 1))&
            &, iond_cartb(3*nion*nion), winvb(ipc*nelorb*nel), winvjb(nelorbjmax*nel)&
            &, winvbarb(ipf*ipc*nelorbh*nel_mat), winvjbarb(ipj*nelorbjh*nel + 1)&
            &, winvjbarszb(nelorbjh*nel + 1)&
            &, prefactorb((indt - istart + 1)*nel), wpseudob(2*lmax)&
            &, legendreb((lmax - 1)*nintpseudo), rmucosb(3*nion*(indt + 1))&
            &, rmusinb(3*nion*(indt + 1)), tmub(nel*indt), ivicb(3*indtmax*nel)&
            &, vjb(size(vj)), vjurb(size(vjur)), duprb(size(dupr))&
            &, vjlb(size(vj)), vjurlb(size(vjur)), duprlb(size(dupr))&
            &, tabpipsav(nel*(indt + 4)), forcedw(3, nion), pulaydw(3, nion))
    if (npar_eagp .gt. 0) then
        allocate (eagp_pfafflb(ipc*ndiff, ndiff))
    else
        allocate (eagp_pfafflb(1, 1))
    end if
    if (yesfast .eq. 0) then
        write (6, *) ' nelcolh inside =', nelorbh*ipf, nelcolh
        allocate (detmatb(ipc*ipf*nelorbh*nelcolh))
        allocate (detmatlb(ipc*ipf*nelorbh*nelcolh))
        allocate (detmat_cb(1), mu_cb(1))
        allocate (detmat_clb(1), mu_clb(1))
    else
        allocate (detmatb(ipc*ipf*nelorbh))
        allocate (detmatlb(ipc*ipf*nelorbh))
        allocate (detmat_cb(ipc*nelorb_c*nelcol_c), mu_cb(ipc*ipf*nelorbh*nelorb_c))
        allocate (detmat_clb(ipc*nelorb_c*nelcol_c), mu_clb(ipc*ipf*nelorbh*nelorb_c))
    end if
    if (npar_eagp .gt. 0) then
        allocate (eagp_pfaffb(ipc*ndiff, ndiff))
    else
        allocate (eagp_pfaffb(1, 1))
    end if
    dim_jasmat = 1
    if (nelorbjh .gt. 0) then
        if (yesfastj .eq. 0) then
            if (yes_sparse) then
                allocate (jasmatb(nnozeroj))
                allocate (jasmatlb(nnozeroj))
                dim_jasmat = nnozeroj
            else
                allocate (jasmatb(ipj*ipj*nelorbjh*nelorbjh))
                allocate (jasmatlb(ipj*ipj*nelorbjh*nelorbjh))
                dim_jasmat = ipj*ipj*nelorbjh*nelorbjh
            end if
            if (iessz) then
                allocate (jasmatszb(nelorbjh*nelorbjh))
                allocate (jasmatszlb(nelorbjh*nelorbjh))
            else
                allocate (jasmatszb(1))
                allocate (jasmatszlb(1))
            end if
            allocate (muj_cb(1))
            allocate (muj_clb(1))
            allocate (jasmat_cb(1))
            allocate (jasmat_clb(1))
            allocate (jasmatsz_cb(1))
            allocate (jasmatsz_clb(1))
        else
            allocate (jasmatb(1))
            allocate (jasmatlb(1))
            allocate (jasmatszb(1))
            allocate (jasmatszlb(1))
            allocate (muj_cb(nelorbj_c*nelorbjh))
            allocate (muj_clb(nelorbj_c*nelorbjh))
            allocate (jasmat_cb(ipj*ipj*nelorbj_c*nelorbj_c))
            allocate (jasmat_clb(ipj*ipj*nelorbj_c*nelorbj_c))
            if (iessz) then
                allocate (jasmatsz_cb(nelorbj_c*nelorbj_c))
                allocate (jasmatsz_clb(nelorbj_c*nelorbj_c))
            else
                allocate (jasmatsz_cb(1))
                allocate (jasmatsz_clb(1))
            end if
        end if
    else
        allocate (jasmatszb(1), jasmatb(1))
        allocate (jasmatszlb(1), jasmatlb(1))
        allocate (muj_cb(1))
        allocate (muj_clb(1))
        allocate (jasmat_cb(1))
        allocate (jasmat_clb(1))
        allocate (jasmatsz_cb(1))
        allocate (jasmatsz_clb(1))
    end if

    cellscaleo = cellscale

    iflagerr = 0

    write (6, *) ' Initial electron coordinates '
    do i = 1, nel
        write (6, *) i, kel(:, i)
    end do

    write (6, *) ' Initial ion  coordinates '
    do i = 1, nion
        write (6, *) i, rion(:, i)
    end do
    psisn(1) = 0.d0

    movedions = 1
    psisn(1) = 0.d0
    singdet = .true.
    !     Reset to zero all output for a more stringent test.
    winvfn = 0.d0
    winvbarfn = 0.d0
    winv = 0.d0
    winvbar = 0.d0
    if (nelorbj .gt. 0) winvj = 0.d0
    jastrowall_ee = 0.d0
    jastrowall_ei = 0.d0
    ainv = 0.d0
    ainvup = 0.d0
    ainvdo = 0.d0
    winvup = 0.d0
    winvdo = 0.d0
    tabpip = 0.d0
    psidetln = 0.d0
    singdet = .true.

!   test_aad=.true.
!   iflagnorm=3
    write (6, *) ' indshell before = ', indshell_tab(1)

    write (6, *) ' Test local energy ', ipc, Lbox, nel
    if (contraction .ne. 0) then
        write (6, *) ' detmat_c input compute  ', sum(abs(detmat_c(1:ipc*nelorb_c**2)))
        write (6, *) ' input mu_c compute', sum(abs(mu_c(1:ipf*ipc*nelorbh, 1:nelorb_c)))
    end if
    if (yesfast .eq. 0) write (6, *) ' input detmat=', sum(abs(detmat(1:ipc*ipf*nelorbh*nelcolh)))
    write (6, *) ' nelorb_at inside =', nelorb_at, 2*nelorbh

    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, contractionj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
            &, nintpseudo, prefactor, rcutoff, parshell                            &
            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
            &, jpseudo, pseudolocal, istart, costz, costz3         &
            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
            &, psidetln, 1, nelorbh, nelorbjh                         &
            &, niesd&
            &, iond_cart, mu_c, detmat_c, projm, yesprojm, nelorb_c, firstmol, nmolfn, yesfast, eloc, logpsi&
            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

    !       call ratioflip
    !     do i =1,nel
    !      kel--> kel_new
    !      cambiare nelup e neldo secondo quanto  deallocare e riallocare con la dimensione giusta
    !      ainvup ainvdo e riallocare con la dimensione giusta

    write (6, *) ' after new ', ipc, nelorbh, indt4, indt4j, size(winv)
    write (6, *) ' kel =', sum(kel(:, 1:nel*(indt + 1)))
    write (6, *) ' winv =', sum(winv(1:ipc*nelorbh*(indt4 + 1)*nel))
    write (6, *) ' winvj =', sum(winvj(1:nelorbjh*(indt4j + 1)*nel))
    write (6, *) ' jasmat =', sum(jasmat(1:dim_jasmat))
    write (6, *) ' winvjbar =', sum(winvjbar(1:nelorbjh*nel))
    write (6, *) ' winvbar =', sum(abs(winvbar(1:ipf*ipc*nelorbh*nel_mat)))
    write (6, *) ' ainv =', sum(abs(ainv(1:ipc*nelup_mat*nelup_mat)))
    write (6, *) ' ainvup =', sum(ainvup(1:ipc*nelup*nelorbh))
    write (6, *) ' ainvdo =', sum(ainvdo(1:ipc*neldo*nelorbh))
    write (6, *) ' winvup =', sum(winvup(1:ipc*nelup*(indt + 4)))
    write (6, *) ' |winvup| =', sum(abs(winvup(1:ipc*nelup*(indt + 4))))
    write (6, *) ' winvdo =', sum(winvup(1:ipc*neldo*(indt + 4)))
    write (6, *) ' |winvdo| =', sum(abs(winvup(1:ipc*neldo*(indt + 4))))
    write (6, *) ' tabpip =', sum(abs(tabpip(1:nel*(indt + 4))))
    write (6, *) ' tabpip lap =', sum(abs(tabpip(nel*indt + 1:nel*(indt + 4))))
    write (6, *) ' tabpip lap up=', sum(abs(tabpip(nel*indt + 1:nel*indt + nelup)))
    write (6, *) ' tabpip pseudo =', sum(abs(tabpip(1:nel*indt)))
    write (6, *) ' tabpip pseudo up =', sum(abs(tabpip(1:nelup)))
    write (6, *) ' jastrowall_ee=', sum(jastrowall_ee(:, :, :, :))
    write (6, *) ' jastrowall_ee VMC =', sum(jastrowall_ee(:, :, 0, :))
    write (6, *) ' jastrowall_ei=', sum(jastrowall_ei(:, :, :))
    write (6, *) ' winvfn=', sum(abs(winvfn(:)))
    write (6, *) ' winvbarfn=', sum(abs(winvbarfn(:)))
    write (6, *) ' tmu =', sum(abs(tmu(:)))

    write (6, *) ' New Local energy =', eloc(1:ipc)
    write (6, *) ' New Log wave function =', logpsi(1:ipc)
    write (6, *) ' New  vpot_reg =', vpotreg(1:2, 1)
    write (6, *) ' New psidetln =', psidetln(1)

    test_aad = .true.
    itestrfn = 2

    timep = cclock()
    singdet = .true.

    aadinp(1) = 1.d0
    aadinp(2) = 0.d0
    !      checking two body jastrow
    write (6, *) ' Input kel =', kelind(1:3, 1)

    ispin = -1
    kel(:, 1) = kelind(:, 1)

    go to 23
    write (6, *) ' Gradient Laplacian two body term '

    if (LBox .gt. 0.d0) then
        call jastrowgrad_pbc(kel, vj, iesd, jastrow, grad1, grad2, ispin, 1.d0)
    else
        call jastrowgrad(kel, vj, iesd, jastrow, grad1, grad2, ispin)
    end if
    rc(:) = kelind(:, 1)

    if (LBox .gt. 0.d0) then
        call CartesianToCrystal(rc, 1)
        !       rc(1:3)=cellpi(1:3)*dsin(rc(1:3)/cellpi(1:3))
        do k = 1, 3
            rc(k) = map(rc(k), cellscale(k))
        end do
    end if

    jee = jastrow_ee(rc, vj, iesd, ispin)

    write (6, *) 'Consistency with jastrow_ee ', jastrow, jee
    grad2num = 0.d0
    do k = 1, 3
        kel(:, 1) = kelind(:, 1)
        kel(k, 1) = kelind(k, 1) + epsder
        if (LBox .gt. 0.d0) then
            call jastrowgrad_pbc(kel, vj, iesd, jastrowp, grad1p, grad2p, ispin, 1.d0)
        else
            call jastrowgrad(kel, vj, iesd, jastrowp, grad1p, grad2p, ispin)
        end if
        kel(:, 1) = kelind(:, 1)
        kel(k, 1) = kelind(k, 1) - epsder
        if (LBox .gt. 0.d0) then
            call jastrowgrad_pbc(kel, vj, iesd, jastrowm, grad1m, grad2m, ispin, 1.d0)
        else
            call jastrowgrad(kel, vj, iesd, jastrowm, grad1m, grad2m, ispin)
        end if
        write (6, *) k, (jastrowp - jastrowm)/2.d0/epsder, grad1(k)
        grad2num = grad2num + (grad1p(k) - grad1m(k))/2.d0/epsder
        kel(:, 1) = kelind(:, 1)
    end do

    write (6, *) ' Consistency Laplacian = ', grad2, grad2num

    !      test gradient and laplacian first electron
    allocate (zout(ipc*nelorb, 0:indt + 4))
    allocate (zoutb(ipc*nelorb, 0:indt + 4))
    allocate (zoutsav(ipc*nelorb, 0:indt + 4))
    allocate (lapnum(ipc*nelorb))

    zout = 0.d0
    zoutb = 0.d0
    zoutsav = 0.d0
    lapnum = 0.d0
    kel(:, 1) = kelind(:, 1)

    call upnewwf(indt, 0, 0, 0, nshell, ioptorb, ioccup, kel&
            &, nel, r, rmu, dupr, zetar, rion, psip, zoutsav(1, 0), nelorb, nion, kion&
            &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
            &, indpar_tab, indorb_tab, indshell_tab, .true.)

    do k = 1, 3
        kel(:, 1) = kelind(:, 1)
        kel(k, 1) = kelind(k, 1) + epsder

        !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
        call upnewwf(indt, 0, 0, 0, nshell, ioptorb, ioccup, kel&
                &, nel, r, rmu, dupr, zetar, rion, psip, zout, nelorb, nion, kion&
                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                &, indpar_tab, indorb_tab, indshell_tab, .true.)

        kel(k, 1) = kelind(k, 1) - epsder

        write (6, *) ' indpar_tab  = ', indpar_tab(1)

        !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
        call upnewwf(indt, 0, 0, 0, nshell, ioptorb, ioccup, kel&
                &, nel, r, rmu, dupr, zetar, rion, psip, zoutb, nelorb, nion, kion&
                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                &, indpar_tab, indorb_tab, indshell_tab, .true.)

        kel(k, 1) = kelind(k, 1) ! restore value

        do j = 1, ipc*nelorb
            if (abs((zout(j, 0) - zoutb(j, 0))/2.d0/epsder - zoutsav(j, indt + k)) .gt. 1d-6) then
                write (6, *) ' ERROR gradient ', k, j

            end if
            write (6, *) ' Numerical/Anal grad orbital  ', k, j, (zout(j, 0) - zoutb(j, 0))/2.d0/epsder, zoutsav(j, indt + k)
            lapnum(j) = lapnum(j) + (zout(j, indt + k) - zoutb(j, indt + k))/2.d0/epsder
        end do
        kel(:, 1) = kelind(:, 1)
    end do

    write (6, *) ' Laplacian '
    do j = 1, ipc*nelorb
        if (abs(lapnum(j) - zoutsav(j, indt + 4)) .gt. 1d-6) then
            write (6, *) ' ERROR  lap '
        end if
        write (6, *) j, lapnum(j), zoutsav(j, indt + 4)
    end do

    !       go to 24

    go to 23

    !go to 24

    do k = 1, 3
        do i = 1, nel
            ind = 20*(drand1() - 0.5d0)
            kelind(k, i) = kelind(k, i) + cellscale(k)*ind
        end do
        !         do i=1,nion
        !         ind=20*(drand1()-0.5d0)
        !         rion(k,i)=rion(k,i)+cellscale(k)*ind
        !         enddo
    end do

    write (6, *) ' Compute periodized energy '

    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
            &, nintpseudo, prefactor, rcutoff, parshell                            &
            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
            &, jpseudo, pseudolocal, istart, costz, costz3         &
            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
            &, psidetln, 1, nelorbh, nelorbjh                         &
            &, niesd&
            &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, eloc, logpsi&
            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

    write (6, *) ' New Local energy =', eloc
    write (6, *) ' New Log wave function =', logpsi

    write (6, *) ' test upnewwf ', epsder
    !         allocate(zout(nelorb,0:indt+4))
    !         allocate(zoutb(nelorb,0:indt+4))
    !         allocate(zoutsav(nelorb,0:indt+4))
    !         zoutsav=1.d0
    !         zoutsav(1:nelorb,indt+1:indt+4)=1.d0

    allocate (zout(nelorb, 1))
    allocate (zoutb(nelorb, 1))
    allocate (zoutsav(nelorb, 1))

    zoutsav = 1.d0

    zoutb = zoutsav

    kelb = 0.d0
    rionb = 0.d0

    !     call upnewwf_b(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
    call upnewwf_b(0, 0, 0, 1, nshell, ioptorb, ioccup, kel&
            &, kelb, nel, r, rb, rmu, rmub, dupr, zetar, rion, rionb, psip, psipb, zout, zoutb&
            &, nelorb, nion, kion, iflagnorm, cnorm, LBox, rmucos, rmucosb, rmusin, rmusinb&
            &, 1d-9, cellscale, cellscaleb)

    zoutb = zoutsav

    !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
    !    &,nel,r,rmu,dupr,zetar,rion,psip,zout,nelorb,nion,kion&
    !    &,iflagnorm,cnorm,LBox,rmucos,rmusin,1d-9)

    !        do j=0,indt+4
    !         write(6,*) ' Component =',j
    !         do i=1,nelorb
    !         write(6,*) i,zout(i,j)
    !         enddo
    !        enddo

    !        cost0=sum(zoutb(:,:)*zout(:,:))

    do k = 1, 3
        kel(:, 1) = kelind(:, 1)
        kel(k, 1) = kelind(k, 1) + epsder

        !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
        call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, kel&
                &, nel, r, rmu, dupr, zetar, rion, psip, zout, nelorb, nion, kion&
                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                &, indpar_tab, indorb_tab, indshell_tab, .true.)

        costp = sum(zoutb(:, :)*zout(:, :))

        kel(k, 1) = kelind(k, 1) - epsder

        !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
        call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, kel&
                &, nel, r, rmu, dupr, zetar, rion, psip, zout, nelorb, nion, kion&
                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                &, indpar_tab, indorb_tab, indshell_tab, .true.)

        costm = sum(zoutb(:, :)*zout(:, :))

        write (6, *) ' Numerical der kel ', k, (costp - costm)/2.d0/epsder, kelb(k, 1)

        kel(:, 1) = kelind(:, 1)

        rionsav = rion

        do j = 1, nion

            rion(k, j) = rionsav(k, j) + epsder

            call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, kel &
                    !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
                    &, nel, r, rmu, dupr, zetar, rion, psip, zout, nelorb, nion, kion&
                    &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                    &, indpar_tab, indorb_tab, indshell_tab, .true.)

            costp = sum(zoutb(:, :)*zout(:, :))

            rion(k, j) = rionsav(k, j) - epsder

            call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, kel &
                    !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
                    &, nel, r, rmu, dupr, zetar, rion, psip, zout, nelorb, nion, kion&
                    &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                    &, indpar_tab, indorb_tab, indshell_tab, .true.)

            costm = sum(zoutb(:, :)*zout(:, :))
            rion(k, j) = rionsav(k, j)

            write (6, *) ' Numerical/Anal  der rion  ', k, j&
                    &, (costp - costm)/2.d0/epsder, rionb(k, j)

        end do

    end do

    if (iespbc) then

        cellscaleo = cellscale

        celldm(1) = cellscaleo(1) + epsder
        celldm(2) = cellscaleo(2)/(cellscaleo(1) + epsder)
        celldm(3) = cellscaleo(3)/(cellscaleo(1) + epsder)

        call InitCell(nion, nel, yes_complex)

        call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, kel &
                !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
                &, nel, r, rmu, dupr, zetar, rion, psip, zout, nelorb, nion, kion&
                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                &, indpar_tab, indorb_tab, indshell_tab, .true.)

        costp = sum(zoutb(:, :)*zout(:, :))

        celldm(1) = cellscaleo(1) - epsder
        celldm(2) = cellscaleo(2)/(cellscaleo(1) - epsder)
        celldm(3) = cellscaleo(3)/(cellscaleo(1) - epsder)

        call InitCell(nion, nel, yes_complex)

        call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, kel &
                !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
                &, nel, r, rmu, dupr, zetar, rion, psip, zout, nelorb, nion, kion&
                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                &, indpar_tab, indorb_tab, indshell_tab, .true.)

        costm = sum(zoutb(:, :)*zout(:, :))

        write (6, *) ' First der cell =', (costp - costm)/2.d0/epsder, cellscaleb(1)

        celldm(1) = cellscaleo(1)
        celldm(2) = (cellscaleo(2) + epsder)/cellscaleo(1)
        celldm(3) = cellscaleo(3)/cellscaleo(1)

        call InitCell(nion, nel, yes_complex)

        !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
        call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, kel&
                &, nel, r, rmu, dupr, zetar, rion, psip, zout, nelorb, nion, kion&
                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                &, indpar_tab, indorb_tab, indshell_tab, .true.)

        costp = sum(zoutb(:, :)*zout(:, :))

        celldm(1) = cellscaleo(1)
        celldm(2) = (cellscaleo(2) - epsder)/cellscaleo(1)
        celldm(3) = cellscaleo(3)/cellscaleo(1)

        call InitCell(nion, nel, yes_complex)

        !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
        call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, kel&
                &, nel, r, rmu, dupr, zetar, rion, psip, zout, nelorb, nion, kion&
                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                &, indpar_tab, indorb_tab, indshell_tab, .true.)

        costm = sum(zoutb(:, :)*zout(:, :))

        write (6, *) ' Second der cell =', (costp - costm)/2.d0/epsder, cellscaleb(2)

        celldm(1) = cellscaleo(1)
        celldm(2) = cellscaleo(2)/cellscaleo(1)
        celldm(3) = (cellscaleo(3) + epsder)/cellscaleo(1)

        call InitCell(nion, nel, yes_complex)

        !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
        call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, kel&
                &, nel, r, rmu, dupr, zetar, rion, psip, zout, nelorb, nion, kion&
                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                &, indpar_tab, indorb_tab, indshell_tab, .true.)

        costp = sum(zoutb(:, :)*zout(:, :))

        celldm(1) = cellscaleo(1)
        celldm(2) = cellscaleo(2)/cellscaleo(1)
        celldm(3) = (cellscaleo(3) - epsder)/cellscaleo(1)

        call InitCell(nion, nel, yes_complex)

        !     call upnewwf(indt,0,indt,0,nshell,ioptorb,ioccup,kel&
        call upnewwf(0, 0, 0, 1, nshell, ioptorb, ioccup, kel&
                &, nel, r, rmu, dupr, zetar, rion, psip, zout, nelorb, nion, kion&
                &, iflagnorm, cnorm, LBox, rmucos, rmusin, 1d-9&
                &, indpar_tab, indorb_tab, indshell_tab, .true.)

        costm = sum(zoutb(:, :)*zout(:, :))

        write (6, *) ' Third  der cell =', (costp - costm)/2.d0/epsder, cellscaleb(3)

    end if

23  continue

    rionsav = rion
    elocb = aadinp
    logpsib = 0.d0

    write (6, *) ' epsder =', epsder
    write (6, *) ' Num. der eloc logpsi / Analityc AD '
    !        go to 25
    write (6, *) ' before  compute', iscramax, iscramaxb, indshell_tab(1), logpsi(1), logpsib(1:ipc), eloc(1), elocb(1:ipc)
    timep = cclock()

    !         computing reverse routine
    call COMPUTE_ELOC_LOGPSI_B(indt, nelorb, nelup, neldo, tabpip, &
            &  tabpipb, tabpipsav, kelind, kelindb, kel, kelb, winv, winvb, indt4, winvup, winvupb, winvdo, &
            &  winvdob, ainv, ainvb, ainvup, ainvupbb, ainvdo, ainvdobb, psip, psipb, ipsip, &
            &   psisn, iesdr, vj, vjb, size(vj), dupr, duprb, size(dupr), zetar, rion, rionb, dist, distb, &
            &  ioccup, ioccdo, ioptorb, nshell, nshelldo, ivic, ivicb, vpot, tmu, tmub, &
            &  nion, r, rb, rmu, rmub, kion, iond, iond_cartb, winvj, winvjb, indt4j, ioccj, kionj, &
            &  vjur, vjurb, size(vjur), nelorbj, ioptorbj, nshellj, winvbar, winvbarb, detmat, detmatb, eagp_pfaffb, winvjbar, &
            & winvjbarb, winvjbarsz, winvjbarszb, jasmat, jasmatb, muj_c, muj_cb, jasmat_c, jasmat_cb&
            &, jasmatsz, jasmatszb, jasmatsz_c, jasmatsz_cb, nelorbj_c, yesfastj, iessz, cnorm, iflagerr, &
            &  npsa, lmax, nintpseudo, prefactor, prefactorb, rcutoff, parshell, nparpshell, &
            &  kindion, pshell, wpseudo, wpseudob, legendre, legendreb, versor, wintpseudo &
            &, jpseudo, pseudolocal, istart, costz, costz3, angle, indtm, lbox, &
            & rmucos, rmucosb, rmusin, rmusinb, kappa, vpotreg, cutreg, psidetln, 1, &
            & nelorbh, nelorbjh, niesd, iond_cart, mu_c, mu_cb, detmat_c, detmat_cb, nelorb_c, &
            & firstmol, nmolfn, yesfast, yeszagp, yesdodet, yeszj, yesforce, eloc, elocb, logpsi, logpsib, nelorbjmax, &
            & neldomax, indtmax, nshelljmax, cellscaleb, sr2b, iflagnorm&
            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab&
            &, adr_nion, ind_nion, adrj_nion, indj_nion, .false.)
    !24         continue

    write (6, *) ' Time compute_b local energy =', cclock() - timep

25  continue
    elocb = 0.d0
    logpsib = aadinp
    timep = cclock()

    call COMPUTE_ELOC_LOGPSI_B(indt, nelorb, nelup, neldo, tabpip, &
            &  tabpipb, tabpipsav, kelind, kelindlb, kel, kelb, winv, winvb, indt4, winvup, winvupb, winvdo, &
            &  winvdob, ainv, ainvb, ainvup, ainvupbb, ainvdo, ainvdobb, psip, psipb, ipsip, &
            &   psisn, iesdr, vj, vjlb, size(vj), dupr, duprlb, size(dupr), zetar, rion, rionlb, dist, distb, &
            &  ioccup, ioccdo, ioptorb, nshell, nshelldo, ivic, ivicb, vpot, tmu, tmub, &
            &  nion, r, rb, rmu, rmub, kion, iond, iond_cartb, winvj, winvjb, indt4j, ioccj, kionj, &
            &  vjur, vjurlb, size(vjur), nelorbj, ioptorbj, nshellj, winvbar, winvbarb, detmat, detmatlb, eagp_pfafflb, winvjbar, &
            & winvjbarb, winvjbarsz, winvjbarszb, jasmat, jasmatlb, muj_c, muj_clb, jasmat_c, jasmat_clb&
            &, jasmatsz, jasmatszlb, jasmatsz_c, jasmatsz_clb, nelorbj_c, yesfastj, iessz, cnorm, iflagerr, &
            &  npsa, lmax, nintpseudo, prefactor, prefactorb, rcutoff, parshell, nparpshell, &
            &  kindion, pshell, wpseudo, wpseudob, legendre, legendreb, versor, wintpseudo &
            &, jpseudo, pseudolocal, istart, costz, costz3, angle, indtm, lbox, &
            & rmucos, rmucosb, rmusin, rmusinb, kappa, vpotreg, cutreg, psidetln, 1, &
            & nelorbh, nelorbjh, niesd, iond_cart, mu_c, mu_clb, detmat_c, detmat_clb, nelorb_c, &
            & firstmol, nmolfn, yesfast, yeszagp, yesdodet, yeszj, yesforce, eloc, elocb, logpsi, logpsib, nelorbjmax, &
            & neldomax, indtmax, nshelljmax, cellscalelb, sr2lb, iflagnorm&
            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab&
            &, adr_nion, ind_nion, adrj_nion, indj_nion, .false.)
    write (6, *) ' Time compute_b wf =', cclock() - timep
    warp = .false.
    press = cellscaleb
    press_pulay = cellscalelb

    call updatedwarp(nel, nion, kelind, rion, rmu, r, kelindb, kelindlb&
            &, rionb, rionlb, press, press_pulay, forcedw, pulaydw&
            &, iespbc, warp, powerwarp, atom_number, warpmat, niong)

    write (6, *) ' Differential warp forces / pulay '
    do i = 1, nion
        do k = 1, 3
            write (6, *) k, i, forcedw(k, i), pulaydw(k, i)
        end do
    end do
    if (iespbc) then
        write (6, *) ' Differential Warp cell derivatives '
        do k = 1, 3
            write (6, *) k, press(k), press_pulay(k)
        end do
    end if

    !        goto 1133 ! only cell derivatives.

    if (yesforce) then

        do k = 1, 3
            write (6, *) ' Component ', k

            write (6, *) ' # el '
            do j = 1, nel
                kel_sav = kelind(k, j)
                kelind(k, j) = kel_sav + epsder

                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)
                kelind(k, j) = kel_sav - epsder

                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                !      elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
                !      logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder

                elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
                logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder

                if (abs(elocp(1) - kelindb(k, j)) .gt. 1d-4) then
                    write (6, *) ' ERROR in el. force local energy el ', k, j, elocp(1), kelindb(k, j)
                end if
                if (abs(logpsip(1) - kelindlb(k, j)) .gt. 1d-4) then
                    write (6, *) ' ERROR in el. force wf el ', k, j, logpsip(1), kelindlb(k, j)
                end if

                write (6, *) j, elocp(1), kelindb(k, j), logpsip(1), kelindlb(k, j)

                kelind(k, j) = kel_sav

            end do

            write (6, *) ' # ion '

            do j = 1, nion

                rion(k, j) = rionsav(k, j) + epsder

                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                rion(k, j) = rionsav(k, j) - epsder

                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                rion(k, j) = rionsav(k, j)

                !      elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
                !      logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
                elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
                logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder
                if (abs(elocp(1) - rionb(k, j)) .gt. 1d-4) then
                    write (6, *) ' ERROR in force  ion local energy ', k, j, elocp(1), rionb(k, j)
                end if
                if (abs(logpsip(1) - rionlb(k, j)) .gt. 1d-4) then
                    write (6, *) ' ERROR in force ion wf ', k, j, logpsip(1), rionlb(k, j)
                end if

                write (6, *) j, elocp(1), rionb(k, j), logpsip(1), rionlb(k, j)

            end do

        end do

    end if

    write (6, *) ' Derivative One/Two body Jastrow ', iesd

    do j = 1, iesd
        kel_sav = vj(j)
        vj(j) = kel_sav + epsder

        call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                &, nintpseudo, prefactor, rcutoff, parshell                            &
                &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                &, jpseudo, pseudolocal, istart, costz, costz3         &
                &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                &, psidetln, 1, nelorbh, nelorbjh                         &
                &, niesd&
                &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)
        vj(j) = kel_sav - epsder

        call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                &, nintpseudo, prefactor, rcutoff, parshell                            &
                &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                &, jpseudo, pseudolocal, istart, costz, costz3         &
                &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                &, psidetln, 1, nelorbh, nelorbjh                         &
                &, niesd&
                &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

        !      elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
        !      logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
        elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
        logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder
        if (abs(elocp(1) - vjb(j)) .gt. 1d-4) then
            write (6, *) ' ERROR in Two body Jas  local energy ', j, elocp(1), vjb(j)
        end if
        if (abs(logpsip(1) - vjlb(j)) .gt. 1d-4) then
            write (6, *) ' ERROR in Two body Jas  wf ', j, logpsip(1), vjlb(j)
        end if
        write (6, *) j, elocp(1), logpsip(1), vjb(j), vjlb(j)

        vj(j) = kel_sav

    end do

    if (yeszagp) then
        write (6, *) ' Z determinant '

        do j = 1, nshell

            kel_sav = dupr(j)

            dupr(j) = kel_sav + epsder

            call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                    &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                    &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                    &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                    &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                    &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                    &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                    &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                    &, nintpseudo, prefactor, rcutoff, parshell                            &
                    &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                    &, jpseudo, pseudolocal, istart, costz, costz3         &
                    &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                    &, psidetln, 1, nelorbh, nelorbjh                         &
                    &, niesd&
                    &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                    &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

            dupr(j) = kel_sav - epsder

            call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                    &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                    &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                    &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                    &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                    &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                    &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                    &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                    &, nintpseudo, prefactor, rcutoff, parshell                            &
                    &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                    &, jpseudo, pseudolocal, istart, costz, costz3         &
                    &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                    &, psidetln, 1, nelorbh, nelorbjh                         &
                    &, niesd&
                    &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                    &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

            !      elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
            !      logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
            elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
            logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder
            if (abs(elocp(1) - duprb(j)) .gt. 1d-4) then
                write (6, *) ' ERROR in Z det local en. ', j, elocp(1), duprb(j)
            end if
            if (abs(logpsip(1) - duprlb(j)) .gt. 1d-4) then
                write (6, *) ' ERROR in Z det  wf ', j, logpsip(1), duprlb(j)
            end if
            write (6, *) j, elocp(1), logpsip(1), duprb(j), duprlb(j)

            dupr(j) = kel_sav

            !    stop

        end do

    end if

    if (yeszj) then

        write (6, *) ' Z Jastrow '

        do j = 1, size(vjur)
            kel_sav = vjur(j)
            vjur(j) = kel_sav + epsder

            call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                    &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                    &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                    &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                    &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                    &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                    &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                    &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                    &, nintpseudo, prefactor, rcutoff, parshell                            &
                    &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                    &, jpseudo, pseudolocal, istart, costz, costz3         &
                    &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                    &, psidetln, 1, nelorbh, nelorbjh                         &
                    &, niesd&
                    &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                    &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)
            vjur(j) = kel_sav - epsder

            call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                    &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                    &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                    &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                    &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                    &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                    &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                    &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                    &, nintpseudo, prefactor, rcutoff, parshell                            &
                    &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                    &, jpseudo, pseudolocal, istart, costz, costz3         &
                    &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                    &, psidetln, 1, nelorbh, nelorbjh                         &
                    &, niesd&
                    &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                    &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

            !      elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
            !      logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
            elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
            logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder
            if (abs(elocp(1) - vjurb(j)) .gt. 1d-4) then
                write (6, *) ' ERROR in Z Jas local en. ', j, elocp(1), vjurb(j)
            end if
            if (abs(logpsip(1) - vjurlb(j)) .gt. 1d-4) then
                write (6, *) ' ERROR in Z Jas  wf ', j, logpsip(1), vjurlb(j)
            end if
            write (6, *) j, elocp(1), logpsip(1), vjurb(j), vjurlb(j)

            vjur(j) = kel_sav

        end do

    end if
    !        go to 1133 ! (cell derivatives)

    if (iespbc) then
        go to 1133
    else
        go to 24 ! Compiler produce strange executable if they  do not  deallocate
        ! and also one can trace memory  conflicts  in case.
    end if
    if (iesfree .ne. 0) write (6, *) ' Matrix jasmat '

    if (yesfastj .eq. 0 .and. iesfree .ne. 0) then

        do j = 1, ipj*nelorbjh
            do k = j, ipj*nelorbjh

                kel_sav = jasmat(ipj*nelorbjh*(k - 1) + j)

                jasmat(ipj*nelorbjh*(k - 1) + j) = kel_sav + epsder

                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)
                jasmat(ipj*nelorbjh*(k - 1) + j) = kel_sav - epsder

                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                !       elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
                !       logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
                elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
                logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder

                if (abs(elocp(1) - jasmatb(ipj*nelorbjh*(k - 1) + j)) .gt. 1d-4) then
                    write (6, *) ' ERROR in Jas local en. ', j, k, elocp(1), jasmatb(ipj*nelorbjh*(k - 1) + j)
                end if
                if (abs(logpsip(1) - jasmatlb(ipj*nelorbjh*(k - 1) + j)) .gt. 1d-4) then
                    write (6, *) ' ERROR in Jas  wf ', j, logpsip(1), jasmatlb(ipj*nelorbjh*(k - 1) + j)
                end if
                jasmat(ipj*nelorbjh*(k - 1) + j) = kel_sav
                write (6, *) j, k, elocp(1), logpsip(1), jasmatb(ipj*nelorbjh*(k - 1) + j), jasmatlb(ipj*nelorbjh*(k - 1) + j)

            end do
        end do
    elseif (iesfree .ne. 0) then
        write (6, *) ' Contracted matrix jasmat_c '

        do j = 1, ipj*nelorbj_c
            do k = j, ipj*nelorbj_c

                kel_sav = jasmat_c(ipj*nelorbj_c*(k - 1) + j)

                jasmat_c(ipj*nelorbj_c*(k - 1) + j) = kel_sav + epsder

                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                jasmat_c(ipj*nelorbj_c*(k - 1) + j) = kel_sav - epsder

                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)
                jasmat_c(ipj*nelorbj_c*(k - 1) + j) = kel_sav
                !       elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
                !       logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
                elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
                logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder

                if (abs(elocp(1) - jasmat_cb(ipj*nelorbj_c*(k - 1) + j)) .gt. 1d-4) then
                    write (6, *) ' ERROR in Jas contr. local en. ', j, k, elocp(1), jasmat_cb(ipj*nelorbj_c*(k - 1) + j)
                end if
                if (abs(logpsip(1) - jasmat_clb(ipj*nelorbj_c*(k - 1) + j)) .gt. 1d-4) then
                    write (6, *) ' ERROR in Jas  contr wf ', j, logpsip(1), jasmat_clb(ipj*nelorbj_c*(k - 1) + j)
                end if

                write (6, *) j, k, elocp(1), logpsip(1) &
                    , jasmat_cb(ipj*nelorbj_c*(k - 1) + j), jasmat_clb(ipj*nelorbj_c*(k - 1) + j)
            end do
        end do

        write (6, *) '  Derivatives contracted muj_c '

        do j = 1, nelorbj_c
            do k = i, nelorbjh

                kel_sav = muj_c(k, j)

                muj_c(k, j) = kel_sav + epsder

                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                muj_c(k, j) = kel_sav - epsder

                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                muj_c(k, j) = kel_sav
                !       elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
                !       logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
                elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
                logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder

                if (abs(elocp(1) - muj_cb(nelorbjh*(j - 1) + k)) .gt. 1d-4) then
                    write (6, *) ' ERROR in contr. muj_c local en. ', j, k, elocp(1), muj_cb(nelorbjh*(j - 1) + k)
                end if
                if (abs(logpsip(1) - muj_clb((j - 1)*nelorbjh + k)) .gt. 1d-4) then
                    write (6, *) ' ERROR in contr muj_c  wf ', j, logpsip(1), muj_clb(nelorbjh*(j - 1) + k)
                end if
                write (6, *) j, k, elocp(1), logpsip(1), muj_cb(nelorbjh*(j - 1) + k), muj_clb(nelorbjh*(j - 1) + k)
            end do
        end do
    end if

    if (iessz) then

        if (contractionj .eq. 0) then
            write (6, *) ' Matrix Jasmat spin '
            do j = 1, nelorbjh
                do k = j, nelorbjh

                    kel_sav = jasmatsz(nelorbjh*(k - 1) + j)

                    jasmatsz(nelorbjh*(k - 1) + j) = kel_sav + epsder

                    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                            &, nintpseudo, prefactor, rcutoff, parshell &
                            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                            &, jpseudo, pseudolocal, istart, costz, costz3         &
                            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                            &, psidetln, 1, nelorbh, nelorbjh                         &
                            &, niesd&
                            &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip &
                            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)
                    jasmatsz(nelorbjh*(k - 1) + j) = kel_sav - epsder

                    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                            &, nintpseudo, prefactor, rcutoff, parshell                            &
                            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                            &, jpseudo, pseudolocal, istart, costz, costz3         &
                            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                            &, psidetln, 1, nelorbh, nelorbjh                         &
                            &, niesd&
                            &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)
                    jasmatsz(nelorbjh*(k - 1) + j) = kel_sav
                    !      elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
                    !      logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
                    elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
                    logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder
                    if (abs(elocp(1) - jasmatszb(nelorbjh*(k - 1) + j)) .gt. 1d-4) then
                        write (6, *) ' ERROR in Jas Sz local en. ', j, k, elocp(1), jasmatszb(nelorbjh*(k - 1) + j)
                    end if
                    if (abs(logpsip(1) - jasmatszlb(nelorbjh*(k - 1) + j)) .gt. 1d-4) then
                        write (6, *) ' ERROR in Jas Sz  wf ', j, logpsip(1), jasmatszlb(nelorbjh*(k - 1) + j)
                    end if
                    write (6, *) j, k, elocp(1), logpsip(1), jasmatszb(nelorbjh*(k - 1) + j), jasmatszlb(nelorbjh*(k - 1) + j)
                end do
            end do
        else

            write (6, *) ' Contracted matrix jasmatsz_c Sz-Sz '

            do j = 1, nelorbj_c
                do k = j, nelorbj_c

                    kel_sav = jasmatsz_c(nelorbj_c*(k - 1) + j)

                    jasmatsz_c(nelorbj_c*(k - 1) + j) = kel_sav + epsder

                    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                            &, nintpseudo, prefactor, rcutoff, parshell                            &
                            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                            &, jpseudo, pseudolocal, istart, costz, costz3         &
                            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                            &, psidetln, 1, nelorbh, nelorbjh                         &
                            &, niesd&
                            &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                    jasmatsz_c(nelorbj_c*(k - 1) + j) = kel_sav - epsder

                    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                            &, nintpseudo, prefactor, rcutoff, parshell                            &
                            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                            &, jpseudo, pseudolocal, istart, costz, costz3         &
                            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                            &, psidetln, 1, nelorbh, nelorbjh                         &
                            &, niesd&
                            &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                    jasmatsz_c(nelorbj_c*(k - 1) + j) = kel_sav
                    !       elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
                    !       logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
                    elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
                    logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder

                    if (abs(elocp(1) - jasmatsz_cb(nelorbj_c*(k - 1) + j)) .gt. 1d-4) then
                        write (6, *) ' ERROR in Jas Sz contr. local en. ', j, k, elocp(1), jasmatsz_cb(nelorbj_c*(k - 1) + j)
                    end if
                    if (abs(logpsip(1) - jasmatsz_clb(nelorbj_c*(k - 1) + j)) .gt. 1d-4) then
                        write (6, *) ' ERROR in Jas Sz  contr wf ', j, logpsip(1), jasmatsz_clb(nelorbj_c*(k - 1) + j)
                    end if

                    write (6, *) j, k, elocp(1), logpsip(1) &
                        , jasmatsz_cb(nelorbj_c*(k - 1) + j), jasmatsz_clb(nelorbj_c*(k - 1) + j)
                end do
            end do

        end if

    end if

    if (yesdodet) then

        if (yesfast .ne. 0) then
            !    go to 23541 ! do not compute mc_c

            write (6, *) ' Matrix mu_c '

            do j = 1, nelorb_c
                do k = 1, ipc*ipf*nelorbh

                    kel_sav = mu_c(k, j)

                    mu_c(k, j) = kel_sav + epsder

                    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                            &, nintpseudo, prefactor, rcutoff, parshell                            &
                            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                            &, jpseudo, pseudolocal, istart, costz, costz3         &
                            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                            &, psidetln, 1, nelorbh, nelorbjh                         &
                            &, niesd&
                            &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                    mu_c(k, j) = kel_sav - epsder

                    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                            &, nintpseudo, prefactor, rcutoff, parshell                            &
                            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                            &, jpseudo, pseudolocal, istart, costz, costz3         &
                            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                            &, psidetln, 1, nelorbh, nelorbjh                         &
                            &, niesd&
                            &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)
                    mu_c(k, j) = kel_sav
                    !      elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
                    !      logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
                    elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
                    logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder
                    if (abs(elocp(1) - mu_cb(ipc*ipf*nelorbh*(j - 1) + k)) .gt. 1d-4) then
                        write (6, *) ' ERROR in mu_c local energy ', k, j, elocp(1), mu_cb(ipc*ipf*nelorbh*(j - 1) + k)
                    end if
                    if (abs(logpsip(1) - mu_clb(ipc*ipf*nelorbh*(j - 1) + k)) .gt. 1d-4) then
                        write (6, *) ' ERROR in mu_c wf ', k, j, logpsip(1), mu_clb(ipc*ipf*nelorbh*(j - 1) + k)
                    end if
                    write (6, *) k, j, elocp(1), logpsip(1) &
                        , mu_cb(ipc*ipf*nelorbh*(j - 1) + k), mu_clb(ipc*ipf*nelorbh*(j - 1) + k)
                end do
            end do

23541       continue
            write (6, *) ' Matrix detmat_c '

            do j = 1, nelcol_c
                do k = 1, ipc*nelorb_c

                    kel_sav = detmat_c(ipc*nelorb_c*(j - 1) + k)
                    detmat_c(ipc*nelorb_c*(j - 1) + k) = kel_sav + epsder
                    if (ipf .eq. 2 .and. j .le. nelorb_c) then
                        if (ipc .eq. 1) then
                            detmat_c(ipc*nelorb_c*(k - 1) + j) = -(kel_sav + epsder)
                            if (j .eq. k) then
                                detmat_c(nelorb_c*(k - 1) + j) = 0.d0
                                detmat_c(nelorb_c*(j - 1) + k) = 0.d0
                            end if
                        else
                            if (mod(k, 2) .eq. 0) then
                                detmat_c(ipc*nelorb_c*(k/2 - 1) + 2*j) = -(kel_sav + epsder)
                                if (j .eq. k/2) then
                                    detmat_c(2*nelorb_c*(k/2 - 1) + 2*j) = 0.d0
                                    detmat_c(2*nelorb_c*(j - 1) + k) = 0.d0
                                end if
                            else
                                detmat_c(ipc*nelorb_c*((k + 1)/2 - 1) + 2*j - 1) = -(kel_sav + epsder)
                                if (j .eq. (k + 1)/2) then
                                    detmat_c(2*nelorb_c*((k + 1)/2 - 1) + 2*j - 1) = 0.d0
                                    detmat_c(2*nelorb_c*(j - 1) + k) = 0.d0
                                end if
                            end if
                        end if
                    end if
                    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                            &, nintpseudo, prefactor, rcutoff, parshell                            &
                            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                            &, jpseudo, pseudolocal, istart, costz, costz3         &
                            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                            &, psidetln, 1, nelorbh, nelorbjh                         &
                            &, niesd&
                            &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                    detmat_c(ipc*nelorb_c*(j - 1) + k) = kel_sav - epsder
                    if (ipf .eq. 2 .and. j .le. nelorb_c) then
                        if (ipc .eq. 1) then
                            detmat_c(ipc*nelorb_c*(k - 1) + j) = -(kel_sav - epsder)
                            if (j .eq. k) then
                                detmat_c(nelorb_c*(k - 1) + j) = 0.d0
                                detmat_c(nelorb_c*(j - 1) + k) = 0.d0
                            end if
                        else
                            if (mod(k, 2) .eq. 0) then
                                detmat_c(ipc*nelorb_c*(k/2 - 1) + 2*j) = -(kel_sav - epsder)
                                if (j .eq. k/2) then
                                    detmat_c(2*nelorb_c*(k/2 - 1) + 2*j) = 0.d0
                                    detmat_c(2*nelorb_c*(j - 1) + k) = 0.d0
                                end if
                            else
                                detmat_c(ipc*nelorb_c*((k + 1)/2 - 1) + 2*j - 1) = -(kel_sav - epsder)
                                if (j .eq. (k + 1)/2) then
                                    detmat_c(2*nelorb_c*((k + 1)/2 - 1) + 2*j - 1) = 0.d0
                                    detmat_c(2*nelorb_c*(j - 1) + k) = 0.d0
                                end if
                            end if
                        end if
                    end if

                    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                            &, nintpseudo, prefactor, rcutoff, parshell                            &
                            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                            &, jpseudo, pseudolocal, istart, costz, costz3         &
                            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                            &, psidetln, 1, nelorbh, nelorbjh                         &
                            &, niesd&
                            &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                    detmat_c(ipc*nelorb_c*(j - 1) + k) = kel_sav
                    if (ipf .eq. 2 .and. j .le. nelorb_c) then
                        if (ipc .eq. 1) then
                            detmat_c(ipc*nelorb_c*(k - 1) + j) = -kel_sav
                            if (j .eq. k) then
                                detmat_c(nelorb_c*(k - 1) + j) = 0.d0
                                detmat_c(nelorb_c*(j - 1) + k) = 0.d0
                            end if
                        else
                            if (mod(k, 2) .eq. 0) then
                                detmat_c(ipc*nelorb_c*(k/2 - 1) + 2*j) = -kel_sav
                                if (j .eq. k/2) then
                                    detmat_c(2*nelorb_c*(k/2 - 1) + 2*j) = 0.d0
                                    detmat_c(2*nelorb_c*(j - 1) + k) = 0.d0
                                end if
                            else
                                detmat_c(ipc*nelorb_c*((k + 1)/2 - 1) + 2*j - 1) = -kel_sav
                                if (j .eq. (k + 1)/2) then
                                    detmat_c(2*nelorb_c*((k + 1)/2 - 1) + 2*j - 1) = 0.d0
                                    detmat_c(2*nelorb_c*(j - 1) + k) = 0.d0
                                end if
                            end if
                        end if
                    end if

                    !      elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
                    !      logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
                    elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
                    logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder

                    if (abs(elocp(1) - detmat_cb(ipc*nelorb_c*(j - 1) + k)) .gt. 1d-4) then
                        write (6, *) ' ERROR in detmat_c local energy ', k, j, elocp(1), detmat_cb(ipc*nelorb_c*(j - 1) + k)
                    end if
                    if (abs(logpsip(1) - detmat_clb(ipc*nelorb_c*(j - 1) + k)) .gt. 1d-4) then
                        write (6, *) ' ERROR in detmat_c wf ', k, j, logpsip(1), detmat_clb(ipc*nelorb_c*(j - 1) + k)
                    end if

                    write (6, *) k, j, elocp(1), logpsip(1) &
                        , detmat_cb(ipc*nelorb_c*(j - 1) + k), detmat_clb(ipc*nelorb_c*(j - 1) + k)
                end do
            end do

        else

            write (6, *) ' Matrix detmat '

            do j = 1, nelcolh
                do k = 1, ipf*ipc*nelorbh
                    kel_sav = detmat(ipc*nelorbh*ipf*(j - 1) + k)
                    detmat(ipc*nelorbh*ipf*(j - 1) + k) = kel_sav + epsder
                    if (ipf .eq. 2 .and. j .le. 2*nelorbh) then
                        if (ipc .eq. 1) then
                            detmat(ipc*nelorbh*ipf*(k - 1) + j) = -(kel_sav + epsder)
                            if (j .eq. k) then
                                detmat(ipc*nelorbh*ipf*(k - 1) + j) = 0.d0
                                detmat(ipc*nelorbh*ipf*(j - 1) + k) = 0.d0
                            end if
                        else
                            if (mod(k, 2) .eq. 0) then
                                detmat(ipc*nelorbh*ipf*(k/2 - 1) + 2*j) = -(kel_sav + epsder)
                                if (k/2 .eq. j) then
                                    detmat(ipc*nelorbh*ipf*(k/2 - 1) + 2*j) = 0.d0
                                    detmat(ipc*nelorbh*ipf*(j - 1) + k) = 0.d0
                                end if
                            else
                                detmat(ipc*nelorbh*ipf*((k + 1)/2 - 1) + 2*j - 1) = -(kel_sav + epsder)
                                if ((k + 1)/2 .eq. j) then
                                    detmat(ipc*nelorbh*ipf*((k + 1)/2 - 1) + 2*j - 1) = 0.d0
                                    detmat(ipc*nelorbh*ipf*(j - 1) + k) = 0.d0
                                end if
                            end if
                        end if
                    end if
                    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                            &, nintpseudo, prefactor, rcutoff, parshell                            &
                            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                            &, jpseudo, pseudolocal, istart, costz, costz3         &
                            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                            &, psidetln, 1, nelorbh, nelorbjh                         &
                            &, niesd&
                            &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                    detmat(ipc*nelorbh*ipf*(j - 1) + k) = kel_sav - epsder
                    if (ipf .eq. 2 .and. j .le. 2*nelorbh) then
                        if (ipc .eq. 1) then
                            detmat(ipc*nelorbh*ipf*(k - 1) + j) = -(kel_sav - epsder)
                            if (j .eq. k) then
                                detmat(ipc*nelorbh*ipf*(k - 1) + j) = 0.d0
                                detmat(ipc*nelorbh*ipf*(j - 1) + k) = 0.d0
                            end if
                        else
                            if (mod(k, 2) .eq. 0) then
                                detmat(ipc*nelorbh*ipf*(k/2 - 1) + 2*j) = -(kel_sav - epsder)
                                if (k/2 .eq. j) then
                                    detmat(ipc*nelorbh*ipf*(k/2 - 1) + 2*j) = 0.d0
                                    detmat(ipc*nelorbh*ipf*(j - 1) + k) = 0.d0
                                end if
                            else
                                detmat(ipc*nelorbh*ipf*((k + 1)/2 - 1) + 2*j - 1) = -(kel_sav - epsder)
                                if ((k + 1)/2 .eq. j) then
                                    detmat(ipc*nelorbh*ipf*((k + 1)/2 - 1) + 2*j - 1) = 0.d0
                                    detmat(ipc*nelorbh*ipf*(j - 1) + k) = 0.d0
                                end if
                            end if
                        end if
                    end if

                    call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                            &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                            &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                            &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                            &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                            &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                            &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                            &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax &
                            &, nintpseudo, prefactor, rcutoff, parshell                            &
                            &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                            &, jpseudo, pseudolocal, istart, costz, costz3         &
                            &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                            &, psidetln, 1, nelorbh, nelorbjh                         &
                            &, niesd&
                            &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                            &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                            &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)
                    detmat(ipc*nelorbh*ipf*(j - 1) + k) = kel_sav
                    if (ipf .eq. 2 .and. j .le. 2*nelorbh) then
                        if (ipc .eq. 1) then
                            detmat(ipc*nelorbh*ipf*(k - 1) + j) = -kel_sav
                            if (j .eq. k) then
                                detmat(ipc*nelorbh*ipf*(k - 1) + j) = 0.d0
                                detmat(ipc*nelorbh*ipf*(j - 1) + k) = 0.d0
                            end if
                        else
                            if (mod(k, 2) .eq. 0) then
                                detmat(ipc*nelorbh*ipf*(k/2 - 1) + 2*j) = -kel_sav
                                if (k/2 .eq. j) then
                                    detmat(ipc*nelorbh*ipf*(k/2 - 1) + 2*j) = 0.d0
                                    detmat(ipc*nelorbh*ipf*(j - 1) + k) = 0.d0
                                end if
                            else
                                detmat(ipc*nelorbh*ipf*((k + 1)/2 - 1) + 2*j - 1) = -kel_sav
                                if ((k + 1)/2 .eq. j) then
                                    detmat(ipc*nelorbh*ipf*((k + 1)/2 - 1) + 2*j - 1) = 0.d0
                                    detmat(ipc*nelorbh*ipf*(j - 1) + k) = 0.d0
                                end if
                            end if
                        end if
                    end if
                    !      elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
                    !      logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
                    elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
                    logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder
                    if (abs(elocp(1) - detmatb(ipc*ipf*nelorbh*(j - 1) + k)) .gt. 1d-4) then
                        write (6, *) ' ERROR in detmat local energy ', k, j, elocp(1), detmatb(ipf*ipc*nelorbh*(j - 1) + k)
                    end if
                    if (abs(logpsip(1) - detmatlb(ipc*ipf*nelorbh*(j - 1) + k)) .gt. 1d-4) then
                        write (6, *) ' ERROR in detmat wf ', k, j, logpsip(1), detmatlb(ipc*ipf*nelorbh*(j - 1) + k)
                    end if
                    write (6, *) k, j, elocp(1), logpsip(1) &
                        , detmatb(ipc*ipf*nelorbh*(j - 1) + k), detmatlb(ipc*ipf*nelorbh*(j - 1) + k)

                    !     stop

                end do
            end do

        end if

    end if

    if (npar_eagp .gt. 0) then
        write (6, *) ' Matrix Ghosts variables '

        do j = 1, ndiff
            do k = 1, ipc*ndiff

                kel_sav = eagp_pfaff(k, j)
                eagp_pfaff(k, j) = kel_sav + epsder
                if (ipc .eq. 1) then
                    eagp_pfaff(j, k) = -(kel_sav + epsder)
                    if (j .eq. k) then
                        eagp_pfaff(j, k) = 0.d0
                    end if
                else
                    if (mod(k, 2) .eq. 0) then
                        eagp_pfaff(2*j, k/2) = -(kel_sav + epsder)
                        if (j .eq. k/2) then
                            eagp_pfaff(k, j) = 0.d0
                        end if
                    else
                        eagp_pfaff(2*j - 1, (k + 1)/2) = -(kel_sav + epsder)
                        if (j .eq. (k + 1)/2) then
                            eagp_pfaff(k, j) = 0.d0
                        end if
                    end if
                end if
                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                eagp_pfaff(k, j) = kel_sav - epsder
                if (ipc .eq. 1) then
                    eagp_pfaff(j, k) = -(kel_sav - epsder)
                    if (j .eq. k) then
                        eagp_pfaff(j, k) = 0.d0
                    end if
                else
                    if (mod(k, 2) .eq. 0) then
                        eagp_pfaff(2*j, k/2) = -(kel_sav - epsder)
                        if (j .eq. k/2) then
                            eagp_pfaff(k, j) = 0.d0
                        end if
                    else
                        eagp_pfaff(2*j - 1, (k + 1)/2) = -(kel_sav - epsder)
                        if (j .eq. (k + 1)/2) then
                            eagp_pfaff(k, j) = 0.d0
                        end if
                    end if
                end if

                call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                        &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                        &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                        &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                        &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                        &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                        &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                        &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                        &, nintpseudo, prefactor, rcutoff, parshell                            &
                        &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                        &, jpseudo, pseudolocal, istart, costz, costz3         &
                        &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                        &, psidetln, 1, nelorbh, nelorbjh                         &
                        &, niesd&
                        &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                        &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscaleo&
                        &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

                eagp_pfaff(k, j) = kel_sav
                if (ipc .eq. 1) then
                    eagp_pfaff(j, k) = -kel_sav
                    if (j .eq. k) then
                        eagp_pfaff(j, k) = 0.d0
                    end if
                else
                    if (mod(k, 2) .eq. 0) then
                        eagp_pfaff(2*j, k/2) = -kel_sav
                        if (j .eq. k/2) then
                            eagp_pfaff(k, j) = 0.d0
                        end if
                    else
                        eagp_pfaff(2*j - 1, (k + 1)/2) = -kel_sav
                        if (j .eq. (k + 1)/2) then
                            eagp_pfaff(k, j) = 0.d0
                        end if
                    end if
                end if

                !      elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
                !      logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
                elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
                logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder

                if (abs(elocp(1) - eagp_pfaffb(k, j)) .gt. 1d-4) then
                    write (6, *) ' ERROR in Ghost det local energy ', k, j, elocp(1), eagp_pfaffb(k, j)
                end if
                if (abs(logpsip(1) - eagp_pfafflb(k, j)) .gt. 1d-4) then
                    write (6, *) ' ERROR in Ghost det  wf ', k, j, logpsip(1), eagp_pfafflb(k, j)
                end if
                write (6, *) k, j, elocp(1), logpsip(1), eagp_pfaffb(k, j), eagp_pfafflb(k, j)
            end do
        end do
    end if

1133 continue

    if (iespbc .and. yesforce) then

        write (6, *) ' Derivative with respect to cell '

        cellscaleo = cellscale

        cellderiv = .false.

        !     cellderiv=.true.

        !        if(yes_tilted) then
        nk_cell = 12
        min_cell = 4
        !        else
        !        nk_cell=3
        !        min_cell=1
        !        endif

        do k = min_cell, nk_cell

            cellscalen = cellscaleo

            cellscalen(k) = cellscaleo(k) + epsder

            call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                    &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                    &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                    &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                    &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                    &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                    &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                    &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                    &, nintpseudo, prefactor, rcutoff, parshell                            &
                    &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                    &, jpseudo, pseudolocal, istart, costz, costz3         &
                    &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                    &, psidetln, 1, nelorbh, nelorbjh                         &
                    &, niesd&
                    &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocp, logpsip&
                    &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscalen&
                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

            cellscalen(k) = cellscaleo(k) - epsder

            call compute_eloc_logpsi(indt, indt4, indt4j, nelorb, nelup, neldo&
                    &, tabpip, kelind, kel, winv, winvup, winvdo, ainv, ainvup, ainvdo, psip &
                    &, ipsip, wconfn, psisn, iesdr, vj, dupr                       &
                    &, zetar, rion, dist, ioccup, ioccdo, ioptorb                             &
                    &, nshell, nshelldo, ivic, alat, plat, vpot, tmu&
                    &, nion, r, rmu, kion, iond, winvj, ioccj, kionj, vjur, nelorbj   &
                    &, ioptorbj, nshellj, winvbar, detmat, winvjbar, winvjbarsz, jasmat       &
                    &, jasmatsz, muj_c, jasmat_c, jasmatsz_c, yesfastj, nelorbj_c, iessz, cnorm, iflagerr, npsa, lmax        &
                    &, nintpseudo, prefactor, rcutoff, parshell                            &
                    &, nparpshell, kindion, pshell, wpseudo, legendre, versor, wintpseudo     &
                    &, jpseudo, pseudolocal, istart, costz, costz3         &
                    &, angle, indtm, LBox, rmucos, rmusin, kappa, vpotreg, cutreg           &
                    &, psidetln, 1, nelorbh, nelorbjh                         &
                    &, niesd&
                    &, iond_cart, mu_c, detmat_c, projm, .false., nelorb_c, firstmol, nmolfn, yesfast, elocm, logpsim&
                    &, nelorbjmax, neldomax, indtmax, nshelljmax, cellscalen&
                    &, indpar_tab, indorb_tab, indshell_tab, indparj_tab, indorbj_tab, indshellj_tab)

            !      elocp(1)=(elocp(1)-elocm(1))/2.d0/epsder
            !      logpsip(1)=(logpsip(1)-logpsim(1))/2.d0/epsder
            elocp(1) = sum(aadinp(:)*(elocp(:) - elocm(:)))/2.d0/epsder
            logpsip(1) = sum(aadinp(:)*(logpsip(:) - logpsim(:)))/2.d0/epsder
            if (abs(elocp(1) - sr2b(k - 3)) .gt. 1d-4) then
                write (6, *) ' ERROR in cell derivative loc. energy=', k - 3, elocp(1), sr2b(k - 3)
            end if
            if (abs(logpsip(1) - sr2lb(k - 3)) .gt. 1d-4) then
                write (6, *) ' ERROR in cell derivative wf=', k - 3, logpsip(1), sr2lb(k - 3)
            end if
            write (6, *) k - 3, elocp(1), sr2b(k - 3), logpsip(1), sr2lb(k - 3)

        end do

    end if

    write (6, *) ' Time compute =', cclock() - timep

    go to 24

    !23     continue

    allocate (zout(indt, 1))
    allocate (zoutb(indt, 1))
    allocate (zoutsav(indt, 1))
    do i = 1, indt
        zoutsav(i, 1) = drand1()
    end do
    zoutb = zoutsav
    distb = 0.d0
    legendreb = 0.d0
    ivicb = 0.d0
    wpseudob = 0.d0
    pseudolocalb = 1.d0

    write (6, *) ' Check pseudoset_b ', epsder

    call PSEUDOSET_B(1, kel, kelb, ivic, ivicb, prefactor, &
            &  zoutb, pseudolocal, pseudolocalb, nion, nintpseudo, dist, distb, rcutoff, &
            &  kindion, rion, rionb, pshell, nparpshell, parshell, lmax, wpseudo, &
            &  wpseudob, legendre, legendreb, versor, wintpseudo, jpseudo, npsa, &
            &  indt, indtm, .false., angle, psip, psipb, lbox, 1d-9, iflagerr)

    write (6, *) ' Numerical /Analityc derivative pseudoset_b '

    zout = 0.d0

    do k = 1, 3
        write (6, *) ' Component ', k

        write (6, *) ' # el '
        kel_sav = kelind(k, j)
        kelind(k, 1) = kel_sav + epsder

        call pseudoset(1, kelind, ivic, zout, pseudolocal, nel &
                &, dist, rion, wpseudo, npsa, indt, indtm &
                &, .false., angle, psip, Lbox, 1d-9, iflagerr, rc)

        costp = sum(zoutsav(1:indt, 1)*zout(1:indt, 1)) + pseudolocal(1)

        kelind(k, 1) = kel_sav - epsder

        call pseudoset(1, kelind, ivic, zout, pseudolocal, nel &
                &, dist, rion, wpseudo, npsa, indt, indtm &
                &, .false., angle, psip, Lbox, 1d-9, iflagerr, rc)

        costm = sum(zoutsav(1:indt, 1)*zout(1:indt, 1)) + pseudolocal(1)

        write (6, *) k, (costp - costm)/2.d0/epsder, kelb(k, 1)

        kelind(k, 1) = kel_sav

        do j = 1, nion
            write (6, *) ' ion =', j

            rion(k, j) = rionsav(k, j) + epsder

            call pseudoset(1, kelind, ivic, zout, pseudolocal, nel &
                    &, dist, rion, wpseudo, npsa, indt, indtm &
                    &, .false., angle, psip, Lbox, 1d-9, iflagerr, rc)

            costp = sum(zoutsav(1:indt, 1)*zout(1:indt, 1)) + pseudolocal(1)

            rion(k, j) = rionsav(k, j) - epsder

            call pseudoset(1, kelind, ivic, zout, pseudolocal, nel &
                    &, dist, rion, wpseudo, npsa, indt, indtm &
                    &, .false., angle, psip, Lbox, 1d-9, iflagerr, rc)

            costm = sum(zoutsav(1:indt, 1)*zout(1:indt, 1)) + pseudolocal(1)

            write (6, *) k, (costp - costm)/2.d0/epsder, rionb(k, j)

            rion(k, j) = rionsav(k, j)

        end do

    end do

24  continue

    write (6, *) 'Warning  deallocating everything ...'

    call Finalizeall
    if (allocated(winvfn)) deallocate (winvfn, winvbarfn)
    if (allocated(zout)) deallocate (zout)
    if (allocated(zoutb)) deallocate (zoutb)
    if (allocated(zoutsav)) deallocate (zoutsav)
    if (allocated(lapnum)) deallocate (lapnum)
    deallocate (eagp_pfafflb, eagp_pfaffb)
    deallocate (detmatb)
    deallocate (detmatlb)
    deallocate (detmat_cb, mu_cb)
    deallocate (detmat_clb, mu_clb)
    deallocate (jasmatb)
    deallocate (jasmatlb)
    deallocate (jasmatszb)
    deallocate (jasmatszlb)
    deallocate (muj_cb)
    deallocate (muj_clb)
    deallocate (jasmat_cb)
    deallocate (jasmat_clb)
    deallocate (jasmatsz_cb)
    deallocate (jasmatsz_clb)
    deallocate (forcedw, pulaydw)

    deallocate (kelb, kelindb, rionb&
            &, tabpipb, winvupb, winvdob, ainvb, ainvupbb, ainvdobb&
            &, psipb, distb, rb, rmub, iond_cartb, winvb, winvjb&
            &, winvbarb, winvjbarb, winvjbarszb, prefactorb, wpseudob&
            &, legendreb, rmucosb, rmusinb, tmub, ivicb, tabpipsav)
    deallocate (vjb, vjlb, duprb, duprlb, vjurb, vjurlb)

    stop

contains

    subroutine Initializeall
        !------------Berry phase-------------------------
        !------------------------------------------------
        ! by E. Coccia (9/11/10)
        use extpot, only: mm_restr, n_x, n_y, n_z, delta, x0, ext_pot, link_atom, write_rwalk
        use splines, only: bscoef
        ! by E. Coccia (23/12/10)
        use van_der_waals, only: vdw

        implicit none
        real*8 drand1, enercont, jacobian, dnrm2, mapping, reweight_fn
        real*8 rion_ref(3), mind(3)
        real, external :: ran
        real*8, external :: atom_weight
        integer ind, ii, jj, nmoltry, read_seeds, read_seeds_mpiio, ntpar
        ! by E. Coccia (15/4/11)
        integer, external :: omp_get_max_threads
        real*8, external :: dlamch
#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(1) ! scalar code
        old_threads = omp_get_max_threads()
#endif
        if (nproc .gt. 1 .and. rank .eq. 0) write (6, *) ' Number of mpi proc =', nproc
#ifdef _OPENMP
        old_threads = omp_get_max_threads()
#else
        old_threads = 1
#endif
        if (rank .eq. 0) write (6, *) ' Number of threads/mpi proc =', old_threads
        new_threads = old_threads

        epsmach = 1000.d0*dlamch('e')
        safemin = 10000.d0*dlamch('s')

        !          Initialize flags error
        yesdft = .false.
        iflagerr = 0
        !          iflagerr is the local processor flag for error.
        !          Each subroutine does not make anything if iflagerr is non zero.
        !          Inside any non parallel subroutine iflag is set to 1 if the
        !          subroutine went in error. In the output, the subroutine "checkiflagerr"
        !          will make an mpi_allreduce of iflagerr summing to iflagerrall.
        !          All processor will stop if iflagerrall>0, namely the run
        !          will continue if and only if all processor have not found any error.
        !
        call get_dir(path)
        if (rank == 0) write (*, *) ' Initial path : ', trim(path)
        !
        ! reading input data
        !
        call read_datasmin

#ifdef PARALLEL
        row_comm = MPI_COMM_WORLD
        col_comm = MPI_COMM_WORLD
        if (yesquantum .and. in1 .ne. 1) then
            if (rank .eq. 0) write (6, *) ' ERROR walker = nproc in quantum case !!!'
            call mpi_finalize(ierr)
            stop
        end if
#endif
        row_id = 0
        col_id = rank
        yescomm = .true.

        !       Default values of communicators/id
#ifdef  PARALLEL
        commopt_mpi = MPI_COMM_WORLD
        rankopt = rank
        nprocopt = nproc

        commrep_mpi = MPI_COMM_WORLD
        rankcolrep = rank
        rankrep = rank
        nprocrep = nproc

        row_comm = MPI_COMM_WORLD
        mcol = nproc

        commsr_mpi = MPI_COMM_WORLD
        commcolsr_mpi = 0
        nprocsr = nproc
        ranksr = rank

        commcov_mpi = MPI_COMM_WORLD
        nproccov = nproc

!       commjas_mpi=MPI_COMM_WORLD
!       nprocjas=nproc
#else
        commopt_mpi = 0
        rankopt = rank
        nprocopt = nproc

        commrep_mpi = 0
        rankcolrep = rank
        rankrep = rank
        nprocrep = nproc

        row_comm = 0
        mcol = nproc

        commsr_mpi = 0
        commcolsr_mpi = 0
        nprocsr = nproc
        ranksr = rank

        commcov_mpi = 0
        nproccov = nproc

        !       commjas_mpi=0
        !       nprocjas=nproc
#endif

        if (yesquantum .and. nbead .eq. nproc .and. nrep_bead .eq. 1) yescomm = .false.
#ifdef PARALLEL
        if (yesquantum .and. nbead .gt. 1) then

            mcol = nproc/(nbead/nrep_bead)

            if (mcol*(nbead/nrep_bead) .ne. nproc .or. (nbead/nrep_bead)*nrep_bead .ne. nbead) then
                if (rank .eq. 0) write (6, *) ' ERROR nproc should be multiple of # beads !!!'
                call mpi_finalize(ierr)
                stop
            end if

            irow = rank/mcol
            jcol = mod(rank, mcol)
            call mpi_comm_split(mpi_comm_world, irow, jcol, row_comm, ierr)
            call mpi_comm_split(mpi_comm_world, jcol, irow, col_comm, ierr)
            call mpi_comm_rank(row_comm, row_id, ierr)
            call mpi_comm_rank(col_comm, col_id, ierr)
            call mpi_barrier(MPI_COMM_WORLD, ierr) ! for unreliable networks

            nproccolrep = nbead
            if (nrep_bead .gt. 1) then
                mcol_rep = mcol/nrep_bead ! i.e. = nproc/nbead per bead
                irow = rank/mcol_rep
                nprocrep = mcol_rep
                jcol = mod(rank, mcol_rep)
                call mpi_comm_split(mpi_comm_world, irow, jcol, commrep_mpi, ierr)
                call mpi_comm_split(mpi_comm_world, jcol, irow, commcolrep_mpi, ierr)
                call mpi_comm_rank(commrep_mpi, rankrep, ierr)
                call mpi_comm_rank(commcolrep_mpi, rankcolrep, ierr)
            else
                mcol_rep = mcol ! i.e. = nproc/(nbead/nrep_bead) per bead group
                commrep_mpi = row_comm
                nprocrep = mcol
                commcolrep_mpi = col_comm
                rankrep = row_id
                rankcolrep = col_id
            end if

!       row_id = column number
!       col_id = row number
            if (.not. yesavopt) then
                nprocopt = mcol
                commopt_mpi = row_comm
                rankopt = row_id
            end if
            if (.not. yesavcov) then
                commcov_mpi = commrep_mpi
                nproccov = mcol_rep
            end if
            if (.not. yesavsr) then
                commsr_mpi = row_comm
                commcolsr_mpi = col_comm
                nprocsr = mcol
                ranksr = row_id
            end if
!       if(.not.yesavjas) then
!         commjas_mpi=row_comm
!         nprocjas=mcol
!       endif

        end if

! Define nprocu
        if (nproc_diag .eq. 0) then
            nprocu = nprocopt
        elseif (nproc_diag .eq. 1) then
            nprocu = 1
        else
            nprocu = nprocopt
        end if

#else
        nprocu = 1
        nproc_diag = 1
#endif

        ! create sub communicator for diagonalization
        ! protect the size of nproc_diag
        if (nproc_diag .le. 0 .or. nproc_diag .gt. nprocrep) then
            ! select nproc_diag
            nproc_diag = min(nprocu, nprocrep)
        end if

        if (rank .eq. 0) write (6, *) "sub_comm_diag uses", nproc_diag, "processors"
        call mpi_sub_comm_create(commrep_mpi, nproc_diag, sub_comm_diag, ierr)

        !        write(6,*) ' Topology '
        !        write(6,*) rank,row_id,col_id
        !        call mpi_finalize(ierr)
        !        stop
        !
        call check_scratch(rank, path, scratchpath)
        ranseedfilename = trim(scratchpath)//'randseed.'//trim(chara)
        !
        iese_eff_c = 0 ! index for complex total energy
        iese_eff = min(iese, 3) ! no more than 3 averaged on time corr fun so far
        if (yes_complex) iese_eff_c = 2*iese_eff
        !
        ! reading pseudo potential if any (npsa>0 read by datasmin)
        call read_pseudo
        !
        if (alat .eq. 0.d0) then
            alat = 1.d0
            if (rank .eq. 0) write (6, *) ' warning alat set to one ', alat
        end if
        if (rank .eq. 0) then
            if (alat2 .ne. 0.d0) then
                write (6, *) ' lattice spacing a1 a2 =', abs(alat)              &
                        &, abs(alat*alat2)
            else
                write (6, *) ' Single mesh  =', abs(alat)
            end if
            write (6, *) ' scratch of determinant each ', nscra
        end if
        ! rsignr=0    fixed node
        ! rsignr=1     VMC ref.
        if (tstepfn .eq. 0.d0) then
            rsignr = -rsignr*(1.d0 + gamma)
        elseif (itestr .ne. -5 .or. itestr4 .lt. -10) then
            rsignr = tstepfn
        end if

        if (abs(itestr) .eq. 1 .or. abs(itestr) .eq. 6 .or. itestr .eq. -2 .or. itestr .eq. -3 .or. itestr4 .lt. -10) then

            if (rsignr .ne. 0 .and. rank .eq. 0)                              &
                    &  write (6, *) ' Warning non standard Fixed node !!!!', rsignr
            if (rsignr .eq. 0 .and. rank .eq. 0)                              &
                    &  write (6, *) ' Standard Fixed node '

        end if

        if ((itestr .eq. -5 .and. itestr4 .ge. -10) .and. rsignr .lt. 0) then
            call checkiflagerr(1, rank, ' Error rsign >= 0 in this case !!!')
        end if

        iseed = abs(iseedr)

        ibinit = iboot

        ! setting main indices for
        ! energy and correlation functions

        nmat = np + 1
        ndim = np + 1
        npm = np + 1
        npmn = nbinmax + 1 - ieskin
        if (npmn .gt. npm) npm = npmn

        if (npm .lt. 6) npm = 6

        nmats = nmat - iesm - iesd
        ! Number of VMC parameters
        ndimp = np - ieskin
        ndimpdim = max(ndimp, 1)

        if ((ieskint .ne. 0. .and. epsder .ne. 0.d0) .or. yespress) then
            nwm = nw + 1
            nws = in1 + 1
        else
            nwm = nw
            nws = in1
        end if
        !  The GLOBAL variables are:
        !  wconfn,jbra,enert,zeta

        ! ############### ALLOCATION ############################
        np3p3 = np3 + 3
        allocate (alphavar(npm), wcorw(nfat + 1), alphab(npm))
        ! COMPLEX DEB - allocation of the total energy of double dimension
        if (.not. yes_complex) then
            allocate (etot(npm), wtot(np3p3))
        else
            allocate (etot(2*npm), wtot(2*np3p3))
        end if
        etot = 0.d0
        if (itestrr .eq. -4) then
            if (idyn .lt. 2 .or. idyn .eq. 7 .or. idyn .eq. 17) then
                allocate (cov(1))
                allocate (sov(npmn, npmn, 5))
            else
                allocate (sov(npmn, npmn, 5))
                allocate (cov(max(ieskin*ieskin, 1)))
            end if
        else
            if (idyn .lt. 2 .or. idyn .eq. 7 .or. idyn .eq. 17) then
                allocate (cov(1))
                allocate (sov(npmn, npmn, 3))
            else
                allocate (sov(npmn, npmn, 4))
                allocate (cov(max(ieskin*ieskin, 1)))
            end if
        end if
        ! ########################################################

        do i = 1, nmat
            alphab(i) = 0.d0
        end do
        alphab(nmat) = 1.d0
        cost = 1.d0
        np3m = (np3 - 1)/nmat + 1

        if (rank .eq. 0) then
            write (6, *) ' Number of corr functions written =', np3m*nmat
            write (6, *) ' iopt =', iopt
        end if

        if ((writescratch .ne. 0) .and. (.not. allocated(bufscra))) then
            allocate (bufscra(max((in1*(2*npdim + 2) + 1)*nweight*nmore_force, 1)))
        elseif (.not. allocated(bufscra)) then
            allocate (bufscra(1))
        end if
        ngn = 0
        ngg = 0

        npp = np + 1
        npf = np

        if (npbra .gt. np) then
            call checkiflagerr(1, rank, ' npbra <= np !!!')
        end if

        nwm2 = 0
        nwdim = 0
        nwfix = 0
        nwinv = 0
        nwnel = 0
        nwrep = 0
        nwdw = 0
        nwfree = 0
        nwsw = 0
        nwkin = 0
        nwup = 0
        nwking = 0
        nprest = np
        nindt = 0

        maxoutput = (dble(in1)*(2*dble(npdim) + 2) + 1)*dble(nweight)*8.d0/1d9

        if (rank .eq. 0) then
            if (writescratch .eq. 0) then
                write (6, *) ' Warning TurboRVB needs ', maxoutput, &
                        &' Gygabyte disc space per processor '
            else
                write (6, *) ' Warning TurboRVB needs ', maxoutput, &
                        &' Gygabyte RAM per processor '
            end if
        end if

        if ((writescratch .ne. 0) .and. (.not. allocated(bufscra))) then
            allocate (bufscra(max((in1*(2*npdim + 2) + 1)*nweight, 1)))
        elseif (.not. allocated(bufscra)) then
            allocate (bufscra(1))
        end if
#ifdef PARALLEL
        call rand_init(iseed)
        !    define  random number with certainly different initial iseed
        allocate (ipsip(2*nproc), psip(2*nproc))
        do i = 1, 2*nproc
            irstart = (1.d0 - drand1())*2**29
            psip(i) = 2*irstart + 1
        end do
        call dsortx(psip, 1, 2*nproc, ipsip)
        !     Disregard the lowest iseed (e.g. drand1=1)
        j = 2
        do while (psip(j) .eq. psip(1))
            j = j + 1
        end do
        iseed = psip(j)
        i = 1
        do while (i .le. rank .and. j .lt. 2*nproc)
            j = j + 1
            iseed = psip(j)
            if (psip(j) .ne. psip(j - 1)) i = i + 1
        end do
        iflagerr = 0
        if (i .ne. rank + 1) iflagerr = 1
        call mpi_allreduce(iflagerr, iflagerrall, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        if (iflagerrall .ne. 0) then
            if (rank .eq. 0) write (6, *) ' ERROR in random seed generation, try with different iseed !!! '
            deallocate (ipsip, psip)
            call mpi_finalize(ierr)
            stop
        end if
        deallocate (ipsip, psip)
#else
        irstart = iseed
        iseed = 2*irstart + 1
#endif
        if (rank .eq. 0) write (6, *) ' initial iseed =', rank, iseed
        !     write(6,*) ' Chosen iseed =',rank,iseed
        call rand_init(iseed)

        if (itestr .eq. -5) then
            skipforce = 1
        else
            skipforce = 3
        end if
        !   It has to be opened before open_files because open_files needs some info
        !   (iespbc) read in fort.10
        if (rank .eq. 0) open (unit=10, file='fort.10', form='formatted', position='REWIND')
        ! reading th w.f.
        call read_fort10(10)

        write (6, *) ' yesfast after read_fort10 ', yesfast
        ! open all the rest files
        write (6, *) ' Input membig, indt4, indt4j =', membig, indt4, indt4j
        call open_files(rank, scratchpath)

        if (defparcutg .and. .not. fncont) then
            if (n_body_on .ne. 0) then
                parcutg = 2
            else
                parcutg = 1
            end if
            if (rank .eq. 0) write (6, *) ' Default value of parcutg= ', parcutg
        end if
        if (rank .eq. 0) then
            if (.not. defparcutg .and. itest .eq. 1 .and. .not. fncont .and. n_body_on .ne. 0 .and. parcutg .le. 1) then
                write (6, *) ' Warning parcutg should be equal to 2 in this case !!! You are doing a test? '
            elseif (.not. defparcutg .and. itest .eq. 1 .and. .not. fncont .and. n_body_on .eq. 0 .and. parcutg .ne. 1) then
                write (6, *) ' Warning parcutg should be equal to 1 in this case !!! You are doing a test? '
            end if
        end if
        if (allfit) then
            ntpar = nmax_ion*nmax_ion*4
        else
            ntpar = (nmax_ion*(nmax_ion + 1))/2*4
        end if
        if (npower .gt. 0) then
            initpar = -2 - powermin
            npar = ntpar*npower
            if (rank .eq. 0) write (6, *) ' Default #parameters in the long range Jastrow', npar
        elseif (initpar .lt. -1) then
            if (mod(npar, ntpar) .ne. 0) then
                write (errmsg, *) 'ERROR Jastrow #parameters npar multiple of', ntpar
                call checkiflagerr(1, rank, errmsg)
            end if
        end if

        if (npowersz .gt. 0) then
            initparinv = -2 - powerminsz
            nparinv = ntpar*npowersz
            if (rank .eq. 0) write (6, *) ' Default #parameters in the long range spin Jastrow', nparinv
        elseif (initparinv .lt. -1) then
            if (mod(nparinv, ntpar) .ne. 0) then
                write (errmsg, *) 'ERROR spin Jastrow #parameters nparinv multiple of', ntpar
                call checkiflagerr(1, rank, errmsg)
            end if
        end if
        !      Recomputing ncg_adr as it can be changed above
        ncg_adr = npar + nparsw + nparinv
        if (ncg_adr .gt. 0) then
            allocate (reducel(ncg_adr, kp0), v_adr(ncg_adr), mat_adr(ncg_adr, ncg_adr)&
                    &, ipip_adr(ncg_adr))
            reducel = 0.d0
            v_adr = 0.d0
            mat_adr = 0.d0
            ipip_adr = 0
            if (rank .eq. 0 .and. iopt .eq. 1) write (25, '(a1,10I7)') '#', ntpar, nmax_ion, npower&
                    &, npowersz, initpar, initparinv, initparsw, npar, nparinv, nparsw
        end if

        if (vdw) then
            call vdw_read()
        end if
        ! by E. Coccia (9/5/11): define QM-link and capping atoms
        ! in the Turbo list, and MM-link atoms in the CPMD list
        ! link_read reads the link.dat from CPMD containing
        ! the classical force field for angle and dihedrals
        ! involving the QM-link atom
        if (link_atom) then
            call link_read()
            ! by E. Coccia (31/5/11): exclusion list for vdW interactions
            call exclusion_list()
            ! by E. Coccia (10/5/11): the capping atom
            ! is not an independent object
            call r_capping()
        end if

        firstmol = 1
        if (contraction .gt. 0) then
            nmolfn = nelorb_c
        else
            nmolfn = nelorbh
        end if

        if (contraction .gt. 0) then
            allocate (allowcontr(nelorb_c, 2))
            allowcontr(:, :) = .true.
        end if
        if (contractionj .gt. 0) then
            allocate (allowcontrj(nelorbj_c, 2))
            allowcontrj(:, :) = .true.
        end if

        if (nws .ne. in1) allocate (cnorm_nw(nshell + nshellj))

        !---------Berry phase-------------------------------

        if (developer .eq. -1) allocate (kelsav(3, nel))
        !       stop for trivial input
        if (nelup .le. 0 .or. nelorb .le. 0 .or. nelup .lt. neldo .or. neldo .lt. 0) then

            if (nelup .eq. 0) write (6, *) ' No electrons up  ', nelup
            if (nelup .lt. neldo) write (6, *) ' Please put  #electrons up > #neldo ', nelup, neldo

            if (nelorb .le. 0) write (6, *) ' The electron basis is empty ', nelorb

            if (nelup .lt. 0) write (6, *) ' Negative # electrons up? ', nelup
            if (neldo .lt. 0) write (6, *) ' Negative # electrons down? ', neldo

            call checkiflagerr(1, rank, "Stop for trivial input!")

        end if

        firstmolt = 1
        if (molecular .ne. 0) then
            firstmolt = nelorb_c - molecular + 1
            !       check if only the molecular diagonal are present
            firstmol = firstmolt
            lastmol = firstmolt - 1 ! It could be that all molec. orb. are not used.
            if (ipc .eq. 2) then
                do i = firstmolt, nelorb_c - ndiff
                    do j = firstmolt, nelorb_c - ndiff
                        if (sum(abs(detmat_c(2*nelorb_c*(i - 1) + 2*j - 1:2*nelorb_c*(i - 1) + 2*j))) .ne. 0&
                                &.and. max(i, j) .gt. lastmol) lastmol = max(i, j)
                    end do
                end do
            else
                do i = firstmolt, nelorb_c - ndiff
                    do j = firstmolt, nelorb_c - ndiff
                        if (detmat_c(nelorb_c*(i - 1) + j) .ne. 0 .and. max(i, j) .gt. lastmol) lastmol = max(i, j)
                    end do
                end do
            end if

            !       Update the last molecular orbitals according to the one to be optimized
            if (itestr .eq. -5) then
                call findmaxmat(iessw, nnozero_c, nozero_c, nelorb_c, jbradet, lastmol)
            end if
        elseif (contraction .ne. 0) then
            lastmol = 0
            if (ipc .eq. 2) then
                do j = 1, nelorb_c
                    do i = 1, nelorb_c
                        if (sum(abs(detmat_c(2*nelorb_c*(i - 1) + 2*j - 1:2*nelorb_c*(i - 1) + 2*j))) .ne. 0 &
                            .and. max(i, j) .gt. lastmol) then
                            lastmol = max(i, j)
                        end if
                    end do
                    do i = nelorb_c + 1, nelcol_c
                        if (sum(abs(detmat_c(2*nelorb_c*(i - 1) + 2*j - 1:2*nelorb_c*(i - 1) + 2*j))) .ne. 0 &
                            .and. j .gt. lastmol) then
                            lastmol = j
                        end if
                    end do
                end do
            else
                do j = 1, nelorb_c
                    do i = 1, nelorb_c
                        if (detmat_c(nelorb_c*(i - 1) + j) .ne. 0 .and. max(i, j) .gt. lastmol) lastmol = max(i, j)
                    end do
                    do i = nelorb_c + 1, nelcol_c
                        if (detmat_c(nelorb_c*(i - 1) + j) .ne. 0 .and. j .gt. lastmol) lastmol = j
                    end do
                end do
            end if
            !       Update the last molecular orbitals according to the one to be optimized
            if (itestr .eq. -5) then
                call findmaxmat(iessw, nnozero_c, nozero_c, nelorb_c, jbradet, lastmol)
            end if
        end if

        if (rank .eq. 0 .and. contraction .ne. 0 .and. lastmol .lt&
                    &. nelorb_c - ndiff .and. molecular .ne. 0)&
                    & write (6, *) ' Warning estimated last relevant molecular orbital =', lastmol

        if (yesfast .eq. -1) then

            forceyes = .true.

            if ((itestr .ne. -5 .or. (iesup .eq. 0 .and. iessw .eq. 0)) .and. molecular .ne. 0) then
                yesfast = 1
                firstmol = firstmolt
                nmolfn = lastmol - firstmolt + 1

                cost = 0.d0
                do i = 1, firstmol - 1
                    do j = 1, firstmol - 1
                        cost = cost + abs(detmat_c((ipc*nelorb_c)*(i - 1) + ipc*(j - 1) + 1))
                        if (ipc .gt. 1) cost = cost + abs(detmat_c((ipc*nelorb_c)*(i - 1) + ipc*j))
                    end do
                end do

                write (6, *) ' After 1 ', yesfast

                if (cost .ne. 0.d0) then
                    yesfast = 2
                    firstmol = 1
                    nmolfn = lastmol
                end if

                if (ireadmin .eq. 1 .and. membig) yesfast = 0

            elseif (contraction .ne. 0) then

                if (ireadmin .eq. 1 .and. membig) then
                    yesfast = 0
                else
                    yesfast = 2
                    firstmol = 1
                    nmolfn = lastmol
                end if
                write (6, *) ' After 2 ', yesfast

            else

                yesfast = 0
                write (6, *) ' After 3 ', yesfast

            end if

            if (nelorbh .le. 2*nmolfn .and. symmagp .and. ireadmin .eq. 0 .and. membig) then
                yesfast = 0 ! the basis is too small to be conv.
            end if
            write (6, *) ' After 4 ', yesfast

        else

            forceyes = .false.

            ! here yesfast in changed to zero if there are no MO.
            ! DO NOT USE yesfast=1 if no MO are present!!

            if (yesfast .eq. 1 .and. molecular .ne. 0) then
                firstmol = firstmolt
                nmolfn = lastmol - firstmolt + 1
            elseif (yesfast .eq. 2 .and. contraction .ne. 0) then
                firstmol = 1
                nmolfn = lastmol
            else
                yesfast = 0
            end if

        end if

        write (6, *) ' yesfast after init =', yesfast

        if (contraction .ne. 0) then
#if defined (_OPENMP) && defined (__NOOMP)
            call omp_set_num_threads(old_threads) ! restore the previous threads
#endif
            if (yes_complex) then
                call zgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, (1.d0, 0.d0)&
                        &, mu_c, ipf*nelorbh, detmat_c, nelorb_c, (0.d0, 0.d0), projm, ipf*nelorbh, nprocu&
                        &, rankopt, commopt_mpi)
            else
                call dgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, 1.d0, mu_c, ipf*nelorbh&
                        &, detmat_c, nelorb_c, 0.d0, projm, ipf*nelorbh, nprocu, rankopt, commopt_mpi)
            end if
#if defined (_OPENMP) && defined (__NOOMP)
            call omp_set_num_threads(1) ! restore the scalar
#endif
        end if

        ax = 0.d0
        ay = 0.d0
        az = 0.d0
        nx = 0
        ny = 0
        nz = 0
        nmolmax = 0
        nmolmin = 0
        nmol = 0
        noproj = .true.
        yesmin = 0
        weight_loc = -1.d0
        power = 1.d0

        if (yesfast .eq. 1 .and. molecular .eq. 0) then
            call checkiflagerr(1, rank, ' ERROR yesfast=1 only with molecular orbitals!, Do not &
                    & define yesfast if you have no idea what this variable means !!! ')
        end if

        yeseloc = .true.
        !       No Hessian optimization and ieskin=0
        if (itestrr .ne. -4 .and. ieskin .eq. 0) yeseloc = .false.
        if (rank .eq. 0 .and. epsder .eq. 0 .and. ieskint .ne. 0) then
            if (yeseloc) then
                write (6, *) ' Warning calculation of derivatives local energy in &
                        &forces!!! '
            else
                write (6, *) ' Warning no calculation of derivatives local energy in &
                        &forces!!! '
            end if
        end if

        dofast = .false. ! Initial definition false in any case.
        ! if molopt=0  optimize  in the usual way th

        if (molopt .ne. 0) then

            call read_datasmin_mol

            !        Replacing value of nmol not assumed to change from input
            if (symmagp) then
                if (nmol .ne. molecular - (ndiff)) then
                    nmol = molecular - (ndiff)
                    if (rank .eq. 0) write (6, *) ' Warning replacing nmol =', nmol
                end if
            else
                if (nmol .ne. (molecular - (ndiff))/2) then
                    nmol = (molecular - (ndiff))/2
                    if (rank .eq. 0) write (6, *) ' Warning replacing nmol =', nmol
                end if
            end if

            if (nmol .eq. 0) then
                call checkiflagerr(1, rank, ' ERROR you cannot run |molopt|>0  without &
                        &molecular orbitals in your fort.10, use convertfort10mol.x for &
                        & converting it !!! ')
            end if

            if (.not. symmagp) then
                lastmol = firstmolt + 2*nmolmatw - 1
            else
                lastmol = firstmolt + nmolmatw - 1
            end if

            if (forceyes) then
                if ((yesmin .ne. 0 .and. ireadmin .ne. 1)) then
                    yesfast = 1 ! the fastest code
                    firstmol = nelorb_c - molecular + 1
                    nmolfn = lastmol - firstmol + 1
                    if (2*nmolfn .gt. nelorbh) yesfast = 0
                else
                    if (membig) yesfast = 0
                end if

            else ! forceyes=.false.
                if (yesfast .ne. 0) then
                    if (yesfast .eq. 1) then
                        firstmol = nelorb_c - molecular + 1
                    elseif (yesfast .eq. 2) then
                        if (rank .eq. 0) write (6, *) ' Warning it should be faster with yesfast= 1 !!! '
                        firstmol = 1
                    end if
                    nmolfn = lastmol - firstmol + 1
                end if
            end if ! endif forceyes

            if (yesmin .eq. 1) then
                nmoltry = 2*nelorb_c - molecular
            else
                nmoltry = nelorb_c
            end if

        end if ! endif molopt>0

        if (yesfast .eq. 0 .and. contraction .ne. 0) then
            yesdetmat = .false.
            yesdetmatc = .true.
            if (allocated(detmat)) deallocate (detmat)
            allocate (detmat(ipc*ipf*nelorbh*nelcolh))
            detmat = 0.d0
            if (iscramax .le. ipf*nelorbh*nelcol_c) then
                iscramax = ipf*nelorbh*nelcol_c
                deallocate (psip)
                allocate (psip(iscramax))
                psip = 0.d0
            end if
            call scontract_mat_det(nelorbh, nelorbh, nelcolh, nelorb_c&
                    &, nelcol_c, detmat, detmat_c, mu_c, psip)
        end if

        if (yesmin .ne. 0) then
            ! if(.not.symmagp.and.yesfast.eq.0.and.ngen.ge.nweight) then
            !   call checkiflagerr(1,rank,' ERROR !!! asymmetric AGP does not work with yesfast=0')
            ! endif

            if (yesmin .ne. 0 .and. rank .eq. 0) write (6, *) ' Projection scheme !!! '

            if (detc_proj .and. itestr .eq. -5) then
                if (molecular .eq. 0) then
                    write (errmsg, *) ' ERROR molopt=-1,5 works only&
                            &   with molecular orbitals !!!'
                    call checkiflagerr(1, rank, errmsg)
                end if
                allocate (detmat_proj(nelorb_c*max(nelcol_c, nel)))
                if (rank .eq. 0) write (6, *) ' Warning molecular orbitals optimization &
                        &  with contracted coefficients '
                detmat_proj = 0.d0
            end if

            if (iopt .eq. 1) then
                !       Input detmat  output NEW detmat rank mol, mu_c detmat_c same Z

                ! initialize projmat or orthogonalize the orbit

                if (ndiff .ne. 0) then
                    !        check consistency molecular orbitals
                    flag = .true.
                    do i = 1, ndiff
                        do j = 1, ndiff
                            if (j .ne. i) then
                                if (detmat_c(nelorb_c*(nelorb_c + i - 1) + nelorb_c - ndiff + j) .ne. 0.d0) flag = .false.
                            else
                                if (detmat_c(nelorb_c*(nelorb_c + i - 1) + j + nelorb_c - ndiff) .eq. 0.d0) flag = .false.
                            end if
                        end do
                    end do

                    if (.not. flag) then
                        call checkiflagerr(1, rank, ' The unpaired molecular orbitals should be &
                                & appropriately ordered in input ')
                    end if

                end if ! endif nelup.ne.neldo

                if (symmagp) then
                    nmollimit = molecular - ndiff
                else
                    nmollimit = (molecular - ndiff)/2
                end if

                if (nmolmatw .gt. nmollimit) then
                    write (errmsg, *) ' ERROR not enough molecular orbital in &
                            & input fort.10, nmol/nmolmax  should be less than  ', nmollimit

                    call checkiflagerr(1, rank, errmsg)
                end if

                if (ireadmin .gt. 0 .or. yesfast .eq. 0 .or. yesmin .eq. 1) then ! first big if

                    !               if((ireadmin.eq.1.and.membig).or.&
                    !                    &(yesfast.eq.0.and..not.detc_proj).or&
                    !                    &..not.dofast.or.abs(molopt).gt.1) then ! second big if ireadmin
                    !                  projectagp=0  ! write long output dmrg
                    !                  call convertmol

                    !                  if(detc_proj.and.ireadmin.eq.1) then
                    !                     dofast=.false.
                    !                     yesmin=0
                    !                     molopt=0
                    !                  endif

                    !                  !      if(detc_proj)  detc_proj=.false.
                    !               else
                    !       do ii=1,nelorb_c
                    !       write(6,*) ' Orbital ',ii,sum(abs(mu_c(:,ii)))
                    !       enddo
                    if (ireadmin .eq. 0 .or. membig .or. .not. detc_proj) then !  III big if

                        call convertmol_fast

                        if (allocated(detmat_sav)) then
                            detmat_c = detmat_sav
                            deallocate (detmat_sav)
                        end if

                        !      Preparing projection matrix project_c
                        if (detc_proj) then
                            detmat_c = 0.d0
                            call convertmol_c
                        end if ! endif detc_proj
                    else ! referred to the III big if
                        yesmin = 0
                        !       detc_proj=.false.
                        if (detc_proj) then
                            detc_proj = .false.
                            call convertmol_fast
                            !  below defined slowest more general output for avoiding errors
                            nmolfn = nelorb_c
                            firstmol = 1
                        end if
                    end if ! endif III big if ireadmin
                    !     write(6,*) ' mu_c after convertmol_fast  ',sum(mu_c(:,:))
                    !     do ii=1,nelorb_c
                    !     write(6,*) ' Orbital ',ii,sum(abs(mu_c(:,ii)))
                    !     enddo
                    !      write(6,*) ' Check ortho projmat '
                    !       do i=1,nmolmat
                    !         do j=i,nmolmat
                    !         write(6,*) i,j,sum(projmat(:,i)*projmat(:,j+nmolmat))
                    !         write(6,*) i,j,sum(projmat(:,i)*projmat(:,j+2*nmolmat))&
                    !    &,sum(projmat(:,i+nmolmat)*projmat(:,j+3*nmolmat))
                    !         enddo
                    !       enddo
                    !      stop

                    !               endif ! endif second big if ireadmin
                elseif (ireadmin .gt. 0) then ! referred to the first big if
                    !     prepare projectmat
                    allocate (detmat_sav(nelorb_c*nelcol_c))
                    detmat_sav = detmat_c

                    call convertmol_fast

                    !     restore original detmat_c

                    detmat_c = detmat_sav
                    deallocate (detmat_sav)

                    !     Now project

                    call convertmol_fast

                    !     Now is no longer needed the projection

                    yesmin = 0
                    molopt = 0

                end if ! first big if ireadmin

#if defined (_OPENMP) && defined (__NOOMP)
                call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

                !       In this way readalles makes the average of detmat_c and not mu_c
                !       mu_c is calculated below
                if (contraction .ne. 0) then

                    do jj = 1, iesup_c
                        do kk = 1, multranspip(jj)
                            iy = (transpip(kk)%col(jj) - 1)/nelorbh + 1
                            ix = transpip(kk)%col(jj) - (iy - 1)*nelorbh
                            mu_c(ix, iy) = dup_c(jj)
                        end do
                    end do

                    if (yes_complex) then

                        call zgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, (1.d0, 0.d0), mu_c, ipf*nelorbh&
                                &, detmat_c, nelorb_c, (0.d0, 0.d0), projm, ipf*nelorbh, nprocu, rankopt, commopt_mpi)

                        !    we still need for allio the matrix detmat
                        if (yesdetmatc) call scontract_mat_det(nelorbh, nelorbh, nelcolh, nelorb_c&
                                &, nelcol_c, detmat, detmat_c, mu_c, psip)

                    else

                        call dgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, 1.d0, mu_c, ipf*nelorbh&
                                &, detmat_c, nelorb_c, 0.d0, projm, ipf*nelorbh, nprocu, rankopt, commopt_mpi)

                        !    we still need for allio the matrix detmat
                        if (yesdetmatc) call scontract_mat_det(nelorbh, nelorbh, nelcolh, nelorb_c&
                                &, nelcol_c, detmat, detmat_c, mu_c, psip)

                    end if

                end if

                if (contractionj .ne. 0) then

                    do jj = 1, npar3body_c
                        do kk = 1, multranspipj(jj)
                            iy = (transpipj(kk)%col(jj) - 1)/nelorbjh + 1
                            ix = transpipj(kk)%col(jj) - (iy - 1)*nelorbjh
                            muj_c(ix, iy) = vju_c(jj)
                        end do
                    end do

!                       if(ipj.eq.2) then
!                           call scontract_genj(nelorbjh, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!                       else
!                           call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh&
!                                   &, nelorbj_c, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!                       endif

                    if (iessz) call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh, nelorbj_c&
                            &, nelorbj_c, jasmatsz, jasmatsz_c, muj_c, psip)
                end if

#if defined (_OPENMP) && defined (__NOOMP)
                call omp_set_num_threads(1) ! restore the scalar
#endif

                !       updating scale used later on
                if (yesdetmat) then
                    if (yes_complex) then
                        do ii = 1, nnozero
                            scale(2*ii - 1) = detmat(2*nozero(ii) - 1)
                            scale(2*ii) = detmat(2*nozero(ii))
                        end do
                    else
                        do ii = 1, nnozero
                            scale(ii) = detmat(nozero(ii))
                        end do
                    end if
                end if
                if (contractionj .ne. 0) then
                    !       updating scalej used later on
                    do ii = 1, nnozeroj
                        scalej(ii) = jasmat(nozeroj(ii))
                        !       write(6,*) ii,scalej(ii)
                    end do

                    !       updating scalejsz used later on
                    if (iessz) then
                        do ii = 1, nnozeroj
                            scalejsz(ii) = jasmatsz(nozeroj(ii))
                        end do
                    end if
                end if

            elseif (yesmin .ne. 0) then
                !         if(.not.allocated(projmat)) then

                if (symmagp) then
                    if (detc_proj) then
                        allocate (projmat_c(nelorb_at, 2*nmolmat))
                    end if
                else
                    if (detc_proj) then
                        allocate (projmat_c(nelorb_at, 4*nmolmat))
                    end if
                end if
                !         endif  allocated projmat
            end if
        end if ! endif molopt ne 0

        if (itest .ne. 2) then
            if (ndiff .ne. 0 .and. molecular .ne. 0) then
                nmolfn = nelorb_c - firstmol + 1
                lastmol = nelorb_c
                if (rank .eq. 0) write (6, *) ' Warning changing nmolfn for LRDMC/DMC =', nmolfn
            end if
            if (contraction .ne. 0 .and. firstmol + nmolfn - 1 .gt. nelorb_c) then
                nmolfn = nelorb_c - firstmol + 1
            end if
            if (rank .eq. 0) write (6, *) ' firstmol nmolfn after all ', firstmol, nmolfn
        end if

        !       The projection is the identity if there are enough moleculars

        !  write(6,*) ' projm used =',rank,sum(abs(projm(:))),sum(abs(mu_c(:,:))),sum(abs(detmat_c(:)))

        if ((nelorbh .le. nmolmax .or. iessw .eq. 0) .and. rank .eq. 0 .and. yesmin .ne. 0) then
            write (6, *) ' Warning using fast algorithm, but not necessary !!! '
        end if

        !       write(6,*) ' yesfast =',yesfast
        !       write(6,*) ' symmagp =',symmagp

        !      if(rank.eq.0) write(6,*) 'Warning:  Optimization of Z with no symm. AGP does not work with yesfast>0 '
        !
        !        yesfast=0
        !        endif

        if (rank .eq. 0) then
            write (6, *) ' Default chosen yesfast ', yesfast
            if (yesfast .ne. 0) then
                write (6, *) ' Warning nmolfn, Speeding factor =  ', nmolfn&
                        &, nelorbh/dble(2*nmolfn)
            end if
        end if

        if (molecular .gt. 0 .and. .not. symmagp .and. molopt .gt. 0 .and. yesfast .eq. 0&
                &.and. (ngen .ge. nweight)) then
            call checkiflagerr(1, rank, ' ERROR The program with non symmetric AGP does not &
                    & work with molecular orbitals optimization yesfast=0 ')
        end if

        !       For the time being detmat and nozero have to be allocated in order
        !       to be compatible with readalles when yesmin>0 even when yesfast>0
        if (.not. yesdetmat .and. .not. yesdetmatc) then

            !         in some cases it is also not necessary detmat and nozero

            deallocate (detmat, nozero)
            if (yes_complex) then
                allocate (detmat(2*nelorbh), nozero(1))
            else
                allocate (detmat(nelorbh), nozero(1))
            end if
            detmat = 0.d0
            nozero = 0
            if (itestr .ne. -5) then
                deallocate (nozerodet)
                allocate (nozerodet(1))
            end if
        end if

        if (LBox .gt. 0) then
            call InitEwald(nion, zetar, nel, nws)
            if (rank .eq. 0) then
                write (*, *) ' Number of G vectors ', n_gvec
                write (*, *) ' Ewald Self Energy ', eself
            end if
            kmax = n_gvec
            kmax2 = 2*kmax
        else
            n_gvec = 1
            kmax = 0
            kmax2 = 0
            if (iesbra) then
                allocate (sum_q_cos_gr(1, 1), sum_q_sin_gr(1, 1))
                sum_q_cos_gr = 0.d0
                sum_q_sin_gr = 0.d0
            end if
        end if

        if (rank .eq. 0) then
            write (6, *) '*********************************************'
            write (6, *) '*********************************************'
            if (itestr .eq. -5) then

                if (itestr4 .ge. -10) then
                    write (6, *) 'VMC with ENERGY MINIMIZATION'
                else
                    write (6, *) 'FN   with ENERGY MINIMIZATION'
                end if

                if (itestrr .eq. -4) then
                    write (6, *) 'STOCHASTIC RECONFIGURATION WITH HESSIAN'
                    write (6, *) ' Beta used =', beta
                    write (6, *) ' Cut off used for sr when kl=6,7 =', parr
                elseif (itestrr .eq. -5) then
                    write (6, *) 'SIMPLE STOCHASTIC RECONFIGURATION'
                end if
                write (6, *) ' Stochastic acceleration =', tpar
                if (itestr3 .eq. -9) write (6, *)                                   &
                        &' No optimization Z with contracted '
            elseif (itestr .eq. 2) then
                write (6, *) 'VARIATIONAL MONTE CARLO'
            elseif (itestr .eq. 1 .or. itestr .eq. -2 .or. itestr .eq. -3) then
                write (6, *) 'DIFFUSION MONTE CARLO'
                if (rejweight) then
                    write (6, *) 'weights updated according to acceptance/rejection in diffusion'
                else
                    write (6, *) 'weights always updated'
                end if
            elseif (itestr .eq. -6) then
                write (6, *) 'LATTICE REGULARIZED DIFFUSION MONTE CARLO'
            end if
            write (6, *) '*********************************************'
            write (6, *) '*********************************************'
        end if

        if (itestrr .eq. -4 .or. itestrr .eq. -5) then

            if (contraction .ne. 0 .or. contractionj .ne. 0) then
                if (iesup .ne. 0 .and. .not. yeszagp) then

                    if (rank .eq. 0) write (6, *) 'zeta exponents in determinant excluded'
                end if
                if (iesm .ne. 0 .and. .not. yeszj) then
                    if (rank .eq. 0) write (6, *) 'zeta exponents in Jastrow  excluded'
                end if
            end if

        end if ! itestrr.eq.-4.or.itestrr.eq.-5

        inddsw = iesfree + iesinv + iesm + iesd + 1
        stodim = iesm + iesd + iesup_c + nnozero_c + iesinv + iesfree + 3*nion + 3

        stodim = max(stodim, nmat)

        ! allocation cccccccccccccccccccccccccccccccccccccccccccccccccccccc
        allocate (dup(max(iesup, 1)))
        dup = 0.d0
        if (itestr .eq. -5) then ! only for the optimization this memory is required
            if (stodim .gt. iscramax) then
                if (rank .eq. 0) then
                    write (6, *) ' Warning changing iscramax for the master ', stodim
                    deallocate (psip)
                    iscramax = stodim
                    allocate (psip(iscramax))
                    psip = 0.d0
                end if
            end if
            allocate (ef(ieskindim, 3, nbindim), efenergy(2, nbindim)             &
                    &, err(npm), force(npm), efp(ndimpdim, 2, nbindim))
            ef = 0.d0
            efenergy = 0.d0
            err = 0.d0
            force = 0.d0
            efp = 0.d0
            allocate (fk(nbindim, npdim), fkav(npm), reduce(ncgdim, npdim)&
                    &, okav(npm), skdiag(npm))
            allocate (velion(3, ieskindim), efpress(3, nbindim))
            efpress = 0.d0
            velion = 0.d0
            !     Initialize randomly according to the given temperature
            if (temp .ne. 0.d0) then
                if (rank .eq. 0) write (6, *) ' Initializing velocities at temp (a.u.)', temp*ris(2)
                do ii = 1, ieskindim
                    !      velion(1,ii)=dsqrt(-2.d0*dlog(1.d0-drand1())*temp)*dcos(2.d0*pi*drand1())
                    !      velion(2,ii)=dsqrt(-2.d0*dlog(1.d0-drand1())*temp)*dcos(2.d0*pi*drand1())
                    velion(3, ii) = dsqrt(-2.d0*dlog(1.d0 - drand1())*temp)*dcos(2.d0*pi*drand1())
                end do
            end if
            fk = 0.d0
            fkav = 0.d0
            okav = 0.d0
            skdiag = 0.d0
            reduce = 0.d0
        end if !endif itestr=-5

        ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        stepcg = 0
        !    copy  scale in detmat

        if (yesdetmat) then
            if (yes_complex) then
                call dscalzero(2*ipf*nelorbh*nelcolh, 0.d0, detmat, 1)
                do i = 1, nnozero
                    call upsim_complex(detmat, ipf*nelorbh, nozero(i), scale(2*i - 1), symmagp, ipf)
                end do
            else
                call dscalzero(ipf*nelorbh*nelcolh, 0.d0, detmat, 1)
                do i = 1, nnozero
                    call upsim(detmat, ipf*nelorbh, nozero(i), scale(i), symmagp, ipf)
                end do
            end if
        end if

        !         No longer used this huge memory only in testpsi
        if (developer .eq. 0) then
            deallocate (scale)
            allocate (scale(1))
        end if

#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif

        if (allocated(work)) deallocate (work)

        if (yesfast .eq. 0 .and. rank .eq. 0) then ! main test diagonalization
            deallocate (psip) ! to save memory
            ! calculation eigenvalues
            if (contraction .eq. 0) then
                if (.not. yes_complex) then
                    allocate (work(ipf*ipf*nelorbh*nelorbh + 3*ipf*nelorbh))
                    call dcopy(ipf*ipf*nelorbh*nelorbh, detmat, 1, work, 1)
                    allocate (eig(2*ipf*nelorbh))
                    eig = 0.d0
                    if (symmagp .and. ipf .eq. 1) then
                        call dsyev('N', 'U', nelorbh, work, nelorbh, eig                   &
                                &, work(nelorbh*nelorbh + 1), 3*nelorbh, info)
                    else
                        call dgeev('N', 'N', ipf*nelorbh, work, ipf*nelorbh, eig, eig(ipf*nelorbh + 1), vl, 1, vr, 1&
                                &, work(ipf*ipf*nelorbh*nelorbh + 1), 3*ipf*nelorbh, info)
                        do i = 1, ipf*nelorbh
                            eig(i) = sqrt(eig(i)**2 + eig(i + nelorbh*ipf)**2)
                        end do
                    end if
                else ! complex eigenvalues
                    allocate (work(2*ipf*ipf*nelorbh*nelorbh + 8*ipf*nelorbh))
                    call dcopy(2*ipf*ipf*nelorbh*nelorbh, detmat, 1, work, 1)
                    allocate (eig(2*ipf*nelorbh))
                    call zgeev('N', 'N', ipf*nelorbh, work, ipf*nelorbh, eig, vl, 1, vr, 1&
                            &, work(2*ipf*ipf*nelorbh*nelorbh + 1), 3*ipf*nelorbh&
                            &, work(2*ipf*ipf*nelorbh*nelorbh + 6*ipf*nelorbh + 1), info)
                    do i = 1, 2*ipf*nelorbh, 2
                        eig(i) = sqrt(eig(i)**2 + eig(i + 1)**2)
                    end do
                end if
            else ! contracted orbitals
                if (.not. yes_complex) then
                    allocate (work(nelorb_c*nelorb_c + 5*nelorb_c))
                    call dcopy(nelorb_c*nelorb_c, detmat_c, 1, work, 1)
                    if (symmagp .and. ipf .eq. 1) then
                        allocate (eig(nelorb_c))
                        eig = 0.d0
                        call dsyev('N', 'U', nelorb_c, work, nelorb_c, eig &
                                &, work(nelorb_c*nelorb_c + 1), 3*nelorb_c, info)
                    else
                        allocate (eig(nelorb_c))
                        eig = 0.d0
                        call dgesvd('N', 'N', nelorb_c, nelorb_c, work, nelorb_c, eig, vl, 1, vr, 1&
                                &, work(nelorb_c*nelorb_c + 1), 5*nelorb_c, info)
                    end if
                else ! complex eigenvalues
                    allocate (work(2*nelorb_c*nelorb_c + 11*nelorb_c))
                    call dcopy(2*nelorb_c*nelorb_c, detmat_c, 1, work, 1)
                    allocate (eig(nelorb_c))
                    eig = 0.d0
                    call zgesvd('N', 'N', nelorb_c, nelorb_c, work, nelorb_c, eig, vl, 1, vr, 1&
                            &, work(2*nelorb_c*nelorb_c + 1), 3*nelorb_c, work(2*nelorb_c*nelorb_c + 6*nelorb_c + 1), info)
                end if
            end if

#if defined (_OPENMP) && defined (__NOOMP)
            call omp_set_num_threads(1) ! restore the scalar code
#endif

            !            if(rank.eq.0) then
            !              write(6,*) 'input detmat'
            !               do i=1,nelorb_c
            !                  do j=1,nelorb_c
            !                     write(6,*) i,j,detmat_c(nelorb_c*(i-1)+j),detmat_c(nelorb_c*(j-1)+i)
            !                  enddo
            !               enddo
            !           endif
            !          if(rank.eq.0) then
            !         write(6,*) ' Input detmat '
            !         do i=1,nnozero
            !         write(6,*) i,nozero(i),detmat(nozero(i)),scale(i)
            !         enddo
            !         do i=1,nelorbh
            !          do j=i,nelorbh
            !          write(6,*) i,j,detmat( nelorbh*(i-1)+j),detmat(nelorbh*(j-1)+i)
            !          enddo
            !         enddo
            !         endif

            write (6, *) '%%%%%%%%%%%%%%% DETERMINANTAL GEMINAL %%%%%%%%%%%%%%'
            write (6, *) ' Eigenvalues matrix lambda ', info, rank

            !        write(6,*) ' Input detmat '
            !        do i=1,nnozero
            !        write(6,*) i,nozero(i),detmat(nozero(i)),scale(i)
            !        enddo

            if (contraction .eq. 0) then
                maxdimeig = ipf*nelorbh
            else
                maxdimeig = nelorb_c
            end if

            call print_eigenvalues(rank, eig, maxdimeig, irankdet)

            if (allocated(work)) deallocate (work)
            deallocate (eig)
            allocate (psip(iscramax))
            psip = 0.d0

        end if ! endif test diagonalization

#ifdef PARALLEL
        call mpi_bcast(irankdet, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

        if (yesfast .eq. 0) then ! main if test diagonalization

            if (irankdet .lt. neldo) then
                write (errmsg, *) ' Singular matrix !!! ', irankdet
                call checkiflagerr(1, rank, errmsg)
            elseif (irankdet .eq. neldo) then
                if (rank .eq. 0) write (6, *) ' No resonance simple det !!! ', irankdet
            else
                if (rank .eq. 0) write (6, *) ' RVB is starting !!! ', irankdet
            end if

        end if ! if yesfast=0

        if (contractionj .eq. 0 .and. .not. yes_sparse) then
            call dscalzero(ipj*ipj*nelorbjh*nelorbjh, 0.d0, jasmat, 1)
            do i = 1, nnozeroj
                call upsim(jasmat, ipj*nelorbjh, nozeroj(i), scalej(i), .true., 1)
            end do
        end if
        if (developer .eq. 0) then
            deallocate (scalej)
            allocate (scalej(1))
        end if

        if (iessz) then
            call dscalzero(nelorbjh*nelorbjh, 0.d0, jasmatsz, 1)
            do i = 1, nnozeroj
                call upsim(jasmatsz, nelorbjh, nozeroj(i), scalejsz(i), .true., 1)
            end do
            if (developer .eq. 0) then
                deallocate (scalejsz)
                allocate (scalejsz(1))
            end if
        end if

        !             write the iessw vector that satisfy the constraint
        !             in the chosen basis of non zero elements of lambda

        !         if(ireadmin.eq.0) then
        if (contractionj .ne. 0) then
            call constrbra(iesfree, nnozeroj_c, jbraj, nozeroj_c, jasmat_c   &
                    &, ddw, 1, 1)
            if (iessz)                                                     &
                    &    call constrbra(iesinv, nnozeroj_c, jbrajsz, nozeroj_c            &
                    &, jasmatsz_c, ddwsz, 1, 1)
        else
            if (yes_sparse) then
                call constrbra_sparse(iesfree, nnozeroj, jbraj, jasmat&
                        &, ddw, 1, 1)
            else
                call constrbra(iesfree, nnozeroj, jbraj, nozeroj, jasmat&
                        &, ddw, 1, 1)
            end if
            if (iessz)                                                     &
                    &    call constrbra(iesinv, nnozeroj, jbrajsz, nozeroj                &
                    &, jasmatsz, ddwsz, 1, 1)
        end if

        call constrbr(iesup, iesup_c, jbraiesup, dup_c                  &
                &, dup, 1, 1)
        call constrbr(iesm, npar3body_c, jbraiesm, vju_c                &
                &, vju, 1, 1)

        if (contraction .ne. 0) then
            call constrbra(iessw, nnozero_c, jbradet, nozero_c, psip      &
                    &, dsw, 1, 1)
        else
            call constrbra(iessw, nnozero, jbradet, nozero, detmat            &
                    &, dsw, 1, 1)
        end if

        if (rank .eq. 0) then
            write (6, *) '%%%%%%%%%%%%%%%%%% 2 BODY JASTROW %%%%%%%%%%%%%%%%'
            write (6, *) 'Number 2 body Jastrow parameters ', niesd
            if (niesd .ge. 1) then
                do i = 1, niesd
                    write (6, *) i, vj(i)
                end do
            else
                write (6, *) '2 body Jastrow not present'
            end if
            if (iesdr .le. -5) then
                write (6, *) ' initial costz, costz3, zeta_Q_Caffarel '
                do i = 1, nion
                    write (6, *) i, costz(i), costz3(i), zetaq(i)
                end do
            end if
            if (iesfree .ne. 0 .or. iesm .ne. 0 .or. iesfreesz .ne. 0) then
                write (6, *) '%%%%%%%%%%%%%%%%%% 3 BODY JASTROW %%%%%%%%%%%%%%%%'
                !        if(iesm.ne.0)                                                  &
                !    &   write(6,*) 'Independent orbital parameters in the basis'
                !        do i=1,iesm
                !        write(6,*) i,vju(i)
                !        enddo
                !        if(iesfree.ne.0)                                               &
                !    &   write(6,*) 'Independent geminal Jastrow coefficients'
                !        do i=1,iesfree
                !        write(6,*) i,ddw(i)
                !        enddo
                !        if(iesinv.ne.0) then
                !        write(6,*) 'Independent geminal Jastrow Sz coefficients'
                !        do i=1,iesfreesz
                !        write(6,*) i,ddwsz(i)
                !        enddo
                !        endif
            end if
            write (6, *) 'Total independent parameters in Jastrow sector', &
                    &              inddsw - 1
            write (6, *)

        end if ! endif rank=0

        iesupskip = iesupr + npar3body
        iskip = iesupskip + 1

        !----------------------------------------------
        ! Definition of main indices related to energy,
        ! corr functions and variational paramers
        !----------------------------------------------

        ! first the energy index
        if (iese .ne. 0) then
            if (.not. yes_complex) then
                nwinv = nwinv + in1*iese
                nwm2 = nwm2 + in1*iese
                nwdim = nwdim + in1*iese
                nwfix = nwfix + in1*iese
                nwnel = nwnel + in1*iese
                nwrep = nwrep + in1*iese
                nwdw = nwdw + in1*iese
                nwfree = nwfree + iese*in1
                nwsw = nwsw + iese*in1
                nwkin = nwkin + iese*in1
                nwup = nwup + iese*in1
                nwking = nwking + iese*in1
                nprest = nprest - iese
                nindt = nindt + iese
            else
                nwinv = nwinv + 2*in1*iese
                nwm2 = nwm2 + 2*in1*iese
                nwdim = nwdim + 2*in1*iese
                nwfix = nwfix + 2*in1*iese
                nwnel = nwnel + 2*in1*iese ! this variable is apparently unused
                nwrep = nwrep + 2*in1*iese
                nwdw = nwdw + 2*in1*iese ! this variable is apparently unused
                nwfree = nwfree + 2*iese*in1
                nwsw = nwsw + 2*iese*in1
                nwkin = nwkin + 2*iese*in1
                nwup = nwup + 2*iese*in1
                nwking = nwking + 2*iese*in1
                !
                nprest = nprest - iese
                nindt = nindt + iese ! updated later on
            end if
        end if

        ! then indices related to the w.f. parameters
        if (iesinv .ne. 0) then
            nwm2 = nwm2 + in1*iesinv
            nwdim = nwdim + in1*iesinv
            nwrep = nwrep + in1*iesinv
            nwfix = nwfix + in1*iesinv
            nwnel = nwnel + in1*iesinv
            nwdw = nwdw + in1*iesinv
            nwfree = nwfree + iesinv*in1
            nwsw = nwsw + iesinv*in1
            nwkin = nwkin + iesinv*in1
            nwup = nwup + iesinv*in1
            nwking = nwking + iesinv*in1
            nprest = nprest - iesinv
            nindt = nindt + iesinv
        end if

        if (iesm .ne. 0) then
            nwdim = nwdim + abs(iesm)*in1
            nwrep = nwrep + in1*abs(iesm)
            nwfix = nwfix + abs(iesm)*in1
            nwnel = nwnel + in1*abs(iesm)
            nwdw = nwdw + in1*abs(iesm)
            nwfree = nwfree + abs(iesm)*in1
            nwsw = nwsw + abs(iesm)*in1
            nwkin = nwkin + abs(iesm)*in1
            nwup = nwup + abs(iesm)*in1
            nwking = nwking + abs(iesm)*in1
            nprest = nprest - abs(iesm)
            nindt = nindt + abs(iesm)
        end if

        if (iesd .ne. 0) then
            nwrep = nwrep + in1*abs(iesd)
            nwfix = nwfix + abs(iesd)*in1
            nwnel = nwnel + in1*abs(iesd)
            nwdw = nwdw + in1*abs(iesd)
            nwfree = nwfree + abs(iesd)*in1
            nwsw = nwsw + abs(iesd)*in1
            nwkin = nwkin + abs(iesd)*in1
            nwup = nwup + abs(iesd)*in1
            nwking = nwking + abs(iesd)*in1
            nprest = nprest - abs(iesd)
            nindt = nindt + abs(iesd)
        end if

        if (iesfree .ne. 0) then
            nwrep = nwrep + in1*iesfree
            nwnel = nwnel + in1*iesfree
            nwfix = nwfix + iesfree*in1
            nwkin = nwkin + iesfree*in1
            nwup = nwup + iesfree*in1
            nwking = nwking + iesfree*in1
            nwsw = nwsw + iesfree*in1
            nwdw = nwdw + in1*iesfree
            nprest = nprest - iesfree
            nindt = nindt + iesfree
        end if

        if (iessw .ne. 0) then
            nwrep = nwrep + in1*iessw
            nwnel = nwnel + in1*iessw
            nwkin = nwkin + iessw*in1
            nwfix = nwfix + iessw*in1
            nwup = nwup + iessw*in1
            nwking = nwking + iessw*in1
            nwdw = nwdw + in1*iessw
            nprest = nprest - iessw
            nindt = nindt + iessw
        end if

        if (iesup .ne. 0) then
            nwkin = nwkin + iesup*in1
            nwking = nwking + iesup*in1
            nwrep = nwrep + in1*iesup
            nwfix = nwfix + iesup*in1
            nwnel = nwnel + in1*iesup
            nwdw = nwdw + in1*iesup
            nprest = nprest - iesup
            nindt = nindt + iesup
        end if

        if (iesking .ne. 0) then
            nwkin = nwkin + iesking*in1
            nwrep = nwrep + in1*iesking
            nwfix = nwfix + iesking*in1
            nwnel = nwnel + in1*iesking
            nwdw = nwdw + in1*iesking
            nprest = nprest - iesking
            nindt = nindt + iesking
        end if

        if (ieskin .ne. 0) then
            nwrep = nwrep + in1*ieskin*skipforce
            nwfix = nwfix + ieskin*in1*skipforce
            nwnel = nwnel + in1*ieskin*skipforce
            nwdw = nwdw + in1*ieskin*skipforce
            nprest = nprest - ieskin*skipforce
            nindt = nindt + ieskin*skipforce
        end if

        ! variance and Berry phase
        if (.not. yes_complex) then
            indfix = nindt + 1
            if (isfix .ne. 0) then
                nwrep = nwrep + in1*isfix
                nprest = nprest - isfix
                nindt = nindt + isfix
            end if
            indberry = nindt
            repf = nwrep/in1 + 1
            kinf = nwkin/in1 + 1
        else
            nindt = nindt*2
            indfix = nindt + 1
            if (isfix .ne. 0) then
                nwrep = nwrep + 2*in1*isfix
                nprest = nprest - isfix
                nindt = nindt + 2*isfix
            end if
            indberry = nindt
            repf = nwrep/in1 + 1
            kinf = nwkin/in1 + 1
        end if

        ! COMPLEX DEB TO modify
        if (itestr .ne. -5) then ! if the minimization is not assumed
            nrep = iesm + iesd + iesfree + iessw + iesinv + iesup + iesking
            nprest = nprest - nrep
            nindt = nindt + nrep
        else
            nrep = 0
        end if

        !---------------------------
        ! End of indices definitions
        !---------------------------

        if (rank .eq. 0) then
            write (6, *) ' nwfix given =', nwfix, nwrep
            write (6, *) ' Number of parameters in SR =', nindt
            write (6, *) ' Hartree Atomic Units '
        end if

        if (nprest .lt. 0) then
            write (errmsg, *) ' ERROR INCREASE NP !!! by ', -nprest, '=', np - nprest
            call checkiflagerr(1, rank, errmsg)
        elseif (nprest .gt. 0) then
            write (errmsg, *) ' ERROR DECREASE NP !!! by ', nprest, '=', np - nprest
            call checkiflagerr(1, rank, errmsg)
        end if

        !---------------------------------------------
        ! Definition of dimension andaddresses of the
        ! main matrices
        !---------------------------------------------

        nelnw = nel*nw*(indt + 1)
        nelnw0 = nel*nw
        nrnel = (indt + 1)*nel
        nelkel = 3*nrnel
        dimee = nel*nel*(indt4j + 1)
        dimei = nel*nion
        nrest = 0

        Lzeff = nel*indteff
        Lztab = nel*indt
        if (Lztab .eq. 0) Lztab = 1
        Lztabr = Lztab
        if (yes_complex) Lztab = Lztab*2
        Ltab = nel*(indt + ip4)
        Ltabb = max(nel*indt, 1)

        nwnp = in1*nmat

        ! save the starting ionic configurations
        call dscalzero(ieskindim, 0.d0, dek, 1)
        call dscalzero(ieskingdim, 0.d0, dekg, 1)

        if (iespbc .and. itestr .eq. -5) then
            if (ieskint .eq. ieskinr_pos + 2) then
                dek(ieskinr_pos + 1 - iesking) = cellscale(2)
                dek(ieskinr_pos + 2 - iesking) = cellscale(3)
                if (rank .eq. 0) write (6, *) ' PBC Updating cellscale b,c=', ieskinr_pos &
                        &, dek(ieskinr_pos + 1 - iesking), dek(ieskinr_pos + 2 - iesking)
            elseif (ieskint .eq. ieskinr_pos + 1) then
                dek(ieskinr_pos + 1 - iesking) = cellscale(1)
                if (rank .eq. 0) write (6, *) ' PBC Updating cellscale a=', ieskinr_pos   &
                        &, dek(ieskinr_pos + 1 - iesking)
            elseif (ieskint .eq. ieskinr_pos + 3) then
                dek(ieskinr_pos + 1 - iesking) = cellscale(1)
                dek(ieskinr_pos + 2 - iesking) = cellscale(2)
                dek(ieskinr_pos + 3 - iesking) = cellscale(3)
                if (rank .eq. 0) write (6, *) ' PBC Updating cellscale a,b,c=', ieskinr_pos  &
                        &, dek(ieskinr_pos + 1 - iesking), dek(ieskinr_pos + 2 - iesking), dek(ieskinr_pos + 3 - iesking)
            end if
        end if

        do j = 1, nw
            jbra(j) = j
        end do

        ! doubled indices for complex allocation
        ! only determinantal part
        if (yes_complex) then
            nel2wt = 2*nelorb*nel*(indt4 + 1) ! address for winv
            nel2up = 2*nelup*nelup ! address for ainv
            nel2upt = 2*nelup*(indt + ip4) ! address for winvup
            nel2dot = 2*neldo*(indt + ip4) ! address for winvdo
            nel2bar = 4*nelup*nelorbh ! address for winvbar
        else
            nel2wt = nelorb*nel*(indt4 + 1)
            nel2up = nelup*nelup
            nel2upt = nelup*(indt + ip4)
            nel2dot = neldo*(indt + ip4)
            nel2bar = 2*nelup*nelorbh
        end if

        ! indices for the jastrow
        nel2wtj = nelorbj*nel*(indt4j + 1)
        nel2jbar = (nelup + neldo)*nelorbjh*ipj
        if (iessz) then
            nel2jbarsz = nel2jbar
        else
            nel2jbarsz = 1
        end if

        nelnion = (nelup + neldo)*nion

        indtj_nw = Lztab*in1 + 1
        indtjr_nw = Lztabr*in1 + 1
        indtabj_nw = Ltab*in1 + 1
        indtabbj_nw = Ltabb*in1 + 1
        indkj_nw = nel*(indt + 1)*in1 + 1
        indksj_nw = nel*in1 + 1
        indksij_nw = nelnion*in1 + 1
        indwwj_nw = nel2wt*in1 + 1
        indwwjj_nw = nel2wtj*in1 + 1
        indwupj_nw = nel2upt*in1 + 1
        indwdoj_nw = nel2dot*in1 + 1
        indaupj_nw = nel2up*in1 + 1
        indbar_nw = nel2bar*in1 + 1

        indjbar_nw = nel2jbar*in1 + 1
        if (iessz) then
            indjbarsz_nw = indjbar_nw
        else
            indjbarsz_nw = 1
        end if

        j_nw = in1 + 1 + istm
        j_nws = in1 + 1

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !cccccc parallelization: each process goes from
        !cccccc ist to ien walkers
        ist = rank*(nw/nproc) + 1
        istm = ist - 1
        ien = (rank + 1)*(nw/nproc)
        id1 = rank*(nw/nproc)
#ifdef PARALLEL
        skipreshuff = 15 + Ltab + Ltabb + Lztab + Lztabr + nel2upt + nel2dot                       &
             &+ 8*nel + nelkel + npf
        if (yesivic) skipreshuff = skipreshuff + 3*indt*nel
        if (npsa .gt. 0) skipreshuff = skipreshuff + nel + 9*nel

        skip = skipreshuff*in1
        if (iscramax .lt. skip .and. iesbra) then
            iscramax = skip
            deallocate (psip)
            allocate (psip(skip))
            if (rank .eq. 0) write (6, *) 'iscramax changed!', skip
        end if
#endif
        if (iscramax .lt. nparshellmax + nelorbpp) iscramax = nparshellmax + nelorbpp
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !if(rank.eq.0) write(6,*) ' after 4 ',nelorb,nelorbh

        ! Initialization
        table = 0.d0
        tabler = 0.d0

        if (itestr .eq. -2) then
            ! heat bath after all electron swap
            ! identity splitted into nel particle contributions
            identity = 1.d0/dble(nel)/tbra
        elseif (itestr .eq. -3) then
            ! heat bath after single particle diffusion
            ! non local operator splitted into a product of nel factors < e^{-tbra/nel V_nl} >
            ! tbra --> tbra/nel
            identity = 1.d0/tbra
        end if

        if (itestr .eq. -2 .or. itestr .eq. -3) then
            ! fill table and tabler with the identity in the last entry per particle
            ! this value will never be touched during the simulation
            do kk = 1, nws
                indtj = nel*indt*(kk - 1)
                if (yes_complex) then
                    do k = 1, 2*nel
                        table(k + 2*(indt - 1)*nel + 2*indtj) = dcmplx(identity)
                    end do
                    do k = 1, nel
                        tabler(k + (indt - 1)*nel + indtj) = identity
                    end do
                else
                    do k = 1, nel
                        table(k + (indt - 1)*nel + indtj) = identity
                        tabler(k + (indt - 1)*nel + indtj) = identity
                    end do
                end if
            end do
        end if

        bcost = 1.d0
        ccost = 1.d0
        !     Initialize to zero the global variables
        kel = 0.d0
        wconfn = 0.d0

        if (.not. iesrandoma .and. itest .ne. 2) then
            call findrionfref(nion, alat, rion, rion_ref, mind(1))
            call findrionfref(nion, alat, rion(2, 1), rion_ref(2), mind(2))
            call findrionfref(nion, alat, rion(3, 1), rion_ref(3), mind(3))
            if (rank .eq. 0) write (6, *) ' Minimum ion-mesh distance =', dsqrt(sum(mind(:)**2))
        end if

        call initconf(nel, nelup, psiln, psisn, kel, wconfn, rion, dist&
                &, indt, nion, zetar, in1, rank, ierr, LBox, alat, rion_ref, iesrandoma, itest)

        imin = 0
        itry = 0

        psisn = 0.d0
        countscra = 0
        spsi = 1.d0
        vcut = 0.d0
        yescut = .false.
        flagmu = .false.
        wconfsav(:) = 1.d0

        ntry = 0
        nmatb = 0
        ttry = 0.d0
        ngentry = 0
        rata = 0.d0
        signflip = 0.d0
        nontr = 0.d0
        naccpseudo = 0.d0
        tmes = temp*ris(2)
        epscutu = epscut
        epstlu = epstl
        ngs = 0
        nmp = 0
        ng = 0
        icdiff = 0
        wbra = 1.d0
        wbran = 1.d0
        wbra1 = 0.d0
        wbra2 = 0.d0
        varmin = 0.d0
        enermin = 0.d0
        acclarge = 0.d0
        indvic = 1
        parbest = 0.d0
        iflagnorm = 3
        ! iflagnorm=3 compute the normalization of orbitals and distances
        ! iflagnorm=2 compute only distances (rmu = el-nu , r=|el-nu|)
        coeff = 0.d0
        iesconv = 0
        jmax = 0
        fmax = 0.d0

        if (mod(nbra, nel) .ne. 0 .and. fncont) then
            nbra = ((nbra - 1)/nel + 1)*nel
            if (rank .eq. 0)                                                &
                    &     write (6, *) ' Warning changing nbra multiple of nel =', nbra
        end if
        nbram = nbra - 1

        call dscalzero(nelorbh*(indt + 1 + ip4), 0.d0, psinew, 1)
        call dscalzero(iscramax, 0.d0, psip, 1)

        do j = 1, nwnp
            econf(j) = 1.d0
            econfh(j) = 0.d0
        end do
        do j = 1, in1
            econfh(in1*(npp - 1) + j) = 1.d0
        end do

        if (iesd .gt. 2 .or. (iesd .gt. 1 .and. iesdr .gt. -5)) then
            if (rank .eq. 0)                                                 &
                    &    write (6, *) ' Warning the hessian  does not work !!!'
        end if

        iend = 0
        costpassed = 1.d0

        do i = 1, nfat
            wcorw(i) = 1.d0
        end do
        indcor = nfat
        wcort = 1.d0

        psiav = 0.d0
        psisav = 0.d0
        psirav = 0.d0
        countav = 0.d0
        counttot = 0.d0
        countcut = 0.d0
        countreg = 0.d0

        if (idyn .gt. 0) call dscal(ieskin, ris(2), scalpar(np - ieskin + 1), 1)
        dt = 0.d0
        ii = 1
        !           The first one non zero
        do while (dt .eq. 0.d0 .and. ii .le. ieskin)
            dt = scalpar(np - ieskin + ii)
            ii = ii + 1
        end do
        if (rank .eq. 0 .and. idyn .gt. 0) write (6, *) ' Warning Chosen Dt /component=', dt/ris(2), ii - 1
        if (dt .eq. 0 .and. idyn .gt. 0) then
            if (rank .eq. 0) write (6, *) ' ERROR no component to do dynamics  Dt=0!!!!', dt
#ifdef PARALLEL
            call mpi_finalize(ierr)
#endif
            stop
        end if

        ! iopt=0   continue the run starting from a given conf
        ! iopt=1   begin the calculation from scratch
        ! iopt=2   begin with old conf but restart averages
        if (iopt .eq. 0 .or. iopt .eq. 2) then

            if (rank .eq. 0) call read_fort11_begin
            call checkiflagerr(iflagerr, rank, 'ERROR in read_fort11_begin')
#ifdef PARALLEL
            call mpi_bcast(ngenc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nwr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nbrar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(npow, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(etryr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nrest, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nwr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(tbrar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(npr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(npbrar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(ireadr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nmpr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nfatr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(np3r, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(npsovr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(wcort, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nweightr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(ngs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(iendr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(indcor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(parbest, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
                 &, ierr)
            call mpi_bcast(iflagerr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(enermin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
                 &, ierr)
            call mpi_bcast(varmin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD     &
                 &, ierr)
            call mpi_bcast(iesconv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(fmax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(jmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(inext, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(pippoc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(pippo, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nprocr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(nbinread, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(stepcg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(ncgread, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(epscutur, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD   &
                 &, ierr)
            call mpi_bcast(epstlur, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
                 &, ierr)
            call mpi_bcast(tmes, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD       &
                 &, ierr)
#endif

            if (iopt .eq. 0 .and. (nprocr .eq. nproc)) then
                if (io_level .eq. 1) iflagerr = read_seeds(cstr(trim(ranseedfilename)))
#ifdef PARALLEL
                if (io_level .eq. 2) iflagerr = read_seeds_mpiio(cstr(trim(scratchpath)//'randseed.all'), rank)
#endif
                call checkiflagerr(iflagerr, rank, "Fail in reading random seeds!")
            elseif (iopt .eq. 0) then
                if (rank .eq. 0) then
                    write (6, *) ' Warning continuing with different number of processors '
                    write (6, *) ' Warning the random number are initialized again '
                end if
            end if

            if (itestr .eq. -5) then

                if (iopt .ne. 2) then
                    epscutu = epscutur
                    epstlu = epstlur
                else
                    epscutu = epscut
                    epstlu = epstl
                end if

                if (rank .eq. 0) write (6, *) ' Used epscut,epstl =', epscutu, epstlu
                if (iopt .eq. 2) then
                    if (writescratch .eq. 0) rewind (19)
                    if (abs(kl) .eq. 9) rewind (20)
                else
                    irec = 0
                    if (writescratch .eq. 0) rewind (19)
                    if (abs(kl) .eq. 9) rewind (20)
                    if (rank .eq. 0) then
                        write (6, *) ' record read =', irec
                    end if
                    if (irec .ne. iendr - inext + nweightr) then
                        if (rank .eq. 0) write (6, *) ' Too few number of records read !!! ', irec
                        !#ifdef PARALLEL
                        !         write(6,*) ' Processor =',rank+1,nweightr,iendr,inext
                        !        call mpi_finalize(ierr)
                        !#endif
                        !         stop
                        iendr = inext - nweightr
                        ngs = iendr
                        pippo = 1
                        pippoc = 1
#ifdef PARALLEL
                        jj = inext
                        call mpi_allreduce(jj, inext, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
                        jj = ngs
                        call mpi_allreduce(jj, ngs, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
                        jj = iendr
                        call mpi_allreduce(jj, iendr, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif

                        if (rank .eq. 0) write (6, *) ' Warning changed iendr,ngs,inext'&
                                &, iendr, ngs, inext

                    end if

                end if
            end if

            nmatr = npr + 1

            if (iopt .eq. 0) then
                !         nfat=nfatr
                !         npsov=npsovr
                iend = iendr
                !         np3r=np3
                !         np=npr
                !         npbra=npbrar
                ! fine iopt.eq.0
            end if
            !---------------------------------------------------

            if (rank .eq. 0) call read_fort11_end
            call checkiflagerr(iflagerr, rank, 'Error in reading fort.11')

            if (nwr .eq. nw) then
#ifdef PARALLEL
                if (yesquantum) then
                    if (rankrep .eq. 0) then
                        if (io_level .eq. 1) then
                            rewind (8)
                            if (contraction .ne. 0) then
                                read (8) rion, velion, alphavar, detmat_c, projm, mu_c&
                                  &, jasmat, jasmatsz, jasmat_c, jasmatsz_c, vjur, dupr
                            else
                                read (8) rion, velion, alphavar, detmat, projm, mu_c&
                                  &, jasmat, jasmatsz, jasmat_c, jasmatsz_c, vjur, dupr
                            end if
                            if (allocated(muj_c)) read (8) muj_c
                            if (detc_proj) read (8) detmat_proj, projmat_c
                            if (itestr .eq. -5) read (8) reduce
                        elseif (io_level .eq. 2) then
                            call MPI_File_read_all(quantcont%fp, rion, size(rion), MPI_DOUBLE_PRECISION, status, ierr)
                            call MPI_File_read_all(quantcont%fp, velion, size(velion), MPI_DOUBLE_PRECISION, status, ierr)
                            call MPI_File_read_all(quantcont%fp, alphavar, size(alphavar), MPI_DOUBLE_PRECISION, status, ierr)
                            if (contraction .ne. 0) then
                                call MPI_File_read_all(quantcont%fp, detmat_c, size(detmat_c), MPI_DOUBLE_PRECISION, status, ierr)
                            else
                                call MPI_File_read_all(quantcont%fp, detmat, size(detmat), MPI_DOUBLE_PRECISION, status, ierr)
                            end if
                            call MPI_File_read_all(quantcont%fp, projm, size(projm), MPI_DOUBLE_PRECISION, status, ierr)
                            call MPI_File_read_all(quantcont%fp, mu_c, size(mu_c), MPI_DOUBLE_PRECISION, status, ierr)
                            call MPI_File_read_all(quantcont%fp, jasmat, size(jasmat), MPI_DOUBLE_PRECISION, status, ierr)
                            call MPI_File_read_all(quantcont%fp, jasmat_c, size(jasmat_c), MPI_DOUBLE_PRECISION, status, ierr)
                            call MPI_File_read_all(quantcont%fp, jasmatsz, size(jasmatsz), MPI_DOUBLE_PRECISION, status, ierr)
                            call MPI_File_read_all(quantcont%fp, jasmatsz_c, size(jasmatsz_c), MPI_DOUBLE_PRECISION, status, ierr)
                            call MPI_File_read_all(quantcont%fp, vjur, size(vjur), MPI_DOUBLE_PRECISION, status, ierr)
                            call MPI_File_read_all(quantcont%fp, dupr, size(dupr), MPI_DOUBLE_PRECISION, status, ierr)
                            if (allocated(muj_c))& !  read(8) muj_c
                             &call MPI_File_read_all(quantcont%fp, muj_c, size(muj_c), MPI_DOUBLE_PRECISION, status, ierr)
                            if (detc_proj) then ! read(8) detmat_proj,projmat_c
                                call MPI_File_read_all(quantcont%fp, detmat_proj, size(detmat_proj), MPI_DOUBLE_PRECISION &
                                                       , status, ierr)
                                call MPI_File_read_all(quantcont%fp, projmat_c, size(projmat_c), MPI_DOUBLE_PRECISION &
                                                       , status, ierr)
                            end if
                            if (itestr .eq. -5) & !  read(8) reduce
                             &call MPI_File_read_all(quantcont%fp, reduce, size(reduce), MPI_DOUBLE_PRECISION, status, ierr)
                        end if
                    end if ! endif rankrep

                    if (yescomm) then
                        call bcast_real(rion, size(rion), 0, commrep_mpi)
                        call bcast_real(velion, size(velion), 0, commrep_mpi)
                        call bcast_real(alphavar, size(alphavar), 0, commrep_mpi)
                        if (contraction .ne. 0) then
                            call bcast_real(detmat_c, size(detmat_c), 0, commrep_mpi)
                        else
                            call bcast_real(detmat, size(detmat), 0, commrep_mpi)
                        end if
                        call bcast_real(projm, size(projm), 0, commrep_mpi)
                        call bcast_real(mu_c, size(mu_c), 0, commrep_mpi)
                        call bcast_real(jasmat, size(jasmat), 0, commrep_mpi)
                        call bcast_real(jasmatsz, size(jasmatsz), 0, commrep_mpi)
                        call bcast_real(jasmat_c, size(jasmat_c), 0, commrep_mpi)
                        call bcast_real(jasmatsz_c, size(jasmatsz_c), 0, commrep_mpi)
                        if (allocated(muj_c)) call bcast_real(muj_c, size(muj_c), 0, commrep_mpi)
                        call bcast_real(vjur, size(vjur), 0, commrep_mpi)
                        call bcast_real(dupr, size(dupr), 0, commrep_mpi)
                        if (detc_proj) then
                            call bcast_real(detmat_proj, size(detmat_proj), 0, commrep_mpi)
                            call bcast_real(projmat_c, size(projmat_c), 0, commrep_mpi)
                        end if
                        if (itestr .eq. -5) call bcast_real(reduce, size(reduce), 0, commrep_mpi)
                    end if
                    if (yesdetmatc) then
                        call scontract_mat_det(nelorbh, nelorbh, nelcolh, nelorb_c&
                       &, nelcol_c, detmat, detmat_c, mu_c, psip)
                    end if
                end if

                if (io_level .eq. 1) then
                    rewind (9)
                    read (9) (kel(1:3, (ii - 1)*nrnel + 1:(ii - 1)*nrnel + nel), ii=1, in1), angle
                else if (io_level .eq. 2) then
                    if (rank .eq. 0) write (6, *) " Reading electronic configurations from ", trim(kelcont%name)
                    do ii = 1, in1
                        call MPI_File_read_all(kelcont%fp, kel(1, (ii - 1)*nrnel + 1), nel3,&
                   & MPI_DOUBLE_PRECISION, status, ierr)
                    end do
                    call MPI_File_read_all(kelcont%fp, angle, size(angle), MPI_DOUBLE_PRECISION, status, ierr)
                end if
                !write(6,*) "check Read kel 1", rank, kel(1:3,1)
#else
                read (11) kel, angle
#endif
            end if

            call checkiflagerr(iflagerr, rank, 'Error in main')
            if (rank .eq. 0) write (6, *) ' All files are read correctly ...'
#ifdef PARALLEL
            if (contraction .gt. 0) &
                 &call mpi_bcast(allowcontr, 2*nelorb_c, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#ifdef UNREL
            call bcast_real(wconfn, nw, 0, MPI_COMM_WORLD)
            call bcast_real(wcorw, nfatr, 0, MPI_COMM_WORLD)
            if (.not. yesquantum) then
                call bcast_real(dek, ieskindim, 0, MPI_COMM_WORLD)
                call bcast_real(dekg, ieskingdim, 0, MPI_COMM_WORLD)
            end if
            if (detc_proj .and. .not. yesquantum) then
                call bcast_real(detmat_proj, size(detmat_proj), 0, MPI_COMM_WORLD)
                call bcast_real(projmat_c, size(projmat_c), 0, MPI_COMM_WORLD)
            end if
            if (idyn .eq. 17) call bcast_real(gp, size(gp), 0, MPI_COMM_WORLD)
#else
            call mpi_bcast(wconfn, nw, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD    &
           &, ierr)
            call mpi_bcast(wcorw, nfatr, MPI_DOUBLE_PRECISION, 0                 &
           &, MPI_COMM_WORLD, ierr)
            if (.not. yesquantum) then
                call mpi_bcast(dek, ieskindim, MPI_DOUBLE_PRECISION, 0               &
               &, MPI_COMM_WORLD, ierr)
                call mpi_bcast(dekg, ieskingdim, MPI_DOUBLE_PRECISION, 0             &
               &, MPI_COMM_WORLD, ierr)
            end if
            if (detc_proj .and. .not. yesquantum) then
                call mpi_bcast(detmat_proj, size(detmat_proj), MPI_DOUBLE_PRECISION, 0&
                &, MPI_COMM_WORLD, ierr)
                call mpi_bcast(projmat_c, size(projmat_c), MPI_DOUBLE_PRECISION, 0&
                &, MPI_COMM_WORLD, ierr)
            end if
            if (idyn .eq. 17) call mpi_bcast(gp, size(gp), MPI_DOUBLE_PRECISION, 0&
                &, MPI_COMM_WORLD, ierr)
#endif
            if (itestr .eq. -5) then
#ifdef UNREL
!   For unreliable  networks.
                if (.not. yesquantum) then
                    call bcast_real(velion, ieskindim*3, 0, MPI_COMM_WORLD)
                    call bcast_real(reduce, ncgdim*npdim, 0, MPI_COMM_WORLD)
                end if
                call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp barrier
#else
                if (.not. yesquantum) then
                    call mpi_bcast(velion, ieskindim*3, MPI_DOUBLE_PRECISION            &
                   &, 0, MPI_COMM_WORLD, ierr)
                    call mpi_bcast(reduce, ncgdim*npdim, MPI_DOUBLE_PRECISION           &
                   &, 0, MPI_COMM_WORLD, ierr)
                end if
#endif
                if (rank .eq. 0 .and. idyn .gt. 2) write (6, *) ' Temp velocities read ='&
                     &, dnrm2(ieskin, velion(3, 1), 3)**2/ieskin*ris(2)
            end if
#endif

            alphab = 0.d0

            if (nfat .ne. nfatr) then
                do i = 1, nfat
                    wcorw(i) = 1.d0
                end do
            end if

            ng = ngenc
            if (iopt .eq. 2) then
                ng = 0
                ngs = 0
            end if

            if (rank .eq. 0) then
                write (6, *) ' Number of QMC iterations  in this run', ngen
                if (itestr .eq. -5) write (6, *) ' Number of previous optimization steps', ng
                write (6, *) ' Number of previous QMC iterations ', ngs
            end if

            ! ngs number of optimization steps, ng number of QMC iterations.
            ! rewind the old file to continue to write if iopt.eq.0

            ! fine if iopt=0 or iopt=2
        end if

        call read_alphavar

        if (rank .eq. 0 .and. ieskin .ne. 0) then
            write (6, *) ' Warning rescaling acc forces ion  in Hartree'
        end if

        if (idyn .eq. 0) delta0 = 0.d0

        if (idyn .eq. 4) then
            if (delta0 .lt. 4.d0/3.d0*dt*normcorr) then
                !          if(delta0.lt.2.d0*dt) then
                !          delta0=2.d0*dt
                delta0 = 4.d0/3.d0*dt*normcorr
                if (rank .eq. 0) write (6, *) ' Warning  delta_0 >= 4/3 dt !!! Changed to ' &
                        &, delta0
            end if
        end if

        if (idyn .eq. 6) then
            if (delta0 .lt. dt*normcorr) then
                !          if(delta0.lt.2.d0*dt) then
                !          delta0=2.d0*dt
                delta0 = dt*normcorr
                if (rank .eq. 0) write (6, *) ' Warning  delta_0 >= dt !!! Changed to ' &
                        &, delta0
            end if
        end if

        if (idyn .eq. 3) then
            if (delta0 .lt. dt*normcorr) then
                delta0 = dt*normcorr
                if (rank .eq. 0) write (6, *) ' Warning  delta_0 >  dt !!! Changed to ', delta0/dt, 'dt'
            end if
        end if

        if (rank .eq. 0) then
            write (6, *) ' Initial scaling '
            do i = 1, np
                if (scalpar(i) .ge. 0.d0) write (6, *) i, scalpar(i)
            end do
        end if

        srforce = 0.d0
        srforcew = 0.d0
        wforce = 0.d0
        wforcew = 0.d0

        !          walkers with weight wconf=1
        !   if etry=0 set etry to the classical energy per site
        nacc = 0.d0
        do j = 1, nws
            naccm(j) = 1
        end do
        sumdiff = 0.d0

        time_ratiovar = 0.d0
        time_uptabtot = 0.d0
        time_meas = 0.d0
        timemc = 0.d0
        timeopt = 0.d0
        timescra = 0.d0
        timewf = 0.d0
        timepip = 0.d0
        if (itest .eq. 2) then
            lambda = dabs(-etry)
        else
            lambda = -etry
            lambda = lambda*(1.d0 - rsignr)

            if (rank .eq. 0) write (6, *) ' lambda chosen =', lambda
        end if

        ngg = ng
        ngn = ngs

        do j = ist, ien
            wconfw(j - istm) = wconfn(j)
        end do

        rweight = nweight

        if (rank .eq. 0) then
            write (6, *) ' nweight before starting ', nweight
        end if

        !       write(6,*) 'before  nweight nweightr # proc  ='
        !    1,nweight,nweightr,rank+1

        ndone = nweight - ibinit
        lbin = ndone/nbinr

        if (lbin .eq. 0) then
            write (errmsg, *) ' ERROR # bins larger than nweight,&
                    & reduce nbinr in input !!! <=', ndone
            call checkiflagerr(1, rank, errmsg)
        end if
        flag = .false.
        if (nbinr*lbin .ne. ndone .and. itestr .eq. -5) flag = .true.
        ndone = nbinr*lbin
        nweight = ndone + ibinit
        if (rank .eq. 0 .and. flag)                                                 &
                &write (6, *) ' Warning nweight changed  !!!', nweight

        i = iend
        if (iopt .ne. 0 .or. nweightr .ne. nweight) then
            inext = iend + nweight
            pippoc = 1
            pippo = 1
        end if

        ! Be sure that all walkers starts with the wf initialized correctly.
#ifdef  PARALLEL
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

        if (iesm + abs(iesfree) + abs(iesinv) + abs(iessw) + iesup + ieskin + iesking .gt. 0) then
            iesmeas = .true.
        else
            iesmeas = .false.
        end if

        !       allocating the minimum memory for ipsip

        if (itestr .ne. -5) then
            iscraipsip = nelup + 9*nion
            if (iesbra) iscraipsip = max(iscraipsip, 2*nw)
            deallocate (ipsip)
            allocate (ipsip(iscraipsip))
            ipsip = 0
        end if

#if defined (_OPENMP) && defined (__NOOMP)
        call omp_set_num_threads(old_threads) ! restore the previous threads
#endif
        nelorbjmax = max(nelorbj, 1)
        neldomax = max(neldo, 1)
        indtmax = max(indt, 1)
        nshelljmax = max(nshellj, 1)

        isdistp = 20*(indt + 1)
        !  ALLOCATE all required for reverse mode.

        iscramaxb = isdistp
        iscramaxb = max(iscramaxb, (nion*(nion + 1))/2)

        !  ALLOCATE all required for reverse mode.
        if (iessz) then
            iscramaxb = max(iscramaxb, 2*max(nel, indt + 5) + 2*max(nelorbjh, 1) + 20*(indt + 1)&
                    & + max(nelorbjh, 1)*(2*indt + 11))
        else
            iscramaxb = max(iscramaxb, 2*max(nel, indt + 5) + max(ipj*nelorbjh, 1) + 20*(indt + 1)&
                    & + max(nelorbjh, 1)*(2*indt + 11))
        end if
        !  task9
        if (iessz) then
            isdistp = 27*(indt + 1)*nshelljmax + nelorbjmax*(indt + 6) + 1 + (7 + indt)*nelorbjmax + 2*max(nel, indt + 5)
        else
            isdistp = 27*(indt + 1)*nshelljmax + nelorbjmax*(indt + 6) + 1 + (ipj + 5 + indt)*nelorbjmax + 2*max(nel, indt + 5)
        end if
        iscramaxb = max(iscramaxb, isdistp)

        !   NB in TASK8_B reverse mode psip is not used

        iscramaxb = max(iscramaxb, nparshellmax) ! task2 nparshellmax defined  after read_pseudo -->

        iscramaxb = max(iscramaxb, 2*ipc*nelup_mat*nelorb_c) ! task5 upper bound

        iscramaxb = max(iscramaxb, 2*nelorbj_c*nel) ! task5 Jastrow
        iscramaxb = max(iscramaxb, (nion*(nion - 1))/2, ipc*nelup*nelup)
        iscramaxb = max(iscramaxb, ipc*(indt + 5)*nelorbh + nelorbh*(indt + 6) + 27*(indt + 1)*nshell + nelorbh) ! task10

        !   NB in TASK8_B reverse mode psip is not used

        iscramaxb = max(iscramaxb, nparshellmax) ! task2 nparshellmax defined  after read_pseudo -->

        iscramaxb = max(iscramaxb, 2*ipc*nelup_mat*nelorb_c) ! task5 upper bound
        iscramaxb = max(iscramaxb, 2*nelorbj_c*nel) ! task5 Jastrow
        iscramaxb = max(iscramaxb, (nion*(nion - 1))/2, ipc*nelup*nelup)
        iscramaxb = max(iscramaxb, ipc*(indt + 5)*nelorbh + nelorbh*(indt + 5) + 20*(indt + 1) + nelorbh) ! task10

        iscramaxold = iscramax
        iscramax = iscramaxb

        write (6, *) ' iscramax before =', iscramax

        if (yesfast .ne. 0) then
            iscramax = max(iscramax, nmolfn*(nelup + nelorbh))
        else
            iscramax = max(iscramax, ipf*ipf*nelorbh*nelorbh)
        end if
        iscramax = max(iscramax, nelorbjh*nelorbjh)

        if (iscramaxold .lt. iscramax) then
            deallocate (psip)
            allocate (psip(iscramax))
            if (rank .eq. 0) write (6, *) ' Warning increasing allocation scratch for&
                    & psip', iscramax
            iscramaxb = iscramax
        end if

        if (epsder .eq. 0 .and. ieskint .ne. 0) then
            if (rank .eq. 0) write (6, *) ' DIFFERENTIAL WARP ALGORITHM FOR FORCES  '
        end if

        !        write(6,*) ' Initial tabs AGP '
        !        do ii=1,nshell
        !        write(6,*) ii,indpar_tab(ii),indshell_tab(ii),indorb_tab(ii)
        !        enddo

        !        write(6,*) ' Initial tabs Jastrow  '
        !        do ii=1,nshellj
        !        write(6,*) ii,indparj_tab(ii),indshellj_tab(ii),indorbj_tab(ii)
        !        enddo

        !         stop

        if (rank .eq. 0 .and. .not. yespulay) &
                &write (6, *) ' Warning calculation of pulay suppressed '

        if (ieskin .ne. 0) then
            ieskinion = ieskin - (ieskint - ieskinr_pos)
            if (rank .eq. 0) write (6, *) ' Number of components ions =', ieskinion
        end if
        if (itest .ne. 2) then
            if (yes_complex) then
                allocate (winvfn(2*nmolfn*(indt + ip4)*nel*nws), winvbarfn(nmolfn*2*nel_mat*nws))
                nel2wtfn = 2*nmolfn*nel*(indt + ip4)
                nel2barfn = 4*nelup*nmolfn
            else
                allocate (winvfn(nmolfn*(indt + ip4)*nel*nws), winvbarfn(nmolfn*nel_mat*nws))
                nel2wtfn = nmolfn*nel*(indt + ip4)
                nel2barfn = 2*nelup*nmolfn
            end if
            !        indbarfn=nel2barfn*(j-1)+1

            Ltot = 20 + nelorb*(indt4 + 1) + 2*nion + indt + nelorbj*(indt4j + 1) + 2*nel* &
                    &(indt4j + 1) + nmolfn*(indt + ip4 + 1) + Ltab + nel2upt + nel2dot + nel2up + nelorbjh
            if (iessz) Ltot = Ltot + nelorbjh
            if (iespbc) Ltot = Ltot + 2*n_gvec

        else
            nel2wtfn = 0
            nel2barfn = 0
            if (yes_complex) then
                allocate (winvfn(2), winvbarfn(2))
            else
                allocate (winvfn(1), winvbarfn(1))
            end if
        end if
        indbarfn_nw = nel2barfn*in1 + 1
        indwwjfn_nw = nel2wtfn*in1 + 1
        winvbarfn = 0.d0
        winvfn = 0.d0

        ! by E. Coccia (8/11/10): read_cube and spline interpolation
        if (ext_pot) then
            call extpot_read
            ! by E. Coccia (4/1/11): spline interpolation for the forces
            ! geometry optimization (idyn.ne.0)
            ! only for AD (epsder.eq.0)
            if (idyn .ne. 0 .and. epsder .eq. 0) then
                call forces_interpolate()
            end if
        end if
        ! by E. Coccia (10/12/11): read restr.dat
        if (mm_restr) then
            call restr_read()
        end if

        !by E. Coccia (20/12/11): writing electronic random walk
        if (rank .eq. 0 .and. write_rwalk) open (671, file='rwalk.xyz')

        !by E. Coccia (8/4/14)
        ! Opening clients sockets for NEB calculations
#ifdef NEB
        !Opening the socket
        if (rank .eq. 0) then
            stop_var = 0
            open (740, file='client.info')
            read (740, *) iclient
            read (740, *) ncycles
            read (740, *) stop_var

            write (*, *)
            write (*, *) '********************************'
            write (*, *) 'NEB algorithm with client/server'
            write (*, *) 'This job is the client no.;', iclient
            write (*, *) '********************************'
            write (*, *)
            close (740)
        end if
        allocate (r_sock(nion, 3))
        allocate (f_sock(nion, 3))

        !Initialization for f_sock and r_sock
        f_sock = 0.d0
        do isock = 1, nion
            do ksock = 1, 3
                r_sock(isock, ksock) = rion(ksock, isock)
            end do
        end do

        neb_converged = .false.

#endif

        !if(rank.eq.0) write(6,*) ' after 6 ',nelorb,nelorbh

        yeswrite = .false.
        !       definition jbrasymiesup

        if (symiesup) then
            allocate (jbrasymiesup(iesupr_c))
            jbrasymiesup = 0
            ind = 0
            do jj = 1, iesupind
                ind = ind + 1
                ii = jbraiesup_sav(ind)
                do j = 1, abs(ii)
                    if (jbraiesup_sav(ind + j) .gt. 0) then
                        jbrasymiesup(jbraiesup_sav(ind + j)) = jj
                    else
                        jbrasymiesup(-jbraiesup_sav(ind + j)) = -jj
                    end if
                end do
                ind = ind + abs(ii)
            end do
        end if

        acc_dyn = .true.

        if (scalepulay .le. 0) then
            scalepulay = -scalepulay
            scalepulayw = scalepulay
            if (warp .and. epsder .eq. 0.d0) then
                scalepulay = 1.d0
            end if
        else
            scalepulayw = 1.d0
        end if
        if (min(scalepulayw, scalepulay) .ne. 1.d0 .and. rank .eq. 0) &
                &write (6, *) ' Warning biased Pulay scheme !!!  ', min(scalepulay, scalepulayw)
        dt4 = 1.d0
        if (yesquantum) then
            temp = temp*nbead
            allocate (rionall(3, nion, nbead))
            allocate (fbead(3, nion))
            allocate (mass_ion(3, nion))
            dt4 = 4.d0/temp
            if (rank .eq. 0) write (6, *) ' Mass particles (unit m_e) '
            do ii = 1, 3
                do jj = 1, nion
                    mass_ion(ii, jj) = atom_weight(nint(atom_number(jj)))/mass_unit*scale_mass
                    if (rank .eq. 0 .and. ii .eq. 1) write (6, *) jj, mass_ion(ii, jj)
                    !         Use the same mass in the classical dynamics
                end do
            end do
            if (yesturboq) then
                allocate (kdyn(nbead, nbead), kdyn_eig(nbead))
                kdyn = 0.d0
                do ii = 1, nbead
                    kdyn(ii, mod(ii, nbead) + 1) = -1.d0
                    kdyn(mod(ii, nbead) + 1, ii) = -1.d0
                    kdyn(ii, ii) = 2.d0
                end do
                !          Only for monoatomic species.
                kdyn = kdyn*temp**2*0.5d0*mass_ion(1, 1) ! The factor 1/2 comes from the Rydberg units
                !          setting all the quantum masses equal to the first one
                do ii = 1, 3
                    do jj = 1, nion
                        ind = ii + (jj - 1)*3
                        if (ind .ne. 1) then
                            scalpar(np - ieskin + ind) = scalpar(np - ieskin + 1)*mass_ion(1, 1)/mass_ion(ii, jj)
                        end if
                    end do
                end do

                if (size(psip) .lt. 3*nbead) then
                    iscramax = 3*nbead
                    deallocate (psip)
                    allocate (psip(iscramax))
                    psip = 0.d0
                end if

                call dsyev('V', 'L', nbead, kdyn, nbead, kdyn_eig&
                        &, psip, 3*nbead, info)
                    !!         Identity test
                !          kdyn=0.d0
                !          do ii=1,nbead
                !          kdyn(ii,ii)=1.d0
                !          kdyn_eig(ii)=0.d0
                !          enddo
#ifdef     PARALLEL
!          All should have exactly the same matrix to avoid roundoff
                call bcast_real(kdyn_eig, nbead, 0, MPI_COMM_WORLD)
                call bcast_real(kdyn, size(kdyn), 0, MPI_COMM_WORLD)
#endif
                if (rank .eq. 0) then
                    write (6, *) ' Eigenvalues elastic quantum term '
                    do ii = 1, nbead
                        write (6, *) ii, kdyn_eig(ii)
                    end do
                end if
            end if
        end if
        dimfk = nbin
        perbin = 1
        ndimj = iesinv + iesm + iesd + iesfree
        ndimjp = ndimj + 1
        if (ieskint .ne. 0) then
            allocate (rion_write(3, nion))
            rion_write = rion
        end if
        celldm_write(1:3) = celldm(1:3)
        rs_write = rs
        if (rank .eq. 0) write (6, *) ' Warning  epsder = ', epsder
        if (yes_sparse .and. rank .eq. 0 .and. .not. iessz .and. contractionj .eq. 0 .and. nelorbjh .ne. 0) then
            write (6, *) ' Warning using SPARSE matrix algorithm for Jastrow '
        end if
    end subroutine Initializeall
    subroutine Finalizeall
        implicit none
        real*8 drand1, enercont, jacobian, dnrm2, mapping, reweight_fn
#ifdef __NOOMP
        integer, external :: omp_get_max_threads
        old_threads = omp_get_max_threads()
        call omp_set_num_threads(1) ! scalar code
#endif

        !     at the end write the final  configuration file
        !     after the branching

        iend = i

#ifdef PARALLEL
!  sum the acceptances for each process
        call mpi_reduce(nacc, naccmpi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0 &
   &, MPI_COMM_WORLD, ierr)
        call mpi_reduce(naccpseudo, naccpseudompi, 1                    &
   &, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(nontr, cost, 1                                  &
   &, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(acclarge, acclargempi, 1                        &
   &, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
! communicate to the master the configurations needed for the
!   final write
        call mpi_gather(kel(1, nel*(indt + 1)*(ist - 1) + 1), 3*nel*(indt + 1)*in1 &
      &, MPI_DOUBLE_PRECISION, kel(1, nel*(indt + 1)*(ist - 1) + 1)               &
      &, 3*nel*(indt + 1)*in1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        call mpi_gather(wconfn(ist), in1                                  &
      &, MPI_DOUBLE_PRECISION, wconfn(ist)                                 &
      &, in1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        if (rank .eq. 0) then
            nacc = naccmpi
            naccpseudo = naccpseudompi
            nontr = cost
            acclarge = acclargempi
        end if
        call mpi_reduce(psiav, cost, 1                                    &
     &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (rank .eq. 0) psiav = cost

        call mpi_reduce(counttot, cost, 1                                 &
     &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (rank .eq. 0) counttot = cost

        call mpi_reduce(countcut, cost, 1                                 &
     &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (rank .eq. 0) countcut = cost

        call mpi_reduce(countreg, cost, 1                                 &
     &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (rank .eq. 0) countreg = cost

        call mpi_reduce(psisav, cost, 1                                   &
     &  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (rank .eq. 0) psisav = cost
#endif

        !      writes at the end the new wavefunctions found
        !       replace the changed orbitals

123     format(1000000e15.7)

        if (rank .eq. 0) then
            close (10)
            close (11)
            !---------------------- Berry phase --------------------

            ! this part is executed only in a berry phase calculation
            !         if( isberry == 1 ) then
            !           close(120)
            !         endif

            !-------------------------------------------------------
            close (12)
            if (npsa .ne. 0) close (8)
            if (itestr .eq. -5 .and. rank .eq. 0) then
                close (16)
                if (ieskin .ne. 0) close (17)
                if (idyn .gt. 0) close (22)
                if (iespbc) close (23)
            end if
        end if
        if (itestr .eq. -5) then
            if (writescratch .eq. 0) close (unit=indopen3, status='DELETE')
            if (abs(kl) .eq. 9) close (unit=indopen, status='DELETE')
        end if
        !       if(itestr.eq.-5.and.ndimp.ge.1.and.rank.eq.0) close(18)
#ifdef PARALLEL
        if (iread .eq. 2) close (indopen2)
        if (iread .ge. 6 .and. rank .eq. 0) close (15)
#else
        close (12)
        if (iread .eq. 2 .or. iread .ge. 6) close (15)
#endif

        if (itestr .eq. -5 .and. rank .eq. 0) close (13)

        if (unreliable .eq. 0) then
            call deallocate_all ! sometimes does not work on
            ! fucking infiniband network

            write (6, *) ' After deallocate_all '
            stop
            if (allocated(kelsav)) deallocate (kelsav)
        end if
#ifdef PARALLEL
        call mpi_finalize(ierr)
#endif

    end subroutine Finalizeall

    subroutine read_alphavar
        implicit none
        !       put the initial parameter wavefunction in alphavar

        if (np .gt. 0) then

            if (iopt .eq. 1 .or. .not. yesquantum) then

                dek(1:ieskinr_pos) = 0.d0

                indc = 0
                if (iesinv .ne. 0) then
                    do k = indc + 1, indc + iesinv
                        alphavar(k) = ddwsz(k - indc)
                    end do
                end if
                indc = indc + iesinv
                if (iesm .ne. 0) then
                    do k = indc + 1, indc + iesm
                        alphavar(k) = vju(k - indc)
                    end do
                end if
                indc = indc + iesm

                if (iesd .ne. 0) then
                    do k = indc + 1, indc + iesd
                        alphavar(k) = vj(k - indc)
                    end do
                end if
                indc = indc + iesd

                if (iesfree .ne. 0) then
                    do k = indc + 1, indc + iesfree
                        alphavar(k) = ddw(k - indc)
                    end do
                end if
                indc = indc + iesfree

                if (iessw .ne. 0) then
                    do k = indc + 1, indc + iessw
                        alphavar(k) = dsw(k - indc)
                    end do
                end if
                indc = indc + iessw

                if (iesup .ne. 0) then
                    do k = indc + 1, indc + iesup
                        alphavar(k) = dup(k - indc)
                    end do
                end if

                indc = indc + iesup

                if (iesking .ne. 0) then
                    do k = indc + 1, indc + iesking
                        alphavar(k) = dekg(k - indc)
                    end do
                end if
                indc = indc + iesking

                if (ieskin .ne. 0) then
                    ! The ion mapping coordinates
                    do k = indc + 1, indc + ieskin
                        alphavar(k) = dek(k - indc)
                    end do
                end if

                if (iespbc .and. itestr .eq. -5) then
                    !       Restoring consistency with fort.10
                    if (ieskint .eq. ieskinr_pos + 2) then
                        dek(ieskinr_pos + 1 - iesking) = cellscale(2)
                        dek(ieskinr_pos + 2 - iesking) = cellscale(3)
                    elseif (ieskint .eq. ieskinr_pos + 1) then
                        dek(ieskinr_pos + 1 - iesking) = cellscale(1)
                    elseif (ieskint .eq. ieskinr_pos + 3) then
                        dek(ieskinr_pos + 1 - iesking) = cellscale(1)
                        dek(ieskinr_pos + 2 - iesking) = cellscale(2)
                        dek(ieskinr_pos + 3 - iesking) = cellscale(3)
                    end if
                end if

            else
                !    Read parameters from alphavar
                indc = 0
                if (iesinv .ne. 0) then
                    do k = indc + 1, indc + iesinv
                        ddwsz(k - indc) = alphavar(k)
                    end do
                end if
                indc = indc + iesinv
                if (iesm .ne. 0) then
                    do k = indc + 1, indc + iesm
                        vju(k - indc) = alphavar(k)
                    end do
                end if
                indc = indc + iesm

                if (iesd .ne. 0) then
                    do k = indc + 1, indc + iesd
                        vj(k - indc) = alphavar(k)
                    end do
                end if
                indc = indc + iesd

                if (iesfree .ne. 0) then
                    do k = indc + 1, indc + iesfree
                        ddw(k - indc) = alphavar(k)
                    end do
                end if
                indc = indc + iesfree

                if (iessw .ne. 0) then
                    do k = indc + 1, indc + iessw
                        dsw(k - indc) = alphavar(k)
                    end do
                end if
                indc = indc + iessw

                if (iesup .ne. 0) then
                    do k = indc + 1, indc + iesup
                        dup(k - indc) = alphavar(k)
                    end do
                end if

                indc = indc + iesup

                if (iesking .ne. 0) then
                    do k = indc + 1, indc + iesking
                        dekg(k - indc) = 0.d0
                        alphavar(k) = 0.d0
                    end do
                end if
                indc = indc + iesking

                if (iespbc .and. itestr .eq. -5) then
                    !       Restoring consistency with fort.10
                    if (ieskint .eq. ieskinr_pos + 2) then
                        cellscale(2) = dek(ieskinr_pos + 1 - iesking)
                        cellscale(3) = dek(ieskinr_pos + 2 - iesking)
                        cellscale(1) = omega/cellscale(2)/cellscale(3)
                    elseif (ieskint .eq. ieskinr_pos + 1) then
                        cellscale(1) = dek(ieskinr_pos + 1 - iesking)
                        cellscale(2:3) = 1.d0
                    elseif (ieskint .eq. ieskinr_pos + 3) then
                        cellscale(1) = dek(ieskinr_pos + 1 - iesking)
                        cellscale(2) = dek(ieskinr_pos + 2 - iesking)
                        cellscale(3) = dek(ieskinr_pos + 3 - iesking)
                    end if
                end if

                if (ieskin .ne. 0) then
                    ! The ion mapping coordinates
                    do k = indc + 1, indc + ieskin
                        if (k - indc .le. ieskinr_pos - iesking) then
                            dek(k - indc) = 0.d0
                            alphavar(k) = 0.d0
                        else
                            dek(k - indc) = alphavar(k)
                        end if
                    end do
                end if

            end if ! not yesquantum or iopt=1

        end if ! np.gt.0

    end subroutine read_alphavar

end
