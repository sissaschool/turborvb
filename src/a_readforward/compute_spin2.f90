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

!
! Spin2: a module for spin^2 correlation function used by readforword
! written by Ye Luo
! last update: 20 MAR 2012
! developer notes in doc/S^2
!

!#define DEBUG

!------------------------------------------------------------------------------
module Spin2
    !------------------------------------------------------------------------------
    use allio, only: nelup, neldo, nel, yesfast, nelup_mat, nel_mat
    use allio, only: nelorb, nelorbh, nelorb_c, indt4, firstmol
    use allio, only: nelorbj, nelorbjh, iesdr, indt4j, iessz, vj
    use allio, only: LBox
    use allio, only: rank
    use constants, only: ipj, ipf, ipc
    use cell, only: CartesianToCrystal, cellscale, car2cry, map
    implicit none
    logical ifspin2 !(global) label for spin^2 correlation
    integer nspin2 !(global) the number of output values, =1
    real*8 spin2_local !(global) the only output value
    logical ifjasz !(local)  ture if jastrows are spin sensitive
    real(8), dimension(:, :), allocatable :: Lspin, ratiospin
contains
    !------------------------------------------------------------------------------
    subroutine prepare_jastrow_spin(kel, winvj, winvjbar, winvjbarn, winvjbarsz&
            &, JasSpin, aux, Jasupup, Jasupdo, Jasdoup, Jasdodo, dvet)
        !------------------------------------------------------------------------------
        !This subroutine prepare the JasSpin matrix for spin-flip
        !JasSpin, a nel*nel matrix, contains the spin sensitive part of all jastrows
        implicit none
        real*8 winvj(max(nelorbj, 1), 0:indt4j, nel), winvjbarsz(*), winvjbar(*), winvjbarn(*)
        real*8 kel(3, *)
        real*8 JasSpin(nel, nel), aux(nel), dvet(max((ipj - 1)*nel, 1))
        real*8 Jasupup(max((ipj - 1)*nelup, 1), max((ipj - 1)*neldo, 1))
        real*8 Jasupdo(max((ipj - 1)*nelup, 1), max((ipj - 1)*neldo, 1))
        real*8 Jasdoup(max((ipj - 1)*nelup, 1), max((ipj - 1)*neldo, 1))
        real*8 Jasdodo(max((ipj - 1)*nelup, 1), max((ipj - 1)*neldo, 1))

        integer nelorbj5
        integer iesd
        integer j, k, l, kk
        real*8 rc(3)
        logical ifjas2sz
        !external functions
        integer, external :: iesdr2iesd
        real*8, external :: jastrow_ee

        nelorbj5 = nelorbj*(indt4j + 1)
        iesd = iesdr2iesd(iesdr)
        !       check all spin contaminated according to jastrow_ee
        if (iesd .ne. 1 .and. iesd .ne. -4 .and. iesd .ne. 6 .and. iesd .ne. 8 .and. iesd .ne. 9&
                &.and. iesd .ne. 5 .and. iesd .ne. -2 .and. iesd .ne. 3 .and. iesd .ne. 0) then
            ifjas2sz = .true.
        else
            ifjas2sz = .false.
        end if
        if (ifjas2sz .or. iessz .or. ipj .eq. 2) then
            ifjasz = .true.
        else
            ifjasz = .false.
        end if
        if (ipj .eq. 2) then
            !       Compute stuff for k< nelup
            call dgemm('T', 'N', nelup, nel, nelorbj, 1.d0, winvj, nelorbj5, winvjbar(nelorbj + 1), 2*nelorbjh&
                    &, 0.d0, JasSpin, nel)
            call dgemm('T', 'N', nelup, nel, nelorbj, -1.d0, winvj, nelorbj5, winvjbar, 2*nelorbjh&
                    &, 1.d0, JasSpin, nel)
            !       Compute stuff for k> nelup
            call dgemm('T', 'N', neldo, nel, nelorbj, 1.d0, winvj(1, 0, nelup + 1), nelorbj5, winvjbar(nelorbj + 1)&
                    &, 2*nelorbjh, 0.d0, JasSpin(nelup + 1, 1), nel)
            call dgemm('T', 'N', neldo, nel, nelorbj, -1.d0, winvj(1, 0, nelup + 1), nelorbj5, winvjbar&
                    &, 2*nelorbjh, 1.d0, JasSpin(nelup + 1, 1), nel)
            dvet = 0.d0
            do k = 1, nel
                do j = 1, nel
                    if (j .ne. k) dvet(k) = dvet(k) + JasSpin(k, j)
                end do
            end do
            !       Jasdodo
            call dgemm('T', 'N', nelup, neldo, nelorbj, 1.d0, winvj, nelorbj5&
                    &, winvjbar(2*nelorbj*nelup + nelorbj + 1), 2*nelorbjh, 0.d0, Jasdodo, nelup)
            !       Jasdoup
            call dgemm('T', 'N', nelup, neldo, nelorbj, 1.d0, winvj, nelorbj5&
                    &, winvjbarn(2*nelorbj*nelup + nelorbj + 1), 2*nelorbjh, 0.d0, Jasdoup, nelup)
            !       Jasupup
            call dgemm('T', 'N', nelup, neldo, nelorbj, 1.d0, winvj, nelorbj5, winvjbarn(2*nelorbj*nelup + 1)&
                    &, 2*nelorbjh, 0.d0, Jasupup, nelup)
            !       Jasupdo
            call dgemm('T', 'N', nelup, neldo, nelorbj, 1.d0, winvj, nelorbj5, winvjbar(2*nelorbj*nelup + 1)&
                    &, 2*nelorbjh, 0.d0, Jasupdo, nelup)
            Jasspin = 0.d0 ! Reinitializing Jasspin
        end if
        if (ifjas2sz) then
            !compute e-e two-body jastrow matrix nel*nel for parallel and anti-parallel
            !use jastrow_ee taken iesd and sum up
            do l = 1, nel
                JasSpin(l, l) = 0.d0
                do k = l + 1, nel
                    !compute rc and jastrow_ee
                    rc(:) = kel(:, l) - kel(:, k)
                    !************ PERIODIC WORLD ***********
                    if (LBox .gt. 0.d0) then
                        call CartesianToCrystal(rc, 1)
                        do kk = 1, 3
                            rc(kk) = map(rc(kk), cellscale(kk))
                        end do
                    end if
                    !************ PERIODIC WORLD ***********
                    !for parallel
                    JasSpin(k, l) = jastrow_ee(rc, vj, iesd, 1)
                    !for anti-para
                    JasSpin(k, l) = JasSpin(k, l) - jastrow_ee(rc, vj, iesd, -1)
                    !for symetric part
                    JasSpin(l, k) = JasSpin(k, l)
                end do !k
            end do !l
            !3 body jastrow spin part needs winvj winvjbarsz
            if (iessz) then
                call dgemm('T', 'N', nel, nel, nelorbjh, 1.d0, winvjbarsz, nelorbjh, winvj(1, 0, 1), nelorbj5, 0.5d0, JasSpin, nel)
            else
                JasSpin = 0.5*JasSpin
            end if
        else
            !3 body jastrow spin part needs winvj winvjbarsz
            if (iessz) then
                call dgemm('T', 'N', nel, nel, nelorbjh, 1.d0, winvjbarsz, nelorbjh, winvj(1, 0, 1), nelorbj5, 0.d0, JasSpin, nel)
            end if
        end if !ifjas2sz
        if (ifjasz) then
            aux = 0.d0
            do l = 1, nel
                !               remove self correlation. optimal
                JasSpin(l, l) = 0.d0
                do k = 1, nelup
                    aux(l) = aux(l) + JasSpin(k, l)
                end do !k
                do k = nelup + 1, nel
                    aux(l) = aux(l) - JasSpin(k, l)
                end do !k
                !               remove self correlation
                !                if(l>nelup)then
                !                    aux(l)=aux(l)+JasSpin(l,l)
                !                else
                !                    aux(l)=aux(l)-JasSpin(l,l)
                !                endif
            end do !l
        end if !ifjasz

        !------------------------------------------------------------------------------
    end subroutine prepare_jastrow_spin
    !------------------------------------------------------------------------------

    !It calculates the ratio for every possible flip of spin up and spin down
    !couples of a given configuration (the fastest index is for the up electron)
    !The matrices use the name given by Tomonori in its notes with the suffix
    !spin.
    !------------------------------------------------------------------------------
    subroutine preparepfaff(winv, winvbar, ainv, ratiopfaff)
        !------------------------------------------------------------------------------
        use constants, only: ipc, ipf, zone, zzero
        !Lspin to be defined
        implicit none
        real*8 Ainv(ipc*nelup_mat, nelup_mat)
        real*8 winv(ipc*nelorbh, 0:indt4, nel), winvbar(ipf*ipc*nelorbh, nel_mat)
        real*8 ratiopfaff(ipc*nelup, neldo)
        real*8, dimension(:, :), allocatable :: Vspin, Wspin, Dspin, aux, winvbaraux
        integer j, k, nelorb5
        allocate (Vspin(ipc*nel_mat, nel), Wspin(ipc*nel_mat, nel), Dspin(ipc*nelup, neldo), aux(nelup*ipc, nelorbh), &
                  winvbaraux(ipc*nelorbh, nel_mat))

        nelorb5 = nelorbh*(indt4 + 1)
        winvbaraux(1:ipc*nelorbh, 1:nel_mat) = winvbar(1:ipc*nelorbh, 1:nel_mat) &
                                               - winvbar(1 + ipc*nelorbh:2*ipc*nelorbh, 1:nel_mat)

        Vspin = 0.d0
        if (ipc .eq. 1) then

            !call dgemm('T','N', nel,neldo, nelorbh, 1.d0, winvbar(ipc*nelorbh+1,1), nelorbh&
            !           , winv(1,0,nelup+1), nelorb5, 0.d0, Vspin(1n,nel)
            !definition of Vspin
            call dgemm('T', 'N', nel_mat, nelup, nelorbh, 1.d0, winvbaraux, nelorbh, winv, nelorb5, 0.d0, Vspin, nel_mat)
            call dgemm('T', 'N', nel_mat, neldo, nelorbh, -1.d0, winvbaraux, nelorbh, winv(1, 0, nelup + 1), &
                       nelorb5, 1.d0, Vspin(1, nelup + 1), nel_mat)
            !definition of Wspin
            call dgemm('N', 'N', nel_mat, nel, nel_mat, 1.d0, ainv, nel_mat, Vspin, nel_mat, 0.d0, Wspin, nel_mat)
            !definition of Dspin
            call dgemm('T', 'N', nelup, neldo, nel_mat, 1.d0, Vspin, nel_mat, Wspin(1, 1 + nelup), nel_mat, 0.d0, Dspin, nelup)
            !definition of theta spin (for what we need we can define Dspin=Dspin+thetaspin (of the notes)
            call dgemm('T', 'N', nelup, nelorbh, nelorbh, 1.d0, winv, nelorb5, Lspin, nelorbh, 0.d0, aux, nelup)
            call dgemm('N', 'N', nelup, neldo, nelorbh, 1.d0, aux, nelup, winv(1, 0, 1 + nelup), nelorb5, 1.d0, Dspin, nelup)
            do k = 1, neldo
                do j = 1, nelup
                    ratiopfaff(j, k) = -ainv(j, k + nelup)*Dspin(j, k) + 1.d0 + Wspin(j, j) + Wspin(k + nelup, k + nelup) + &
                                       Wspin(j, j)*Wspin(k + nelup, k + nelup) - Wspin(k + nelup, j)*Wspin(j, k + nelup)
                end do
            end do
            !        write (6,*)  "prepare pfaff", ratiopfaff
        else
            !definition of Vspin
            call zgemm('T', 'N', nel_mat, nelup, nelorbh, zone, winvbaraux, nelorbh, winv, nelorb5, zzero, Vspin, nel_mat)
            call zgemm('T', 'N', nel_mat, neldo, nelorbh, -zone, winvbaraux, nelorbh, winv(1, 0, nelup + 1), &
                       nelorb5, zone, Vspin(1, nelup + 1), nel_mat)
            !         call zgemm('T','N', nel,nel, nelorbh, zone,  winv, nelorb5, winvbaraux, nelorbh,zzero, Vspin,nel)
            !definition of Wspin
            call zgemm('N', 'N', nel_mat, nel, nel_mat, zone, ainv, nel_mat, Vspin, nel_mat, zzero, Wspin, nel_mat)
            !definition of Dspin
            call zgemm('T', 'N', nelup, neldo, nel_mat, zone, Vspin, nel_mat, Wspin(1, 1 + nelup), nel_mat, zzero, Dspin, nelup)
            !definition of theta spin (for what we need we can define Dspin=Dspin+thetaspin (of the notes)
            call zgemm('T', 'N', nelup, nelorbh, nelorbh, zone, winv, nelorb5, Lspin, nelorbh, zzero, aux, nelup)
            call zgemm('N', 'N', nelup, neldo, nelorbh, zone, aux, nelup, winv(1, 0, 1 + nelup), nelorb5, zone, Dspin, nelup)
            do k = 1, neldo
                do j = 1, nelup
                    ratiopfaff(2*j - 1, k) = -ainv(2*j - 1, k + nelup)*Dspin(2*j - 1, k) + ainv(2*j, k + nelup)*Dspin(2*j, k) + &
                                             1.d0 + Wspin(2*j - 1, j) + Wspin(2*(k + nelup) - 1, k + nelup) &
                                             - Wspin(2*(k + nelup) - 1, j)*Wspin(2*j - 1, k + nelup) &
                                             + Wspin(2*(k + nelup), j)*Wspin(2*j, k + nelup) &
                                             - Wspin(2*j, j)*Wspin(2*(k + nelup), k + nelup) &
                                             + Wspin(2*j - 1, j)*Wspin(2*(k + nelup) - 1, k + nelup)
                    ratiopfaff(2*j, k) = -ainv(2*j, k + nelup)*Dspin(2*j - 1, k) - ainv(2*j - 1, k + nelup)*Dspin(2*j, k) + &
                                         Wspin(2*j, j) + Wspin(2*(k + nelup), k + nelup) &
                                         - Wspin(2*(k + nelup), j)*Wspin(2*j - 1, k + nelup) &
                                         - Wspin(2*(k + nelup) - 1, j)*Wspin(2*j, k + nelup) &
                                         + Wspin(2*j - 1, j)*Wspin(2*(k + nelup), k + nelup) &
                                         + Wspin(2*j, j)*Wspin(2*(k + nelup) - 1, k + nelup)
                end do
            end do
        end if
        deallocate (Vspin, Wspin, Dspin, aux, winvbaraux)

        !------------------------------------------------------------------------------
    end subroutine preparepfaff
    !------------------------------------------------------------------------------

    !Initialization of the matrix Lspin for the uncontracted case
    !------------------------------------------------------------------------------
    subroutine inits2pfaff(detmat_c)
        !------------------------------------------------------------------------------
        use constants, only: ipc, ipf, zone, zzero
        implicit none
        integer :: k, j
        real(8) :: detmat_c(ipc*2*nelorbh, 2*nelorbh)
        !   compute once for all the lambda matrix

        allocate (Lspin(ipc*nelorbh, nelorbh))
        do k = 1, nelorbh
            do j = 1, nelorbh*ipc
                Lspin(j, k) = detmat_c(j, k + nelorbh) + detmat_c(j + ipc*nelorbh, k) &
                              - detmat_c(j + ipc*nelorbh, k + nelorbh) - detmat_c(j, k)
            end do
        end do

        !      write (6,*) "Initializing Spin2 pfaffian uncontracted case", rank
        !------------------------------------------------------------------------------
    end subroutine inits2pfaff
    !------------------------------------------------------------------------------

    !Initialization of the matrix Lspin for the contracted case
    !------------------------------------------------------------------------------
    subroutine inits2pfaff_c(detmat_c, mu_c)
        !------------------------------------------------------------------------------
        use constants, only: ipc, ipf, zone, zzero
        !Lspin to be defined
        implicit none
        real(8) :: detmat_c(ipc*nelorb_c, nelorb_c), mu_c(2*ipc*nelorbh, nelorb_c)
        real(8), dimension(:, :), allocatable :: auxgen1, auxgen2

        allocate (Lspin(ipc*nelorbh, nelorbh))
        allocate (auxgen1(ipc*2*nelorbh, nelorb_c), auxgen2(ipc*2*nelorbh, 2*nelorbh))
        if (ipc .eq. 1) then
            call dgemm('N', 'N', nelorbh*2, nelorb_c, nelorb_c, 1.d0, mu_c, nelorbh*2&
                    &, detmat_c, nelorb_c, 0.d0, auxgen1, nelorbh*2)
            call dgemm('N', 'T', nelorbh*2, nelorbh*2, nelorb_c, 1.d0, auxgen1, nelorbh*2&
                    &, mu_c, nelorbh*2, 0.d0, auxgen2, nelorbh*2)
        else
            call zgemm('N', 'N', nelorbh*2, nelorb_c, nelorb_c, zone, mu_c, nelorbh*2&
                    &, detmat_c, nelorb_c, zzero, auxgen1, nelorbh*2)
            call zgemm('N', 'T', nelorbh*2, nelorbh*2, nelorb_c, zone, auxgen1, nelorbh*2&
                    &, mu_c, nelorbh*2, zzero, auxgen2, nelorbh*2)
        end if
        Lspin(1:ipc*nelorbh, 1:nelorbh) = auxgen2(1 + ipc*nelorbh:ipc*nelorbh*2, 1:nelorbh) + &
                                          auxgen2(1:ipc*nelorbh, 1 + nelorbh:nelorbh*2) - auxgen2(1:ipc*nelorbh, 1:nelorbh) &
                                          - auxgen2(1 + ipc*nelorbh:2*ipc*nelorbh, 1 + nelorbh:nelorbh*2)
        deallocate (auxgen1, auxgen2)
        if (rank .eq. 0) write (6, *) "Initializing Spin2 pfaffian uncontracted case"
        !------------------------------------------------------------------------------
    end subroutine inits2pfaff_c
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    subroutine compute_spin2(kel, Ainv, winv, winvj, winvbar, winvjbar, winvjbarsz, detmat, projm, mu_c, nmol)
        !------------------------------------------------------------------------------
        !This subroutine explores the spin-swaps of all anti-parallel e-e pairs
        !and compute S^2
        use constants, only: ipc, ipj, ipf, zone, zzero
        use allio, only: jasmat, jasmat_c, muj_c, contractionj, nelorbj_c
        implicit none
        real*8 kel(3, *)
        real*8 Ainv(ipc*nelup_mat, nelup_mat)
        !CCC Perche' winvbar e' dichiarata ipf*ipc*nelorb,2*nel??
        real*8 winv(ipc*nelorb, 0:indt4, nel), winvbar(ipf*ipc*nelorb, nel_mat)
        real*8 winvj(max(nelorbj, 1), 0:indt4j, nel), winvjbarsz(*), winvjbar(*)
        real*8 detmat(ipc*nelorbh, nelorbh)
        real*8 projm(*), mu_c(*)
        integer nmol
        integer ierr
        real(8), dimension(:, :), allocatable :: winvbar_ext, det_update_vec1, det_update_vec2
        real(8), dimension(:, :), allocatable :: Amatrix, UpUp, DownDown, DownUp, Buu, Bud
        real(8), dimension(:, :), allocatable :: JasSpin, tempMat
        real(8), dimension(:, :), allocatable :: Jasupup, Jasupdo, Jasdoup, Jasdodo
        real(8), dimension(:, :), allocatable :: ratiopfaff
        real(8), dimension(:), allocatable :: aux, dvet, winvjbarn, psip

        real*8 det_update_matrix(ipc*2, 2)
        real*8 Sz, spin2_jastrow, spin2_det
        integer i, j, k, l !index only used for loops
        integer nelorb5, nelorbj5, firstmmu, nelorbjcp

        if (ipf .eq. 1) then
            allocate (winvbar_ext(ipc*nelorb, nelup), det_update_vec1(ipc*nelup, 2), det_update_vec2(ipc*nelup, 2))
            allocate (Amatrix(ipc*nelup, nelup), UpUp(ipc*nelup, nelup), DownDown(ipc*neldo, nelup)&
                    &, DownUp(ipc*neldo, nelup), Buu(ipc*nelup, nelup), Bud(ipc*nelup, neldo))
        else
            allocate (ratiopfaff(ipc*nelup, neldo))
        end if
        allocate (JasSpin(nel, nel), aux(nel), tempMat(nmol, nelup))

        if (ipj .eq. 2) then
            allocate (dvet(nel))
            allocate (Jasupup(nelup, neldo), Jasupdo(nelup, neldo)&
                    &, Jasdoup(nelup, neldo), Jasdodo(nelup, neldo))
            allocate (winvjbarn(2*nelorbj*nel))
            allocate (psip(2*nelup*nelorbj_c))
        else
            allocate (dvet(1), Jasupup(1, 1), Jasupdo(1, 1), Jasdoup(1, 1), Jasdodo(1, 1))
            allocate (winvjbarn(1))
            allocate (psip(1))
        end if
        psip = 0.d0
        dvet = 0.d0
        Jasupup = 0.d0
        Jasupdo = 0.d0
        Jasdoup = 0.d0
        Jasdodo = 0.d0

        !initialize variables
        spin2_local = 0.d0
        nelorb5 = nelorb*(indt4 + 1)
        nelorbj5 = nelorbj*(indt4j + 1)
        Sz = (nelup - neldo)*0.5d0

        !prepare Amatrix(actually UpDown), UpUp, DownDown, DownUp  in advance
        if (ipf .eq. 1) then
            if (ipc .eq. 1) then
                call dgemm('T', 'N', nelup, nelup, nelorbh, 1.d0, winv, nelorb5, winvbar(1, nelup + 1) &
                           , nelorbh, 0.d0, Amatrix, nelup)
                call dgemm('T', 'N', neldo, nelup, nelorbh, 1.d0, winv(1, 0, nelup + 1), nelorb5, winvbar(1, nelup + 1) &
                           , nelorbh, 0.d0, DownDown, max(neldo, 1))
            else
                call zgemm('T', 'N', nelup, nelup, nelorbh, zone, winv, nelorb5, winvbar(1, nelup + 1) &
                           , nelorbh, zzero, Amatrix, nelup)
                call zgemm('T', 'N', neldo, nelup, nelorbh, zone, winv(1, 0, nelup + 1), nelorb5, winvbar(1, nelup + 1) &
                           , nelorbh, zzero, DownDown, max(neldo, 1))
            end if
            if (yesfast == 0) then
                if (ipc .eq. 1) then
                    call dgemm('N', 'N', nelorbh, nelup, nelorbh, 1.d0, detmat, nelorbh, winv &
                               , nelorb5, 0.d0, winvbar_ext, nelorbh)
                else
                    call zgemm('N', 'N', nelorbh, nelup, nelorbh, zone, detmat, nelorbh, winv &
                               , nelorb5, zzero, winvbar_ext, nelorbh)
                end if
            else
                firstmmu = ipc*(firstmol - 1)*nelorbh + 1
                if (ipc .eq. 1) then
                    call dgemm('T', 'N', nmol, nelup, nelorbh, 1.d0, mu_c(firstmmu), nelorbh &
                               , winv, nelorb5, 0.d0, tempMat, nmol)
                    call dgemm('N', 'N', nelorbh, nelup, nmol, 1.d0, projm(firstmmu), nelorbh &
                               , tempMat, nmol, 0.d0, winvbar_ext, nelorbh)
                else
                    call zgemm('T', 'N', nmol, nelup, nelorbh, zone, mu_c(firstmmu), nelorbh &
                               , winv, nelorb5, zzero, tempMat, nmol)
                    call zgemm('N', 'N', nelorbh, nelup, nmol, zone, projm(firstmmu), nelorbh &
                               , tempMat, nmol, zzero, winvbar_ext, nelorbh)
                end if
            end if
            if (ipc .eq. 1) then
                call dgemm('T', 'N', nelup, nelup, nelorbh, 1.d0, winv, nelorb5 &
                           , winvbar_ext, nelorbh, 0.d0, UpUp, nelup)
                call dgemm('T', 'N', neldo, nelup, nelorbh, 1.d0, winv(1, 0, nelup + 1), nelorb5 &
                           , winvbar_ext, nelorbh, 0.d0, DownUp, max(neldo, 1))
                !prepare Buu and Bud matrix to reach N^3 scaling
                call dgemm('N', 'N', nelup, nelup, nelup, 1.d0, Ainv, nelup, UpUp, nelup, 0.d0, Buu, nelup)
                call dgemm('N', 'N', nelup, neldo, nelup, 1.d0, Ainv, nelup, Amatrix, nelup, 0.d0, Bud, nelup)
            else
                call zgemm('T', 'N', nelup, nelup, nelorbh, zone, winv, nelorb5 &
                           , winvbar_ext, nelorbh, zzero, UpUp, nelup)
                call zgemm('T', 'N', neldo, nelup, nelorbh, zone, winv(1, 0, nelup + 1), nelorb5 &
                           , winvbar_ext, nelorbh, zzero, DownUp, max(neldo, 1))
                !prepare Buu and Bud matrix to reach N^3 scaling
                call zgemm('N', 'N', nelup, nelup, nelup, zone, Ainv, nelup, UpUp, nelup, zzero, Buu, nelup)
                call zgemm('N', 'N', nelup, neldo, nelup, zone, Ainv, nelup, Amatrix, nelup, zzero, Bud, nelup)
            end if
        else
            call preparepfaff(winv, winvbar, ainv, ratiopfaff)
        end if

#ifdef DEBUG
        write (6, *) "winvbar"
        write (6, *) winvbar(:, 1:nelup)
        write (6, *) "winvbar_ext"
        write (6, *) winvbar_ext(:, 1:nelup)
        write (6, *) "Amatrix"
        write (6, *) Amatrix
        write (6, *) "Ainv"
        write (6, *) Ainv
        write (6, *) "UpUp"
        write (6, *) UpUp
        write (6, *) "DownDown"
        write (6, *) DownDown
        write (6, *) "DownUp"
        write (6, *) DownUp
#endif

        if (ipj .eq. 2) then
            !   computed winvjbarn (spin flipped)
            if (contractionj .eq. 0) then
                call dgemm('N', 'N', 2*nelorbj, nelup, nelorbj, 1.d0, jasmat(2*nelorbj*nelorbj + 1)&
                        &, 2*nelorbj, winvj(1, 0, 1), nelorbj5, 0.d0, winvjbarn, 2*nelorbj)
                call dgemm('N', 'N', 2*nelorbj, neldo, nelorbj, 1.d0, jasmat, 2*nelorbj&
                        &, winvj(1, 0, nelup + 1), nelorbj5, 0.d0, winvjbarn(1 + 2*nelorbj*nelup), 2*nelorbj)

            else

                nelorbjcp = nelup*nelorbj_c + 1

                !  up-up
                call dgemm('T', 'N', nelorbj_c, nelup, nelorbjh, 1.d0, muj_c, nelorbj&
                        &, winvj, nelorbj5, 0.d0, psip(nelorbjcp), nelorbj_c)
                call dgemm('N', 'N', nelorbj_c, nelup, nelorbj_c, 1.d0&
                        &, jasmat_c(2*nelorbj_c*nelorbj_c + 1), 2*nelorbj_c&
                        &, psip(nelorbjcp), nelorbj_c, 0.d0, psip, nelorbj_c)
                call dgemm('N', 'N', nelorbjh, nelup, nelorbj_c, 1.d0, muj_c, nelorbjh&
                        &, psip, nelorbj_c, 0.d0, winvjbarn, 2*nelorbjh)
                !  down-up
                call dgemm('N', 'N', nelorbj_c, nelup, nelorbj_c, 1.d0&
                        &, jasmat_c(2*nelorbj_c*nelorbj_c + 1 + nelorbj_c), 2*nelorbj_c&
                        &, psip(nelorbjcp), nelorbj_c, 0.d0, psip, nelorbj_c)
                call dgemm('N', 'N', nelorbjh, nelup, nelorbj_c, 1.d0, muj_c, nelorbjh&
                        &, psip, nelorbj_c, 0.d0, winvjbarn(1 + nelorbjh), 2*nelorbjh)

                !  up-down
                call dgemm('T', 'N', nelorbj_c, neldo, nelorbj, 1.d0, muj_c, nelorbj&
                        &, winvj(1, 0, 1 + nelup), nelorbj5, 0.d0, psip(nelorbjcp), nelorbj_c)
                call dgemm('N', 'N', nelorbj_c, neldo, nelorbj_c, 1.d0, jasmat_c, 2*nelorbj_c&
                        &, psip(nelorbjcp), nelorbj_c, 0.d0, psip, nelorbj_c)
                call dgemm('N', 'N', nelorbj, neldo, nelorbj_c, 1.d0, muj_c, nelorbj&
                        &, psip, nelorbj_c, 0.d0, winvjbarn(1 + 2*nelorbj*nelup), 2*nelorbj)

                !  down-down
                call dgemm('N', 'N', nelorbj_c, neldo, nelorbj_c, 1.d0&
                        &, jasmat_c(1 + nelorbj_c), 2*nelorbj_c, psip(nelorbjcp), nelorbj_c, 0.d0&
                        &, psip, nelorbj_c)
                call dgemm('N', 'N', nelorbj, neldo, nelorbj_c, 1.d0, muj_c, nelorbj&
                        &, psip, nelorbj_c, 0.d0, winvjbarn(1 + 2*nelup*nelorbj + nelorbj), 2*nelorbj)

            end if

        end if

        !prepare a matrix and a auxiliary vec for the spin flip in jastrow
        call prepare_jastrow_spin(kel, winvj, winvjbar, winvjbarn, winvjbarsz&
                &, JasSpin, aux, Jasupup, Jasupdo, Jasdoup, Jasdodo, dvet)

        do l = 1, neldo
            if (ipf .eq. 1) then
                det_update_vec1 = 0.d0
                det_update_vec1(ipc*(l - 1) + 1, 2) = 1.d0
            end if
            do k = 1, nelup
                !           determinant part starts
                if (ipf .eq. 2) then
                    !   CCC Here put you pfaffian ratio r_kl.
                    !   SSS Immagino che spin2_det sia la somma su k ed l di r_kl
                    if (ipc .eq. 1) then
                        spin2_det = -ratiopfaff(k, l)
                    else
                        spin2_det = -ratiopfaff(2*k - 1, l)
                    end if
                else
                    if (ipc .eq. 2) then
                        call det_update(l, k, nelup, neldo, det_update_vec1, det_update_vec2, DownDown, Amatrix, Ainv, Buu&
                                &, Bud, Downup, UpUp, det_update_matrix, spin2_det)
                    else
                        call dcopy(nelup, DownDown(l, 1), neldo, det_update_vec1, 1)
                        call daxpy(nelup, -1.d0, Amatrix(k, 1), nelup, det_update_vec1, 1)
                        call dcopy(nelup, Ainv(1, k), 1, det_update_vec2, 1)
                        do i = 1, nelup
                            det_update_vec2(i, 2) = Buu(i, k) - Bud(i, l) &
                                                    + Ainv(i, k)*(DownUp(l, k) + Amatrix(k, l) - DownDown(l, l) - UpUp(k, k))
                        end do
                        !construct the 2*2 matrix and compute its deteminant which is the difference after the interchange
                        det_update_matrix = 0.d0
                        det_update_matrix(1, 1) = 1.d0
                        det_update_matrix(2, 2) = 1.d0
#ifdef DEBUG
                        write (6, *) "det_update_vec1"
                        write (6, *) det_update_vec1
                        write (6, *) "det_update_vec2"
                        write (6, *) det_update_vec2
#endif
                        call dgemm('T', 'N', 2, 2, nelup, 1.d0, det_update_vec1, nelup &
                                   , det_update_vec2, nelup, 1.d0, det_update_matrix, 2)
                        spin2_det = det_update_matrix(1, 1)*det_update_matrix(2, 2) &
                                    - det_update_matrix(2, 1)*det_update_matrix(1, 2)
#ifdef DEBUG
                        write (6, *) "det", spin2_det
#endif
                        !           determinant part ends
                        !           jastrow part starts

                    end if

                end if

                if (ifjasz) then
                    spin2_jastrow = -2.d0*(aux(k) + JasSpin(nelup + l, k)) + 2.0*(aux(nelup + l) - JasSpin(k, nelup + l))
#ifdef DEBUG
                    write (6, *) "JasSpin"
                    do i = 1, nel
                        write (6, *) (JasSpin(i, j), j=1, nel)
                    end do
                    write (6, *) "jastrow", spin2_jastrow
                    call test_jastrow_spin_flip(JasSpin, k, l + nelup)
#endif
                    if (ipj .eq. 2) then
                        !  Shirakawa eq.
                        spin2_jastrow = spin2_jastrow + dvet(k) - dvet(l + nelup) &
                                        + Jasdoup(k, l) - Jasdodo(k, l) - Jasupup(k, l) + Jasupdo(k, l)
                    end if
                else
                    spin2_jastrow = 0.d0
                end if
                !           jastrow parts ends
                spin2_local = spin2_local + spin2_det*dexp(spin2_jastrow)
                if (.not. allocated(ratiospin)) write (6, *) "DEBUG HERE IT IS THE BUG"
                ratiospin(k, l) = -spin2_det*dexp(spin2_jastrow)
            end do !k
        end do !l

        spin2_local = Sz*Sz + nel/2.d0 - spin2_local

#ifdef DEBUG
        write (6, *) "spin2", spin2_local
#endif
        if (ipf .eq. 2) then
            deallocate (ratiopfaff)
        else
            deallocate (winvbar_ext, det_update_vec1, det_update_vec2)
            deallocate (Amatrix, UpUp, DownDown, DownUp, Buu, Bud)
        end if
        deallocate (JasSpin, aux, tempMat)
        deallocate (dvet, Jasupup, Jasdodo, Jasupdo, Jasdoup, winvjbarn, psip)

        !------------------------------------------------------------------------------
    end subroutine compute_spin2
    !------------------------------------------------------------------------------

#ifdef DEBUG
!------------------------------------------------------------------------------
    subroutine test_jastrow_spin_flip(JasSpin, k, l)
!------------------------------------------------------------------------------
        real*8 JasSpin(nel, nel)
        integer k, l !assume k<l

        integer i, j, spinprod
        real*8 ref, cost

        write (6, *) "switch", k, l
        ref = 0.d0
        call up3bodypsi_sz(nel, nelup, ref, 1, JasSpin)
        cost = 0.d0
        do i = 1, nel - 1
            do j = i + 1, nel
                spinprod = 1
                if (i <= nelup .and. j > nelup) spinprod = -1
                if (i == k) then
                    if (j /= l) spinprod = -spinprod
                elseif (i == l) then
                    spinprod = -spinprod
                elseif (j == k .or. j == l) then
                    spinprod = -spinprod
                end if
                write (6, *) i, j, " spinprod ", spinprod
                cost = cost + JasSpin(j, i)*spinprod
            end do
        end do
        write (6, *) "ref = ", ref, "cost = ", cost
        write (6, *) "jastrow_spin_flip refence difference = ", cost - ref
!------------------------------------------------------------------------------
    end subroutine test_jastrow_spin_flip
!------------------------------------------------------------------------------
#endif
end module Spin2

subroutine det_update(l, k, nelup, neldo, det_update_vec1, det_update_vec2, DownDown, Amatrix, Ainv, Buu&
        &, Bud, Downup, UpUp, det_update_matrix, spin2_det)
    use constants, only: zzero, zone, zmone
    implicit none
    integer nelup, neldo, l, k, i
    complex*16 det_update_vec1(nelup, 2), det_update_vec2(nelup, 2), Amatrix(nelup, nelup), Ainv(nelup, nelup)&
            &, UpUp(nelup, nelup), DownDown(neldo, nelup), DownUp(neldo, nelup), Buu(nelup, nelup), Bud(nelup, neldo)&
            &, det_update_matrix(2, 2)
    real*8 spin2_det

    !           determinant part starts
    call zcopy(nelup, DownDown(l, 1), neldo, det_update_vec1, 1)
    call zaxpy(nelup, zmone, Amatrix(k, 1), nelup, det_update_vec1, 1)
    call zcopy(nelup, Ainv(1, k), 1, det_update_vec2, 1)
    do i = 1, nelup
        det_update_vec2(i, 2) = Buu(i, k) - Bud(i, l) + Ainv(i, k)*(DownUp(l, k) + Amatrix(k, l) - DownDown(l, l) - UpUp(k, k))
    end do
    !construct the 2*2 matrix and compute its deteminant which is the difference after the interchange
    det_update_matrix = zzero
    det_update_matrix(1, 1) = zone
    det_update_matrix(2, 2) = zone
#ifdef DEBUG
    write (6, *) "det_update_vec1"
    write (6, *) det_update_vec1
    write (6, *) "det_update_vec2"
    write (6, *) det_update_vec2
#endif
    call zgemm('T', 'N', 2, 2, nelup, zone, det_update_vec1, nelup, det_update_vec2, nelup, zone, det_update_matrix, 2)
    spin2_det = det_update_matrix(1, 1)*det_update_matrix(2, 2) - det_update_matrix(2, 1)*det_update_matrix(1, 2)
#ifdef DEBUG
    write (6, *) "det", spin2_det
#endif
    !           determinant part ends
    !           jastrow part starts
end subroutine det_update
