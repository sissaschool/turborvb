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

function ioptorbcontr(ioptorb, LBox, j)
    implicit none
    integer ioptorb, ioptorbcontr, j, iflagerr
    real*8 LBox
    logical ieso

    ! contracted --> uncontracted orbitals
    ! Lbox=3 --> Crystal basis set, it needs open boundary conditions to work
    if (LBox .le. 0.d0 .or. abs(LBox) .eq. 2.d0 .or. Lbox .eq. 3.d0) then
        ieso = .true.
    else
        ieso = .false.
    end if

    select case (ioptorb)

    case (1:299) ! non contracted orbitals
        ioptorbcontr = ioptorb
    case (1000:1299)
        ioptorbcontr = ioptorb
    case (2000:2299)
        ioptorbcontr = ioptorb

    case (300) ! gaussian contraction
        if (ieso) then
            ioptorbcontr = 16
        else
            ioptorbcontr = 161
        end if

    case (301) ! Slater contraction
        ioptorbcontr = 34

    case (302) ! Slater contraction r^2
        ioptorbcontr = 10

    case (303) ! gaussian contraction r^2
        ioptorbcontr = 17

    case (307) ! mixed orbitals first Exp other gaussians
        if (j .eq. 1) then
            ioptorbcontr = 34
        else
            ioptorbcontr = 16
        end if

    case (308) ! mixed orbitals first Exp -gaussians-gaussian x r^3
        if (j .eq. 1) then
            ioptorbcontr = 34 ! the first is Exp x (1 + z r)   (no cusp)
        elseif (mod(j, 2) .eq. 0) then
            ioptorbcontr = 16 ! Gaussian
        else
            ioptorbcontr = 60 ! Gaussian x r^3
        end if

    case (309) ! mixed orbitals first Exp -gaussians-gaussian x r^2
        if (j .eq. 1) then
            ioptorbcontr = 34 ! the first is Exp x (1 + z r)   (no cusp)
        elseif (mod(j, 2) .eq. 0) then
            ioptorbcontr = 16 ! Gaussian
        else
            ioptorbcontr = 17 ! Gaussian x r^2
        end if

    case (310) ! mixed orbitals s  Gaussian + r^3 Gaussian
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 16
        else
            ioptorbcontr = 60
        end if

    case (311) ! mixed orbitals (1+z r) *exp(-z r )+r^2 x exp
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 34
        else
            ioptorbcontr = 10
        end if

    case (312) ! mixed orbitals Exp x (1+z r) +r^3 x Exp
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 34
        else
            ioptorbcontr = 12
        end if

    case (313) ! mixed orbitals r^3 x Exp+ Exp x (1+z r) (reversed 312).
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 12
        else
            ioptorbcontr = 34
        end if

    case (314) ! mixed orbitals s   Gaussian + r^2 Gaussian
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 16
        else
            ioptorbcontr = 17
        end if

    case (315) ! mixed orbitals s  r^2 Gaussian+Gaussian
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 17
        else
            ioptorbcontr = 16
        end if

    case (316)
        ioptorbcontr = 80

    case (320) ! Smoothed STO + Gaussian
        if (j .eq. 1) then
            ioptorbcontr = 28
        else
            ioptorbcontr = 16
        end if

    case (321) ! New regolarized STO contraction
        ioptorbcontr = 57

    case (322) ! mixed orbitals first NewSTO -gaussians-gaussian x r^3
        if (j .eq. 1) then
            ioptorbcontr = 57 ! the first is Exp * P(x)   (no cusp)
        elseif (mod(j, 2) .eq. 0) then
            ioptorbcontr = 16 ! Gaussian
        else
            ioptorbcontr = 60 ! Gaussian x r^3
        end if

    case (323) ! mixed orbitals first 1s(NewSTO) - 3s(STO)- 4s(STO)
        if (j .eq. 1) then
            ioptorbcontr = 57 ! the first is Exp * P(x)   (no cusp)
        elseif (mod(j, 2) .eq. 0) then
            ioptorbcontr = 10 ! r^2 x Exp
        else
            ioptorbcontr = 12 ! r^3 x Exp
        end if

    case (324) ! mixed orbitals first 1s(NewSTO) - 2s( (1+x)e^-x ) - 3s(STO)- 4s(STO)
        if (mod(j, 4) .eq. 1) then
            ioptorbcontr = 57 ! P(zr) * Exp(-zr)   (1s, no cusp)
        elseif (mod(j, 4) .eq. 2) then
            ioptorbcontr = 34 ! (1+zx) * Exp(-zr)  (2s, no cusp)
        elseif (mod(j, 4) .eq. 3) then
            ioptorbcontr = 10 ! r^2 * Exp
        elseif (mod(j, 4) .eq. 0) then
            ioptorbcontr = 12 ! r^3 * Exp
        end if

    case (325) ! mixed orbitals Exp +   Gaussian + Gaussian x r^2+ Gaussian x r^3
        if (j .eq. 1) then
            ioptorbcontr = 34
        else
            if (mod(j, 3) .eq. 2) then
                ioptorbcontr = 16
            elseif (mod(j, 3) .eq. 0) then
                ioptorbcontr = 17
            else
                ioptorbcontr = 60
            end if
        end if

    case (3000) ! Jastrow gaussian contraction
        ioptorbcontr = 100

    case (3001) ! Jastrow gaussian mixed contraction
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 100 ! Simple Gaussian
        else
            ioptorbcontr = 152 ! Gaussian x r^3
        end if

    case (3002) !  Gaussian + Gaussian x r^2
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 100
        else
            ioptorbcontr = 131
        end if

    case (3011) ! same as 311 for Jastrow
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 34
        else
            ioptorbcontr = 10
        end if

    case (3012) ! same as 312 for Jastrow
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 34
        else
            ioptorbcontr = 12
        end if

    case (3024) ! mixed orbitals first 1s(NewSTO) - 2s( (1+x)e^-x ) - 3s(STO)- 4s(STO)
        if (mod(j, 4) .eq. 1) then
            ioptorbcontr = 57 ! P(zr) * Exp(-zr)   (1s, no cusp)
        elseif (mod(j, 4) .eq. 2) then
            ioptorbcontr = 34 ! (1+zx) * Exp(-zr)  (2s, no cusp)
        elseif (mod(j, 4) .eq. 3) then
            ioptorbcontr = 10 ! r^2 * Exp
        elseif (mod(j, 4) .eq. 0) then
            ioptorbcontr = 12 ! r^3 * Exp
        end if

        !     P contraction
    case (400) ! gaussian contraction, works also for PBC
        ioptorbcontr = 36
    case (401) ! Exp  contraction
        ioptorbcontr = 20 ! Exp  contraction
    case (402) ! Exp x r  contraction
        ioptorbcontr = 22
    case (407) ! mixed orbitals p first Exp other gaussians
        if (j .eq. 1) then
            ioptorbcontr = 20
        else
            ioptorbcontr = 36
        end if
    case (408) ! mixed orbitals p first Exp- gaussian -gaussian x r
        if (j .eq. 1) then
            ioptorbcontr = 20 ! Exp
        elseif (mod(j, 2) .eq. 0) then
            ioptorbcontr = 36 ! Gaussian
        else
            ioptorbcontr = 62 ! Gaussian x r
        end if
    case (410) !  mixed orbitals p  gaussian -gaussian x r
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 36
        else
            ioptorbcontr = 62
        end if
    case (411) ! mixed orbital  p   Exp - Exp x r
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 20
        else
            ioptorbcontr = 22
        end if
    case (412) ! mixed orbital Slater  p
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 22
        else
            ioptorbcontr = 20
        end if
    case (415) !  mixed orbitals p  gaussian -gaussian x r
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 62
        else
            ioptorbcontr = 36
        end if
    case (416)
        ioptorbcontr = 82

    case (4000) ! Jastrow gaussian contraction
        ioptorbcontr = 103
    case (4001) ! Jastrow gaussian x r  contraction
        ioptorbcontr = 150
    case (4002) !  mixed basis gaussian - gaussian x r
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 103 ! Jastrow gaussian contraction
        else
            ioptorbcontr = 150 ! Jastrow gaussian x r  contraction
        end if
    case (4003) !  mixed basis gaussian-gaussian x r^2
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 103 ! Jastrow gaussian contraction
        else
            ioptorbcontr = 1101 ! Jastrow gaussian x r^2  contraction
        end if
    case (4010) ! Slater  Exp  for Jastrow
        ioptorbcontr = 20
    case (4011) ! mixed orbital p  Exp  for Jastrow   (same as 411)
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 20
        else
            ioptorbcontr = 22
        end if
    case (4012) ! mixed orbital p   Exp  for Jastrow   (same as 412)
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 22
        else
            ioptorbcontr = 20
        end if
        !  D orbitals
    case (500) ! gaussian contraction
        if (ieso) then
            ioptorbcontr = 37
        else
            !         ioptorbcontr=66
            ioptorbcontr = 68 ! minimal rmucos orbital
        end if
    case (504) ! gaussian contraction old  version
        ioptorbcontr = 37
    case (501) ! Exp  contraction
        ioptorbcontr = 30
    case (502) ! Exp x r contraction
        ioptorbcontr = 33
    case (503) ! mixed orbitals Slater d
        if (mod(j, 2) .eq. 1) then
            ioptorbcontr = 30
        else
            ioptorbcontr = 33
        end if
    case (507) ! mixed orbitals d Exp -gaussian
        if (j .eq. 1) then
            ioptorbcontr = 30 ! Exp
        else
            if (ieso) then
                ioptorbcontr = 37 ! Gaussian
            else
                !       ioptorbcontr=66  ! Gaussian
                ioptorbcontr = 68 ! minimal rmucos orbital
            end if
        end if

    case (508)
        if (j .eq. 1) then ! mixed orbitals d Exp -gaussian- gaussian x r
            ioptorbcontr = 30 ! Exp
        elseif (mod(j, 2) .eq. 0) then
            ioptorbcontr = 37 ! Gaussian
        else
            ioptorbcontr = 64 ! Gaussian x r
        end if

    case (510)
        if (mod(j, 2) .eq. 1) then ! mixed orbitals d -gaussian- gaussian x r
            if (ieso) then
                ioptorbcontr = 37
            else
                !       ioptorbcontr=66
                ioptorbcontr = 68 ! minimal rmucos orbital
            end if
        else
            ioptorbcontr = 64
        end if

    case (516)
        ioptorbcontr = 84

    case (5000) ! Jastrow gaussian contraction
        ioptorbcontr = 147 ! gaussian

        ! F orbitals
        !     case(3200)        ! constant orbital contracted
        !     ioptorbcontr=200
    case (600) ! Gaussian contracted minimal rmucos power
        if (.not. ieso) then
            ioptorbcontr = 58
        else
            ioptorbcontr = 48
        end if

    case (601) ! Slater contracted f orbitals
        ioptorbcontr = 70

    case (602) ! Gaussian contracted previous version
        ioptorbcontr = 48

    case (607) ! mixed f- orbitals Exp -gaussian
        if (j .eq. 1) then
            ioptorbcontr = 70 ! Exp
        else
            if (.not. ieso) then
                ioptorbcontr = 58
            else
                ioptorbcontr = 48
            end if
        end if

    case (616)
        ioptorbcontr = 86

    case (6000) ! Jastrow gaussian contraction
        ioptorbcontr = 154 ! gaussian

        ! G orbitals
    case (700) ! contracted Gaussian for g orbitals
        if (.not. ieso) then
            ioptorbcontr = 53
        else
            ioptorbcontr = 51
        end if

    case (701) ! contracted Slater for g orbitals
        ioptorbcontr = 55

    case (702)
        ioptorbcontr = 51

    case (707) ! mixed g- orbitals Exp -gaussian
        if (j .eq. 1) then
            ioptorbcontr = 55 ! Exp
        else
            if (.not. ieso) then
                ioptorbcontr = 53
            else
                ioptorbcontr = 51
            end if
        end if

    case (716)
        ioptorbcontr = 88

        ! H orbitals
    case (800) ! contracted Gaussian for H orbitals
        ioptorbcontr = 72
        ! I orbitals
    case (900)
        ioptorbcontr = 73
    case (1000000) ! molecular orbital
        ioptorbcontr = 999999
    case (900000)
        ioptorbcontr = 999998
    case default
        ioptorbcontr = ioptorb
        !     The  program will stop in ioptorbder or makefun makefun_pbc
    end select
    return
end function ioptorbcontr

function multioptorb(ioptorb)
    implicit none
    integer multioptorb, ioptorb

    select case (ioptorb)
    case (900000, 1000000)
        multioptorb = 1 ! hybrid or molecular orbitals
    case (10, 12, 16:17, 34, 57, 60, 80, 100, 131, 152, 161, 200, 300:399, 3000:3999)
        multioptorb = 1 ! s
    case (20, 22, 36, 62, 82, 103, 150, 400:499, 4000:4999)
        multioptorb = 3 ! p
    case (30, 33, 37, 64, 68, 147, 500:599, 5000:5999)
        multioptorb = 5 ! d
    case (48, 58, 70, 86, 154, 600:699)
        multioptorb = 7 ! f
    case (51, 53, 55, 88, 700:799)
        multioptorb = 9 ! g
    case (72, 800)
        multioptorb = 11 ! h
    case (73, 900)
        multioptorb = 13 ! i
    case(90:99)
        multioptorb = (ioptorb - 90 + 2) * (ioptorb - 90 + 1) / 2
    case default
        multioptorb = 0
    end select
    return
end function multioptorb

function check_multioptorb(multi, ioptorb, switch)
    implicit none
    integer :: check_multioptorb, multi, ioptorb
    character :: switch
    integer, external :: multioptorb

    check_multioptorb = 0
    if (multioptorb(ioptorb) .ne. 0) then
        if (multi .ne. multioptorb(ioptorb)) then
            if (switch .eq. 'E') then
                write (6, *) ' ERROR: the multiplicity of orbital ', ioptorb
                write (6, *) '              should be corrected as ', multioptorb(ioptorb)
            elseif (switch .eq. 'W') then
                multi = multioptorb(ioptorb)
                write (6, *) ' Warning the multiplicity of orbital ', ioptorb, &
                        &' has been replaced by ', multi
            else
                write (6, *) 'Non-exist switch in check_multioptorb!'
                stop
            end if
            check_multioptorb = 1
        end if
    else
        write (6, *) ' Warning the default multiplicity of ', ioptorb, &
                &' is not defined! Please update function multioptorb!'
    end if

    return
end function check_multioptorb
