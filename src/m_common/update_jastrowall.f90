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

! -----------------------------------------
! for updating gradient and laplacian part
! -----------------------------------------
subroutine upjastrowall(nel, jastrowall, psip)
    implicit none
    integer j, k, nel
    real*8 psip(nel, *), jastrowall(nel, *)

    do j = 1, nel
        do k = 1, nel
            if (k .ne. j) jastrowall(k, j) = jastrowall(k, j) + psip(j, k)
        end do
    end do
    return
end subroutine upjastrowall

subroutine upjastrowall_sz(nel, nelup, jastrowall, psip)
    implicit none
    integer j, k, nel, nelup
    real*8 psip(nel, *), jastrowall(nel, *)

    do j = 1, nelup
        do k = 1, nelup
            if (k .ne. j) jastrowall(k, j) = jastrowall(k, j) + psip(j, k)
        end do
        do k = nelup + 1, nel
            jastrowall(k, j) = jastrowall(k, j) - psip(j, k)
        end do
    end do
    do j = nelup + 1, nel
        do k = 1, nelup
            jastrowall(k, j) = jastrowall(k, j) - psip(j, k)
        end do
        do k = nelup + 1, nel
            if (k .ne. j) jastrowall(k, j) = jastrowall(k, j) + psip(j, k)
        end do
    end do
    return
end subroutine upjastrowall_sz

! -----------------------------------------
! for updating (1:indt) part of jastrowall
! pseudo and lattice positions in LRDMC
! -----------------------------------------
subroutine upjastrowallfat(nel, jastrowall, psipmu)
    implicit none
    integer nel, i, j
    real*8 jastrowall(nel, *), psipmu(nel, *)

    do i = 1, nel
        do j = 1, nel
            if (j .ne. i) jastrowall(j, i) = jastrowall(j, i) + psipmu(i, j)
        end do
    end do

    return
end subroutine upjastrowallfat

subroutine upjastrowallfat_sz(nel, nelup, jastrowall, psipmu)
    implicit none
    integer nel, nelup, i, j
    real*8 jastrowall(nel, *), psipmu(nel, *)

    do i = 1, nelup
        do j = 1, nelup
            if (j .ne. i) jastrowall(j, i) = jastrowall(j, i) + psipmu(i, j)
        end do
        do j = nelup + 1, nel
            jastrowall(j, i) = jastrowall(j, i) - psipmu(i, j)
        end do
    end do

    do i = nelup + 1, nel
        do j = 1, nelup
            jastrowall(j, i) = jastrowall(j, i) - psipmu(i, j)
        end do
        do j = nelup + 1, nel
            if (j .ne. i) jastrowall(j, i) = jastrowall(j, i) + psipmu(i, j)
        end do
    end do

    return
end subroutine upjastrowallfat_sz

! ------------------------------------------
! for updating the part related to the wave
! function in jastrowall
! ------------------------------------------
subroutine upjastrowallpsi(nel, jastrowall, psip)
    implicit none
    integer j, k, nel
    real*8 psip(nel, *), jastrowall(nel, *)

    do j = 1, nel
        do k = 1, nel
            if (j .ne. k) jastrowall(j, k) = jastrowall(j, k) + psip(k, j)
        end do
    end do

    return
end subroutine upjastrowallpsi

subroutine upjastrowallpsi_sz(nel, nelup, jastrowall, psip)
    implicit none
    integer j, k, nel, nelup
    real*8 psip(nel, *), jastrowall(nel, *)

    do j = 1, nelup
        do k = j + 1, nelup
            jastrowall(k, j) = jastrowall(k, j) + psip(k, j)
        end do
        do k = nelup + 1, nel
            jastrowall(k, j) = jastrowall(k, j) - psip(k, j)
        end do
    end do
    do j = nelup + 1, nel
        do k = j + 1, nel
            jastrowall(k, j) = jastrowall(k, j) + psip(k, j)
        end do
    end do

    do j = 1, nel
        do k = j + 1, nel
            jastrowall(j, k) = jastrowall(k, j)
        end do
    end do

    return
end subroutine upjastrowallpsi_sz
