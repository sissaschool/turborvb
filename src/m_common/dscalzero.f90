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

subroutine dscalzero(n, zero, vet, m)
    implicit none
    integer n, m, i
    real*8 zero, vet(m*(n - 1) + 1)
    do i = 1, n
        vet(m*(i - 1) + 1) = zero
    end do
    return
end
subroutine dscalzero_(n, zero, vet, m)
    implicit none
    integer n, m, i
    real*8 zero, vet(m*(n - 1) + 1)
#ifdef  _OFFLOAD
!$omp target teams distribute parallel do
#endif
    do i = 1, n
        vet(m*(i - 1) + 1) = zero
    end do
#ifdef  _OFFLOAD
!$omp end target teams distribute parallel do
#endif
    return
end
subroutine dscalzero__(n, zero, vet, m)
    use constants, only: yes_ontarget
    implicit none
    integer n, m, i
    real*8 zero, vet(m*(n - 1) + 1)
!    if(n.le.16384.and.m.eq.1) then
!    do i = 1, n
!    vet(m*(i-1)+1) = zero
!    enddo
!#ifdef  _OFFLOAD
!!$omp target update to (vet)  if(yes_ontarget)
!#endif
!    else
#ifdef  _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
    do i = 1, n
        vet(m*(i - 1) + 1) = zero
    end do
#ifdef  _OFFLOAD
!$omp end target teams distribute parallel do
#endif
    return
end
subroutine iscalzero_(n, zero, vet, m)
    implicit none
    integer n, m, i
    integer zero, vet(m*(n - 1) + 1)
#ifdef  _OFFLOAD
!$omp target teams distribute parallel do
#endif
    do i = 1, n
        vet(m*(i - 1) + 1) = zero
    end do
#ifdef  _OFFLOAD
!$omp end target teams distribute parallel do
#endif
    return
end
subroutine zscalzero_(n, zero, vet, m)
    implicit none
    integer n, m, i
    complex*16 zero, vet(m*(n - 1) + 1)
#ifdef  _OFFLOAD
!$omp target teams distribute parallel do
#endif
    do i = 1, n
        vet(m*(i - 1) + 1) = zero
    end do
#ifdef  _OFFLOAD
!$omp end target teams distribute parallel do
#endif
    return
end
