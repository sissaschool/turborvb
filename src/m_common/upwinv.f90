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

subroutine upwinv(nel, jel, indt, nelorb, winv, v, psi)
    use constants, only: yes_ontarget
    implicit none

    ! argument parameters
    integer, intent(in) :: nel, indt, jel, nelorb
    real*8, intent(in) :: psi(indt, nel), v(nel)
    real*8, intent(inout) :: winv(nel, indt)

    ! local variables
    integer i, j

#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2) if(yes_ontarget)
#endif
    do j = 1, nel
        !              if(j.ne.jel) then
        do i = 1, indt
            winv(j, i) = winv(j, i) + v(j)*psi(i, j)
        end do
        !              else
        !              do i=1,indt
        !              winv(j,i)=psi(i,j)*(1.d0+v(jel))
        !              enddo
        !              endif
    end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
    do i = 1, indt
        winv(jel, i) = psi(i, jel)*(1.d0 + v(jel))
    end do
    return
end

subroutine upwinv_complex(nel, jel, indt, nelorb, winv, v, psi)
    use constants, only: yes_ontarget
    implicit none

    ! argument parameters
    integer, intent(in) :: nel, indt, jel, nelorb
    complex*16, intent(in) :: psi(indt, nel), v(nel)
    complex*16, intent(inout) :: winv(nel, indt)

    ! local variables
    integer i, j
#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2) if(yes_ontarget)
#endif
    do j = 1, nel
        !              if(j.ne.jel) then
        do i = 1, indt
            winv(j, i) = winv(j, i) + v(j)*psi(i, j)
        end do
        !              else
        !              do i=1,indt
        !              winv(j,i)=psi(i,j)*(1.d0+v(jel))
        !              enddo
        !              endif
    end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
#ifdef _OFFLOAD
!$omp target teams distribute parallel do if(yes_ontarget)
#endif
    do i = 1, indt
        winv(jel, i) = psi(i, jel)*(1.d0 + v(jel))
    end do
    return
end
