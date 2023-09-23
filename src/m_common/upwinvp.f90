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

!> subroutine to update the matrix winv, for real(8)
subroutine upwinvp(nel, indt, winv, ainv, ainvn, psi)
    use constants, only: yes_ontarget
    implicit none

    ! argument parameters
    integer, intent(in) :: nel, indt
    real*8, intent(in) :: psi(indt, nel, 2), ainv(nel), ainvn(nel)
    real*8, intent(inout) :: winv(nel, indt)

    ! local variables
    integer :: i, j

#ifdef _OFFLOAD
    if (yes_ontarget) then
!$omp target teams distribute parallel do collapse(2)
        do i = 1, indt
            do j = 1, nel
                winv(j, i) = winv(j, i) + ainv(j)*psi(i, j, 1) + ainvn(j)*psi(i, j, 2)
            end do
        end do
    else
        do i = 1, indt
            do j = 1, nel
                winv(j, i) = winv(j, i) + ainv(j)*psi(i, j, 1) + ainvn(j)*psi(i, j, 2)
            end do
        end do
    end if
#else
    do i = 1, indt
        do j = 1, nel
            winv(j, i) = winv(j, i) + ainv(j)*psi(i, j, 1) + ainvn(j)*psi(i, j, 2)
        end do
    end do
#endif
    return
end

!> subroutine to update the matrix winv, for complex(16)
subroutine upwinvp_complex(nel, indt, winv, ainv, ainvn, psi)
    use constants, only: yes_ontarget
    implicit none

    ! argument parameters
    integer, intent(in) :: nel, indt
    complex*16, intent(in) :: psi(indt, nel, 2), ainv(nel), ainvn(nel)
    complex*16, intent(inout) :: winv(nel, indt)

    ! local variables
    integer i, j

#ifdef _OFFLOAD
!$omp target teams distribute parallel do collapse(2) if(yes_ontarget)
#endif
    do i = 1, indt
        do j = 1, nel
            winv(j, i) = winv(j, i) + ainv(j)*psi(i, j, 1) + ainvn(j)*psi(i, j, 2)
        end do
    end do
#ifdef _OFFLOAD
!$omp end target teams distribute parallel do
#endif
    return
end

!> subroutine to update the matrix winv of Pfaffian, for real(8)
subroutine upwinvp_pfaff(nelc, nelup, neldo, nmol, nmolipf, nmolshift, indt, winvup, winvdo, psi, ainv)
    use constants, only: yes_ontarget
    implicit none

    ! argument parameters
    integer, intent(in) :: nelc, nelup, neldo, nmol, nmolipf, nmolshift, indt
    real*8, intent(in) :: psi(nmolipf, indt), ainv(nmol)
    real*8, intent(inout) :: winvup(nelup, indt), winvdo(neldo, indt)

    ! local variables
    integer i, j, nelcdo
    real*8 csum

    if (yes_ontarget) then
    if (nelc .le. nelup) then
#ifdef _OFFLOAD
!!$omp target teams distribute  parallel do private(csum)
!$omp target teams distribute parallel do reduction(+:winvup)
#else
!!$omp parallel do default(shared) private(i,j,csum) reduction(+:csum)
!$omp parallel do default(shared) private(i,j) reduction(+:winvup)
#endif
        do i = 1, indt
            do j = 1, nmolipf
                !csum = csum + psi(j, i) * ainv(j)
                winvup(nelc, i) = sum(psi(1:nmolipf, i)*ainv(1:nmolipf))
            end do
        end do
#ifdef _OFFLOAD
!$omp end target teams distribute  parallel do
#else
!$omp end parallel do
#endif

    else
        nelcdo = nelc - nelup
#ifdef _OFFLOAD
!!$omp target teams distribute  parallel do private(csum)
!$omp target teams distribute parallel do reduction(+:winvup)
#else
!!$omp parallel do default(shared) private(i,j,csum) reduction(+:csum)
!$omp parallel do default(shared) private(i,j) reduction(+:winvup)
#endif
        do i = 1, indt
            csum = 0.d0
            do j = 1, nmolipf
                winvdo(nelc - nelup, i) = sum(psi(1:nmolipf, i)*ainv(nmolshift + 1:nmol))
            end do
        end do
#ifdef _OFFLOAD
!$omp end target teams distribute  parallel do
#else
!$omp end parallel do
#endif
    end if
    else
    if (nelc .le. nelup) then
!$omp parallel do default(shared) private(i,j,csum)
        do i = 1, indt
            csum = 0.d0
            do j = 1, nmolipf
                csum = csum + psi(j, i)*ainv(j)
                !          winvup(nelc,i)=sum(psi(1:nmolipf,i)*ainv(1:nmolipf))
            end do
            winvup(nelc, i) = csum
        end do
    else
        nelcdo = nelc - nelup
!$omp parallel do default(shared) private(i,j,csum)
        do i = 1, indt
            csum = 0.d0
            do j = 1, nmolipf
                csum = csum + psi(j, i)*ainv(nmolshift + j)
            end do
            winvdo(nelcdo, i) = csum
            !          winvdo(nelc-nelup,i)=sum(psi(1:nmolipf,i)*ainv(nmolshift+1:nmol))
        end do
    end if
    end if
    return
end

!> subroutine to update the matrix winv of Pfaffian, for complex(16)
subroutine upwinvp_pfaff_complex(nelc, nelup, neldo, nmol, nmolipf, nmolshift, indt, winvup, winvdo, psi, ainv)
    use constants, only: yes_ontarget
    implicit none

    ! argument parameters
    integer, intent(in) :: nelc, nelup, neldo, nmol, nmolipf, nmolshift, indt
    complex*16, intent(in) :: psi(nmolipf, indt), ainv(nmol)
    complex*16, intent(inout) :: winvup(nelup, indt), winvdo(neldo, indt)

    ! local variables
    integer i, j, nelcdo
    complex*16 csum
    if (yes_ontarget) then
        if (nelc .le. nelup) then
#ifdef _OFFLOAD
!$omp target teams distribute  parallel do  private(csum)
#else
!$omp parallel do default(shared) private(i,j,csum)
#endif
            do i = 1, indt
                csum = (0.d0, 0.d0)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                do j = 1, nmolipf
                    csum = csum + psi(j, i)*ainv(j)
                    !          winvup(nelc,i)=sum(psi(1:nmolipf,i)*ainv(1:nmolipf))
                end do
                winvup(nelc, i) = csum
            end do
#ifdef _OFFLOAD
!$omp end target teams distribute  parallel do
#else
!$omp end parallel do
#endif
        else
            nelcdo = nelc - nelup
#ifdef _OFFLOAD
!$omp target teams distribute  parallel do  private(csum)
#else
!$omp parallel do default(shared) private(i,j,csum)
#endif
            do i = 1, indt
                csum = (0.d0, 0.d0)
#ifdef _OFFLOAD
!$omp parallel do reduction(+:csum)
#endif
                do j = 1, nmolipf
                    csum = csum + psi(j, i)*ainv(nmolshift + j)
                end do
                winvdo(nelcdo, i) = csum
                !          winvdo(nelc-nelup,i)=sum(psi(1:nmolipf,i)*ainv(nmolshift+1:nmol))
            end do
#ifdef _OFFLOAD
!$omp end target teams distribute  parallel do
#else
!$omp end parallel do
#endif
        end if
    else
        if (nelc .le. nelup) then
!$omp parallel do default(shared) private(i,j,csum)
            do i = 1, indt
                csum = (0.d0, 0.d0)
                do j = 1, nmolipf
                    csum = csum + psi(j, i)*ainv(j)
                    !          winvup(nelc,i)=sum(psi(1:nmolipf,i)*ainv(1:nmolipf))
                end do
                winvup(nelc, i) = csum
            end do
        else
!$omp parallel do default(shared) private(i,j,csum)
            do i = 1, indt
                csum = (0.d0, 0.d0)
                do j = 1, nmolipf
                    csum = csum + psi(j, i)*ainv(nmolshift + j)
                end do
                winvdo(nelcdo, i) = csum
                !          winvdo(nelc-nelup,i)=sum(psi(1:nmolipf,i)*ainv(nmolshift+1:nmol))
            end do
        end if
    end if
    return
end
