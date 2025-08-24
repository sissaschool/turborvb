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

!> @brief Update wave function inverse matrix for real wave functions with two terms
!>
!> This subroutine updates the wave function inverse matrix WINV using
!> a rank-2 update formula with two sets of coefficients and wave functions.
!> The update follows the Sherman-Morrison-Woodbury formula for rank-2 updates.
!>
!> Parameters
!> ----------
!> NEL : integer, in
!>     Number of electrons.
!> INDT : integer, in
!>     Number of basis functions (dimension of orbital space).
!> WINV : real*8 array, inout
!>     Wave function inverse matrix (NEL × INDT).
!>     On entry: current inverse matrix.
!>     On exit: updated inverse matrix.
!> AINV : real*8 array, in
!>     First vector of coefficients for the update (NEL).
!> AINVN : real*8 array, in
!>     Second vector of coefficients for the update (NEL).
!> PSI : real*8 array, in
!>     Wave function matrix with two components (INDT × NEL × 2).
!>
!> Notes
!> -----
!> - Uses OpenMP offloading for parallel computation on accelerators.
!> - The update follows: WINV_new = WINV_old + ainv * psi(:,:,1)^T + ainvn * psi(:,:,2)^T.
!> - Optimized for quantum Monte Carlo applications with rank-2 updates.
!>
!> Algorithm
!> ---------
!> WINV(j,i) += ainv(j) * psi(i,j,1) + ainvn(j) * psi(i,j,2)
!>
!> Example
!> -------
!> Used in variational Monte Carlo for efficient wave function updates with two terms.
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

!> @brief Update wave function inverse matrix for complex wave functions with two terms
!>
!> This subroutine updates the complex wave function inverse matrix WINV using
!> a rank-2 update formula with two sets of complex coefficients and wave functions.
!> The update follows the Sherman-Morrison-Woodbury formula for complex rank-2 updates.
!>
!> Parameters
!> ----------
!> NEL : integer, in
!>     Number of electrons.
!> INDT : integer, in
!>     Number of basis functions (dimension of orbital space).
!> WINV : complex*16 array, inout
!>     Complex wave function inverse matrix (NEL × INDT).
!>     On entry: current inverse matrix.
!>     On exit: updated inverse matrix.
!> AINV : complex*16 array, in
!>     First vector of complex coefficients for the update (NEL).
!> AINVN : complex*16 array, in
!>     Second vector of complex coefficients for the update (NEL).
!> PSI : complex*16 array, in
!>     Complex wave function matrix with two components (INDT × NEL × 2).
!>
!> Notes
!> -----
!> - Uses OpenMP offloading for parallel computation on accelerators.
!> - The update follows: WINV_new = WINV_old + ainv * psi(:,:,1)^T + ainvn * psi(:,:,2)^T.
!> - Complex arithmetic is used throughout the computation.
!> - Optimized for quantum Monte Carlo applications with complex rank-2 updates.
!>
!> Algorithm
!> ---------
!> WINV(j,i) += ainv(j) * psi(i,j,1) + ainvn(j) * psi(i,j,2)
!>
!> Example
!> -------
!> Used in variational Monte Carlo for efficient complex wave function updates with two terms.
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

!> @brief Update wave function inverse matrix for real Pfaffian wave functions
!>
!> This subroutine updates the wave function inverse matrices WINVUP and WINVDO
!> for Pfaffian wave functions. It handles spin-up and spin-down electrons
!> separately and updates the appropriate matrix based on the electron index.
!>
!> Parameters
!> ----------
!> NELC : integer, in
!>     Index of the current electron being updated.
!> NELUP : integer, in
!>     Number of spin-up electrons.
!> NELDO : integer, in
!>     Number of spin-down electrons.
!> NMOL : integer, in
!>     Total number of molecular orbitals.
!> NMOLIPF : integer, in
!>     Number of molecular orbitals in Pfaffian.
!> NMOLSHIFT : integer, in
!>     Shift index for molecular orbitals.
!> INDT : integer, in
!>     Number of basis functions (dimension of orbital space).
!> WINVUP : real*8 array, inout
!>     Wave function inverse matrix for spin-up electrons (NELUP × INDT).
!> WINVDO : real*8 array, inout
!>     Wave function inverse matrix for spin-down electrons (NELDO × INDT).
!> PSI : real*8 array, in
!>     Pfaffian wave function matrix (NMOLIPF × INDT).
!> AINV : real*8 array, in
!>     Vector of coefficients for the update (NMOL).
!>
!> Notes
!> -----
!> - Uses OpenMP offloading for parallel computation on accelerators.
!> - Handles spin-up and spin-down electrons separately.
!> - For NELC ≤ NELUP: updates WINVUP using ainv(1:NMOLIPF).
!> - For NELC > NELUP: updates WINVDO using ainv(NMOLSHIFT+1:NMOL).
!> - Optimized for Pfaffian wave function calculations.
!>
!> Algorithm
!> ---------
!> If NELC ≤ NELUP:
!>   WINVUP(nelc,i) = sum(psi(1:nmolipf,i) * ainv(1:nmolipf))
!> Else:
!>   WINVDO(nelc-nelup,i) = sum(psi(1:nmolipf,i) * ainv(nmolshift+1:nmol))
!>
!> Example
!> -------
!> Used in variational Monte Carlo for Pfaffian wave function updates.
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
!$omp target teams distribute parallel do reduction(+:winvup)
#else
!$omp parallel do default(shared) private(i,j) reduction(+:winvup)
#endif
        do i = 1, indt
            do j = 1, nmolipf
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
!$omp target teams distribute parallel do reduction(+:winvup)
#else
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
        end do
    end if
    end if
    return
end

!> @brief Update wave function inverse matrix for complex Pfaffian wave functions
!>
!> This subroutine updates the complex wave function inverse matrices WINVUP and WINVDO
!> for complex Pfaffian wave functions. It handles spin-up and spin-down electrons
!> separately and updates the appropriate matrix based on the electron index.
!>
!> Parameters
!> ----------
!> NELC : integer, in
!>     Index of the current electron being updated.
!> NELUP : integer, in
!>     Number of spin-up electrons.
!> NELDO : integer, in
!>     Number of spin-down electrons.
!> NMOL : integer, in
!>     Total number of molecular orbitals.
!> NMOLIPF : integer, in
!>     Number of molecular orbitals in Pfaffian.
!> NMOLSHIFT : integer, in
!>     Shift index for molecular orbitals.
!> INDT : integer, in
!>     Number of basis functions (dimension of orbital space).
!> WINVUP : complex*16 array, inout
!>     Complex wave function inverse matrix for spin-up electrons (NELUP × INDT).
!> WINVDO : complex*16 array, inout
!>     Complex wave function inverse matrix for spin-down electrons (NELDO × INDT).
!> PSI : complex*16 array, in
!>     Complex Pfaffian wave function matrix (NMOLIPF × INDT).
!> AINV : complex*16 array, in
!>     Vector of complex coefficients for the update (NMOL).
!>
!> Notes
!> -----
!> - Uses OpenMP offloading for parallel computation on accelerators.
!> - Handles spin-up and spin-down electrons separately.
!> - For NELC ≤ NELUP: updates WINVUP using ainv(1:NMOLIPF).
!> - For NELC > NELUP: updates WINVDO using ainv(NMOLSHIFT+1:NMOL).
!> - Complex arithmetic is used throughout the computation.
!> - Optimized for complex Pfaffian wave function calculations.
!>
!> Algorithm
!> ---------
!> If NELC ≤ NELUP:
!>   WINVUP(nelc,i) = sum(psi(1:nmolipf,i) * ainv(1:nmolipf))
!> Else:
!>   WINVDO(nelc-nelup,i) = sum(psi(1:nmolipf,i) * ainv(nmolshift+1:nmol))
!>
!> Example
!> -------
!> Used in variational Monte Carlo for complex Pfaffian wave function updates.
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
            end do
        end if
    end if
    return
end
