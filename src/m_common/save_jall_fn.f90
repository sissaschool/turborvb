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

!> @brief Jastrow factor support vector update for DMC and LR-DMC
!>
!> This subroutine updates the support vector jastrowall_ee used in
!> Diffusion Monte Carlo (DMC) and Lattice Regularized DMC (LR-DMC)
!> calculations. The updated vector is used by the fast update routine
!> "uptabtot_new.f90" to speed up Jastrow factor calculations.
!>
!> The subroutine handles both 3/4-body Jastrow factors and wave function
!> components, with support for spin-dependent calculations and OpenMP
!> offloading for GPU acceleration.
!>
!> @param[in] yesfn Flag indicating whether to update gradients and Laplacians
!> @param[in,out] jastrowall_ee Support vector for Jastrow factor calculations
!> @param[in] winvjbar Jastrow factor basis coefficients for up electrons
!> @param[in] winvjbarsz Jastrow factor basis coefficients for spin-dependent terms
!> @param[in] winvj Jastrow factor basis vectors
!> @param[in,out] psip Working array for matrix operations
!>
!> @details
!> The subroutine performs the following operations:
!> 1. Initializes working array psip to zero
!> 2. If yesfn is true, updates gradients and Laplacians:
!>    - Processes 3/4-body terms (k = 1 to ip4)
!>    - Processes lattice terms (k = 1 to indt)
!>    - Uses matrix multiplication to compute Jastrow contributions
!>    - Handles both up and down electron components
!>    - Updates spin-dependent terms if iessz is true
!> 3. Updates wave function components:
!>    - Computes Jastrow factor contributions to wave function
!>    - Handles both up and down electron components
!>    - Updates spin-dependent terms if iessz is true
!>
!> @note The 2-body part of jastrowall_ee is computed in "compute_fast.f90"
!> @note The subroutine supports OpenMP offloading for GPU acceleration
!> @note Handles both paired (ipj=1) and unpaired (ipj=2) electron configurations
!> @note Spin-dependent calculations are performed when iessz is true
!> @note Used in DMC and LR-DMC calculations for efficient Jastrow updates
subroutine save_jall(yesfn, jastrowall_ee, winvjbar, winvjbarsz, winvj, psip)

    use constants, only: ip4, ipj
    use allio, only: nel, nelup, neldo, indt, indt4j, nelorbj, nelorbjh, iessz

    implicit none
    ! input
    real(8), intent(inout) :: jastrowall_ee(nelup + neldo, nelup + neldo, 0:indt4j), psip(nel, nel)
    real(8), intent(in) :: winvj(max(nelorbjh, 1), 0:indt4j, *), winvjbar(*), winvjbarsz(*)
    integer :: nelused
    ! local
    integer :: i, k, nelorbj5
    logical yesfn

    if (nelorbj .le. 0) return

    nelorbj5 = nelorbj*(indt4j + 1)
    nelused = nelup + mod(ipj, 2)*neldo
    !!!!!!!!!!
    ! 3/4body
    !!!!!!!!!!

    ! in the following I assume that indt4j.ge.indt+ip4,
    ! therefore the basis vector winvj must be already
    ! completely filled.

    ! clearing working space
    psip = 0.d0
#ifdef _OFFLOAD
!$omp target update to (psip)
#endif
    if (yesfn) then
        ! update gradients and laplacians
        do k = 1, ip4
            call dgemm_('T', 'N', nelused, nel, nelorbjh, 1.d0, winvj(1, indt + k, 1) &
                    &, nelorbj5, winvjbar, ipj*nelorbjh, 0.d0, psip, nel)
            if (ipj == 2) then
                call dgemm_('T', 'N', neldo, nel, nelorbjh, 1.d0, winvj(1, indt + k, 1 + nelup) &
                        &, nelorbj5, winvjbar(1 + nelorbjh), 2*nelorbjh, 0.d0, psip(1 + nelup, 1), nel)
            end if
#ifdef _OFFLOAD
!$omp target update from(psip)
#endif
            call upjastrowall(nel, jastrowall_ee(1, 1, indt + k), psip)
            if (iessz) then
                call dgemm('T', 'N', nel, nel, nelorbjh, 1.d0, winvj(1, indt + k, 1) &
                        &, nelorbj5, winvjbarsz, nelorbjh, 0.d0, psip, nel)
                call upjastrowall_sz(nel, nelup, jastrowall_ee(1, 1, indt + k), psip)
            end if
        end do
        ! finally the 1:indt positions in the LRDMC lattice
        do k = 1, indt
            call dgemm_('T', 'N', nelused, nel, nelorbjh, 1.d0, winvj(1, k, 1)      &
                    &, nelorbj5, winvjbar, ipj*nelorbjh, 0.d0, psip, nel)
            if (ipj == 2) then
                call dgemm_('T', 'N', neldo, nel, nelorbjh, 1.d0, winvj(1, k, 1 + nelup)      &
                        &, nelorbj5, winvjbar(1 + nelorbjh), 2*nelorbjh, 0.d0, psip(1 + nelup, 1), nel)
            end if

#ifdef _OFFLOAD
!$omp target update from(psip)
#endif
            call upjastrowallfat(nel, jastrowall_ee(1, 1, k), psip)
            if (iessz) then
                call dgemm('T', 'N', nel, nel, nelorbjh, 1.d0, winvj(1, k, 1) &
                        &, nelorbj5, winvjbarsz, nelorbjh, 0.d0, psip, nel)
                call upjastrowallfat_sz(nel, nelup, jastrowall_ee(1, 1, k), psip)
            end if
        end do
    end if

    ! then update jastrowall_ee for the wave function
    !call dgemm('T','N',nel,nelused,ipj*nelorbjh,1.d0,winvjbar   &
    !    &,ipj*nelorbjh,winvj,nelorbj5,0.d0,psip,nel)
    call dgemm_('T', 'N', nel, nelused, nelorbjh, 1.d0, winvjbar   &
            &, ipj*nelorbjh, winvj, nelorbj5, 0.d0, psip, nel)
    if (ipj == 2) then
        call dgemm_('T', 'N', nel, neldo, nelorbjh, 1.d0, winvjbar(1 + nelorbjh)   &
                &, 2*nelorbjh, winvj(1, 0, 1 + nelup), nelorbj5, 0.d0, psip(1, 1 + nelup), nel)
    end if

    !  write(*,*) " psip in save_jastrowall:", psip(1:nel,1:nel)
#ifdef _OFFLOAD
!$omp target update from(psip)
#endif
    call upjastrowallpsi(nel, jastrowall_ee, psip)
    if (iessz) then
        call dgemm('T', 'N', nel, nel, nelorbjh, 1.d0, winvjbarsz &
                &, nelorbjh, winvj, nelorbj5, 0.d0, psip, nel)
        call upjastrowallpsi_sz(nel, nelup, jastrowall_ee, psip)
    end if

    return

end subroutine save_jall
