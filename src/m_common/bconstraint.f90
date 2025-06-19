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

!> @file bconstraint.f90
!> @brief Module containing constraint handling for wave function optimization
!> @author TurboRVB group
!> @date 2022
!> @version 1.0

!> @brief Translates changes back to the original unconstrained representation
!> @details This subroutine handles the transformation of parameter changes from a 
!>          constrained optimization space back to the original unconstrained 
!>          representation used in the wave function. It supports both real and 
!>          complex wave functions, as well as different types of constraints.
!> 
!>          The subroutine performs the following operations:
!>          1. Maps constrained parameter changes to their original positions
!>          2. Handles symmetry constraints for AGP (Antisymmetrized Geminal Power) matrices
!>          3. Updates either the derivative matrix or the parameter matrix based on yes_update flag
!>          4. Supports both real (ipc=1) and complex (ipc=2) wave functions
!>          5. Handles EAGP (Extended AGP) Pfaffian matrices for complex wave functions
!> 
!>          The algorithm uses a mapping array jbradet to translate between constrained
!>          and unconstrained parameter spaces, where positive values indicate direct
!>          mapping and negative values indicate sign-flipped mapping.
!> 
!> @param[in] iessw Switch flag for constraint handling (0: no constraints, >0: apply constraints)
!> @param[in,out] derl Derivative matrix or parameter matrix to be updated
!> @param[in] nelorb Number of orbitals
!> @param[in] n Number of constrained parameters
!> @param[in] nozero Array of non-zero element indices
!> @param[in,out] psip Temporary workspace array
!> @param[in] econf Configuration energy array
!> @param[in] nw Number of walkers
!> @param[in] jbradet Mapping array from constrained to unconstrained parameters
!> @param[in] symmagp Logical flag for symmetric AGP matrix
!> @param[in] yes_update Logical flag to determine update mode
!> 
!> @note If iessw = 0, the subroutine returns immediately without any operations
!> 
!> @warning The subroutine modifies the derl array in-place based on the constraint mapping
!> 
!> @par Algorithm Details:
!> The subroutine implements a two-step process:
!> 1. <b>Parameter Mapping:</b> Uses jbradet array to map constrained parameters back to original positions
!> 2. <b>Matrix Update:</b> Updates either derivative matrix (upsim/upsimp) or parameter matrix based on yes_update
!> 
!> @par Real vs Complex Wave Functions:
!> - <b>Real (ipc=1):</b> Direct parameter mapping with sign handling
!> - <b>Complex (ipc=2):</b> Complex parameter mapping with real/imaginary part handling
!> 
!> @par EAGP Support:
!> For complex wave functions, the subroutine also handles EAGP Pfaffian matrices,
!> which require special antisymmetric structure maintenance.
!> 
!> @see constants module for ipc, ipf, deps parameters
!> @see cell module for cellscale, phase2pi parameters  
!> @see allio module for rank, yes_correct, nnozero_eagp, eagp_pfaff, ndiff parameters
subroutine bconstraint(iessw, derl, nelorb, n, nozero              &
        &, psip, econf, nw, jbradet, symmagp, yes_update)
    use constants, only: ipc, ipf, deps
    use cell, only: cellscale, phase2pi
    use allio, only: rank, yes_correct, nnozero_eagp, eagp_pfaff, ndiff
    implicit none
    
    !> @param[in] iessw Switch flag for constraint handling (0: no constraints, >0: apply constraints)
    integer, intent(in) :: iessw
    
    !> @param[in] nozero Array of non-zero element indices
    integer, intent(in) :: nozero(*)
    
    !> @param[in] nelorb Number of orbitals
    integer, intent(in) :: nelorb
    
    !> @param[in] n Number of constrained parameters
    integer, intent(in) :: n
    
    !> @param[in] nw Number of walkers
    integer, intent(in) :: nw
    
    !> @param[in] jbradet Mapping array from constrained to unconstrained parameters
    integer, intent(in) :: jbradet(*)
    
    !> @param[in] symmagp Logical flag for symmetric AGP matrix
    logical, intent(in) :: symmagp
    
    !> @param[in] yes_update Logical flag to determine update mode
    logical, intent(in) :: yes_update
    
    !> @param[in] econf Configuration energy array
    real*8, intent(in) :: econf(nw, *)
    
    !> @param[in,out] psip Temporary workspace array
    real*8, intent(inout) :: psip(*)
    
    !> @param[in,out] derl Derivative matrix or parameter matrix to be updated
    real*8, intent(inout) :: derl(*)
    
    ! Local variables
    integer :: i, ix, iy, j, ind, iesswread
    
    ! Early return if no constraints are applied
    if (iessw .eq. 0) return

    ! Handle real wave functions (ipc = 1)
    if (ipc .eq. 1) then
        ! Initialize workspace array to zero
        call dscalzero(n + nnozero_eagp, 0.d0, psip, 1)
        
        ! Map constrained parameters back to original positions
        ! Positive jbradet values: direct mapping
        ! Negative jbradet values: sign-flipped mapping
        do i = 1, n + nnozero_eagp
            j = jbradet(i)
            if (j .gt. 0) then
                psip(i) = econf(1, j)
            elseif (j .lt. 0) then
                psip(i) = -econf(1, -j)
            end if
        end do

        ! Update mode: add changes to existing matrix
        if (yes_update) then
            ! Update determinant matrix elements
            do i = 1, n
                if (psip(i) .ne. 0.d0)                                         &
                        &  call upsimp(derl, nelorb, nozero(i), psip(i), symmagp, ipf)
            end do
            
            ! Update EAGP Pfaffian matrix elements (antisymmetric structure)
            do i = 1, nnozero_eagp
                ind = i + n
                ! Convert linear index to 2D matrix indices
                iy = (nozero(ind) - 1)/ndiff + 1
                ix = nozero(ind) - (iy - 1)*ndiff
                j = jbradet(ind)
                if (j .gt. 0) then
                    eagp_pfaff(ix, iy) = eagp_pfaff(ix, iy) + psip(ind)
                    eagp_pfaff(iy, ix) = -eagp_pfaff(ix, iy)  ! Antisymmetric
                elseif (j .lt. 0) then
                    eagp_pfaff(ix, iy) = eagp_pfaff(ix, iy) - psip(ind)
                    eagp_pfaff(iy, ix) = -eagp_pfaff(ix, iy)  ! Antisymmetric
                end if
            end do

        ! Set mode: replace matrix elements
        else
            ! Set determinant matrix elements
            do i = 1, n
                if (abs(jbradet(i)) .ne. 0)                                    &
                        &  call upsim(derl, nelorb, nozero(i), psip(i), symmagp, ipf)
            end do
            
            ! Set EAGP Pfaffian matrix elements
            do i = 1, nnozero_eagp
                ind = i + n
                iy = (nozero(ind) - 1)/ndiff + 1
                ix = nozero(ind) - (iy - 1)*ndiff
                j = jbradet(ind)
                if (j .gt. 0) then
                    eagp_pfaff(ix, iy) = psip(ind)
                    eagp_pfaff(iy, ix) = -psip(ind)  ! Antisymmetric
                elseif (j .lt. 0) then
                    eagp_pfaff(ix, iy) = -psip(ind)
                    eagp_pfaff(iy, ix) = psip(ind)   ! Antisymmetric
                end if
            end do
        end if

    ! Handle complex wave functions (ipc = 2)
    else
        ! Initialize workspace array to zero (complex: 2x larger)
        call dscalzero(2*n + 2*nnozero_eagp, 0.d0, psip, 1)
        
        ! Handle symmetric AGP with correction
        if (symmagp .and. yes_correct) then
            ! Map complex parameters with symmetry constraints
            do i = 1, n + nnozero_eagp
                j = jbradet(i)
                if (j .gt. 0) then
                    ! Real and imaginary parts for symmetric case
                    psip(2*i - 1) = econf(1, 4*j - 3)  ! Real part
                    psip(2*i) = econf(1, 4*j - 1)      ! Imaginary part
                elseif (j .lt. 0) then
                    psip(2*i - 1) = -econf(1, -4*j - 3)  ! Sign-flipped real part
                    psip(2*i) = -econf(1, -4*j - 1)      ! Sign-flipped imaginary part
                end if
            end do

        ! Standard complex parameter mapping
        else
            do i = 1, n + nnozero_eagp
                j = jbradet(i)
                if (j .gt. 0) then
                    psip(2*i - 1) = econf(1, 2*j - 1)  ! Real part
                    psip(2*i) = econf(1, 2*j)          ! Imaginary part
                elseif (j .lt. 0) then
                    psip(2*i - 1) = -econf(1, -2*j - 1)  ! Sign-flipped real part
                    psip(2*i) = -econf(1, -2*j)          ! Sign-flipped imaginary part
                end if
            end do
        end if

        ! Update mode for complex wave functions
        if (yes_update) then
            ! Update complex determinant matrix elements
            do i = 1, n
                if (sum(abs(psip(2*i - 1:2*i))) .ne. 0.d0)                  &
                        & call upsimp_complex(derl, nelorb, nozero(i), psip(2*i - 1), symmagp, ipf)
            end do
            
            ! Update complex EAGP Pfaffian matrix elements
            do i = 1, nnozero_eagp
                ind = i + n
                iy = (nozero(ind) - 1)/ndiff + 1
                ix = nozero(ind) - (iy - 1)*ndiff
                j = jbradet(ind)
                if (j .gt. 0) then
                    ! Add complex values maintaining antisymmetry
                    eagp_pfaff(2*ix - 1, iy) = eagp_pfaff(2*ix - 1, iy) + psip(2*ind - 1)
                    eagp_pfaff(2*ix, iy) = eagp_pfaff(2*ix, iy) + psip(2*ind)
                    eagp_pfaff(2*iy - 1, ix) = -eagp_pfaff(2*ix - 1, iy)  ! Antisymmetric
                    eagp_pfaff(2*iy, ix) = -eagp_pfaff(2*ix, iy)          ! Antisymmetric
                elseif (j .lt. 0) then
                    ! Subtract complex values maintaining antisymmetry
                    eagp_pfaff(2*ix - 1, iy) = eagp_pfaff(2*ix - 1, iy) - psip(2*ind - 1)
                    eagp_pfaff(2*ix, iy) = eagp_pfaff(2*ix, iy) - psip(2*ind)
                    eagp_pfaff(2*iy - 1, ix) = -eagp_pfaff(2*ix - 1, iy)  ! Antisymmetric
                    eagp_pfaff(2*iy, ix) = -eagp_pfaff(2*ix, iy)          ! Antisymmetric
                end if
            end do

        ! Set mode for complex wave functions
        else
            ! Set complex determinant matrix elements
            do i = 1, n
                if (abs(jbradet(i)) .ne. 0)&
                        & call upsim_complex(derl, nelorb, nozero(i), psip(2*i - 1), symmagp, ipf)
            end do

            ! Set complex EAGP Pfaffian matrix elements
            do i = 1, nnozero_eagp
                ind = i + n
                iy = (nozero(ind) - 1)/ndiff + 1
                ix = nozero(ind) - (iy - 1)*ndiff
                j = jbradet(ind)
                if (j .gt. 0) then
                    ! Set complex values maintaining antisymmetry
                    eagp_pfaff(2*ix - 1, iy) = psip(2*ind - 1)
                    eagp_pfaff(2*ix, iy) = psip(2*ind)
                    eagp_pfaff(2*iy - 1, ix) = -eagp_pfaff(2*ix - 1, iy)  ! Antisymmetric
                    eagp_pfaff(2*iy, ix) = -eagp_pfaff(2*ix, iy)          ! Antisymmetric
                elseif (j .lt. 0) then
                    ! Set sign-flipped complex values maintaining antisymmetry
                    eagp_pfaff(2*ix - 1, iy) = -psip(2*ind - 1)
                    eagp_pfaff(2*ix, iy) = -psip(2*ind)
                    eagp_pfaff(2*iy - 1, ix) = -eagp_pfaff(2*ix - 1, iy)  ! Antisymmetric
                    eagp_pfaff(2*iy, ix) = -eagp_pfaff(2*ix, iy)          ! Antisymmetric
                end if
            end do
        end if ! endif yesupdate

    end if ! endif ipc=1

    return
end
