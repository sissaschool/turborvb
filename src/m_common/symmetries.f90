! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2010-2011 Quantum ESPRESSO group
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

module symmetries

    ! isymm contains all symmetry matrices in crystal lattice
    ! nrot is the total number of Bravais lattice symmetries
    ! sname contains the names of the symmetries
    ! allowed_sym check whether the symmetries of BL are allowed by the point group of the supercell
    integer :: nrot, nrot_full
    integer :: isymm(3, 3, 48)
    character :: sname(64)*45

    ! full name of the rotational part of each symmetry operation
    character :: s0name(64)*45

    data s0name/ &
        'identity                                     ', &
        '180 deg rotation - cart. axis [0,0,1]        ', &
        '180 deg rotation - cart. axis [0,1,0]        ', &
        '180 deg rotation - cart. axis [1,0,0]        ', &
        '180 deg rotation - cart. axis [1,1,0]        ', &
        '180 deg rotation - cart. axis [1,-1,0]       ', &
        ' 90 deg rotation - cart. axis [0,0,-1]       ', &
        ' 90 deg rotation - cart. axis [0,0,1]        ', &
        '180 deg rotation - cart. axis [1,0,1]        ', &
        '180 deg rotation - cart. axis [-1,0,1]       ', &
        ' 90 deg rotation - cart. axis [0,1,0]        ', &
        ' 90 deg rotation - cart. axis [0,-1,0]       ', &
        '180 deg rotation - cart. axis [0,1,1]        ', &
        '180 deg rotation - cart. axis [0,1,-1]       ', &
        ' 90 deg rotation - cart. axis [-1,0,0]       ', &
        ' 90 deg rotation - cart. axis [1,0,0]        ', &
        '120 deg rotation - cart. axis [-1,-1,-1]     ', &
        '120 deg rotation - cart. axis [-1,1,1]       ', &
        '120 deg rotation - cart. axis [1,1,-1]       ', &
        '120 deg rotation - cart. axis [1,-1,1]       ', &
        '120 deg rotation - cart. axis [1,1,1]        ', &
        '120 deg rotation - cart. axis [-1,1,-1]      ', &
        '120 deg rotation - cart. axis [1,-1,-1]      ', &
        '120 deg rotation - cart. axis [-1,-1,1]      ', &
        ' 60 deg rotation - cryst. axis [0,0,1]       ', &
        ' 60 deg rotation - cryst. axis [0,0,-1]      ', &
        '120 deg rotation - cryst. axis [0,0,1]       ', &
        '120 deg rotation - cryst. axis [0,0,-1]      ', &
        '180 deg rotation - cryst. axis [1,-1,0]      ', &
        '180 deg rotation - cryst. axis [2,1,0]       ', &
        '180 deg rotation - cryst. axis [0,1,0]       ', &
        '180 deg rotation - cryst. axis [1,1,0]       ', &
        'inversion                                    ', &
        'inv. 180 deg rotation - cart. axis [0,0,1]   ', &
        'inv. 180 deg rotation - cart. axis [0,1,0]   ', &
        'inv. 180 deg rotation - cart. axis [1,0,0]   ', &
        'inv. 180 deg rotation - cart. axis [1,1,0]   ', &
        'inv. 180 deg rotation - cart. axis [1,-1,0]  ', &
        'inv.  90 deg rotation - cart. axis [0,0,-1]  ', &
        'inv.  90 deg rotation - cart. axis [0,0,1]   ', &
        'inv. 180 deg rotation - cart. axis [1,0,1]   ', &
        'inv. 180 deg rotation - cart. axis [-1,0,1]  ', &
        'inv.  90 deg rotation - cart. axis [0,1,0]   ', &
        'inv.  90 deg rotation - cart. axis [0,-1,0]  ', &
        'inv. 180 deg rotation - cart. axis [0,1,1]   ', &
        'inv. 180 deg rotation - cart. axis [0,1,-1]  ', &
        'inv.  90 deg rotation - cart. axis [-1,0,0]  ', &
        'inv.  90 deg rotation - cart. axis [1,0,0]   ', &
        'inv. 120 deg rotation - cart. axis [-1,-1,-1]', &
        'inv. 120 deg rotation - cart. axis [-1,1,1]  ', &
        'inv. 120 deg rotation - cart. axis [1,1,-1]  ', &
        'inv. 120 deg rotation - cart. axis [1,-1,1]  ', &
        'inv. 120 deg rotation - cart. axis [1,1,1]   ', &
        'inv. 120 deg rotation - cart. axis [-1,1,-1] ', &
        'inv. 120 deg rotation - cart. axis [1,-1,-1] ', &
        'inv. 120 deg rotation - cart. axis [-1,-1,1] ', &
        'inv.  60 deg rotation - cryst. axis [0,0,1]  ', &
        'inv.  60 deg rotation - cryst. axis [0,0,-1] ', &
        'inv. 120 deg rotation - cryst. axis [0,0,1]  ', &
        'inv. 120 deg rotation - cryst. axis [0,0,-1] ', &
        'inv. 180 deg rotation - cryst. axis [1,-1,0] ', &
        'inv. 180 deg rotation - cryst. axis [2,1,0]  ', &
        'inv. 180 deg rotation - cryst. axis [0,1,0]  ', &
        'inv. 180 deg rotation - cryst. axis [1,1,0]  '/

    public :: isymm, nrot, sname, nrot_full
    public :: set_sym_bl, purge_isymm, transform_point, invmat
    private :: s0name

contains

    !   subroutine find_symmetries( nion, rion, atom_number, iespbc, nx, ny, nz )
    !
    !     implicit none
    !
    !     integer :: nion
    !     logical :: iespbc
    !     real(8) :: at(3,3),rion(3,nion),atom_number(nion),cellscale(3)
    !     integer,optional :: nx,ny,nz
    !
    !     integer i,j
    !
    !     allowed_sym(1:48) = .true.
    !
    !     call set_sym_bl()
    !     nrot_full = nrot
    !     if( present(nx) .and. present(ny) .and. present(nz) ) then
    !        call set_sym_cell( nion, rion, atom_number, iespbc, nx, ny, nz )
    !     else
    !        call set_sym_cell( nion, rion, atom_number, iespbc )
    !     endif
    !
    !     return
    !
    !   end subroutine find_symmetries

    ! Copyright (C) 2010-2011 Quantum ESPRESSO group
    ! This file is distributed under the terms of the
    ! GNU General Public License. See the file `License'
    ! in the root directory of the present distribution,
    ! or http://www.gnu.org/copyleft/gpl.txt .

    subroutine set_sym_bl(at)

        !
        ! Provides symmetry operations for all bravais lattices
        ! Tests first the 24 proper rotations for the cubic lattice;
        ! then the 8 rotations specific for the hexagonal axis (special axis c);
        ! then inversion is added
        !

        implicit none
        ! sin3 = sin(pi/3), cos3 = cos(pi/3), msin3 = -sin(pi/3), mcos3 = -cos(pi/3)
        !
        real(8), dimension(3, 3), intent(in) :: at
        real(8), parameter :: sin3 = 0.866025403784438597d0, cos3 = 0.5d0, &
                              msin3 = -0.866025403784438597d0, mcos3 = -0.5d0
        real(8) :: s0(3, 3, 32), overlap(3, 3), rat(3), rot(3, 3), value
        ! s0: the s matrices in cartesian axis
        ! overlap: inverse overlap matrix between direct lattice
        ! rat: the rotated of a direct vector ( cartesian )
        ! rot: the rotated of a direct vector ( crystal axis )
        ! value: component of the s matrix in axis basis
        integer :: jpol, kpol, mpol, irot, i
        ! counters over the polarizations and the rotations

        data s0/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, &
            -1.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, 0.d0, 0.d0, 1.d0, &
            -1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, -1.d0, &
            1.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, 0.d0, 0.d0, -1.d0, &
            0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, -1.d0, &
            0.d0, -1.d0, 0.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, -1.d0, &
            0.d0, -1.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, &
            0.d0, 1.d0, 0.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, &
            0.d0, 0.d0, 1.d0, 0.d0, -1.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
            0.d0, 0.d0, -1.d0, 0.d0, -1.d0, 0.d0, -1.d0, 0.d0, 0.d0, &
            0.d0, 0.d0, -1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
            0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, -1.d0, 0.d0, 0.d0, &
            -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, &
            -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, -1.d0, 0.d0, &
            1.d0, 0.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, 1.d0, 0.d0, &
            1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, -1.d0, 0.d0, &
            0.d0, 0.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, &
            0.d0, 0.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, &
            0.d0, 0.d0, -1.d0, 1.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, &
            0.d0, 0.d0, 1.d0, -1.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, &
            0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 1.d0, 0.d0, 0.d0, &
            0.d0, -1.d0, 0.d0, 0.d0, 0.d0, -1.d0, 1.d0, 0.d0, 0.d0, &
            0.d0, -1.d0, 0.d0, 0.d0, 0.d0, 1.d0, -1.d0, 0.d0, 0.d0, &
            0.d0, 1.d0, 0.d0, 0.d0, 0.d0, -1.d0, -1.d0, 0.d0, 0.d0, &
            cos3, sin3, 0.d0, msin3, cos3, 0.d0, 0.d0, 0.d0, 1.d0, &
            cos3, msin3, 0.d0, sin3, cos3, 0.d0, 0.d0, 0.d0, 1.d0, &
            mcos3, sin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0, 1.d0, &
            mcos3, msin3, 0.d0, sin3, mcos3, 0.d0, 0.d0, 0.d0, 1.d0, &
            cos3, msin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
            cos3, sin3, 0.d0, sin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
            mcos3, msin3, 0.d0, msin3, cos3, 0.d0, 0.d0, 0.d0, -1.d0, &
            mcos3, sin3, 0.d0, sin3, cos3, 0.d0, 0.d0, 0.d0, -1.d0/

        !
        ! the routine computes the BL symmetries using
        ! the lattice vectors "s2r(3,3)". Since so far the cell
        ! is always orthorombic, the versors "at(3,3)" are always the same.
        !
        !    compute the overlap matrix for crystal axis
        do jpol = 1, 3
            do kpol = 1, 3
                rot(kpol, jpol) = at(1, kpol)*at(1, jpol) + &
                                  at(2, kpol)*at(2, jpol) + &
                                  at(3, kpol)*at(3, jpol)
            end do
        end do
        !
        !    then its inverse (rot is used as work space)
        !
        call invmat(3, rot, overlap, value)
        nrot = 1
        do irot = 1, 32
            !
            !   for each possible symmetry
            !
            do jpol = 1, 3
                do mpol = 1, 3
                    !
                    !   compute, in cartesian coordinates the rotated vector
                    !
                    rat(mpol) = s0(mpol, 1, irot)*at(1, jpol) + &
                                s0(mpol, 2, irot)*at(2, jpol) + &
                                s0(mpol, 3, irot)*at(3, jpol)
                end do

                do kpol = 1, 3
                    !
                    !   the rotated vector is projected on the direct lattice
                    !
                    rot(kpol, jpol) = at(1, kpol)*rat(1) + &
                                      at(2, kpol)*rat(2) + &
                                      at(3, kpol)*rat(3)
                end do
            end do
            !
            !  and the inverse of the overlap matrix is applied
            !
            do jpol = 1, 3
                do kpol = 1, 3
                    value = overlap(jpol, 1)*rot(1, kpol) + &
                            overlap(jpol, 2)*rot(2, kpol) + &
                            overlap(jpol, 3)*rot(3, kpol)
                    if (abs(dble(nint(value)) - value) > 1.0d-6) then
                        !
                        ! if a noninteger is obtained, this implies that this operation
                        ! is not a symmetry operation for the given lattice
                        !
                        go to 10
                    end if
                    isymm(kpol, jpol, nrot) = nint(value)
                    sname(nrot) = s0name(irot)
                end do
            end do
            nrot = nrot + 1
10          continue
        end do
        nrot = nrot - 1
        !
        !     set the inversion symmetry ( Bravais lattices have always inversion
        !     symmetry )
        !
        do irot = 1, nrot
            do kpol = 1, 3
                do jpol = 1, 3
                    isymm(kpol, jpol, irot + nrot) = -isymm(kpol, jpol, irot)
                    sname(irot + nrot) = s0name(irot + 32)
                end do
            end do
        end do

        nrot = 2*nrot

        return
        !
    end subroutine set_sym_bl

    !   subroutine set_sym_cell ( nion, rion, atom_number, iespbc, nx, ny, nz )
    !
    !     implicit none
    !     !
    !     ! input variables
    !     !
    !     integer,intent(in) :: nion
    !     real(8),intent(in) :: rion(3,nion),atom_number(nion)
    !     logical,intent(in) :: iespbc
    !     integer,optional   :: nx,ny,nz
    !     !
    !     ! local variables
    !     !
    !     integer :: irot,nrot_,na,nb,ind,i,j,isymm_(3,3,48),ierr
    !     logical,dimension(:,:),allocatable :: sym_map
    !     logical :: found_sym
    !     real(8) :: pos(3), pos_rot(3), diff(3)
    !     character :: sname_(64)*45
    !
    !     allowed_sym(1:nrot) = .false.
    !
    !     allocate(sym_map(nion,nrot))
    !     sym_map(:,:) = .false.
    !
    !     do na = 1,nion
    !        do irot = 1,nrot
    !           pos_rot(:) = matmul( isymm(:,:,irot),rion(:,na) )
    !           found_sym = .true.
    !           nb = 1
    !           do while( nb.le.nion .and. found_sym )
    !              pos(:) = rion(:,nb)
    !              diff(:) = pos_rot(:) - pos(:)
    !              if(iespbc) call applyPBC(diff,1)
    !              if( (sum(abs(diff(:))) .lt. 1.d-5) .and. &
    !                   (nint(atom_number(na)) .eq. nint(atom_number(nb)))  ) then
    !                 sym_map(na,irot) = .true.
    !                 found_sym  = .false.
    !              endif
    !              nb = nb + 1
    !           enddo
    !        enddo
    !     enddo
    !     ! identity
    !     do irot = 1,nrot
    !        sym_map(1,irot) = .true.
    !     enddo
    !     !
    !     ! fill vector of symmetries allowed by the supercell,
    !     ! in general less than the ones allowed by the Bravais lattice.
    !     !
    !     do irot = 1,nrot
    !        if( all(sym_map(:,irot)) ) allowed_sym(irot) = .true.
    !     enddo
    !     !
    !     ! check that the grid is compatible with the S rotation
    !     !
    !     if( present(nx) .and. present(ny) .and. present(nz) ) then
    !        do irot = 1, nrot
    !           if ( mod (isymm(2, 1, irot) * nx, ny) /= 0 .or. &
    !                mod (isymm(3, 1, irot) * nx, nz) /= 0 .or. &
    !                mod (isymm(1, 2, irot) * ny, nx) /= 0 .or. &
    !                mod (isymm(3, 2, irot) * ny, nz) /= 0 .or. &
    !                mod (isymm(1, 3, irot) * nz, nx) /= 0 .or. &
    !                mod (isymm(2, 3, irot) * nz, ny) /= 0 ) then
    !              allowed_sym(irot) = .false.
    !           endif
    !        enddo
    !     endif
    !     !
    !     ! fill symmetry vector isymm with new allowed symmetries
    !     ! update the number of allowed symmetries
    !     !
    !     nrot_ = nrot
    !     nrot  = 0
    !     do irot = 1,nrot_
    !        if(allowed_sym(irot)) then
    !           nrot = nrot + 1
    !           isymm_(:,:,nrot) = isymm(:,:,irot)
    !           sname_(nrot) = sname(irot)
    !        endif
    !     enddo
    !     isymm(:,:,:) = 0
    !     isymm = isymm_
    !     sname(:) = " "
    !     sname(:) = sname_(:)
    !
    !     deallocate(sym_map)
    !
    !     return
    !
    !   end subroutine set_sym_cell

    subroutine purge_isymm(iespbc, nsym, isymm, sname, t_rev, rion, nion, zeta, cellscale, nx, ny, nz)

        use constants, only: deps
        implicit none

        integer :: j, nion, nsym, nsym_acc, is, i2, npip(3), isymm(3, 3, *), t_rev(*)
        integer, optional :: nx, ny, nz
        logical :: found, accept, iespbc
        real(8) :: rionmap(3), rion(3, nion), dist(3), zeta(*), cellscale(3)
        character :: sname(64)*45
        logical, dimension(:), allocatable :: yessite

        npip = 0
        allocate (yessite(nion))
        nsym_acc = 0

        nrot_full = nrot

        do is = 1, nrot_full
            accept = .true.
            yessite = .false.
            do j = 1, nion
                rionmap(1:3) = matmul(isymm(:, :, is), rion(:, j))
                i2 = 1
                found = .true.
                do while (i2 .le. nion .and. found)
                    dist(1:3) = rionmap(1:3) - rion(1:3, i2)
                    if (iespbc) then
                        call makeimage(dist, cellscale, deps)
                    end if
                    if (sum(abs(dist(:))) .lt. deps .and. nint(zeta(i2)) .eq. nint(zeta(j))) then
                        found = .false.
                        yessite(i2) = .true.
                    end if
                    i2 = i2 + 1
                end do
            end do
            if (found) then
                accept = .false.
            else
                accept = .true.
                do j = 1, nion
                    if (.not. yessite(j)) accept = .false.
                end do
            end if
            !
            ! check that the grid is compatible with the S rotation
            !
            if (present(nx) .and. present(ny) .and. present(nz)) then
                if (mod(isymm(2, 1, is)*nx, ny) /= 0 .or. &
                    mod(isymm(3, 1, is)*nx, nz) /= 0 .or. &
                    mod(isymm(1, 2, is)*ny, nx) /= 0 .or. &
                    mod(isymm(3, 2, is)*ny, nz) /= 0 .or. &
                    mod(isymm(1, 3, is)*nz, nx) /= 0 .or. &
                    mod(isymm(2, 3, is)*nz, ny) /= 0) then
                    accept = .false.
                end if
            end if

            if (accept) then
                nsym_acc = nsym_acc + 1
                isymm(:, :, nsym_acc) = isymm(:, :, is)
                sname(nsym_acc) = sname(is)
                t_rev(nsym_acc) = t_rev(is)
            end if
        end do

        nrot = nsym_acc
        deallocate (yessite)

        return

    end subroutine purge_isymm

    subroutine transform_point(s, i, j, k, nx, ny, nz, ri, rj, rk)

        !    This routine computes the rotated of the point i,j,k throught
        !    the symmetry (s,f). Then it computes the equivalent point
        !    on the original real space mesh
        !    This routine is used in the DFT code (routine updenorb_new) for
        !    charge density symmetrisation in case of k-points sampling.

        implicit none
        !
        !    first the dummy variables
        !
        integer, intent(in) :: s(3, 3), i, j, k, nx, ny, nz
        integer, intent(inout) :: ri, rj, rk
        ! input: the rotation matrix
        ! input: the fractionary translation
        ! \
        !   input: the point to rotate
        ! /
        ! \
        !   input: the dimension of the mesh
        ! /
        !\
        !  output: the rotated point
        !/
        !
        !  local variable
        !
        ! the rotation matrix in scaled crystallographic coordinates

        ri = s(1, 1)*(i - 1) + s(2, 1)*(j - 1) + s(3, 1)*(k - 1)
        ri = mod(ri, nx) + 1
        if (ri .lt. 1) ri = ri + nx
        rj = s(1, 2)*(i - 1) + s(2, 2)*(j - 1) + s(3, 2)*(k - 1)
        rj = mod(rj, ny) + 1
        if (rj .lt. 1) rj = rj + ny
        rk = s(1, 3)*(i - 1) + s(2, 3)*(j - 1) + s(3, 3)*(k - 1)
        rk = mod(rk, nz) + 1
        if (rk .lt. 1) rk = rk + nz

        return

    end subroutine transform_point

    ! adapted from QuantumESPRESSO
    subroutine invmat(n, a, a_inv, da)
        !
        ! computes the inverse "a_inv" of matrix "a", both dimensioned (n,n)
        ! matrix "a" is unchanged on output - LAPACK
        !
        implicit none
        integer :: n
        real(8), dimension(n, n) :: a, a_inv
        real(8) :: da
        integer :: info, lda, lwork, ipiv(n)
        ! info=0: inversion was successful
        ! lda   : leading dimension (the same as n)
        ! ipiv  : work space for pivoting (assumed of length lwork=n)
        real(8) :: work(n)
        ! more work space
        lda = n
        lwork = n
        a_inv(:, :) = a(:, :)
        !
        call DGETRF(n, n, a_inv, lda, ipiv, info)
        if (info .ne. 0) call errore('invmat', 'error in DGETRF', abs(info))
        call DGETRI(n, a_inv, lda, ipiv, work, lwork, info)
        if (info .ne. 0) call errore('invmat', 'error in DGETRI', abs(info))
        !
        if (n == 3) then
            da = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) + &
                 a(1, 2)*(a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3)) + &
                 a(1, 3)*(a(2, 1)*a(3, 2) - a(3, 1)*a(2, 2))
            if (abs(da) < 1.d-10) call errore(' invmat ', ' singular matrix ', 1)
        else
            da = 0.d0
        end if
        return
    end subroutine invmat

end module symmetries
