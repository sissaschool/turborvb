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

program real_to_complex

    ! This program take in input a real wave function fort.10 and
    ! build the correspondent complex fort.10_complex setting to zero (or
    ! to a random number) the imaginary parts of the complex coefficients.

    use allio

    implicit none
    integer :: nshell_o, indpar, indpar_read, is, i, j, ufort10, ndetmatR, &
               ndupcR, cfl, nnozeroR, ind, nelorb_cnew, ix1, iy1, ii, count_mol
    real(8), dimension(:), allocatable :: detmatRead, dupcRead
    real*8 detmat_sav
    integer, dimension(:), allocatable :: ioccup_cread, ioptorb_cread, mult_cread, nparam_cread, nozero_cnew, kion_cread
    character(30) :: filling
    logical :: molyesRead, match
    !
    character(100) name_tool
    character(20) str
    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then
        ! input the name of the file exactly as it is in /doc
        name_tool = 'real_to_complex'
        call help_online(name_tool)
        stop
    end if
    !
    ! read the input wave function in fort.10
    !
    ufort10 = 10
    call default_allocate
    open (unit=ufort10, file='fort.10', status='old', form='formatted', err=101)
    call read_fort10(10)
    if (yes_complex) then
        go to 102
    else
        yes_complex = .true.
        forcecomplex = .true.
    end if
    !
    ! set basis set coefficients complex where needed
    !
    ndupcR = iesupr_c
    allocate (dupcRead(ndupcR))
    dupcRead(1:ndupcR) = dup_c(1:ndupcR)
    deallocate (dup_c)
    nshell_o = nshell_c
    if (molecular .ne. 0 .and. symmagp) then
        nmol = molecular - ndiff
        iesup_c = iesup_c + 2*nmol*ipf*nelorbh
        allocate (ioccup_cread(occ_c + nmol), nparam_cread(nshell_c), &
                  mult_cread(nshell_c), ioptorb_cread(nshell_c), kion_cread(nshell_c))
        nparam_cread = nparam_c
        ioccup_cread = ioccup_c
        mult_cread = mult_c
        ioptorb_cread = ioptorb_c
        kion_cread = kion_c
        nshell_c = nshell_c + nmol

        deallocate (mult_c, nparam_c, ioptorb_c, ioccup_c, kion_c)
        occ_c = occ_c + nmol
        allocate (ioccup_c(occ_c), mult_c(nshell_c), nparam_c(nshell_c), ioptorb_c(nshell_c), kion_c(nshell_c))
        ioccup_c(1:occ_c - nmol) = ioccup_cread(1:occ_c - nmol)
        ioccup_c(occ_c - nmol + 1:occ_c) = 1
        mult_c(1:nshell_o) = mult_cread(1:nshell_o)
        nparam_c(1:nshell_o) = nparam_cread(1:nshell_o)
        ioptorb_c(1:nshell_o) = ioptorb_cread(1:nshell_o)
        kion_c(1:nshell_o) = kion_cread(1:nshell_o)
        kion_c(nshell_o + 1:nshell_c) = 1
        deallocate (ioptorb_cread, mult_cread, nparam_cread, ioccup_cread, kion_cread)
        nelorb_cnew = nelorb_c + nmol
        allocate (dup_c(2*ndupcR + (nmol + ndiff)*ipf*nelorbh*4), nozero_cnew(nnozero_c))

        do i = 1, nnozero_c - molecular
            iy = (nozero_c(i) - 1)/nelorb_c + 1
            ix = nozero_c(i) - (iy - 1)*nelorb_c
            if (iy .gt. nelorb_c) iy = iy + nmol
            nozero_cnew(i) = (iy - 1)*nelorb_cnew + ix
        end do
        ! First molecular orbital ME
        iy1 = (nozero_c(nnozero_c - molecular + 1) - 1)/nelorb_c + 1
        ix1 = nozero_c(nnozero_c - molecular + 1) - (iy1 - 1)*nelorb_c
        indnn = 1
        match = .true.
        do while (match)
            ii = jbradetn(indnn)
            if (ii .ne. -1) then
                indnn = indnn - 2*ii + 1
            else
                ix = jbradetn(indnn + 1)
                iy = jbradetn(indnn + 2)
                if (ix .eq. ix1 .and. iy .eq. iy1) then
                    match = .false.
                else
                    indnn = indnn + 3
                end if
            end if
        end do
        indnn = indnn - 1

        do i = nshell_o + 1, nshell_c
            mult_c(i) = 1
            ioptorb_c(i) = 1000000
            nparam_c(i) = nparam_c(nshell_o)
        end do
        ind = nnozero_c - molecular
        do i = 1, nmol
            ix = nelorb_cnew - 2*nmol - ndiff + 2*i
            iy = ix - 1
            ind = ind + 1
            nozero_cnew(ind) = nelorb_cnew*(iy - 1) + ix
            indnn = indnn + 1
            jbradetn(indnn) = -1
            jbradetn(indnn + 1) = ix
            jbradetn(indnn + 2) = iy
            indnn = indnn + 2
        end do
        do i = nmol + 1, nmol + ndiff
            ind = ind + 1
            j = nelorb_cnew - ndiff - nmol + i
            nozero_cnew(ind) = nelorb_cnew*(nelorb_cnew + i - 1 - nmol) + j
            indnn = indnn + 1
            jbradetn(indnn) = -1
            jbradetn(indnn + 1) = j
            jbradetn(indnn + 2) = nelorb_cnew + i - nmol
            indnn = indnn + 2
        end do

    else
        nelorb_cnew = nelorb_c
        allocate (dup_c(2*ndupcR))
    end if
    count_mol = 0
    dup_c = 0.d0
    indpar = 0
    indpar_read = 0
    !  molyes is defined only in read_fort10_fast and not in read_fort10
    if (molecular .gt. 0) then
        molyesread = .true.
    else
        molyesread = .false.
    end if
    molyes = .false.
    do is = 1, nshell_o
        if (ioptorb_c(is) .lt. 900000) then
            if (nparam_c(is) .eq. 1) then
                dup_c(2*indpar + 1) = dupcRead(indpar_read + 1)
                dup_c(2*indpar + 2) = 0.d0
            else
                do j = 1, 2*nparam_c(is), 2
                    dup_c(2*indpar + j) = dupcRead(indpar_read + j/2 + 1)
                    dup_c(2*indpar + j + 1) = 0.d0
                end do
            end if
        elseif (ioptorb_c(is) .eq. 900000) then
            do j = 1, nparam_c(is), 2
                dup_c(2*indpar + j) = dupcRead(indpar_read + j/2 + 1)
                dup_c(2*indpar + j + 1) = 0.d0
            end do
            do j = nparam_c(is) + 1, 2*nparam_c(is), 2
                dup_c(2*indpar + j) = dupcRead(indpar_read + j/2 + 1)
                dup_c(2*indpar + j + 1) = 0.d0
            end do
        elseif (ioptorb_c(is) .eq. 1000000) then
            count_mol = count_mol + 1
            do j = 1, nparam_c(is), 2
                dup_c(2*indpar + j) = dupcRead(indpar_read + j/2 + 1)
                dup_c(2*indpar + j + 1) = 0.d0
            end do
            do j = nparam_c(is) + 1, 2*nparam_c(is), 2
                dup_c(2*indpar + j) = dupcRead(indpar_read + j/2 + 1)
                dup_c(2*indpar + j + 1) = 0.d0
            end do
            if (symmagp .and. count_mol .le. nmol) then
                indpar = indpar + nparam_c(is)
                do j = 1, nparam_c(is), 2
                    dup_c(2*indpar + j) = dupcRead(indpar_read + j/2 + 1)
                    dup_c(2*indpar + j + 1) = 0.d0
                end do
                do j = nparam_c(is) + 1, 2*nparam_c(is), 2
                    dup_c(2*indpar + j) = dupcRead(indpar_read + j/2 + 1)
                    dup_c(2*indpar + j + 1) = 0.d0
                end do
            end if
        end if
        ! check if shells are correctly ordered
        if (molyes .and. ioptorb_c(is) .lt. 900000) go to 103
        if (ioptorb_c(is) .ge. 1000000) molyes = .true.
        indpar = indpar + nparam_c(is)
        indpar_read = indpar_read + nparam_c(is)
    end do

    if (molyesRead .neqv. molyes) go to 104
    !
    ! transform geminal matrix to complex

    !
    call set_new_detmat(contraction)
    do is = 1, nnozeroR
        if (contraction .ne. 0) then
            if (molecular .ne. 0 .and. symmagp) then
                if (2*nozero_cnew(is) .gt. size(detmat_c)) then
                    write (6, *) ' ERROR ', 2*nozero_cnew(is), size(detmat_c), 2*ndetmatR
                    stop
                end if
                detmat_c(2*nozero_cnew(is) - 1) = detmatRead(nozero_c(is))
                detmat_c(2*nozero_cnew(is)) = 0.d0
            else
                detmat_c(2*nozero_c(is) - 1) = detmatRead(nozero_c(is))
                detmat_c(2*nozero_c(is)) = 0.d0
            end if
        else
            detmat(2*nozero(is) - 1) = detmatRead(nozero(is))
            detmat(2*nozero(is)) = 0.d0
        end if
    end do
    if (molecular .ne. 0 .and. symmagp) then
        nozero_c = nozero_cnew
        nelorb_c = nelorb_cnew
        deallocate (nozero_cnew)
    end if

    do i = nnozero_c - molecular + 1, nnozero_c
        iy = (nozero_c(i) - 1)/nelorb_c + 1
        ix = nozero_c(i) - (iy - 1)*nelorb_c
    end do
#if defined DEBUG
    write (6, *) ' Check basis set '
    do i = 1, iesupr_c
        write (6, *) dupcRead(i), dup_c(2*(i - 1) + 1), dup_c(2*(i - 1) + 2)
    end do
    write (6, *) ' Check detmat '
    if (contraction .ne. 0) then
        do i = 1, nnozero_c
            write (6, *) detmatRead(nozero_c(i)), detmat_c(2*nozero_c(i) - 1), detmat_c(2*nozero_c(i))
        end do
    else
        do i = 1, nnozero
            write (6, *) detmatRead(nozero(i)), detmat(2*nozero(i) - 1), detmat(2*nozero(i))
        end do
    end if
#endif
    !
    ! write the output complex wave function in fort.10_complex
    !
    close (ufort10)
    open (unit=ufort10, file='fort.10_complex', status='unknown', form='formatted', err=105)
    nelcol_c = nelorb_cnew
    ipc = 2
    call write_fort10(ufort10)
    close (ufort10)

    deallocate (dupcRead, detmatRead)
    call deallocate_all

    stop

    !!!!!!! ERRORS !!!!!!!

101 write (6, *) ' Input fort.10 not found!!'
    stop
102 write (6, *) ' Input fort.10 is already complex!!'
    stop
103 write (6, *) ' You cannot use atomic orbital after molecular!!'
    stop
104 write (6, *) ' Some error in reading determinant shells!!'
    stop
105 write (6, *) ' Error in writing output fort.10_complex!!'
    stop

contains

    subroutine set_new_detmat(contraction)
        implicit none
        integer, intent(in) :: contraction
        if (contraction .ne. 0) then
            ndetmatR = nelorb_c*max(nelcol_c, nel)
            allocate (detmatRead(ndetmatR))
            detmatRead(:) = detmat_c(:)
            deallocate (detmat_c)
            if (molecular .ne. 0 .and. symmagp) then
                ndetmatR = nelorb_cnew*max(nelorb_cnew + ndiff, nel)
            end if
            allocate (detmat_c(2*ndetmatR))
            detmat_c = 0.d0
            nnozeroR = nnozero_c
        else
            ndetmatR = ipf*nelorbh*max(nelcolh, nel)
            allocate (detmatRead(ndetmatR))
            detmatRead(:) = detmat(:)
            deallocate (detmat)
            allocate (detmat(2*ndetmatR))
            detmat = 0.d0
            nnozeroR = nnozero
        end if
        return
    end subroutine set_new_detmat

end program real_to_complex
