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

module exact_diagonalization

    !use constants
    implicit none

    integer :: param_number, points_number
    real(8), dimension(:), allocatable :: norm_ed
    real(8), dimension(:), allocatable :: param0_ed
    real(8), dimension(:, :), allocatable :: points_ed
    logical, dimension(:), allocatable :: sel_ed
    real(8), external :: pseudofun

    integer :: bs_number
    real(8), dimension(:), allocatable :: bs_coeff
    real(8), dimension(:), allocatable :: dx, vecscratch
    real(8), dimension(:, :), allocatable :: overlap
    real(8), dimension(:), allocatable :: proj
    integer, dimension(:), allocatable :: point2param

    ! INPUT ###################################################

    ! PARAMETERS
    ! param_number: number of parameters
    ! param0_ed: initial set of parameters
    ! norm_ed: normalization of basis set functions
    ! sel_ed: constraint on the parameter freedom
    ! true --> free to move
    ! false --> fixed

    ! MAIN FUNCTION
    ! pseudofun defined elsewhere

    ! DATA SET
    ! points_number: number of data points to fit
    ! points_ed: data to fit
    ! points_ed(1,i): x_i
    ! points_ed(2,i): y_i

    private
    public :: param0_ed, points_ed, sel_ed, &
              alloc_datavec_ed, dealloc_datavec_ed, matrix_setup, dealloc_matrix_setup, &
              overlap_matrix, update_param_ed, orthonormalization, norm_ed

contains

    ! basis set functions
    ! reference

    ! compute the overlap matrix

    ! orthonormalization of the basis set

    subroutine alloc_datavec_ed(param_number_read, points_number_read)
        integer, intent(in) :: param_number_read, points_number_read
        param_number = param_number_read
        points_number = points_number_read
        allocate (param0_ed(param_number), sel_ed(param_number))
        allocate (points_ed(2, points_number))
        allocate (norm_ed(param_number/3))
    end subroutine alloc_datavec_ed

    subroutine dealloc_datavec_ed
        if (allocated(param0_ed)) deallocate (param0_ed)
        if (allocated(norm_ed)) deallocate (norm_ed)
        if (allocated(sel_ed)) deallocate (sel_ed)
        if (allocated(points_ed)) deallocate (points_ed)
    end subroutine dealloc_datavec_ed

    subroutine dealloc_matrix_setup
        if (allocated(bs_coeff)) deallocate (bs_coeff)
        if (allocated(point2param)) deallocate (point2param)
        if (allocated(dx)) deallocate (dx)
        if (allocated(vecscratch)) deallocate (vecscratch)
        if (allocated(overlap)) deallocate (overlap)
        if (allocated(proj)) deallocate (proj)
    end subroutine dealloc_matrix_setup

    subroutine matrix_setup(weights)

        implicit none
        integer :: i, j
        real(8) :: p(3), norm
        real(8), intent(in) :: weights(*)

        ! count the number of basis set elements
        j = 0
        do i = 1, param_number
            if (sel_ed(i)) then
                j = j + 1
                if (mod(i, 3) .ne. 1) then
                    write (6, *) 'you are not allowed to change the power or the gaussian exponent by ED'
                    write (6, *) 'please change the sel_edection options'
                    stop
                end if
            end if
        end do
        bs_number = j

        allocate (bs_coeff(bs_number), point2param(bs_number))

        ! initialize the linear combination coefficients
        j = 0
        do i = 1, param_number
            if (sel_ed(i)) then
                j = j + 1
                bs_coeff(j) = param0_ed(i)
                point2param(j) = i
            end if
        end do

        ! set integration grid dx for the interval [0,points_ed(points_number)]
        ! DEFINE THE WEIGHTED METRICS
        allocate (dx(points_number))

        ! old integration
        !    dx(1)=(points_ed(1,2)-points_ed(1,1))/2.d0+points_ed(1,1)
        !    dx(1)=(points_ed(1,2)-points_ed(1,1))/2.d0
        !    dx(1)=points_ed(1,2)/2.d0
        !    dx(1)=0.d0
        !    dx(2)=0.d0
        !    dx(points_number)=(points_ed(1,points_number)-points_ed(1,points_number-1))/2.d0
        !    do i=2,points_number-1
        !       dx(i)=0.5d0*(points_ed(1,i+1)-points_ed(1,i-1))
        !    enddo

        do i = 1, points_number
            dx(i) = points_ed(1, i)
        end do

        do i = 1, points_number
            dx(i) = dx(i)*weights(i)
        end do

        ! for consistency the same as vecsigma in amoeba
        ! pin the value at the origin
        !    dx(1)=1.d0/abs(points_ed(1,1)-points_ed(1,2))/dble(points_number)/1000.d0
        !    do i=2,points_number-1
        !       dx(i)=1.d0/abs(points_ed(1,i)-points_ed(1,i+1))/dble(points_number)
        !    enddo
        !    dx(points_number)=dx(points_number-1)

        !    norm=1.d0
        !    do i=1,points_number
        !       norm=norm+dx(i)
        !    enddo
        !    norm=norm*abs(points_ed(1,points_number)-points_ed(1,1))/dble(points_number)**2
        !    do i=1,points_number
        !       dx(i)=dx(i)*sigma/norm
        !       dx(i)=1.d0/dx(i)**2
        !    enddo

        allocate (vecscratch(param_number/3))

        ! reset the data to fit based on the fixed functions
        ! new y_i = y_i - \sum_j alpha_j f_j (x_i)

        i = 0
        do while (i .lt. param_number)
            i = i + 1
            if (.not. sel_ed(i) .and. mod(i, 3) .eq. 1) then
                p(1) = param0_ed(i)*norm_ed(1 + (i - 1)/3)
                p(2) = param0_ed(i + 1)
                p(3) = param0_ed(i + 2)

                do j = 1, points_number
                    points_ed(2, j) = points_ed(2, j) - pseudofun(1, points_ed(1, j), p, vecscratch)
                end do

                i = i + 2

            end if
        end do

        allocate (overlap(bs_number, bs_number), proj(bs_number))

    end subroutine matrix_setup

    subroutine update_param_ed

        implicit none
        integer :: i, j

        ! set linear combination coefficients
        j = 0
        do i = 1, param_number
            if (sel_ed(i)) then
                j = j + 1
                param0_ed(i) = bs_coeff(j)
            end if
        end do

    end subroutine update_param_ed

    subroutine overlap_matrix

        implicit none
        integer :: i, j
        real(8) :: integral, p1(3), p2(3)

        do i = 1, bs_number
            p1(1) = norm_ed(1 + (point2param(i) - 1)/3)
            p1(2) = param0_ed(point2param(i) + 1)
            p1(3) = param0_ed(point2param(i) + 2)

            do j = 1, i
                p2(1) = norm_ed(1 + (point2param(j) - 1)/3)
                p2(2) = param0_ed(point2param(j) + 1)
                p2(3) = param0_ed(point2param(j) + 2)

                !          write(*,*) i,j,p1,p2

                call oned_integral(p1, p2, integral)

                overlap(i, j) = integral

                !          if(i.eq.1.and.j.eq.1) overlap(1,1)=0.5d0
                !          if(i.eq.2.and.j.eq.1) overlap(2,1)=0.0d0
                !          if(i.eq.2.and.j.eq.2) overlap(2,2)=1.5d0

                !          if(i.eq.j) write(6,*) i,i,overlap(i,i)

                if (i .ne. j) overlap(j, i) = overlap(i, j)

            end do
        end do

        do i = 1, bs_number
            p1(1) = norm_ed(1 + (point2param(i) - 1)/3)
            p1(2) = param0_ed(point2param(i) + 1)
            p1(3) = param0_ed(point2param(i) + 2)
            call oned_integral2(p1, integral)
            proj(i) = integral
            !       if(i.eq.1) proj(1)=0.5d0
            !       if(i.eq.2) proj(2)=1.5d0
        end do

    end subroutine overlap_matrix

    subroutine oned_integral(p1, p2, integral)

        implicit none

        real(8) :: integral, p1(3), p2(3)
        integer :: i

        integral = 0.d0
        do i = 1, points_number
            integral = integral + pseudofun(1, points_ed(1, i), p1, vecscratch)*pseudofun(1, points_ed(1, i), p2, vecscratch)*dx(i)
        end do

    end subroutine oned_integral

    subroutine oned_integral2(p1, integral)

        implicit none

        real(8) :: integral, p1(3)
        integer :: i

        integral = 0.d0
        do i = 1, points_number
            integral = integral + pseudofun(1, points_ed(1, i), p1, vecscratch)*points_ed(2, i)*dx(i)
        end do

    end subroutine oned_integral2

    subroutine orthonormalization(diag_type, eps)
        ! it uses either the cholesky factorization and diagonalization or
        ! the Sorella and Cavazzoni algorithm for overlap matrix preconditioning
        ! (still experimental in this particular case...)
        ! for the moment cholesky recommended

        implicit none

        integer :: nelorb_c, n, lwork, mine, lda, info, i, j, neig
        real(8) :: abstol, eps, condnumber, cost, dlamch, norm
        real(8), dimension(:, :), allocatable :: umatl
        real(8), dimension(:), allocatable :: work, eig, eigmat
        integer, dimension(:), allocatable :: iwork, ifail
        real(8), dimension(:, :), allocatable :: mat_in, overs
        character(60) :: diag_type

        nelorb_c = bs_number
        lda = bs_number
        abstol = 2.d0*dlamch('S')
        n = nelorb_c
        mine = 1

        allocate (overs(lda, nelorb_c))

        ! calculation normalized orbitals and eigenvectors of the overlap matrix
        ! maximum accuracy
        ! by stable pre-diagonalization of the matrix overs

        do j = 1, nelorb_c
            do i = 1, lda
                overs(i, j) = overlap(i, j)
                !          write(6,*) i,j,overs(i,j)
            end do
        end do

        if (trim(diag_type) .eq. 'cholesky') then

            !       allocate(umatl(nelorb_c,nelorb_c))
            lwork = 3*n
            allocate (work(lwork), iwork(n))

            ! cholesky factorization
            call DPOTRF('U', n, overs, n, info)
            if (info .ne. 0) write (6, *) 'problems in dpotrf'

            ! compute 1-norm of the matrix overs
            norm = 0.d0
            do j = 1, n
                norm = max(norm, sum(abs(overs(:, j))))
            end do

            call DPOCON('U', n, overs, n, norm, condnumber, work, iwork, info)
            if (info .ne. 0) write (6, *) 'problems in dpocon'

            write (6, *) '1-norm of overlap and condition number', norm
            write (6, *) 'condition number', condnumber

            ! compute inverse
            call DPOTRI('U', n, overs, n, info)
            if (info .ne. 0) write (6, *) 'problems in dpotri'

            ! full matrix
            do i = 1, n
                do j = 1, i - 1
                    overs(i, j) = overs(j, i)
                end do
            end do

            ! check
            !       call dgemm('N','N',n,n,n,1.d0,overs,n,overlap,n,0.d0,umatl,n)
            !       do i=1,n
            !          do j=1,n
            !             write(6,*) i,j,umatl(i,j)
            !          enddo
            !       enddo

            call dgemv('T', n, n, 1.d0, overs, n, proj, 1, 0.d0, bs_coeff, 1)

            ! other way directly calling dpotrs after dpotrf
            ! equivalent to the above
            !       do i=1,n
            !          bs_coeff(i)=proj(i)
            !       enddo

            !       call DPOTRS('U',n,1,overs,n,bs_coeff,n,info)

            !       if(info.lt.0) write(6,*) 'problems on dpotrs'

            deallocate (work, iwork)
            !       deallocate(umatl)

        elseif (trim(diag_type) .eq. 'sorella_cavazzoni') then
            ! perform Sorella Cavazzoni algorithm

            lwork = 30*n
            allocate (umatl(nelorb_c, nelorb_c), eig(nelorb_c), eigmat(nelorb_c))
            allocate (work(lwork), iwork(5*nelorb_c), ifail(nelorb_c))

            ! diagonalize the overlap matrix
            call dsyevx('V', 'A', 'L', n, overs, lda, 0.d0, 0.d0, 1, 1, abstol &
                        , neig, eig, umatl, nelorb_c, work, lwork, iwork, ifail, info)

            write (6, *) ' Lowest/Max  eigenvalue overlap mat =', eig(1), eig(n)
            condnumber = abs(eig(1)/eig(n))
            do i = 2, n
                cost = abs(eig(i)/eig(n))
                if (cost .lt. condnumber) condnumber = cost
            end do
            write (6, *) ' Condition number basis set =', condnumber

            do i = 1, n
                if (eig(i)/eig(n) .gt. abs(eps)) then ! the condition number criterium
                    eigmat(i) = dsqrt(1.d0/eig(i))
                else
                    mine = mine + 1
                    eigmat(i) = 0.d0
                end if
            end do

            if (info .gt. 0) write (6, *) ' info > 0 in dsyevx !!! ', info

            !       we assume here that the garbage eigenvectors are the ones
            !       close to zero eigenvalue.
            do i = 1, info
                if (eigmat(ifail(i)) .ne. 0.d0) then
                    if (ifail(i) .ge. mine) mine = ifail(i) + 1
                    eigmat(1:ifail(i)) = 0.d0
                end if
            end do

            if (mine .ne. 1) write (6, *) ' disregarded coll. =', mine - 1

            ! first transformation  umatl
            do i = 1, n
                umatl(:, i) = umatl(:, i)*eigmat(i)
            end do

            ! second diagonalization
            ! now compute the overlap matrix in the new basis

            allocate (mat_in(nelorb_c, nelorb_c))

            ! destroy previously computed overlap
            call dgemm('N', 'N', n, n, n, 1.d0, overlap, lda, umatl, n, 0.d0, mat_in, n)
            call dgemm('T', 'N', n, n, n, 1.d0, umatl, n, mat_in, lda, 0.d0, overlap, lda)

            do i = 1, n
                if (eigmat(i) .eq. 0) overlap(i, i) = -1000.d0 - dble(i)
            end do

            ! diagonalize again the newly rotated overlap matrix
            call dsyevx('V', 'A', 'L', n, overlap, lda, 0.d0, 0.d0, 1, 1, abstol &
                        , neig, eig, mat_in, n, work, lwork, iwork, ifail, info)

            ! check again eigenvalues
            do i = 1, n
                if (eig(i) .gt. eps) then
                    mat_in(:, i) = mat_in(:, i)/dsqrt(eig(i))
                else
                    if (eig(i) .gt. -1.d0) write (6, *) ' Further singular eigenvalue ', i
                    mat_in(:, i) = 0.d0
                end if
            end do
            overs(1:n, 1:n) = umatl(1:n, 1:n)

            ! Final dgemm: total rotation matrix
            call dgemm('N', 'N', n, n, n, 1.d0, overs, lda, mat_in, n, 0.d0, umatl, n)

            deallocate (mat_in)

            ! end second diagonalization

            ! solution for the coefficients of the linear combination
            call dgemv('T', n, n, 1.d0, umatl, lda, proj, 1, 0.d0, eigmat, 1)
            call dgemv('N', n, n, 1.d0, umatl, lda, eigmat, 1, 0.d0, bs_coeff, 1)

            !    write(6,*) bs_coeff

            deallocate (umatl, eig, eigmat)
            deallocate (work, iwork, ifail)

        else

            write (6, *) 'I do not know the diagonalization algorithm', trim(diag_type)
            write (6, *) 'Choose between "cholesky" and "sorella_cavazzoni"'
            stop

        end if

        deallocate (overs)

    end subroutine orthonormalization

end module exact_diagonalization
