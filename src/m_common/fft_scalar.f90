! Copyright (C) 2022 TurboRVB group based on code by
! Copyright (C) 2001-2008 Quantum-ESPRESSO group
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

!--------------------------------------------------------------------------!
! FFT scalar drivers Module - contains machine-dependent routines for: !
! FFTW, FFTW3, ESSL, LINUX_ESSL, SCSL, SUNPERF, NEC ASL and ACML libraries !
! (both 3d for serial execution and 1d+2d FFTs for parallel execution, !
! excepted NEC ASL, 3d only, no parallel execution) !
! Written by Carlo Cavazzoni, modified by P. Giannozzi, contributions !
! by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc, !
! Nicolas Lacorne, Filippo Spiga - Last update Aug 2008 !
!--------------------------------------------------------------------------!

!=----------------------------------------------------------------------=!
module fft_scalar
    !=----------------------------------------------------------------------=!
    use kinds

    implicit none
    save

    private
    public :: cft_1z, cft_2xy, cft_b, cfft3d, cfft3ds
    public :: good_fft_dimension, allowed, good_fft_order

    ! ... Local Parameter

    ! ndims Number of different FFT tables that the module
    ! could keep into memory without reinitialization
    ! nfftx Max allowed fft dimension

    integer, parameter :: ndims = 3, nfftx = 65537

    ! Workspace that is statically allocated is defined here
    ! in order to avoid multiple copies of the same workspace
    ! lwork: Dimension of the work space array (if any)
    !=----------------------------------------------------------------------=!
contains
    !=----------------------------------------------------------------------=!
    !
    !=----------------------------------------------------------------------=!
    !
    !
    !
    ! FFT along "z"
    !
    !
    !
    !=----------------------------------------------------------------------=!
    !
    subroutine cft_1z(c, nsl, nz, ldz, isign)
        ! driver routine for nsl 1d complex fft's of length nz
        ! ldz >= nz is the distance between sequences to be transformed
        ! (ldz>nz is used on some architectures to reduce memory conflicts)
        ! input : c(ldz*nsl) (complex)
        ! output : c(ldz*nsl) (complex - NOTA BENE: transform is in-place!)
        ! isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
        ! Up to "ndims" initializations (for different combinations of input
        ! parameters nz, nsl, ldz) are stored and re-used if available
        integer, intent(IN) :: isign
        integer, intent(IN) :: nsl, nz, ldz
        complex(DP) :: c(:)
        real(DP) :: tscale
        integer :: i, err, idir, ip
        integer, save :: zdims(3, ndims) = -1
        integer, save :: icurrent = 1
        logical :: done
        ! ... Machine-Dependent parameters, work arrays and tables of factors
        ! ltabl Dimension of the tables of factors calculated at the
        ! initialization stage
        integer*8, save :: fw_planz(ndims) = 0
        integer*8, save :: bw_planz(ndims) = 0
        ! Pointers to the "C" structures containing FFT factors ( PLAN )
        ! integer*8 is defined in include/f_defs.h
        ! for 32bit executables, integer*8 is integer(4)
        ! for 64bit executables, integer*8 is integer(8)
        if (nsl < 0) then
            call errore(" fft_scalar: cft_1z ", " nsl out of range ", nsl)
        end if
        !
        ! Here initialize table only if necessary
        !
        do ip = 1, ndims
            ! first check if there is already a table initialized
            ! for this combination of parameters
            done = (nz == zdims(1, ip))
            if (done) exit
        end do
        if (.not. done) then
            ! no table exist for these parameters
            ! initialize a new one
            ! WRITE( stdout, fmt="('DEBUG cft_1z, reinitializing tables ', I3)" ) icurrent
            if (fw_planz(icurrent) /= 0) call DESTROY_PLAN_1D(fw_planz(icurrent))
            if (bw_planz(icurrent) /= 0) call DESTROY_PLAN_1D(bw_planz(icurrent))
            idir = -1; call CREATE_PLAN_1D(fw_planz(icurrent), nz, idir)
            idir = 1; call CREATE_PLAN_1D(bw_planz(icurrent), nz, idir)
            zdims(1, icurrent) = nz; zdims(2, icurrent) = nsl; zdims(3, icurrent) = ldz
            ip = icurrent
            icurrent = mod(icurrent, ndims) + 1
        end if
        !
        ! Now perform the FFTs using machine specific drivers
        !
        if (isign < 0) then
            call FFT_Z_STICK(fw_planz(ip), c(1), ldz, nsl)
            c(1:ldz*nsl) = c(1:ldz*nsl)/nz
        else if (isign > 0) then
            call FFT_Z_STICK(bw_planz(ip), c(1), ldz, nsl)
        end if
        return
    end subroutine cft_1z
    !
    !
    !=----------------------------------------------------------------------=!
    !
    !
    !
    ! FFT along "x" and "y" direction
    !
    !
    !
    !=----------------------------------------------------------------------=!
    !
    !
    subroutine cft_2xy(r, nzl, nx, ny, ldx, ldy, isign)
        ! driver routine for nzl 2d complex fft's of lengths nx and ny
        ! input : r(ldx*ldy) complex, transform is in-place
        ! ldx >= nx, ldy >= ny are the physical dimensions of the equivalent
        ! 2d array: r2d(ldx, ldy) (x first dimension, y second dimension)
        ! (ldx>nx, ldy>ny used on some architectures to reduce memory conflicts)
        ! isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
        ! Up to "ndims" initializations (for different combinations of input
        ! parameters nx,ny,nzl,ldx) are stored and re-used if available
        implicit none
        integer, intent(IN) :: isign, ldx, ldy, nx, ny, nzl
        complex(DP) :: r(:)
        integer :: i, k, j, err, idir, ip, kk
        real(DP) :: tscale
        integer, save :: icurrent = 1
        integer, save :: dims(4, ndims) = -1
        logical :: done
        integer, parameter :: stdout = 6
        integer*8, save :: fw_planx(ndims) = 0, fw_plany(ndims) = 0
        integer*8, save :: bw_planx(ndims) = 0, bw_plany(ndims) = 0
        !
        ! Here initialize table only if necessary
        !
        do ip = 1, ndims
            ! first check if there is already a table initialized
            ! for this combination of parameters
            done = (ny == dims(1, ip)) .and. (nx == dims(3, ip))
            if (done) exit
        end do
        if (.not. done) then
            ! no table exist for these parameters
            ! initialize a new one
            ! WRITE( stdout, fmt="('DEBUG cft_2xy, reinitializing tables ', I3)" ) icurrent
            if (fw_plany(icurrent) /= 0) call DESTROY_PLAN_1D(fw_plany(icurrent))
            if (bw_plany(icurrent) /= 0) call DESTROY_PLAN_1D(bw_plany(icurrent))
            idir = -1; call CREATE_PLAN_1D(fw_plany(icurrent), ny, idir)
            idir = 1; call CREATE_PLAN_1D(bw_plany(icurrent), ny, idir)
            if (fw_planx(icurrent) /= 0) call DESTROY_PLAN_1D(fw_planx(icurrent))
            if (bw_planx(icurrent) /= 0) call DESTROY_PLAN_1D(bw_planx(icurrent))
            idir = -1; call CREATE_PLAN_1D(fw_planx(icurrent), nx, idir)
            idir = 1; call CREATE_PLAN_1D(bw_planx(icurrent), nx, idir)
            dims(1, icurrent) = ny; dims(2, icurrent) = ldx
            dims(3, icurrent) = nx; dims(4, icurrent) = nzl
            ip = icurrent
            icurrent = mod(icurrent, ndims) + 1
        end if
        !
        ! Now perform the FFTs using machine specific drivers
        !
        if (isign < 0) then
            call FFT_X_STICK(fw_planx(ip), r(1), nx, ny, nzl, ldx, ldy)
            do i = 1, nx
                do k = 1, nzl
                    j = i + ldx*ldy*(k - 1)
                    call FFT_Y_STICK(fw_plany(ip), r(j), ny, ldx)
                end do
            end do
            tscale = 1.0_dp/(nx*ny)
            call ZDSCAL(ldx*ldy*nzl, tscale, r(1), 1)
        else if (isign > 0) then
            do i = 1, nx
                do k = 1, nzl
                    j = i + ldx*ldy*(k - 1)
                    call FFT_Y_STICK(bw_plany(ip), r(j), ny, ldx)
                end do
            end do
            call FFT_X_STICK(bw_planx(ip), r(1), nx, ny, nzl, ldx, ldy)
        end if
        return
    end subroutine cft_2xy
    !
    !=----------------------------------------------------------------------=!
    !
    !
    !
    ! 3D scalar FFTs
    !
    !
    !
    !=----------------------------------------------------------------------=!
    !
    subroutine cfft3d(f, nx, ny, nz, ldx, ldy, ldz, isign)
        ! driver routine for 3d complex fft of lengths nx, ny, nz
        ! input : f(ldx*ldy*ldz) complex, transform is in-place
        ! ldx >= nx, ldy >= ny, ldz >= nz are the physical dimensions
        ! of the equivalent 3d array: f3d(ldx,ldy,ldz)
        ! (ldx>nx, ldy>ny, ldz>nz may be used on some architectures
        ! to reduce memory conflicts - not implemented for FFTW)
        ! isign > 0 : f(G) => f(R) ; isign < 0 : f(R) => f(G)
        !
        ! Up to "ndims" initializations (for different combinations of input
        ! parameters nx,ny,nz) are stored and re-used if available
        implicit none
        integer, intent(IN) :: nx, ny, nz, ldx, ldy, ldz, isign
        complex(DP) :: f(:)
        integer :: i, k, j, err, idir, ip
        real(DP) :: tscale
        integer, save :: icurrent = 1
        integer, save :: dims(3, ndims) = -1
        integer*8, save :: fw_plan(ndims) = 0
        integer*8, save :: bw_plan(ndims) = 0
        if (nx < 1) &
            call errore('cfft3', ' nx is less than 1 ', 1)
        if (ny < 1) &
            call errore('cfft3', ' ny is less than 1 ', 1)
        if (nz < 1) &
            call errore('cfft3', ' nz is less than 1 ', 1)
        !
        ! Here initialize table only if necessary
        !
        ip = -1
        do i = 1, ndims
            ! first check if there is already a table initialized
            ! for this combination of parameters
            if ((nx == dims(1, i)) .and. &
                (ny == dims(2, i)) .and. &
                (nz == dims(3, i))) then
                ip = i
                exit
            end if
        end do
        if (ip == -1) then
            ! no table exist for these parameters
            ! initialize a new one
            if (nx /= ldx .or. ny /= ldy .or. nz /= ldz) &
                call errore('cfft3', 'not implemented', 1)
            if (fw_plan(icurrent) /= 0) call DESTROY_PLAN_3D(fw_plan(icurrent))
            if (bw_plan(icurrent) /= 0) call DESTROY_PLAN_3D(bw_plan(icurrent))
            idir = -1; call CREATE_PLAN_3D(fw_plan(icurrent), nx, ny, nz, idir)
            idir = 1; call CREATE_PLAN_3D(bw_plan(icurrent), nx, ny, nz, idir)
            dims(1, icurrent) = nx; dims(2, icurrent) = ny; dims(3, icurrent) = nz
            ip = icurrent
            icurrent = mod(icurrent, ndims) + 1
        end if
        !
        ! Now perform the 3D FFT using the machine specific driver
        !
        if (isign < 0) then
            call FFTW_INPLACE_DRV_3D(fw_plan(ip), 1, f(1), 1, 1)
            tscale = 1.0_dp/dble(nx*ny*nz)
            call ZDSCAL(nx*ny*nz, tscale, f(1), 1)
        else if (isign > 0) then
            call FFTW_INPLACE_DRV_3D(bw_plan(ip), 1, f(1), 1, 1)
        end if
        return
    end subroutine cfft3d
    !
    !=----------------------------------------------------------------------=!
    !
    !
    !
    ! 3D scalar FFTs, but using sticks!
    !
    !
    !
    !=----------------------------------------------------------------------=!
    !
    subroutine cfft3ds(f, nx, ny, nz, ldx, ldy, ldz, isign, &
                       do_fft_x, do_fft_y)
        !
        ! driver routine for 3d complex "reduced" fft - see cfft3d
        ! The 3D fft are computed only on lines and planes which have
        ! non zero elements. These lines and planes are defined by
        ! the two integer vectors do_fft_x(ldy*nz) and do_fft_y(nz)
        ! (1 = perform fft, 0 = do not perform fft)
        ! The routine is implemented for essl and fftw library only
        !
        !----------------------------------------------------------------------
        !
        implicit none
        integer :: nx, ny, nz, ldx, ldy, ldz, isign
        !
        ! logical dimensions of the fft
        ! physical dimensions of the f array
        ! sign of the transformation
        complex(DP) :: f(ldx*ldy*ldz)
        integer :: do_fft_x(:), do_fft_y(:)
        !
        integer :: m, incx1, incx2
        integer :: i, k, j, err, idir, ip, ii, jj
        real(DP) :: tscale
        integer, save :: icurrent = 1
        integer, save :: dims(3, ndims) = -1
        integer*8, save :: fw_plan(3, ndims) = 0
        integer*8, save :: bw_plan(3, ndims) = 0
        tscale = 1.0_dp
        !
        ! ESSL sign convention for fft's is the opposite of the "usual" one
        !
        ! WRITE( stdout, fmt="('DEBUG cfft3ds :',6I6)") nx, ny, nz, ldx, ldy, ldz
        ! WRITE( stdout, fmt="('DEBUG cfft3ds :',24I2)") do_fft_x
        ! WRITE( stdout, fmt="('DEBUG cfft3ds :',24I2)") do_fft_y
        if (ny /= ldy) &
            call errore(' cfft3ds ', ' wrong dimensions: ny /= ldy ', 1)
        ip = -1
        do i = 1, ndims
            ! first check if there is already a table initialized
            ! for this combination of parameters
            if ((nx == dims(1, i)) .and. (ny == dims(2, i)) .and. &
                (nz == dims(3, i))) then
                ip = i
                exit
            end if
        end do
        if (ip == -1) then
            ! no table exist for these parameters
            ! initialize a new one
            if (fw_plan(1, icurrent) /= 0) call DESTROY_PLAN_1D(fw_plan(1, icurrent))
            if (bw_plan(1, icurrent) /= 0) call DESTROY_PLAN_1D(bw_plan(1, icurrent))
            if (fw_plan(2, icurrent) /= 0) call DESTROY_PLAN_1D(fw_plan(2, icurrent))
            if (bw_plan(2, icurrent) /= 0) call DESTROY_PLAN_1D(bw_plan(2, icurrent))
            if (fw_plan(3, icurrent) /= 0) call DESTROY_PLAN_1D(fw_plan(3, icurrent))
            if (bw_plan(3, icurrent) /= 0) call DESTROY_PLAN_1D(bw_plan(3, icurrent))
            idir = -1; call CREATE_PLAN_1D(fw_plan(1, icurrent), nx, idir)
            idir = 1; call CREATE_PLAN_1D(bw_plan(1, icurrent), nx, idir)
            idir = -1; call CREATE_PLAN_1D(fw_plan(2, icurrent), ny, idir)
            idir = 1; call CREATE_PLAN_1D(bw_plan(2, icurrent), ny, idir)
            idir = -1; call CREATE_PLAN_1D(fw_plan(3, icurrent), nz, idir)
            idir = 1; call CREATE_PLAN_1D(bw_plan(3, icurrent), nz, idir)
            dims(1, icurrent) = nx; dims(2, icurrent) = ny; dims(3, icurrent) = nz
            ip = icurrent
            icurrent = mod(icurrent, ndims) + 1
        end if
        if (isign > 0) then
            !
            ! i - direction ...
            !
            incx1 = 1; incx2 = ldx; m = 1
            do k = 1, nz
                do j = 1, ny
                    jj = j + (k - 1)*ldy
                    ii = 1 + ldx*(jj - 1)
                    if (do_fft_x(jj) == 1) then
                        call FFTW_INPLACE_DRV_1D(bw_plan(1, ip), m, f(ii), incx1, incx2)
                    end if
                end do
            end do
            !
            ! ... j-direction ...
            !
            incx1 = ldx; incx2 = 1; m = nx
            do k = 1, nz
                ii = 1 + ldx*ldy*(k - 1)
                if (do_fft_y(k) == 1) then
                    call FFTW_INPLACE_DRV_1D(bw_plan(2, ip), m, f(ii), incx1, incx2)
                end if
            end do
            !
            ! ... k-direction
            !
            incx1 = ldx*ldy; incx2 = 1; m = ldx*ny
            call FFTW_INPLACE_DRV_1D(bw_plan(3, ip), m, f(1), incx1, incx2)
        else
            !
            ! ... k-direction
            !
            incx1 = ldx*ny; incx2 = 1; m = ldx*ny
            call FFTW_INPLACE_DRV_1D(fw_plan(3, ip), m, f(1), incx1, incx2)
            !
            ! ... j-direction ...
            !
            incx1 = ldx; incx2 = 1; m = nx
            do k = 1, nz
                ii = 1 + ldx*ldy*(k - 1)
                if (do_fft_y(k) == 1) then
                    call FFTW_INPLACE_DRV_1D(fw_plan(2, ip), m, f(ii), incx1, incx2)
                end if
            end do
            !
            ! i - direction ...
            !
            incx1 = 1; incx2 = ldx; m = 1
            do k = 1, nz
                do j = 1, ny
                    jj = j + (k - 1)*ldy
                    ii = 1 + ldx*(jj - 1)
                    if (do_fft_x(jj) == 1) then
                        call FFTW_INPLACE_DRV_1D(fw_plan(1, ip), m, f(ii), incx1, incx2)
                    end if
                end do
            end do
            call DSCAL(2*ldx*ldy*nz, 1.0_dp/(nx*ny*nz), f(1), 1)
        end if
        return
    end subroutine cfft3ds
    !
    !=----------------------------------------------------------------------=!
    !
    !
    !
    ! 3D parallel FFT on sub-grids
    !
    !
    !
    !=----------------------------------------------------------------------=!
    !
    subroutine cft_b(f, nx, ny, nz, ldx, ldy, ldz, imin3, imax3, sgn)
        ! driver routine for 3d complex fft's on box grid, parallel case
        ! fft along xy is done only on planes that correspond to dense grid
        ! planes on the current processor, i.e. planes with imin3 <= nz <= imax3
        ! implemented for essl, fftw, scsl, complib, only for sgn=1 (f(R) => f(G))
        ! (beware: here the "essl" convention for the sign of the fft is used!)
        !
        implicit none
        integer nx, ny, nz, ldx, ldy, ldz, imin3, imax3, sgn
        complex(8) :: f(:)
        integer isign, naux, ibid, nplanes, nstart, k
        real(DP) :: tscale
        integer :: ip, i
        integer, save :: icurrent = 1
        integer, save :: dims(4, ndims) = -1
        integer*8, save :: bw_planz(ndims) = 0
        integer*8, save :: bw_planxy(ndims) = 0
        isign = -sgn
        tscale = 1.0_dp
        if (isign > 0) then
            call errore('cft_b', 'not implemented', isign)
        end if
        !
        ! 2d fft on xy planes - only needed planes are transformed
        ! note that all others are left in an unusable state
        !
        nplanes = imax3 - imin3 + 1
        nstart = (imin3 - 1)*ldx*ldy + 1
        !
        ! Here initialize table only if necessary
        !
        ip = -1
        do i = 1, ndims
            ! first check if there is already a table initialized
            ! for this combination of parameters
            if ((nx == dims(1, i)) .and. (ny == dims(2, i)) .and. &
                (nz == dims(3, i)) .and. (nplanes == dims(4, i))) then
                ip = i
                exit
            end if
        end do
        if (ip == -1) then
            ! no table exist for these parameters
            ! initialize a new one
            if (bw_planz(icurrent) /= 0) &
                call DESTROY_PLAN_1D(bw_planz(icurrent))
            call CREATE_PLAN_1D(bw_planz(icurrent), nz, 1)
            if (bw_planxy(icurrent) /= 0) &
                call DESTROY_PLAN_2D(bw_planxy(icurrent))
            call CREATE_PLAN_2D(bw_planxy(icurrent), nx, ny, 1)
            !
            dims(1, icurrent) = nx; dims(2, icurrent) = ny
            dims(3, icurrent) = nz; dims(4, icurrent) = nplanes
            ip = icurrent
            icurrent = mod(icurrent, ndims) + 1
        end if
        call FFTW_INPLACE_DRV_1D(bw_planz(ip), ldx*ldy, f(1), ldx*ldy, 1)
        call FFTW_INPLACE_DRV_2D(bw_planxy(ip), nplanes, f(nstart), 1, ldx*ldy)
        return
    end subroutine cft_b
    !
    !=----------------------------------------------------------------------=!
    !
    !
    !
    ! FFT support Functions/Subroutines
    !
    !
    !
    !=----------------------------------------------------------------------=!
    !
    !
    integer function good_fft_dimension(n)
        !
        ! Determines the optimal maximum dimensions of fft arrays
        ! Useful on some machines to avoid memory conflicts
        !
        use kinds, only: DP
        implicit none
        integer :: n, nx
        real(DP) :: log2n
        !
        ! this is the default: max dimension = fft dimension
        nx = n
        !
        !
        good_fft_dimension = nx
        return
    end function good_fft_dimension
    !=----------------------------------------------------------------------=!
    function allowed(nr)
        ! find if the fft dimension is a good one
        ! a "bad one" is either not implemented (as on IBM with ESSL)
        ! or implemented but with awful performances (most other cases)
        use kinds
        implicit none
        integer :: nr
        logical :: allowed
        integer :: pwr(5)
        integer :: mr, i, fac, p, maxpwr
        integer :: factors(5) = (/2, 3, 5, 7, 11/)
        ! find the factors of the fft dimension
        mr = nr
        pwr = 0
        factors_loop: do i = 1, 5
            fac = factors(i)
            maxpwr = nint(log(dble(mr))/log(dble(fac))) + 1
            do p = 1, maxpwr
                if (mr == 1) exit factors_loop
                if (mod(mr, fac) == 0) then
                    mr = mr/fac
                    pwr(i) = pwr(i) + 1
                end if
            end do
        end do factors_loop
        if (nr /= (mr*2**pwr(1)*3**pwr(2)*5**pwr(3)*7**pwr(4)*11**pwr(5))) &
            call errore(' allowed ', ' what ?!? ', 1)
        if (mr /= 1) then
            ! fft dimension contains factors > 11 : no good in any case
            allowed = .false.
        else
            ! fftw and all other cases: no factors 7 and 11
            allowed = ((pwr(4) == 0) .and. (pwr(5) == 0))
        end if
        return
    end function allowed
    !=----------------------------------------------------------------------=!
    integer function good_fft_order(nr, np)
        !
        ! This function find a "good" fft order value grather or equal to "nr"
        !
        ! nr (input) tentative order n of a fft
        !
        ! np (optional input) if present restrict the search of the order
        ! in the ensamble of multiples of np
        !
        ! Output: the same if n is a good number
        ! the closest higher number that is good
        ! an fft order is not good if not implemented (as on IBM with ESSL)
        ! or implemented but with awful performances (most other cases)
        !
        implicit none
        integer, intent(IN) :: nr
        integer, optional, intent(IN) :: np
        integer :: new
        new = nr
        if (present(np)) then
            do while (((.not. allowed(new)) .or. (mod(new, np) /= 0)) .and. (new <= nfftx))
                new = new + 1
            end do
        else
            do while ((.not. allowed(new)) .and. (new <= nfftx))
                new = new + 1
            end do
        end if
        if (new > nfftx) &
            call errore(' good_fft_order ', ' fft order too large ', new)
        good_fft_order = new
        return
    end function good_fft_order
    !=----------------------------------------------------------------------=!
end module fft_scalar
!=----------------------------------------------------------------------=!
