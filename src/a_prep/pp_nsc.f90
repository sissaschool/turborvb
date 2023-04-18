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

!-----------------------
subroutine pp_nsc(task)
!-----------------------
    use allio, only: rank, qp
    implicit none
    integer, intent(in) :: task

    select case (task)

    case (1) ! band structure calculation
        if (rank .eq. 0) call plot_bands

    case (2) ! DOS calculation
        if (rank .eq. 0) call evaluate_dos

        !case(3) ! Fermi surface calculation (TO BE IMPLEMENTED)
        !if(rank.eq.0) call evaluate_fermis

    case (5) ! evaluate all
        if (rank .eq. 0) then
            call plot_bands
            call evaluate_dos
        end if

    case (0)
        call error(' PP_nsc ', ' No tasks required!! ', -1, rank)

    case default
        call error(' PP_nsc ', ' Task type not recognized!! ', -1, rank)

    end select

    return

end subroutine pp_nsc

!---------------------
subroutine plot_bands
!---------------------

    ! Routine for plotting band structures from NSC calculations.
    ! kp_type must be 2 or 3 in order to run this routine.
    ! In the case of a kpoints path (kp_type=3), one already knows the number of k-points and the
    ! number of lines in the path, which are "nk" and "nk1-1" respectively.

    use setup, only: eigmol_sav, eigmoldo_sav, yeslsda, bands
    use kpoints_mod, only: nk, xkp
    use constants, only: ipc, energy_unit ! Ha to eV conversion
    implicit none
    logical, dimension(:), allocatable :: high_symmetry, in_range, in_rangedo
    integer i, j, k
    real(8), dimension(:), allocatable :: kp_path
    real(8) tmp_norm, k1(3), k2(3), dkl, dkl_sav, emin, emax
    character(len=80) :: filename

    ! initialization
    write (6, *) '---------------------------'
    write (6, *) ' Evaluating band structure '
    write (6, *) '---------------------------'

    filename = 'band_structure.dat'
    allocate (high_symmetry(nk), in_range(bands))
    in_range = .true.
    if (yeslsda .or. ipc .eq. 2) then
        allocate (in_rangedo(bands))
        in_rangedo = .true.
    end if
    allocate (kp_path(nk))
    high_symmetry = .false.
    kp_path = 0.d0
    !
    ! find high symmetry points in the 2D plot
    !
    do k = 1, nk
        if (k .eq. 1 .or. k .eq. nk) then
            high_symmetry(k) = .true.
        else
            k1(:) = xkp(:, k) - xkp(:, k - 1)
            k2(:) = xkp(:, k + 1) - xkp(:, k)
            tmp_norm = (k1(1)*k2(1) + k1(2)*k2(2) + k1(3)*k2(3))/ &
                       dsqrt(k1(1)*k1(1) + k1(2)*k1(2) + k1(3)*k1(3))/ &
                       dsqrt(k2(1)*k2(1) + k2(2)*k2(2) + k2(3)*k2(3))
            ! if changing direction in the 1BZ, then it's a high symmetry point
            high_symmetry(k) = (abs(tmp_norm - 1.d0) .gt. 1.0d-4)
        end if
        if (sum(xkp(:, k)**2) .lt. 1.0d-9) high_symmetry(k) = .true.
        ! distance between 2 kps.
        if (k .eq. 2) dkl_sav = dsqrt(k1(1)**2 + k1(2)**2 + k1(3)**2)
    end do
    !
    ! create 2D path
    !
    do k = 2, nk

        dkl = dsqrt((xkp(1, k) - xkp(1, k - 1))**2 + &
                    (xkp(2, k) - xkp(2, k - 1))**2 + &
                    (xkp(3, k) - xkp(3, k - 1))**2)

        if (dkl .gt. 5.d0*dkl_sav) then
            kp_path(k) = kp_path(k - 1)
        elseif (dkl .gt. 1.0d-5) then
            kp_path(k) = kp_path(k - 1) + dkl
            dkl_sav = dkl
        else
            kp_path(k) = kp_path(k - 1) + dkl
        end if

    end do
    ! find minimum and maximum values
    emin = eigmol_sav(1, 1)
    emax = eigmol_sav(1, 1)
    do k = 1, nk
        do j = 1, bands
            emin = min(emin, eigmol_sav(j, k))
            if (yeslsda .or. ipc .eq. 2) emin = min(emin, eigmoldo_sav(j, k))
            emax = max(emax, eigmol_sav(j, k))
            if (yeslsda .or. ipc .eq. 2) emax = max(emax, eigmoldo_sav(j, k))
        end do
    end do
    ! find bands values which are in the specified energy range
    in_range = .true.
    in_rangedo = .true.
    do i = 1, bands
        in_range(i) = any(eigmol_sav(i, 1:nk) >= emin .and. eigmol_sav(i, 1:nk) <= emax)
        in_rangedo(i) = any(eigmoldo_sav(i, 1:nk) >= emin .and. eigmoldo_sav(i, 1:nk) <= emax)
    end do

    write (6, *) in_range
    write (6, *) kp_path

    write (6, *) ' Band structures INFORMATION: '
    write (6, *) 'Minimum/maximum bands values:', emin, emax
    write (6, *) 'High symmetry points: '
    do k = 1, nk
        if (high_symmetry(k)) then
            write (6, *) kp_path(k), xkp(:, k)
        end if
    end do
    !
    ! now draw bands
    open (unit=80, file=trim(filename), form='formatted', status='unknown')
    write (80, *) ' # Spin up DFT bands (eV) '
    do i = 1, bands
        if (in_range(i)) then
            if (mod(i, 2) /= 0) then
                write (80, 100) (kp_path(j), energy_unit*eigmol_sav(i, j), j=1, nk)
            else
                write (80, 100) (kp_path(j), energy_unit*eigmol_sav(i, j), j=nk, 1, -1)
            end if
        end if
    end do
    if (yeslsda .or. ipc .eq. 2) then
        ! add two carriage return for gnuplot
        write (80, *)
        write (80, *)
        write (80, *) ' # Spin down DFT bands (eV)'
        do i = 1, bands
            if (in_rangedo(i)) then
                if (mod(i, 2) /= 0) then
                    write (80, 100) (kp_path(j), energy_unit*eigmoldo_sav(i, j), j=1, nk)
                else
                    write (80, 100) (kp_path(j), energy_unit*eigmoldo_sav(i, j), j=nk, 1, -1)
                end if
            end if
        end do
    end if
    !
    close (80)
    deallocate (high_symmetry, kp_path, in_range)
    if (yeslsda .or. ipc .eq. 2) deallocate (in_rangedo)

100 format(2f10.4)

    return

end subroutine plot_bands

!-----------------------
subroutine evaluate_dos
!-----------------------

    !
    ! Routine to calculate the density of states.
    ! Use optocc>0 and a non-zero value of the smearing (epsshell)
    ! in order to compute the density of states.
    !

    use constants, only: ipc, energy_unit
    use kpoints_mod, only: nk, xkp, wkp
    use setup, only: eigmol_sav, eigmoldo_sav, epsshell, optocc, deltaE, &
                     emin, emax, yeslsda, bands
    use compute_efermi, only: smearD

    implicit none
    integer :: i, j, k, ibnd, indk, smear_type
    integer :: ndos ! # of points in the DOS
    real(8), dimension(:), allocatable :: tdos, dos_up, dos_down
    real(qp) :: e, dt
    character(len=80) :: filename

    ! initialization
    write (6, *) '------------------------------'
    write (6, *) ' Evaluating density of states '
    write (6, *) '------------------------------'
    filename = 'density_of_states.dat'
    if (optocc .eq. 0 .or. epsshell .eq. 0.d0) then
        write (6, *) 'Warning: DOS must be updated with a smearing function! Setting default values.'
        optocc = 1
        epsshell = 0.01_qp
    end if
    smear_type = optocc
    dt = epsshell
    ! setting min/max values of the energies
    ! if not chosen from the user
    if (emin .eq. 0.d0 .and. emax .eq. 0.d0) then
        emin = eigmol_sav(1, 1)
        emax = eigmol_sav(1, 1)
        do k = 1, nk
            do j = 1, bands
                emin = min(emin, eigmol_sav(j, k))
                if (yeslsda .or. ipc .eq. 2) emin = min(emin, eigmoldo_sav(j, k))
                emax = max(emax, eigmol_sav(j, k))
                if (yeslsda .or. ipc .eq. 2) emax = max(emax, eigmoldo_sav(j, k))
            end do
        end do
    end if
    !
    ndos = int((emax - emin)/deltaE)
    allocate (dos_up(ndos + 1), dos_down(ndos + 1), tdos(ndos + 1))

    dos_up = 0.d0
    dos_down = 0.d0
    tdos = 0.d0
    e = 0.d0

    open (unit=80, file=trim(filename), form='formatted', status='unknown')
    write (80, *) '#  energy (eV)       Density_of_states (spin up)'
    do i = 0, ndos
        e = emin + deltaE*i
        do indk = 1, nk
            do ibnd = 1, bands
                dos_up(i + 1) = dos_up(i + 1) + wkp(indk)*smearD(eigmol_sav(ibnd, nk), dt, e, smear_type)
            end do
        end do
        write (80, 100) e*energy_unit, dos_up(i + 1)
    end do

    if (yeslsda .or. ipc .eq. 2) then
        write (80, *)
        write (80, *)
        write (80, *) '#  energy (eV)       Density_of_states (spin down)'
        do i = 0, ndos
            e = emin + deltaE*i
            do indk = 1, nk
                do ibnd = 1, bands
                    dos_down(i + 1) = dos_down(i + 1) + wkp(indk)*smearD(eigmoldo_sav(ibnd, indk), dt, e, smear_type)
                end do
            end do
            write (80, 100) e*energy_unit, dos_down(i + 1)
        end do
    end if

    deallocate (dos_up, dos_down, tdos)
    write (6, *)
    close (80)

100 format(2f10.4)

    return

end subroutine evaluate_dos
