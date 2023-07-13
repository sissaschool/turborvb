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
!

!
! It converts an AGP wf in a Pfaffian wave function
!

program convertpfaff
    use allio
    use constants
    use IO_m
    use sub_comm
    implicit none
    integer i, j, k, nelorbin, nelcolin, nelorbout, nelcolout, nunp, nshell_in, indparo, indpar, nelorbh_in, ipf_in, nelorbin2
    real(8), allocatable :: detmatin(:, :), detmatout(:, :), unpaired(:, :), mat_unp(:, :), dup_in(:)
    integer, allocatable :: ioptorb_in(:), nparam_in(:)

    logical :: rotmagn, unp2
    character(20) :: str
    character(100) :: name_tool
    real*8 scale_unp, angle_rot, phi_rot
    real*8, dimension(:, :), allocatable :: surot, suscra, sutry

    ! output version information
    call print_version

    !reading and saving the quantities of the fort.10_in

    open (unit=10, file='fort.10_in', form='formatted', status='unknown')
    call default_allocate
    call read_fort10(10)
    nshell_in = nshell_c
    if (contraction .ne. 0) then
        allocate (ioptorb_in(nshell_c), nparam_in(nshell_c), dup_in(size(dup_c)))
        ioptorb_in(1:nshell_c) = ioptorb_c(1:nshell_c)
        nparam_in(1:nshell_c) = nparam_c(1:nshell_c)
        dup_in = dup_c
    end if
    if (molyes) then
        write (6, *) ' ERROR the code does not work with molecular orbitals '
        stop
    end if
    !Understand if we have to rotate the spin on the plan or not
    ipf_in = ipf

    if (symmagp) then
        rotmagn = .false.
    else
        rotmagn = .true.
    end if

    call get_command_argument(1, str)

    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then
        ! Input the name of the file exactly as it is in /doc
        name_tool = 'convertpfaff'
        call help_online(name_tool)
        stop
    end if

    if (trim(str) .eq. "rotate") then
        rotmagn = .true.
        write (6, *) ' Warning forcing rotation magnetization '
    elseif (trim(str) .eq. "norotate") then
        rotmagn = .false.
        write (6, *) ' Warning forcing no rotation in magnetization '
    end if

    !Activate the creation of the upup or/and downdown part of the pfaffian
    !using two unpaired orbitals
    unp2 = .false.

    nunp = ndiff
    !  write (6,*) nelorb_c
    nelorbh_in = nelorbh
    if (contraction .eq. 0) then
        nelorbin = nelorbh*ipf
        nelcolin = nelcolh
    else
        nelorbin = nelorb_c
        nelcolin = nelcol_c
    end if

    allocate (detmatin(nelorbin*ipc, nelcolin))

    if (contraction .eq. 0) then
        do j = 1, nelcolin
            do i = 1, ipc*nelorbin
                detmatin(i, j) = detmat(i + (j - 1)*(ipc*nelorbin))
            end do
        end do
    else
        do j = 1, nelcolin
            do i = 1, ipc*nelorbin
                detmatin(i, j) = detmat_c(i + (j - 1)*(ipc*nelorbin))
            end do
        end do
    end if
    allocate (mat_unp(ipc*nelorbin, nelorbin))
    mat_unp = 0.d0
    scale_unp = 100.d0
    if (nunp .gt. 0) then
        write (6, *) "The agp has ", nunp, "unpaired"
        !     if (rotmagn) then
        !        write (6,*) "ERROR: AGP with unpaired cannot work with the rotated magnetization"
        !        stop
        write (6, *) ' Scale mean field ( infty=exact)? '
        read (5, *) scale_unp
        !     end if
        allocate (unpaired(ipc*nelorbin, 1:nunp))
        unp2 = .true.
        unpaired(:, 1:nunp) = detmatin(:, nelorbin + 1:nelorbin + nunp)
    end if

    close (10)
    call deallocate_all

    !reading the fort.10_out
    open (unit=10, file='fort.10_out', form='formatted', status='unknown')
    call default_allocate
    call read_fort10(10)

    close (10)
    if (molyes) then
        write (6, *) ' ERROR the code does not work with molecular orbitals '
        stop
    end if
    if (contraction .ne. 0) then
        if (ipf_in .eq. 1) then

            dup_c = 0.d0
            ! put in order dup_c
            indpar = 0
            indparo = 0
            do i = 1, nshell_c
                if (ioptorb_c(i) .eq. 900000) then
                    if (ipc .eq. 2) then
                        do j = 1, nparam_in(i)/2
                            dup_c(2*(indpar + j) - 1) = dup_in(2*(indparo + j) - 1)
                            dup_c(2*(indpar + j + nparam_in(i)/2) - 1) = dup_in(2*(indparo + j) - 1) + nelorbh_in
                        end do
                        do j = 1, nparam_in(i)
                            dup_c(2*(indpar + nparam_c(i)/2) + j) = dup_in(2*(indparo + nparam_in(i)/2) + j)
                            dup_c(2*(indpar + nparam_c(i)/2 + nparam_in(i)/2) + j) = dup_in(2*(indparo + nparam_in(i)/2) + j)
                        end do
                    else
                        do j = 1, nparam_in(i)/2
                            dup_c(indpar + j) = dup_in(indparo + j)
                            dup_c(indpar + j + nparam_in(i)/2) = dup_in(indparo + j) + nelorbh_in
                        end do
                        do j = 1, nparam_in(i)/2
                            dup_c(indpar + nparam_c(i)/2 + j) = dup_in(indparo + j + nparam_in(i)/2)
                            dup_c(indpar + nparam_c(i)/2 + nparam_in(i)/2 + j) = dup_in(indparo + j + nparam_in(i)/2)
                        end do
                    end if
                elseif (nparam_in(i) .gt. 1) then
                    do j = 1, ipc*nparam_in(i)
                        dup_c(ipc*indpar + j) = dup_in(ipc*indparo + j)
                    end do
                else
                    dup_c(ipc*indpar + 1) = dup_in(ipc*indparo + 1)
                end if
                indpar = indpar + nparam_c(i)
                indparo = indparo + nparam_in(i)
            end do
        else
            dup_c = dup_in
        end if
    end if

    if (contraction .eq. 0) then
        nelorbout = nelorbh*2
        nelcolout = nelcolh
    else
        nelorbout = nelorb_c
        nelcolout = nelcol_c
    end if

    if (nelorbout .ne. 2*nelorbin/ipf_in .or. nshell_c .ne. nshell_in) then
        write (6, *) "ERROR AGP not compatible with the pfaffian"
        write (6, *) "Contraction=", contraction
        write (6, *) "# shell input/output=", nshell_in, nshell_c
        write (6, *) "nelorbout=", nelorbout, "nelorbin=", nelorbin
        stop
    end if

    allocate (detmatout(ipc*nelorbout, nelcolout))
    detmatout = 0.d0

    write (6, *) 'nelorbin nelorbout =', nelorbin, nelorbout

    if (ipf_in .eq. 1) then

        if (ipc .eq. 2) then
            do i = 1, nelorbin
                do j = 1, nelorbin
                    detmatout(2*i - 1, j + nelorbin) = -detmatin(2*i - 1, j)
                    detmatout(2*i, j + nelorbin) = -detmatin(2*i, j)
                    detmatout(2*(j + nelorbin) - 1, i) = detmatin(2*i - 1, j)
                    detmatout(2*(j + nelorbin), i) = detmatin(2*i, j)
                end do
            end do
        else
            do i = 1, nelorbin
                do j = 1, nelorbin
                    detmatout(i, j + nelorbin) = -detmatin(i, j)
                    detmatout(j + nelorbin, i) = detmatin(i, j)
                end do
            end do
        end if

    else
        allocate (surot(ipc*2, 2), sutry(ipc*2, 2), suscra(ipc*2, 2))

        if (ipc .eq. 1) then
            write (6, *) ' Input angle (unit  Pi) '
            read (5, *) angle_rot
            angle_rot = angle_rot*pi/2.d0 ! For the spin
            surot(1, 1) = cos(angle_rot)
            surot(2, 2) = surot(1, 1)
            surot(2, 1) = -sin(angle_rot)
            surot(1, 2) = -surot(2, 1)
            write (6, *) ' Matrix rotation'
            do i = 1, 2
                do j = 1, 2
                    write (6, *) i, j, surot(i, j)
                end do
            end do
        else
            write (6, *) ' Input theta, phi (unit  Pi) '
            read (5, *) angle_rot, phi_rot
            angle_rot = angle_rot*pi/2.d0
            phi_rot = phi_rot*pi/2.d0

            call fill_surot(surot, angle_rot, phi_rot)
            write (6, *) ' Matrix rotation'
            do i = 1, 2
                do j = 1, 2
                    write (6, *) i, j, surot(2*i - 1, j), surot(2*i, j)
                end do
            end do

        end if

        suscra = 0.d0
        sutry = 0.d0

    end if

    if (unp2) then
        if (ipc .eq. 2) then
            do k = 1, nunp - ndiff, 2
                do i = 1, nelorbin
                    do j = 1, nelorbin
                        !Here I considered that the product is between two complex numbers
                        if (rotmagn) then
                            mat_unp(i*2 - 1, j) = mat_unp(2*i - 1, j) + unpaired(2*i - 1, k)*unpaired(2*j - 1, k + 1) &
                                                &- unpaired(2*i - 1, k + 1)*unpaired(2*j - 1, k) &
                                                &- (unpaired(2*i, k)*unpaired(2*j, k + 1) &
                                                &- unpaired(2*i, k + 1)*unpaired(2*j, k))
                            mat_unp(i*2, j) = mat_unp(2*i, j) + unpaired(2*i - 1, k)*unpaired(2*j, k + 1) &
                                            &- unpaired(2*i - 1, k + 1)*unpaired(2*j, k) &
                                            &+ unpaired(2*i, k)*unpaired(2*j - 1, k + 1) &
                                            &- unpaired(2*i, k + 1)*unpaired(2*j - 1, k)
                        else
                            detmatout(i*2 - 1, j) = detmatout(2*i - 1, j) &
                                                    &+ scale_unp*(unpaired(2*i - 1, k)*unpaired(2*j - 1, k + 1) &
                                                    &- unpaired(2*i - 1, k + 1)*unpaired(2*j - 1, k) &
                                                    &- (unpaired(2*i, k)*unpaired(2*j, k + 1) &
                                                                    - unpaired(2*i, k + 1)*unpaired(2*j, k)))
                            detmatout(i*2, j) = detmatout(2*i, j) &
                                                &+ scale_unp*(unpaired(2*i - 1, k)*unpaired(2*j, k + 1) &
                                                &- unpaired(2*i - 1, k + 1)*unpaired(2*j, k) &
                                                &+ unpaired(2*i, k)*unpaired(2*j - 1, k + 1) &
                                                &- unpaired(2*i, k + 1)*unpaired(2*j - 1, k))
                        end if
                    end do
                end do
            end do
        else
            write (6, *) ' Read unpaired ', nunp, ndiff
            do k = 1, nunp
                write (6, *) ' Unpaired # =', k
                do i = 1, nelorbin
                    write (6, *) i, unpaired(i, k)
                end do
            end do
            do k = 1, nunp - ndiff, 2
                do i = 1, nelorbin
                    do j = 1, nelorbin
                        if (rotmagn) then
                            mat_unp(i, j) = mat_unp(i, j) &
                                            &+ unpaired(i, k)*unpaired(j, k + 1) &
                                            &- unpaired(i, k + 1)*unpaired(j, k)
                        else
                            detmatout(i, j) = detmatout(i, j) &
                                            &+ scale_unp*(unpaired(i, k)*unpaired(j, k + 1) &
                                            &- unpaired(i, k + 1)*unpaired(j, k))
                        end if
                    end do
                end do
            end do
            if (ndiff .eq. 1) then
                !   There is only one unpaired in the Pfaffian
                if (rotmagn .and. ipf_in .eq. 1) then
                    do i = 1, nelorbin*ipc
                        detmatout(i, nelcolout) = unpaired(i, nunp)/dsqrt(2.d0)
                        detmatout(i + ipc*nelorbin, nelcolout) = unpaired(i, nunp)/dsqrt(2.d0)
                    end do
                else
                    do i = 1, nelorbin*ipc
                        detmatout(i, nelcolout) = unpaired(i, nunp)
                    end do
                end if
            end if
        end if
        mat_unp = mat_unp*scale_unp
        !I'm not sure I really got how this uppfaff works
        if (.not. pfaffup .and. symmagp) then
            write (6, *) ' Warning defining also down-down by symmetry '
            if (ipc .eq. 1 .or. .not. yes_hermite) then
                do i = 1, ipc*nelorbin
                    do j = 1, nelorbin
                        detmatout(i + nelorbin*ipc, j + nelorbin) = detmatout(i, j)
                    end do
                end do
            else
                !       ipc=2 and yes_hermite--> the complex conjugate
                do i = 1, 2*nelorbin, 2
                    do j = 1, nelorbin
                        detmatout(i + 2*nelorbin, j + nelorbin) = detmatout(i, j)
                        detmatout(i + 1 + 2*nelorbin, j + nelorbin) = -detmatout(i + 1, j)
                    end do
                end do
            end if
        end if
    end if

    if (rotmagn) then
        write (6, *) ' Warning rotating magnetic moment // x '
        if (ipf_in .eq. 2) then
            nelorbin2 = nelorbin/2
            if (ipc .eq. 1) then
                do i = 1, nelorbin2
                    do j = 1, nelorbin2
                        suscra(1, 1) = detmatin(i, j)
                        suscra(2, 1) = detmatin(i + nelorbin2, j)
                        suscra(1, 2) = detmatin(i, j + nelorbin2)
                        suscra(2, 2) = detmatin(i + nelorbin2, j + nelorbin2)
                        !        su2 rotation
                        call dgemm('N', 'N', 2, 2, 2, 1.d0, surot, 2, suscra, 2, 0.d0, sutry, 2)
                        call dgemm('N', 'T', 2, 2, 2, 1.d0, sutry, 2, surot, 2, 0.d0, suscra, 2)
                        detmatout(i, j) = suscra(1, 1)
                        detmatout(i + nelorbin2, j) = suscra(2, 1)
                        detmatout(i, j + nelorbin2) = suscra(1, 2)
                        detmatout(i + nelorbin2, j + nelorbin2) = suscra(2, 2)
                    end do
                end do
            else
                do i = 1, nelorbin2
                    do j = 1, nelorbin2
                        suscra(1, 1) = detmatin(2*i - 1, j)
                        suscra(2, 1) = detmatin(2*i, j)
                        suscra(3, 1) = detmatin(2*i - 1 + nelorbin2*2, j)
                        suscra(4, 1) = detmatin(2*i + nelorbin2*2, j)
                        suscra(1, 2) = detmatin(2*i - 1, j + nelorbin2)
                        suscra(2, 2) = detmatin(2*i, j + nelorbin2)
                        suscra(3, 2) = detmatin(2*i - 1 + nelorbin2*2, j + nelorbin2)
                        suscra(4, 2) = detmatin(2*i + nelorbin2*2, j + nelorbin2)
                        !        su2 rotation
                        call zgemm('N', 'N', 2, 2, 2, zone, surot, 2, suscra, 2, zzero, sutry, 2)
                        call zgemm('N', 'T', 2, 2, 2, zone, sutry, 2, surot, 2, zzero, suscra, 2)
                        detmatout(2*i - 1, j) = suscra(1, 1)
                        detmatout(2*i, j) = suscra(2, 1)
                        detmatout(2*i - 1 + nelorbin2*2, j) = suscra(3, 1)
                        detmatout(2*i + nelorbin2*2, j) = suscra(4, 1)
                        detmatout(2*i - 1, j + nelorbin2) = suscra(1, 2)
                        detmatout(2*i, j + nelorbin2) = suscra(2, 2)
                        detmatout(2*i - 1 + nelorbin2*2, j + nelorbin2) = suscra(3, 2)
                        detmatout(2*i + nelorbin2*2, j + nelorbin2) = suscra(4, 2)
                    end do
                end do

            end if

        else
            if (ipc .eq. 2) then
                do i = 1, nelorbin
                    do j = 1, nelorbin
                        detmatout(2*i - 1, j + nelorbin) = -(detmatin(2*i - 1, j) + detmatin(2*j - 1, i))*0.5
                        detmatout(2*i, j + nelorbin) = -(detmatin(2*i, j) + detmatin(2*j, i))*0.5
                        detmatout(2*(j + nelorbin) - 1, i) = (detmatin(2*i - 1, j) + detmatin(2*j - 1, i))*0.5
                        detmatout(2*(j + nelorbin), i) = (detmatin(2*i, j) + detmatin(2*j, i))*0.5
                        detmatout(2*i - 1, j) = (detmatin(2*i - 1, j) - detmatin(2*j - 1, i))*0.5
                        detmatout(2*i, j) = (detmatin(2*i, j) - detmatin(2*j, i))*0.5
                        detmatout(2*(i + nelorbin) - 1, j + nelorbin) = -(detmatin(2*i - 1, j) - detmatin(2*j - 1, i))*0.5
                        detmatout(2*(i + nelorbin), j + nelorbin) = -(detmatin(2*i, j) - detmatin(2*j, i))*0.5
                    end do
                end do
            else
                do i = 1, nelorbin
                    do j = 1, nelorbin
                        detmatout(i, j + nelorbin) = -(detmatin(i, j) + detmatin(j, i))*0.5d0
                        detmatout(j + nelorbin, i) = (detmatin(i, j) + detmatin(j, i))*0.5d0
                        detmatout(i, j) = (detmatin(i, j) - detmatin(j, i))*0.5d0
                        detmatout(i + nelorbin, j + nelorbin) = -(detmatin(i, j) - detmatin(j, i))*0.5d0
                    end do
                end do
            end if
            if (nunp .gt. 0) then
                if (ipc .eq. 2) then
                    do i = 1, nelorbin
                        do j = 1, nelorbin
                            detmatout(2*i - 1:2*i, j) = detmatout(2*i - 1:2*i, j) + 0.5d0*mat_unp(2*i - 1:2*i, j)
                            detmatout(2*i - 1 + nelorbin:2*i + nelorbin, j + nelorbin) &
                                = detmatout(2*i - 1 + nelorbin:2*i + nelorbin, j + nelorbin) + 0.5d0*mat_unp(2*i - 1:2*i, j)
                            detmatout(2*i - 1:2*i, j + nelorbin) = detmatout(2*i - 1:2*i, j + nelorbin)&
                                    & + 0.5d0*mat_unp(2*i - 1:2*i, j)
                            detmatout(2*j - 1 + nelorbin:2*j + nelorbin, i) = -detmatout(2*i - 1:2*i, j + nelorbin)
                        end do
                    end do
                else
                    do i = 1, nelorbin
                        do j = 1, nelorbin
                            detmatout(i, j) = detmatout(i, j) + 0.5d0*mat_unp(i, j)
                            detmatout(i + nelorbin, j + nelorbin) = detmatout(i + nelorbin, j + nelorbin) + 0.5d0*mat_unp(i, j)
                            detmatout(i, j + nelorbin) = detmatout(i, j + nelorbin) + 0.5d0*mat_unp(i, j)
                            detmatout(j + nelorbin, i) = -detmatout(i, j + nelorbin)
                        end do
                    end do
                end if
            end if ! nunp>0
        end if ! ipf=1
    end if ! rotmagn
    !Writing the output matrices

    if (contraction .eq. 0) then
        detmat = 0.d0
        do j = 1, nelcolout
            do i = 1, nelorbout*ipc
                detmat(i + (j - 1)*(nelorbout*ipc)) = detmatout(i, j)
            end do
        end do
    else
        detmat_c = 0.d0
        do j = 1, nelcolout
            do i = 1, nelorbout*ipc
                detmat_c(i + (j - 1)*(nelorbout*ipc)) = detmatout(i, j)
            end do
        end do
    end if

    open (unit=10, file='fort.10_new', form='formatted', status='unknown')
    call write_fort10(10)
    close (10)
    call deallocate_all
    if (allocated(surot)) deallocate (surot, sutry, suscra)
end program convertpfaff
subroutine fill_surot(surot, angle_rot, phi_rot)
    implicit none
    complex*16 surot(2, 2)
    real*8 angle_rot, phi_rot
    !  surot= U_z( 2 phi_rot) U_y( 2 angle_rot)  Landau eq.58.6
    surot(1, 1) = dcos(angle_rot)*exp(dcmplx(0.d0, phi_rot))
    surot(2, 2) = dconjg(surot(1, 1))
    surot(1, 2) = dsin(angle_rot)*exp(dcmplx(0.d0, phi_rot))
    surot(2, 1) = -dconjg(surot(1, 2))
    return
end

!########################################################
!namelist /option/

!#######################################################

