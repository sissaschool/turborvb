#include "c_lister.h"

program makefun_tester

    implicit none

    interface
        subroutine list_directory(path, path_length, files, num_files, prefix, prefix_length) bind(C, name="list_directory")
            use, intrinsic :: iso_c_binding
            character(c_char), intent(in) :: path(*)
            character(c_char), intent(in) :: prefix(*)
            character(c_char), intent(out) :: files(MAX_FILE_LENGTH, MAX_FILES)
            integer(c_int), intent(in), value :: prefix_length
            integer(c_int), intent(in), value :: path_length
            integer(c_int), intent(out) :: num_files
        end subroutine list_directory
    end interface

    integer*4 :: iorb, indt, indpar, indorb, indshell, i0, iflagnorm_unused&
              &, indtm, indtmin, typec, nelskip, num_lines, ii
    real*8 :: cr 
    real*8, dimension(:), allocatable :: dd, zeta, r
    real*8, dimension(:,:), allocatable :: z, rmu, distp

    integer*4, dimension(:), allocatable :: iorbs, multiplicities, npars, failed_test
    character*80 :: dummy_1, dummy_2, dummy_3, dummy_4

    character*12000 :: error_message = "", tmp

    logical :: failed

    call initialize()

    cr = 0.0d0
    iflagnorm_unused = 0
    indtm = 5
    indt = 10

    ! Load parameters form file parameters.csv
    open(unit=10, file='parameters.csv', status='old', action='read')
    ! Skip header
    read(10,*)
    ! Count number of lines
    num_lines = 0
    do
        read(10,*,iostat=iflagnorm_unused)
        if (iflagnorm_unused /= 0) exit
        num_lines = num_lines + 1
    end do
    ! Allocate arrays
    allocate(iorbs(num_lines))
    allocate(multiplicities(num_lines))
    allocate(npars(num_lines))
    allocate(failed_test(num_lines))

    ! -1 means that test not failed
    failed_test = -1

    ! Rewind file
    rewind(10)
    ! Skip header
    read(10,*)
    ! Read data
    do ii = 1, num_lines
        read(10,*) iorbs(ii), dummy_1, dummy_2, dummy_3, dummy_4, multiplicities(ii), npars(ii)
    end do
    close(10)

    do ii = 1, num_lines
        failed = .false.
        write(*,'("Checking orbital index ", I3)') iorbs(ii)
        write(*, *)
    
        if (allocated(dd)) deallocate(dd)
        allocate(dd(npars(ii)))
    
        if (allocated(rmu)) deallocate(rmu)
        allocate(rmu(3,0:indtm))
    
        if (allocated(r)) deallocate(r)
        allocate(r(0:indtm))
    
        if (allocated(z)) deallocate(z)
        allocate(z(multiplicities(ii),0:indt+4))

        if (allocated(distp)) deallocate(distp)
        allocate(distp(0:indt+4,20))
    
        !call create_single_value()
        !call create_pa_value()
        !call create_svgl_value()
        call check_index_movement()
        call check_single_value()
        call check_svgl_value()
        call check_pa_value()

        if (failed) then
            failed_test(ii) = 1
            tmp = trim(error_message)
            write(error_message, '(A,A,I3.3)') trim(tmp), ",", iorbs(ii)
        end if
    
        deallocate(distp)
        deallocate(z)
        deallocate(r)
        deallocate(dd)
        deallocate(rmu)

        write(*, *)
        write(*, '("Result = ", L1)') .not. failed
        write(*, *)
        write(*, '("######################################")')
    end do
    
    allocate(zeta(1))
    allocate(r(1))
    
    if (allocated(iorbs)) deallocate(iorbs)
    if (allocated(multiplicities)) deallocate(multiplicities)
    if (allocated(npars)) deallocate(npars)
    if (allocated(failed_test)) deallocate(failed_test)
    
    if (allocated(zeta)) deallocate(zeta)
    if (allocated(r)) deallocate(r)

    print *, "Failed tests: ", trim(error_message)
    if (len_trim(error_message) > 0) then
        stop 1
    else
        stop 0
    end if

contains

function random_string(length) result(str)

    implicit none

    integer, intent(in) :: length
    character(len=length) :: str
    integer :: i
    real*8 :: r
    character(len=29) :: letters = 'abcdefghijklmnopqrstuvwxyz148'

    do i = 1, length
        ! Generate random number between 1 and 29
        call random_number(r)
        r = r * 29 + 1
        str(i:i) = letters(int(r):int(r))
    end do

end function random_string

subroutine initialize

    implicit none
    logical :: data_dir_exists

    ! Initialize random number generator
    call random_seed()

    ! Check if data directory exists, and create it if not

    inquire(file='data', exist=data_dir_exists)
    if (.not. data_dir_exists) then
        ! Create directory TODO: make it cross-platform
        call system('mkdir data')
    end if

end subroutine initialize

subroutine run_random

    implicit none

    ! Generate random points
    call random_number(rmu)
    rmu = rmu * 2.0d0 - 1.0d0
    
    r = sqrt(sum(rmu**2, dim=1))
    
    ! Randomize parameters
    call random_number(dd)

    i0 = 0
    indtmin = 0
    typec = 0
    indpar = 0
    indorb = 0
    indshell = 0
    z = 0.0d0
    
    call makefun(iorbs(ii),indt,i0,indtmin,indtm,typec,indpar      &
               &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
               &,iflagnorm_unused,cr)

    print *, "indorb = ", indorb
    print *, "indshell = ", indshell
    print *, "indpar = ", indpar

end subroutine run_random

subroutine check_index_movement

    implicit none

    write(*, '(A)') 'Test: index movement'

    failed = .false.

    ! Generate random points
    call random_number(rmu)
    rmu = rmu * 2.0d0 - 1.0d0
    
    r = sqrt(sum(rmu**2, dim=1))
    
    ! Randomize parameters
    call random_number(dd)

    i0 = 0
    indtmin = 0
    typec = 0
    indpar = 0
    indorb = 0
    indshell = 0
    z = 0.0d0
    
    call makefun(iorbs(ii),indt,i0,indtmin,indtm,typec,indpar      &
               &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
               &,iflagnorm_unused,cr)

    write(*,'("indorb   = ",I3)', advance="no") indorb
    if (indorb /= multiplicities(ii)) then
        write(*,'(" FAILED")')
        failed = .true.
    else
        write(*,'(" OK")')
    end if

    write(*,'("indshell = ",I3)', advance="no") indshell
    if (indshell /= multiplicities(ii)) then
        write(*,'(" FAILED")')
        failed = .true.
    else
        write(*,'(" OK")')
    end if

    write(*,'("indpar   = ",I3)', advance="no") indpar
    if (indpar /= npars(ii)) then
        write(*,'(" FAILED")')
        failed = .true.
    else
        write(*,'(" OK")')
    end if

    write(*, *)

end subroutine check_index_movement

subroutine create_pa_value

    implicit none
    character*80 :: filename

    ! Generate random points
    call random_number(rmu)
    rmu = rmu * 4.0d0 - 2.0d0
    
    r = sqrt(sum(rmu**2, dim=1))
    
    ! Randomize parameters
    call random_number(dd)

    i0 = 0
    indtmin = 0
    indpar = 0
    indorb = 0
    indshell = 0
    typec = 1
    z = 0.0d0

    call makefun(iorbs(ii),indt,i0,indtmin,indtm,typec,indpar      &
               &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
               &,iflagnorm_unused,cr)

    ! In a binary file data/sv.dat we store the following data:
    ! 1. iorbs(ii)
    ! 2. dd array
    ! 3. rmu array

    write(filename, '(A,I3.3,A,A,A,A)') 'data/pa_', iorbs(ii), "_", random_string(15), '.dat'
    open(unit=20, file=filename, form='unformatted', status='replace', action='write')
    write(20) dd
    write(20) rmu(:,:)
    write(20) z(:, 0:indtm)
    close(20)

end subroutine create_pa_value

subroutine create_svgl_value

    implicit none
    character*80 :: filename

    ! Generate random points
    call random_number(rmu)
    rmu = rmu * 4.0d0 - 2.0d0
    
    r = sqrt(sum(rmu**2, dim=1))
    
    ! Randomize parameters
    call random_number(dd)

    i0 = 0
    indtmin = 0
    indpar = 0
    indorb = 0
    indshell = 0
    typec = 0
    z = 0.0d0

    call makefun(iorbs(ii),0,i0,indtmin,0,typec,indpar      &
               &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
               &,iflagnorm_unused,cr)

    ! In a binary file data/sv.dat we store the following data:
    ! 1. iorbs(ii)
    ! 2. dd array
    ! 3. rmu array

    write(filename, '(A,I3.3,A,A,A,A)') 'data/svgl_', iorbs(ii), "_", random_string(15), '.dat'
    open(unit=20, file=filename, form='unformatted', status='replace', action='write')
    write(20) dd
    write(20) rmu(:,0)
    write(20) z(:, 0:3)
    close(20)

end subroutine create_svgl_value

subroutine create_single_value

    implicit none
    character*80 :: filename

    ! Generate random points
    call random_number(rmu)
    rmu = rmu * 2.0d0 - 1.0d0
    
    r = sqrt(sum(rmu**2, dim=1))
    
    ! Randomize parameters
    call random_number(dd)

    i0 = 0
    indtmin = 0
    indpar = 0
    indorb = 0
    indshell = 0
    typec = 1
    z = 0.0d0

    call makefun(iorbs(ii),indt,i0,indtmin,0,typec,indpar      &
               &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
               &,iflagnorm_unused,cr)

    ! In a binary file data/sv.dat we store the following data:
    ! 1. iorbs(ii)
    ! 2. dd array
    ! 3. rmu array

    write(filename, '(A,I3.3,A,A,A,A)') 'data/sv_', iorbs(ii), "_", random_string(15), '.dat'
    open(unit=20, file=filename, form='unformatted', status='replace', action='write')
    write(20) dd
    write(20) rmu(:,0)
    write(20) z(:, 0)
    close(20)

end subroutine create_single_value

subroutine check_single_value

    implicit none

    character :: files(MAX_FILE_LENGTH, MAX_FILES)
    character*MAX_FILE_LENGTH :: filename
    integer*4 :: num_files, i, j, iorb_test, ios, count_tests, count_failed
    real*8 :: z_test(multiplicities(ii), 0:0)

    write(*, '(A)') 'Test: single value'

    call list_directory('data', 4, files, num_files, 'sv_', 3)
    
    count_tests = 0
    count_failed = 0
    do i = 1, num_files
        do j = 1, MAX_FILE_LENGTH
            filename(j:j) = files(j,i)
        end do
        
        read(filename(4:6), '(I3)', iostat=ios ) iorb_test
        if (ios /= 0) then
            cycle
        end if
        if (iorb_test /= iorbs(ii)) then
            cycle
        else
        end if

        open(unit=20, file='data/'//trim(filename), form='unformatted', status='old', action='read')
        read(20) dd
        read(20) rmu(:,0)
        read(20) z_test(:, 0)
        close(20)

        i0 = 0
        indtmin = 0
        typec = 1
        indpar = 0
        indorb = 0
        indshell = 0
        r = sqrt(sum(rmu**2, dim=1))
        z = 0.0d0

        call makefun(iorbs(ii),indt,i0,indtmin,0,typec,indpar      &
                   &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
                   &,iflagnorm_unused,cr)

        if (any(abs(z_test - z) > 1.0d-10)) then
            failed = .true.
            count_failed = count_failed + 1
        end if
        count_tests = count_tests + 1

    end do

    write(*, '("tests/failed = ",I3, "/",I3)') count_tests, count_failed
    write(*, *)

end subroutine check_single_value

subroutine check_pa_value

    implicit none

    character :: files(MAX_FILE_LENGTH, MAX_FILES)
    character*MAX_FILE_LENGTH :: filename
    integer*4 :: num_files, i, j, iorb_test, ios, count_tests, count_failed
    real*8 :: z_test(multiplicities(ii), 0:indtm)
    character(*), parameter :: prefix = 'pa_'

    write(*, '(A)') 'Test: calculate values for pseudo average'

    call list_directory('data', 4, files, num_files, prefix, len_trim(prefix))
    
    count_tests = 0
    count_failed = 0
    do i = 1, num_files
        do j = 1, MAX_FILE_LENGTH
            filename(j:j) = files(j,i)
        end do
        
        read(filename(len_trim(prefix)+1:len_trim(prefix)+3), '(I3)', iostat=ios ) iorb_test
        if (ios /= 0) then
            cycle
        end if
        if (iorb_test /= iorbs(ii)) then
            cycle
        else
        end if

        open(unit=20, file='data/'//trim(filename), form='unformatted', status='old', action='read')
        read(20) dd
        read(20) rmu(:,:)
        read(20) z_test(:, 0:indtm)
        close(20)

        i0 = 0
        indtmin = 0
        typec = 1
        indpar = 0
        indorb = 0
        indshell = 0
        r = sqrt(sum(rmu**2, dim=1))
        z = 0.0d0

        call makefun(iorbs(ii),indt,i0,indtmin,indtm,typec,indpar      &
                   &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
                   &,iflagnorm_unused,cr)

        ! Check if values are the same
        ! Check if other values remained untouched
        if (     any(abs(z_test(:,0:indtm) - z(:,0:indtm)) > 1.0d-10) &
            .or. any(abs(z(:,indtm+1:indt)) > 10d-10)) then
            failed = .true.
            count_failed = count_failed + 1
        end if
        count_tests = count_tests + 1

    end do

    write(*, '("tests/failed = ",I3, "/",I3)') count_tests, count_failed
    write(*, *)

end subroutine check_pa_value

subroutine check_svgl_value

    implicit none

    character :: files(MAX_FILE_LENGTH, MAX_FILES)
    character*MAX_FILE_LENGTH :: filename
    integer*4 :: num_files, i, j, iorb_test, ios, count_tests, count_failed
    real*8 :: z_test(multiplicities(ii), 0:3)
    character(*), parameter :: prefix = 'svgl_'

    write(*, '(A)') 'Test: single value gradient and laplacian'

    call list_directory('data', 4, files, num_files, prefix, len_trim(prefix))
    
    count_tests = 0
    count_failed = 0
    do i = 1, num_files
        do j = 1, MAX_FILE_LENGTH
            filename(j:j) = files(j,i)
        end do
        
        read(filename(len_trim(prefix)+1:len_trim(prefix)+3), '(I3)', iostat=ios ) iorb_test
        if (ios /= 0) then
            cycle
        end if
        if (iorb_test /= iorbs(ii)) then
            cycle
        else
        end if

        open(unit=20, file='data/'//trim(filename), form='unformatted', status='old', action='read')
        read(20) dd
        read(20) rmu(:,0)
        read(20) z_test(:, 0:3)
        close(20)

        i0 = 0
        indtmin = 0
        typec = 0
        indpar = 0
        indorb = 0
        indshell = 0
        r = sqrt(sum(rmu**2, dim=1))
        z = 0.0d0

        call makefun(iorbs(ii),0,i0,indtmin,0,typec,indpar      &
                   &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
                   &,iflagnorm_unused,cr)

        if (any(abs(z_test(:,0:3) - z(:,0:3)) > 1.0d-10)) then
            failed = .true.
            count_failed = count_failed + 1
        end if
        count_tests = count_tests + 1

    end do

    write(*, '("tests/failed = ",I3, "/",I3)') count_tests, count_failed
    write(*, *)

end subroutine check_svgl_value

end program makefun_tester
