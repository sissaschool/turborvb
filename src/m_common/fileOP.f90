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

subroutine check_scratch(rank, path, scratchpath)
    use allio, only: oldscra, wherescratch, chara, charaq, rankrep, rankcolrep&
            &, yesquantum, io_level, manyfort10, yesdft
    use io_m
    implicit none
    character(*) :: path, scratchpath
    integer :: iflagerr, rank, idigit(6)

    if (oldscra) then
        scratchpath = './'
        wherescratch = './scratch'
    else
        if (trim(wherescratch) .ne. './' .and. trim(wherescratch) .ne. trim(path)) then
            !     check if path wherescratch exist in each processor
            !     but only if wherescratch is not the current directory that exist by def.
            !     inquire(directory=trim(wherescratch),exist=dir_exist)
#ifndef __KCOMP

            call cd_dir(trim(wherescratch))
            call get_dir(scratchpath)
            call cd_dir(trim(path))

            iflagerr = 0
            if (trim(scratchpath) .eq. trim(path)) iflagerr = 1
            call checkiflagerr(iflagerr, rank, 'ERROR the directory '//trim(wherescratch)// &
                    &' does not exist at least in some node, pls. check wherescratch ')
#endif

        end if

#ifdef __KCOMP
        scratchpath = "./"
#else
        if (io_level .ne. 0) then
            scratchpath = trim(wherescratch)//"turborvb.scratch/"
            call mk_dir(trim(scratchpath))
        else
            scratchpath = "./"
        end if
#endif
        wherescratch = trim(scratchpath)//"tmp"
    end if

    if (rank .eq. 0) write (6, '(a,a)') '  Scratch extension ', trim(wherescratch)
    if (len_trim(scratchpath) .ge. lchlen) then
        if (rank .eq. 0) write (6, *) ' Warning your probably truncated, please &
       &  try to increase lchlen in mod_IO'
    end if

    call convertdec(rank, idigit)
    chara = char(idigit(1))//char(idigit(2))//char(idigit(3))// &
            &char(idigit(4))//char(idigit(5))//char(idigit(6))
    !  charaq defined always as does not cost anything
#ifdef __KCOMP
    call convertdec(rank, idigit)
#else
    call convertdec(rankcolrep, idigit)
#endif
    charaq = char(idigit(1))//char(idigit(2))//char(idigit(3))// &
             char(idigit(4))//char(idigit(5))//char(idigit(6))

end subroutine check_scratch

subroutine open_files(rank, scratchpath)
    use allio, only: iopt, itestr, ieskint, idyn, iespbc, ncg_adr, npower, npowersz&
            &, iread, wherescratch, writescratch, kl, chara, write_cov&
            &, yesquantum, rankrep, charaq, yeswrite10, io_level&
            &, kel, angle, indt, kelcont, commcolrep_mpi, quantcont&
            &, details_SP, details_DP, decoupled_run, molyes
    use allio, only: manyfort10
    use io_m
    use mpiio

    implicit none

    integer rank
    character(lchlen) :: scratchpath
    character(lchlen + 80) :: charascratch
    character(6) :: positi
    character(9) :: message

#ifdef PARALLEL
    integer(kind=MPI_OFFSET_KIND) :: disp_zero = 0
#endif
    if (iopt .eq. 0 .or. iopt .eq. 3) then
        positi = 'APPEND'
        message = 'appended!'
    else
        positi = 'REWIND'
        message = 'rewinded!'
    end if

    if (rank .eq. 0) then
        open (unit=7, file='stop.dat', form='formatted', status='unknown')
        close (unit=7, status='DELETE')
        !   open(unit=10,file='fort.10',form='formatted',position='REWIND')
        open (unit=11, file='fort.11', form='unformatted', position='REWIND')
        !     it contains energy
        open (unit=12, file='fort.12', form='unformatted', position=positi)
        !   they contains ionic forces and positions
        if (itestr .eq. -5) then
            write (6, *) 'Energy files ', message
            open (unit=16, file='forces.dat', form='formatted', position=positi)
            if (write_cov) then
                open (unit=18, file='covmat.dat', form='formatted', position=positi)
            end if
            open (unit=13, file='fort.12.fn', form='unformatted', position=positi)
            if (idyn .gt. 1) then
                write (6, *) 'Velocity file ', message
                open (unit=22, file='velocity.dat', form='formatted', position=positi)
            end if
            if (ieskint .ne. 0) then
                write (6, *) 'forces and position files ', message
                open (unit=17, file='position.dat', form='formatted', position=positi)
                !       if(yesquantum) then
                !         write(6,*) 'Beads position files (PIMD) ',message
                !         open(unit=27,file='beadposition.dat',form='formatted',position=positi)
                !       endif
            end if
            if (iespbc) then
                write (*, *) ' Pressure file ', message
                open (unit=23, file='pressure.dat', form='formatted', position=positi)
            end if
            if (ncg_adr .gt. 0 .or. npower + npowersz .gt. 0) then
                write (*, *) ' parametrization file ', message
                open (unit=25, file='parametrization.dat', form='formatted', position=positi)
            end if

        end if !itestr.eq.-5

        if (iread .ge. 6) then
            open (unit=15, file='fort.12.new', form='unformatted', position=positi)
        end if

    end if !rank.eq.0

    if (io_level .ne. 0) then
#ifdef PARALLEL
        if (io_level .eq. 1) then
            ! POSIX IO
            charascratch = trim(scratchpath)//'kelcont.'//trim(chara)
            open (unit=9, file=charascratch, form='unformatted', position='rewind')
        else if (io_level .eq. 2) then
            ! MPI-IO
            if (iopt .eq. 1) then
                call mpiio_file_open(MPI_COMM_WORLD, trim(scratchpath)//'kelcont.all', MPI_MODE_WRONLY + MPI_MODE_CREATE, kelcont)
            else
                call mpiio_file_open(MPI_COMM_WORLD, trim(scratchpath)//'kelcont.all', MPI_MODE_RDWR, kelcont)
            end if
        end if
        if (iread .eq. 2 .or. iread .eq. 3) then
            if (io_level .eq. 1) then
                ! POSIX IO
                charascratch = trim(scratchpath)//'details.'//trim(chara)
                open (unit=15, file=charascratch, form='unformatted', position=positi)
            else if (io_level .eq. 2) then
                ! MPI-IO
                if (iopt .eq. 0 .or. iopt .eq. 3) then
                    call mpiio_file_open(MPI_COMM_WORLD, trim(scratchpath)//'details_SP.all', &
                                         MPI_MODE_WRONLY + MPI_MODE_APPEND, details_SP)
                    call mpiio_file_open(MPI_COMM_WORLD, trim(scratchpath)//'details_DP.all', &
                                         MPI_MODE_WRONLY + MPI_MODE_APPEND, details_DP)
                else
                    call mpiio_file_open(MPI_COMM_WORLD, trim(scratchpath)//'details_SP.all', &
                                         MPI_MODE_WRONLY + MPI_MODE_CREATE, details_SP)
                    call mpiio_file_open(MPI_COMM_WORLD, trim(scratchpath)//'details_DP.all', &
                                         MPI_MODE_WRONLY + MPI_MODE_CREATE, details_DP)
                    call mpiio_file_set_zero(details_SP)
                    call mpiio_file_set_zero(details_DP)
                end if
            end if
        end if
        if (itestr .eq. -5) then
!   definition name file scratch
            charascratch = trim(wherescratch)//'conf.'//trim(chara)
            if (writescratch .eq. 0) then
                open (unit=19, file=charascratch, form='unformatted', position='rewind')
            end if
            if (abs(kl) .eq. 9) then
                charascratch = trim(wherescratch)//'lan.'//trim(chara)
                open (unit=20, file=charascratch, form='unformatted', position='rewind')
            end if
        end if
#ifdef __KCOMP
        charascratch = trim(scratchpath)//'quantcont.'//trim(chara)
#else
        charascratch = trim(scratchpath)//'quantcont.'//trim(charaq)
#endif
        if (yesquantum) then
            if (rankrep .eq. 0) then
                if (io_level .eq. 1) then
                    open (unit=8, file=charascratch, form='unformatted', position='rewind')
                elseif (io_level .eq. 2) then
                    if (iopt .eq. 1) then
                        call mpiio_file_open(commcolrep_mpi, trim(scratchpath)//'quantcont.all', &
                                             MPI_MODE_WRONLY + MPI_MODE_CREATE, quantcont)
                    else
                        call mpiio_file_open(commcolrep_mpi, trim(scratchpath)//'quantcont.all', &
                                             MPI_MODE_RDWR, quantcont)
                    end if
                    call mpiio_file_create_view(quantcont, 512, MPI_DOUBLE_PRECISION) ! using 512 as block size for simplicity
                    quantcont%disp = 0
                    call mpiio_file_reset_view(quantcont)
                end if
#ifdef __KCOMP
                charascratch = trim(scratchpath)//'fort.10_'//trim(chara)
#else
                charascratch = trim(scratchpath)//'fort.10_'//trim(charaq)
#endif
                if (yeswrite10) open (unit=27, file=charascratch, form='formatted', position='rewind')
#ifdef __KCOMP
            elseif (io_level .eq. 1) then
!   Open in any event an empty file
                open (unit=8, file=charascratch, form='unformatted', position='rewind')
                charascratch = trim(scratchpath)//'fort.10_'//trim(chara)
                if (yeswrite10) open (unit=27, file=charascratch, form='formatted', position='rewind')
#endif
            end if ! endif rankrep==0
        end if !endif yesquantum
        !
        ! open DFT wavefunction for each k-point!
        !
        if (manyfort10 .and. molyes) then
            charascratch = trim(scratchpath)//'fort.10_'//trim(charaq)
            if (rankrep .eq. 0) then
!         write(6,*) '# fort10name: ',trim(charascratch),' rank: ',rank
                open (unit=27, file=charascratch, form='formatted', position='rewind', err=101)
#ifdef __KCOMP
            else
                open (unit=27, file=charascratch, form='formatted', position='rewind', err=101)
#endif
            end if
        end if
#else
        ! serial
        if (iread .eq. 2 .or. iread .eq. 3) then
            open (unit=15, file='fort.12.new', form='unformatted', position=positi)
        end if
        if (itestr .eq. -5) then
            charascratch = trim(wherescratch)//'.conf'
            if (writescratch .eq. 0) then
                open (unit=19, file=charascratch, form='unformatted', position='rewind')
            end if
            if (abs(kl) .eq. 9) then
                charascratch = trim(wherescratch)//'.lan'
                open (unit=20, file=charascratch, form='unformatted', position='rewind')
            end if
        end if
#endif
    end if
    return
101 call error(' open_files ', ' fort.10 for each k-point needed!!! ', 1, rank)
end subroutine open_files

subroutine close_files(rank)
    use allio, only: itestr, ieskint, idyn, iespbc, ncg_adr, iread, writescratch&
            &, kl, write_cov, yesquantum, rankrep, yeswrite10, io_level, kelcont&
            &, quantcont, details_SP, details_DP, decoupled_run, molyes
    use allio, only: manyfort10
    use mpiio
    implicit none
    integer :: rank

    if (rank .eq. 0) then
        ! close(6)
        close (10)
        close (11)
        close (12)

        if (itestr .eq. -5) then
            close (13)
            close (16)
            if (write_cov) close (18)
            if (ieskint .ne. 0) close (17)
            if (idyn .gt. 1) close (22)
            if (iespbc) close (23)
            if (ncg_adr .gt. 0) close (25)
        end if
        if (iread .ge. 6) close (15)
    end if !rank.eq.0
    if (io_level .ne. 0) then
#ifdef PARALLEL
        if (io_level .eq. 1) then
            close (9)
        else if (io_level .eq. 2) then
            call mpiio_file_close(kelcont)
        end if
        if (rankrep .eq. 0) then
            if (yesquantum) then
                if (io_level .eq. 1) then
                    close (8)
                elseif (io_level .eq. 2) then
                    call mpiio_file_close(quantcont)
                end if
            end if
            if ((yesquantum .and. yeswrite10) .or. (manyfort10 .and. molyes)) close (27)
        end if
#endif
        if (iread .eq. 2 .or. iread .eq. 3) then
            if (io_level .eq. 1) then
                close (15)
#ifndef PARALLEL
            end if
#else
        else if (io_level .eq. 2) then
            call mpiio_file_close(details_SP)
            call mpiio_file_close(details_DP)
        end if
#endif
    end if

    if (itestr .eq. -5) then
        if (writescratch .eq. 0) close (unit=19, status='DELETE')
        if (abs(kl) .eq. 9) close (unit=20, status='DELETE')
    end if

    end if
end subroutine close_files
