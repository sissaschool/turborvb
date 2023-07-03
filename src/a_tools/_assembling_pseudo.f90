program assembling_pseudo
    use allio
    implicit none
    real*8 atomic_sub
    real*8 ps_par1, ps_par2, ps_par3, ps_cut, dist_kel(3)
    character(20) extension
    character(220) file_pseudo
    character(100) dir_pseudo
    character(100) name_dir
    integer unit_number, line_to_read, num_sh, i, j, jj, fake_num, ios
    integer, dimension(:), allocatable :: shwork
!   AAA    Lines to be added just after all definitions of variables.
    character(60) name_tool
    character(20) str

    call getarg(1, str)
    if (str .eq. "--help" .or. str .eq. "-help" .or. str .eq. "help") then

!          Input the name of the file exactly as it is in /doc
        name_tool = 'assembling_pseudo'
        call help_online(name_tool)

        stop
    end if
!    AAA   end lines to be added

    name_dir = 'NAME_DIR'

    call default_allocate
    open (unit=10, file='fort.10', status='old', form='formatted')
    call read_fort10(10)

!  ##########################################################
!       write(6,*) ' Input directory pseudo e.g. /home/trunk/pseudo/ '
!                                  !path of the dir for the PP files
!       read(5,'(A70)') dir_pseudo
    dir_pseudo = trim(name_dir)//'/pseudo/'

    write (6, *) ' dir_pseudo =', dir_pseudo
    write (6, *) ' Input file extension (e.g. filippi) '
    !nome_file.extension
    read (5, '(A20)') extension

    write (6, *) ' Extension read =', extension

    open (unit=33, file='pseudo.dat', status='unknown', form='formatted')

    write (33, *) 'ECP'

    do j = 1, nion
        if (zetar(j) .ne. atom_number(j) .and. atom_number(j) .gt. 0) then

            call open_file_pseudo(int(zetar(j)), int(atom_number(j))      &
      &                        , extension, file_pseudo)

!          write(6,*) ' file pseudo after =',file_pseudo

            unit_number = 100
            file_pseudo = trim(dir_pseudo)//trim(file_pseudo)
            open (unit=unit_number, file=file_pseudo, status='old'          &
      &, form='formatted', iostat=ios)

            if (ios .ne. 0) then

                if (trim(extension) .ne. 'distance') then

                    write (*, *) '##################################'
                    write (*, *) 'Error opening the pseudo potential file'
                    write (*, *) 'File  ', trim(file_pseudo), ' not found'
                    write (*, *) 'Program ends'

                    write (*, *) '##################################'
                end if

            else

                read (unit_number, *) pseudoname
                read (unit_number, *) fake_num, ps_cut, num_sh

                write (33, *) j, ps_cut, num_sh

                allocate (shwork(num_sh))
                read (unit_number, *) (shwork(jj), jj=1, num_sh)
                write (33, *) (shwork(jj), jj=1, num_sh)
                line_to_read = sum(shwork)
                deallocate (shwork)

                do jj = 1, line_to_read
                    read (unit_number, *) ps_par1, ps_par2, ps_par3
                    write (33, *) ps_par1, ps_par2, ps_par3
                end do
            end if

            close (unit_number)
        end if
    end do

    if (trim(extension) .eq. 'distance') then
        write (6, *) ' Distance atoms #i,#j,val distance '
        do i = 1, nion
            do j = i, nion
                dist_kel(:) = rion(:, i) - rion(:, j)
                if (iespbc) call ApplyPBC(dist_kel, 1)
                write (6, *) i, j, dsqrt(sum(dist_kel(:)**2))
            end do
        end do
    end if

    close (10)
    stop

end program assembling_pseudo

subroutine open_file_pseudo(z, atom, extension, file_pseudo)
    implicit none
    integer z, atom
    character(6) char_z, char_atom
    character(20) extension
    character(220) file_pseudo

    if (z < 10) then
        write (char_z, FMT="(I1)") z
    else if (z < 100) then
        write (char_z, FMT="(I2)") z
    else if (z < 1000) then
        write (char_z, FMT="(I3)") z
    end if

    if (atom < 10) then
        write (char_atom, FMT="(I1)") atom
    else if (atom < 100) then
        write (char_atom, FMT="(I2)") atom
    else if (atom < 1000) then
        write (char_atom, FMT="(I3)") atom
    end if

    file_pseudo = 'Z'//trim(char_z)//'_atomnumber'//trim(char_atom)
    file_pseudo = trim(file_pseudo)//'.'//trim(extension)

end
